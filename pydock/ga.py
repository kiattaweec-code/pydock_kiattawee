"""Island-model genetic algorithm with niching, plus local refinement.

The chromosome layout is

    [ tx ty tz | qw qx qy qz | torsion_1 .. torsion_T | recOH_1 .. recOH_R ]

which is the same decomposition GOLD uses: rigid-body placement, ligand
torsions, and rotatable receptor hydroxyls as genes of their own.

Search strategy follows GOLD's in shape rather than in detail: several
semi-isolated subpopulations evolve in parallel, a steady-state scheme
replaces one individual at a time, and *niching* stops a single binding mode
from taking over the population before alternatives have been explored.
"""
from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np
from numba import njit

from .geometry import FASTMATH, apply_pose, random_quaternion, rmsd


@njit(**FASTMATH)
def niche_target(child_pose, poses, heavy, fit, lo, n, cutoff, niche_size):
    """Pick the individual the child should replace.

    Members of the child's own niche (poses within `cutoff` RMSD) compete first,
    and only if that niche is already at capacity; otherwise the child competes
    against the whole subpopulation.  Returns the index of the weakest rival.

    Written as one kernel because calling an RMSD routine once per population
    member from Python costs more than every other part of a GA step combined.
    """
    count = 0
    worst_n = -1
    worst_n_fit = 1e30
    worst_a = -1
    worst_a_fit = 1e30
    nh = heavy.shape[0]
    for k in range(n):
        p = lo + k
        s = 0.0
        for t in range(nh):
            i = heavy[t]
            dx = child_pose[i, 0] - poses[p, i, 0]
            dy = child_pose[i, 1] - poses[p, i, 1]
            dz = child_pose[i, 2] - poses[p, i, 2]
            s += dx * dx + dy * dy + dz * dz
        r = np.sqrt(s / nh)
        if fit[p] < worst_a_fit:
            worst_a_fit = fit[p]
            worst_a = p
        if r < cutoff:
            count += 1
            if fit[p] < worst_n_fit:
                worst_n_fit = fit[p]
                worst_n = p
    if count >= niche_size and worst_n >= 0:
        return worst_n
    return worst_a


@dataclass
class GAConfig:
    n_islands: int = 5
    pop_size: int = 100            # individuals per island
    n_evals: int = 60_000          # score evaluations per run
    niche_size: int = 2            # max individuals allowed per niche
    niche_rmsd: float = 1.0        # A, niche radius
    p_crossover: float = 0.90
    p_mutate: float = 0.90
    migrate_every: int = 1_500     # evaluations between migrations
    n_runs: int = 10               # independent GA runs per ligand
    tournament: int = 3
    trans_sigma: float = 0.8       # A
    rot_sigma: float = 0.25        # quaternion perturbation
    tors_sigma: float = 0.45       # rad
    local_opt: bool = True
    local_opt_maxiter: int = 400
    early_stop: int = 12_000       # stop a run after this many evals with no gain
    seed: int = 1
    autoscale: bool = True         # scale the budget with ligand flexibility

    # Lamarckian local search (Solis-Wets), applied to a fraction of children.
    ls_frequency: float = 0.06
    ls_steps: int = 90
    ls_rho: float = 0.45
    ls_rho_min: float = 0.015
    ls_window: int = 4             # successes/failures before rho is adjusted

    # Cavity-constrained placement
    use_cavity: bool = True
    cavity_jump: float = 0.25      # probability a translation mutation relocates


@dataclass
class Solution:
    score: float
    chromosome: np.ndarray
    coords: np.ndarray
    terms: dict = field(default_factory=dict)
    run: int = 0
    rmsd: float | None = None      # to the reference ligand, when one is given


class Docker:
    """Runs the GA for one ligand against one prepared scorer."""

    def __init__(self, scorer, lig, centre, half_widths,
                 config: GAConfig | None = None, cavity_points=None):
        self.scorer = scorer
        self.lig = lig
        self.centre = np.asarray(centre, dtype=np.float64)
        self.half = np.asarray(half_widths, dtype=np.float64)
        self.cfg = config or GAConfig()
        self.cavity_points = (np.ascontiguousarray(cavity_points)
                              if cavity_points is not None
                              and len(cavity_points) and self.cfg.use_cavity
                              else None)

        self.n_tors = int(lig.n_rot)
        self.n_roth = int(scorer.n_rot_h)
        self.n_genes = 7 + self.n_tors + self.n_roth
        self._buf = lig.new_buffer()
        self.n_evals_used = 0
        # Optional progress sink: on_progress(fraction_in_0_1, message).
        self.on_progress = None
        self._run_index = 0
        self._n_runs = 1

        if self.cfg.autoscale:
            # GOLD scales effort with the number of degrees of freedom; so do we.
            dof = self.n_tors + 0.5 * self.n_roth
            self.cfg.n_evals = int(np.clip(8_000 + 4_000 * dof, 20_000, 250_000))
            self.cfg.early_stop = max(6_000, self.cfg.n_evals // 4)

    # ------------------------------------------------------------------ #
    # Chromosome handling
    # ------------------------------------------------------------------ #
    def decode(self, chrom: np.ndarray) -> np.ndarray:
        apply_pose(self.lig.coords0,
                   chrom[7:7 + self.n_tors],
                   self.lig.bond_a, self.lig.bond_b,
                   self.lig.move_ptr, self.lig.move_idx,
                   chrom[3:7], chrom[0:3], self._buf)
        return self._buf

    def evaluate(self, chrom: np.ndarray, breakdown: bool = False) -> float:
        xyz = self.decode(chrom)
        roth = chrom[7 + self.n_tors:] if self.n_roth else None
        self.n_evals_used += 1
        return self.scorer.score(xyz, roth, breakdown=breakdown)

    def random_position(self, rng) -> np.ndarray:
        """A starting centroid: a cavity point when one is known, else the box."""
        if self.cavity_points is not None:
            k = int(rng.integers(self.cavity_points.shape[0]))
            return self.cavity_points[k] + rng.normal(0.0, 0.4, 3)
        return self.centre + (rng.random(3) * 2 - 1) * (self.half * 0.7)

    def random_chromosome(self, rng) -> np.ndarray:
        c = np.empty(self.n_genes)
        c[0:3] = self.random_position(rng)
        c[3:7] = random_quaternion(rng)
        if self.n_tors:
            c[7:7 + self.n_tors] = rng.uniform(-np.pi, np.pi, self.n_tors)
        if self.n_roth:
            c[7 + self.n_tors:] = rng.uniform(-np.pi, np.pi, self.n_roth)
        return c

    def _clip(self, c: np.ndarray) -> np.ndarray:
        lo = self.centre - self.half
        hi = self.centre + self.half
        c[0:3] = np.clip(c[0:3], lo, hi)
        n = np.linalg.norm(c[3:7])
        if n < 1e-8:
            c[3:7] = np.array([1.0, 0.0, 0.0, 0.0])
        else:
            c[3:7] /= n
        if self.n_genes > 7:
            c[7:] = (c[7:] + np.pi) % (2 * np.pi) - np.pi
        return c

    # ------------------------------------------------------------------ #
    # Operators
    # ------------------------------------------------------------------ #
    def crossover(self, a: np.ndarray, b: np.ndarray, rng) -> np.ndarray:
        """Block-wise uniform crossover: position, orientation, each torsion."""
        child = a.copy()
        if rng.random() < 0.5:
            child[0:3] = b[0:3]
        if rng.random() < 0.5:
            child[3:7] = b[3:7]
        for k in range(7, self.n_genes):
            if rng.random() < 0.5:
                child[k] = b[k]
        return child

    def mutate(self, c: np.ndarray, rng) -> np.ndarray:
        cfg = self.cfg
        c = c.copy()
        r = rng.random()
        if r < 0.3:
            if self.cavity_points is not None and rng.random() < cfg.cavity_jump:
                c[0:3] = self.random_position(rng)
            else:
                c[0:3] += rng.normal(0.0, cfg.trans_sigma, 3)
        elif r < 0.6:
            c[3:7] += rng.normal(0.0, cfg.rot_sigma, 4)
        elif self.n_genes > 7:
            k = rng.integers(7, self.n_genes)
            if rng.random() < 0.25:
                c[k] = rng.uniform(-np.pi, np.pi)      # bold jump
            else:
                c[k] += rng.normal(0.0, cfg.tors_sigma)
        else:
            c[0:3] += rng.normal(0.0, cfg.trans_sigma, 3)
        return self._clip(c)

    # ------------------------------------------------------------------ #
    # One GA run
    # ------------------------------------------------------------------ #
    def run_once(self, rng) -> Solution:
        cfg = self.cfg
        n_isl = cfg.n_islands
        pop_per = cfg.pop_size
        total = n_isl * pop_per

        pop = np.empty((total, self.n_genes))
        fit = np.empty(total)                 # higher is better
        poses = np.empty((total, self.lig.n_atoms, 3))

        for i in range(total):
            pop[i] = self._clip(self.random_chromosome(rng))
            fit[i] = -self.evaluate(pop[i])
            poses[i] = self._buf

        best = int(np.argmax(fit))
        best_fit = fit[best]
        best_chrom = pop[best].copy()
        evals = total
        since_gain = 0
        next_report = cfg.n_evals // 20 or cfg.n_evals

        heavy = self.lig.heavy_idx

        while evals < cfg.n_evals and since_gain < cfg.early_stop:
            if self.on_progress is not None and evals >= next_report:
                next_report = evals + max(1, cfg.n_evals // 20)
                within = min(1.0, evals / cfg.n_evals)
                self.on_progress(
                    (self._run_index + within) / max(1, self._n_runs),
                    f"run {self._run_index + 1}/{self._n_runs}: "
                    f"{evals:,}/{cfg.n_evals:,} evaluations, "
                    f"best score {-best_fit:.3f}")

            isl = int(rng.integers(n_isl))
            lo = isl * pop_per
            hi = lo + pop_per

            p1 = self._tournament(fit, lo, hi, rng)
            p2 = self._tournament(fit, lo, hi, rng)
            child = pop[p1].copy()
            if rng.random() < cfg.p_crossover:
                child = self.crossover(pop[p1], pop[p2], rng)
            if rng.random() < cfg.p_mutate:
                child = self.mutate(child, rng)
            child = self._clip(child)

            f = -self.evaluate(child)
            evals += 1

            # --- Lamarckian step: improve the child, write the result back --
            if cfg.ls_frequency > 0 and rng.random() < cfg.ls_frequency:
                child, f, used = self.solis_wets(child, f, rng)
                evals += used

            self.decode(child)
            child_pose = self._buf.copy()

            # --- niching: who does this child compete with? ----------------
            worst = niche_target(child_pose, poses, heavy, fit, lo, pop_per,
                                 cfg.niche_rmsd, cfg.niche_size)
            if f > fit[worst]:
                pop[worst] = child
                fit[worst] = f
                poses[worst] = child_pose

            if f > best_fit:
                best_fit = f
                best_chrom = child.copy()
                since_gain = 0
            else:
                since_gain += 1

            if cfg.migrate_every and evals % cfg.migrate_every == 0 and n_isl > 1:
                self._migrate(pop, fit, poses, rng)

        # --- local refinement of the best pose ----------------------------
        if cfg.local_opt:
            best_chrom, best_fit = self._refine(best_chrom, best_fit)

        score = -best_fit
        self.evaluate(best_chrom, breakdown=True)
        return Solution(score=score, chromosome=best_chrom,
                        coords=self.decode(best_chrom).copy(),
                        terms=self.scorer.terms())

    def solis_wets(self, chrom, fitness, rng):
        """Solis-Wets adaptive random local search.

        Derivative-free, cheap, and biased towards the direction that has been
        working: on a successful step the bias vector moves towards the step
        taken, and the step size doubles after a run of successes and halves
        after a run of failures.  This is the 'Lamarckian' half of the search --
        the improved chromosome replaces the one the GA produced.
        """
        cfg = self.cfg
        n = self.n_genes
        scale = np.empty(n)
        scale[0:3] = 1.0
        scale[3:7] = 0.15
        if n > 7:
            scale[7:] = 0.5

        bias = np.zeros(n)
        rho = cfg.ls_rho
        best = chrom.copy()
        best_f = fitness          # fitness, higher is better
        succ = fail = 0
        used = 0

        for _ in range(cfg.ls_steps):
            if rho < cfg.ls_rho_min:
                break
            dev = rng.normal(0.0, rho, n) * scale
            cand = self._clip(best + bias + dev)
            f = -self.evaluate(cand)
            used += 1
            if f > best_f:
                best = cand
                best_f = f
                bias = 0.4 * bias + 0.6 * dev
                succ += 1
                fail = 0
            else:
                cand = self._clip(best - bias - dev)
                f = -self.evaluate(cand)
                used += 1
                if f > best_f:
                    best = cand
                    best_f = f
                    bias = bias - 0.4 * dev
                    succ += 1
                    fail = 0
                else:
                    bias *= 0.5
                    fail += 1
                    succ = 0
            if succ >= cfg.ls_window:
                rho *= 2.0
                succ = 0
            elif fail >= cfg.ls_window:
                rho *= 0.5
                fail = 0
        return best, best_f, used

    def _tournament(self, fit, lo, hi, rng) -> int:
        idx = rng.integers(lo, hi, self.cfg.tournament)
        return int(idx[np.argmax(fit[idx])])

    def _migrate(self, pop, fit, poses, rng):
        """Copy each island's best individual over a random neighbour."""
        n_isl = self.cfg.n_islands
        pop_per = self.cfg.pop_size
        src = int(rng.integers(n_isl))
        dst = (src + 1 + int(rng.integers(n_isl - 1))) % n_isl
        b = src * pop_per + int(np.argmax(fit[src * pop_per:(src + 1) * pop_per]))
        w = dst * pop_per + int(np.argmin(fit[dst * pop_per:(dst + 1) * pop_per]))
        pop[w] = pop[b]
        fit[w] = fit[b]
        poses[w] = poses[b]

    def _refine(self, chrom, fit):
        """Nelder-Mead polish in chromosome space."""
        from .optimize import nelder_mead

        def f(x):
            return self.evaluate(self._clip(x.copy()))

        res = nelder_mead(f, chrom.copy(), maxiter=self.cfg.local_opt_maxiter,
                          maxfev=self.cfg.local_opt_maxiter * 2,
                          xatol=1e-3, fatol=1e-4, adaptive=True)
        if -res.fun > fit:
            return self._clip(np.asarray(res.x)), -res.fun
        return chrom, fit

    # ------------------------------------------------------------------ #
    def dock(self, n_runs: int | None = None, verbose: bool = False,
             on_progress=None, on_run: bool = None) -> list[Solution]:
        cfg = self.cfg
        n_runs = n_runs or cfg.n_runs
        self.on_progress = on_progress
        self._n_runs = n_runs
        out: list[Solution] = []
        for r in range(n_runs):
            self._run_index = r
            sol = self.run_once(np.random.default_rng(cfg.seed + 1000 * r))
            sol.run = r + 1
            out.append(sol)
            if verbose:
                print(f"  run {r + 1:2d}/{n_runs}  score = {sol.score:9.3f}")
            if on_progress is not None:
                best_so_far = min(s.score for s in out)
                on_progress((r + 1) / n_runs,
                            f"run {r + 1}/{n_runs} finished: score "
                            f"{sol.score:.3f} (best {best_so_far:.3f})")
        out.sort(key=lambda s: s.score)
        return out


# --------------------------------------------------------------------------- #
# Pose clustering
# --------------------------------------------------------------------------- #
def cluster_solutions(solutions: list[Solution], heavy_idx: np.ndarray,
                      cutoff: float = 1.5) -> list[list[Solution]]:
    """Greedy single-linkage clustering of poses by heavy-atom RMSD."""
    clusters: list[list[Solution]] = []
    for s in sorted(solutions, key=lambda x: x.score):
        placed = False
        for cl in clusters:
            if rmsd(s.coords[heavy_idx], cl[0].coords[heavy_idx]) < cutoff:
                cl.append(s)
                placed = True
                break
        if not placed:
            clusters.append([s])
    return clusters
