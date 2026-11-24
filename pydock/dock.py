"""Top-level orchestration: prepare, dock, cluster, write."""
from __future__ import annotations

import os
import time
from dataclasses import dataclass, field

import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem  # noqa: F401

from . import cavity as cavmod
from . import grid as gridmod
from . import ligand as ligmod
from . import receptor as recmod
from . import scoring as scoremod
from .ga import Docker, GAConfig, Solution, cluster_solutions
from .geometry import symmetry_rmsd


@dataclass
class DockingConfig:
    receptor: str
    ligands: list[str]
    # Site definition: give one of these.
    reference_ligand: str | None = None
    centre: tuple[float, float, float] | None = None
    radius: float | None = None
    box_size: tuple[float, float, float] | None = None   # full widths, in A
    site_residues: list[tuple[str, int]] | None = None
    padding: float = 4.0

    scoring: str = "empirical"
    spacing: float = 0.4
    detect_cavity: bool = True
    cavity_spacing: float = 0.6
    cavity_probe: float = 1.4
    cavity_min_psp: int = 4
    rot_h_cutoff: float = 5.5
    n_runs: int = 10
    n_poses: int = 5
    cluster_rmsd: float = 1.5
    keep_waters: bool = False
    rotate_receptor_oh: bool = True
    rotate_ligand_polar_h: bool = True
    exclude_resnames: set[str] = field(default_factory=set)
    output_dir: str = "pydock_out"
    ga: GAConfig = field(default_factory=GAConfig)
    verbose: bool = True
    write_poses_to_disk: bool = True
    # progress(fraction_in_0_1, message) -- for web clients and long runs.
    progress: object = None
    # cancelled() -> bool, polled between GA runs so a client can stop a job.
    cancelled: object = None


@dataclass
class LigandResult:
    name: str
    solutions: list[Solution]
    clusters: list[list[Solution]]
    mol: Chem.Mol
    seconds: float
    rmsd_to_reference: float | None = None
    ligand: object = None          # the prepared Ligand, for pose serialisation
    pose_file: str | None = None   # where the poses were written, if anywhere


def resolve_site(cfg: DockingConfig, rec):
    # An explicit box wins over a reference ligand: when both are given the
    # user has chosen the site by hand and the reference is there to measure
    # the result against, not to place the box.
    ref = None
    if cfg.reference_ligand:
        ref = ligmod.load_ligand(cfg.reference_ligand)
        if cfg.centre is None and not cfg.site_residues:
            return recmod.site_from_ligand(
                ref.coords0[ref.heavy_idx], cfg.padding), ref
    if cfg.centre is not None:
        if cfg.box_size is not None:
            half = np.asarray(cfg.box_size, dtype=float) / 2.0
        else:
            half = np.full(3, float(cfg.radius if cfg.radius is not None else 10.0))
        return (np.asarray(cfg.centre, dtype=float),
                np.maximum(half, 3.0)), ref
    if cfg.site_residues:
        return recmod.site_from_residues(rec, cfg.site_residues, cfg.padding), ref
    raise ValueError("define the site with reference_ligand, centre and size, "
                     "or site_residues")


class JobCancelled(RuntimeError):
    """Raised when a client cancels a running job."""


def run_docking(cfg: DockingConfig) -> list[LigandResult]:
    if cfg.write_poses_to_disk:
        os.makedirs(cfg.output_dir, exist_ok=True)

    def log(*a, **k):
        line = " ".join(str(x) for x in a)
        if cfg.verbose:
            print(line, **k)
        if cfg.progress is not None and line.strip():
            cfg.progress(None, line)

    n_lig = max(1, len(cfg.ligands))

    def report(frac, message, ligand_index=0):
        if cfg.progress is None:
            return
        if frac is None:
            cfg.progress(None, message)
        else:
            cfg.progress((ligand_index + frac) / n_lig, message)

    def check_cancel():
        if cfg.cancelled is not None and cfg.cancelled():
            raise JobCancelled("cancelled by the client")

    t0 = time.time()
    rec = recmod.prepare_receptor(
        cfg.receptor, keep_waters=cfg.keep_waters,
        rotate_oh=cfg.rotate_receptor_oh,
        exclude_resnames=cfg.exclude_resnames)
    n_added = int(np.sum(rec.is_donor_h))
    log(f"receptor  {cfg.receptor}: {rec.n_atoms} atoms, "
        f"{n_added} polar H, {rec.rot_h_axis_b.size} rotatable OH/NH groups "
        f"({time.time() - t0:.1f}s)")

    (centre, half), ref = resolve_site(cfg, rec)
    log(f"site      centre {np.round(centre, 2)}  half-widths {np.round(half, 2)}")

    cav = None
    if cfg.detect_cavity:
        t1 = time.time()
        c_origin, c_shape = gridmod.make_lattice(centre, half, cfg.cavity_spacing)
        cav = cavmod.detect_cavity(
            rec, c_origin, c_shape, cfg.cavity_spacing,
            probe=cfg.cavity_probe, min_psp=cfg.cavity_min_psp,
            seed_point=centre, seed_radius=float(np.max(half)) + 2.0)
        box_vol = float(np.prod(2 * half))
        log(f"cavity    {cav.points.shape[0]} points, {cav.volume:.0f} A^3 "
            f"({100 * cav.volume / box_vol:.0f}% of the box)  "
            f"[{time.time() - t1:.1f}s]")
        if cav.points.shape[0] < 20:
            log("cavity    too small to trust; falling back to the whole box")
            cav = None

    # Expand every input file into its records: an SDF or SMILES file may hold
    # a whole library, and each molecule in it is docked separately.
    ligands: list = []
    dropped: list[tuple[str, str]] = []
    for path in cfg.ligands:
        check_cancel()
        batch = ligmod.load_ligands(
            path, failures=dropped, rotate_polar_h=cfg.rotate_ligand_polar_h)
        if len(batch) > 1:
            log(f"library   {path}: {len(batch)} molecules")
        ligands.extend(batch)
    for name, why in dropped:
        log(f"skipped   {name}: {why}")
    if not ligands:
        raise ValueError("no ligand could be prepared from the given files")
    n_lig = len(ligands)
    log(f"ligands   {n_lig} to dock")

    # The empirical grids depend only on the receptor and the box, so for a
    # library they are built once and shared.  The force-field maps are keyed
    # on ligand atom types, so those are rebuilt whenever the types change.
    shared_grid = None
    results: list[LigandResult] = []

    for lig_i, lig in enumerate(ligands):
        check_cancel()
        t1 = time.time()
        log(f"\nligand    {lig.name} [{lig_i + 1}/{n_lig}]: {lig.n_atoms} atoms, "
            f"{lig.n_rot} rotatable bonds ({lig.n_frozen_rot} heavy)")

        if shared_grid is not None:
            scorer = scoremod.Scorer(kind=shared_grid[0], grid=shared_grid[1],
                                     lig=lig, rec=rec, w_plp=shared_grid[2])
            reused = True
        else:
            scorer = scoremod.build_scorer(cfg.scoring, rec, lig, centre, half,
                                           cfg.spacing)
            reused = False
            if scorer.kind == "empirical" and n_lig > 1:
                shared_grid = (scorer.kind, scorer.grid, scorer.w_plp)
        if cav is not None:
            scorer.restrict_rot_h(cav.points, cfg.rot_h_cutoff)
        log(f"grids     {cfg.scoring}: {len(scorer.grid.map_names)} maps "
            f"{scorer.grid.shape} @ {cfg.spacing} A"
            f"{' (reused)' if reused else ''}, "
            f"{scorer.n_rot_h} receptor OH genes")

        ga_cfg = GAConfig(**{**cfg.ga.__dict__})
        ga_cfg.n_runs = cfg.n_runs
        docker = Docker(scorer, lig, centre, half, ga_cfg,
                        cavity_points=None if cav is None else cav.points)
        log(f"search    {ga_cfg.n_islands} islands x {ga_cfg.pop_size}, "
            f"{ga_cfg.n_evals} evals/run, {cfg.n_runs} runs, "
            f"{docker.n_genes} genes")

        def ga_progress(frac, message, _i=lig_i):
            check_cancel()
            report(frac, message, _i)

        sols = docker.dock(verbose=cfg.verbose, on_progress=ga_progress)
        clusters = cluster_solutions(sols, lig.heavy_idx, cfg.cluster_rmsd)

        rms = None
        if ref is not None and ref.n_atoms == lig.n_atoms:
            target = np.ascontiguousarray(ref.coords0[ref.heavy_idx])
            maps = np.ascontiguousarray(lig.sym_maps)
            for s in sols:
                s.rmsd = float(symmetry_rmsd(
                    np.ascontiguousarray(s.coords[lig.heavy_idx]), target, maps))
            rms = sols[0].rmsd

        dt = time.time() - t1
        res = LigandResult(name=lig.name, solutions=sols, clusters=clusters,
                           mol=lig.mol, seconds=dt, rmsd_to_reference=rms,
                           ligand=lig)
        results.append(res)

        log(f"best      {sols[0].score:.3f}   "
            f"{len(clusters)} distinct binding modes   {dt:.1f}s")
        if rms is not None:
            log(f"rmsd      {rms:.2f} A to the reference pose")
        if cfg.write_poses_to_disk:
            # Index the file name: library records often share a title.
            prefix = f"{lig_i + 1:03d}_" if n_lig > 1 else ""
            out = os.path.join(cfg.output_dir,
                               f"{prefix}{_safe(lig.name)}_poses.sdf")
            write_poses(out, lig, sols[:cfg.n_poses], cfg.scoring)
            res.pose_file = out
            log(f"written   {out}")

    if n_lig > 1:
        ranked = sorted(results, key=lambda r: r.solutions[0].score)
        log("\nscreening ranking (best score per ligand)")
        for rank, r in enumerate(ranked, start=1):
            log(f"  {rank:3d}. {r.name[:34]:<34} {r.solutions[0].score:9.3f}")
        if cfg.write_poses_to_disk:
            _write_ranking(os.path.join(cfg.output_dir, "ranking.csv"), ranked)

    return results


def ligand_efficiency(result: LigandResult) -> float:
    """Best score per heavy atom.

    Raw docking scores grow with molecular size, so a ranking by score alone
    rewards whatever is biggest.  Dividing by heavy-atom count is the standard
    correction and usually tells a different, more useful story in screening.
    """
    n = max(1, len(result.ligand.heavy_idx))
    return result.solutions[0].score / n


def _write_ranking(path: str, ranked: list[LigandResult]):
    """One row per ligand, best pose first -- the screening deliverable."""
    import csv

    terms: list[str] = []
    for r in ranked:
        for k in r.solutions[0].terms:
            if k != "total" and k not in terms:
                terms.append(k)
    with open(path, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["rank", "ligand", "score", "ligand_efficiency", "n_atoms",
                    "n_rotatable", "binding_modes", "rmsd_to_reference",
                    "seconds"] + terms)
        for rank, r in enumerate(ranked, start=1):
            s = r.solutions[0]
            w.writerow([rank, r.name, f"{s.score:.4f}",
                        f"{ligand_efficiency(r):.4f}", r.ligand.n_atoms,
                        r.ligand.n_frozen_rot, len(r.clusters),
                        "" if r.rmsd_to_reference is None
                        else f"{r.rmsd_to_reference:.2f}",
                        f"{r.seconds:.1f}"]
                       + [f"{float(s.terms.get(t, 0)):.4f}" for t in terms])


def poses_to_sdf(lig, solutions: list[Solution], kind: str) -> str:
    """Serialise poses to an SDF string instead of a file (for web responses)."""
    import io

    buf = io.StringIO()
    w = Chem.SDWriter(buf)
    for rank, s in enumerate(solutions, start=1):
        m = Chem.Mol(lig.mol)
        conf = m.GetConformer()
        for i in range(m.GetNumAtoms()):
            conf.SetAtomPosition(i, s.coords[i].tolist())
        m.SetProp("_Name", f"{lig.name}_pose{rank}")
        m.SetProp("pydock_rank", str(rank))
        m.SetProp("pydock_scoring", kind)
        m.SetProp("pydock_score", f"{s.score:.4f}")
        m.SetProp("pydock_run", str(s.run))
        for k, v in s.terms.items():
            if k != "total":
                m.SetProp(f"pydock_{k}", f"{float(v):.4f}")
        w.write(m)
    w.close()
    return buf.getvalue()


def _safe(name: str) -> str:
    return "".join(c if c.isalnum() or c in "-_." else "_" for c in name)[:60] or "ligand"


def write_poses(path: str, lig, solutions: list[Solution], kind: str):
    w = Chem.SDWriter(path)
    for rank, s in enumerate(solutions, start=1):
        m = Chem.Mol(lig.mol)
        conf = m.GetConformer()
        for i in range(m.GetNumAtoms()):
            conf.SetAtomPosition(i, s.coords[i].tolist())
        m.SetProp("_Name", f"{lig.name}_pose{rank}")
        m.SetProp("pydock_rank", str(rank))
        m.SetProp("pydock_scoring", kind)
        m.SetProp("pydock_score", f"{s.score:.4f}")
        m.SetProp("pydock_run", str(s.run))
        for k, v in s.terms.items():
            if k == "total":
                continue
            m.SetProp(f"pydock_{k}", f"{float(v):.4f}")
        w.write(m)
    w.close()
