"""The two scoring functions, behind one interface.

* ``empirical``  -- ChemScore-style free-energy terms (hydrogen bond, metal,
  lipophilic, rotor penalty) plus a piecewise-linear steric term.  This is the
  closer analogue of GOLD's ChemScore / ChemPLP.
* ``forcefield`` -- AutoDock4-style van der Waals, screened electrostatics,
  volume-based desolvation and a directional 12-10 hydrogen bond.

Both read gridded, type-dependent terms by interpolation and evaluate
hydrogen bonds and metal coordination pairwise, so full donor *and* acceptor
directionality is kept and rotatable receptor hydroxyls can move.

Lower is better for both.  ``Scorer.terms()`` returns the breakdown.
"""
from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np
from numba import njit

from . import grid as gridmod
from . import params
from .geometry import FASTMATH, rotate_about_axis

DEG = np.pi / 180.0
HB_CUTOFF = 4.0
METAL_CUTOFF = 3.2


# --------------------------------------------------------------------------- #
# Geometric weights for a hydrogen bond
# --------------------------------------------------------------------------- #
@njit(**FASTMATH)
def _ramp(x, full, zero):
    """1 below `full`, 0 above `zero`, linear in between."""
    if x <= full:
        return 1.0
    if x >= zero:
        return 0.0
    return (zero - x) / (zero - full)


# NOTE: these kernels take whole coordinate arrays plus integer indices rather
# than per-atom slices.  Slicing inside an njit loop allocates an array view
# per call, which dominated the runtime (150 ns/pair) before this rewrite.
@njit(**FASTMATH)
def _hb_geometry(hxyz, h, d, axyz, a, ante, ante_row):
    """Geometry of donor(d)-H(h) ... acceptor(a).

    `hxyz` holds the donor side, `axyz` the acceptor side; `ante` is the
    acceptor side's antecedent table and `ante_row` the row to read from it.
    Returns (r_H..A, geometric weight in [0, 1]).
    """
    dx = hxyz[h, 0] - axyz[a, 0]
    if dx > HB_CUTOFF or dx < -HB_CUTOFF:
        return 99.0, 0.0
    dy = hxyz[h, 1] - axyz[a, 1]
    if dy > HB_CUTOFF or dy < -HB_CUTOFF:
        return 99.0, 0.0
    dz = hxyz[h, 2] - axyz[a, 2]
    if dz > HB_CUTOFF or dz < -HB_CUTOFF:
        return 99.0, 0.0
    d2 = dx * dx + dy * dy + dz * dz
    if d2 > HB_CUTOFF * HB_CUTOFF:
        return 99.0, 0.0
    r = np.sqrt(d2)

    # Angle D-H...A at the hydrogen; ideal is 180 degrees.
    ax = hxyz[d, 0] - hxyz[h, 0]
    ay = hxyz[d, 1] - hxyz[h, 1]
    az = hxyz[d, 2] - hxyz[h, 2]
    bx = axyz[a, 0] - hxyz[h, 0]
    by = axyz[a, 1] - hxyz[h, 1]
    bz = axyz[a, 2] - hxyz[h, 2]
    na = np.sqrt(ax * ax + ay * ay + az * az)
    nb = np.sqrt(bx * bx + by * by + bz * bz)
    if na < 1e-8 or nb < 1e-8:
        return r, 0.0
    c = (ax * bx + ay * by + az * bz) / (na * nb)
    if c > 1.0:
        c = 1.0
    elif c < -1.0:
        c = -1.0
    alpha = 180.0 - np.arccos(c) / DEG
    w_ang = _ramp(alpha, params.CS_HB_DA1, params.CS_HB_DA2)
    if w_ang <= 0.0:
        return r, 0.0

    # Lone-pair directionality at the acceptor.
    n = 0
    mx = 0.0
    my = 0.0
    mz = 0.0
    for k in range(ante.shape[1]):
        j = ante[ante_row, k]
        if j < 0:
            continue
        mx += axyz[j, 0]
        my += axyz[j, 1]
        mz += axyz[j, 2]
        n += 1
    if n == 0:
        return r, w_ang
    mx /= n
    my /= n
    mz /= n
    lx = axyz[a, 0] - mx
    ly = axyz[a, 1] - my
    lz = axyz[a, 2] - mz
    hx = hxyz[h, 0] - axyz[a, 0]
    hy = hxyz[h, 1] - axyz[a, 1]
    hz = hxyz[h, 2] - axyz[a, 2]
    nl = np.sqrt(lx * lx + ly * ly + lz * lz)
    nh = np.sqrt(hx * hx + hy * hy + hz * hz)
    if nl < 1e-8 or nh < 1e-8:
        return r, w_ang
    c = (lx * hx + ly * hy + lz * hz) / (nl * nh)
    if c > 1.0:
        c = 1.0
    elif c < -1.0:
        c = -1.0
    w_acc = _ramp(np.arccos(c) / DEG, 60.0, 100.0)
    return r, w_ang * w_acc


# --------------------------------------------------------------------------- #
# Hydrogen bond / metal sums
# --------------------------------------------------------------------------- #
@njit(cache=True, nogil=True)
def hbond_chemscore(lig_xyz, rec_xyz,
                    lig_dh, lig_dh_parent, lig_acc, lig_acc_ante,
                    rec_acc, rec_acc_ante, rec_dh, rec_dh_parent):
    """Sum of block-function hydrogen bonds, in units of 'ideal H-bonds'."""
    total = 0.0
    for a in range(lig_dh.shape[0]):
        h = lig_dh[a]
        d = lig_dh_parent[h]
        for b in range(rec_acc.shape[0]):
            j = rec_acc[b]
            r, w = _hb_geometry(lig_xyz, h, d, rec_xyz, j, rec_acc_ante, j)
            if w <= 0.0:
                continue
            total += w * _ramp(abs(r - params.CS_HB_R_IDEAL),
                               params.CS_HB_DR1, params.CS_HB_DR2)
    for b in range(rec_dh.shape[0]):
        h = rec_dh[b]
        d = rec_dh_parent[h]
        for a in range(lig_acc.shape[0]):
            i = lig_acc[a]
            r, w = _hb_geometry(rec_xyz, h, d, lig_xyz, i, lig_acc_ante, i)
            if w <= 0.0:
                continue
            total += w * _ramp(abs(r - params.CS_HB_R_IDEAL),
                               params.CS_HB_DR1, params.CS_HB_DR2)
    return total


@njit(cache=True, nogil=True)
def hbond_12_10(lig_xyz, rec_xyz,
                lig_dh, lig_dh_parent, lig_acc, lig_acc_ante,
                rec_acc, rec_acc_ante, rec_dh, rec_dh_parent):
    """Directional 12-10 hydrogen bond energy (kcal/mol)."""
    total = 0.0
    R = params.AD_HB_R
    eps = params.AD_HB_EPS
    for a in range(lig_dh.shape[0]):
        h = lig_dh[a]
        d = lig_dh_parent[h]
        for b in range(rec_acc.shape[0]):
            j = rec_acc[b]
            r, w = _hb_geometry(lig_xyz, h, d, rec_xyz, j, rec_acc_ante, j)
            if w <= 0.0:
                continue
            if r < 0.6:
                r = 0.6
            x = R / r
            total += w * eps * (5.0 * x ** 12 - 6.0 * x ** 10)
    for b in range(rec_dh.shape[0]):
        h = rec_dh[b]
        d = rec_dh_parent[h]
        for a in range(lig_acc.shape[0]):
            i = lig_acc[a]
            r, w = _hb_geometry(rec_xyz, h, d, lig_xyz, i, lig_acc_ante, i)
            if w <= 0.0:
                continue
            if r < 0.6:
                r = 0.6
            x = R / r
            total += w * eps * (5.0 * x ** 12 - 6.0 * x ** 10)
    return total


@njit(cache=True, nogil=True)
def metal_term(lig_xyz, rec_xyz, lig_acc, rec_metal):
    total = 0.0
    for b in range(rec_metal.shape[0]):
        j = rec_metal[b]
        for a in range(lig_acc.shape[0]):
            i = lig_acc[a]
            dx = lig_xyz[i, 0] - rec_xyz[j, 0]
            dy = lig_xyz[i, 1] - rec_xyz[j, 1]
            dz = lig_xyz[i, 2] - rec_xyz[j, 2]
            r = np.sqrt(dx * dx + dy * dy + dz * dz)
            if r > METAL_CUTOFF:
                continue
            total += _ramp(r, params.CS_METAL_R1, params.CS_METAL_R2)
    return total


# --------------------------------------------------------------------------- #
# Internal (ligand strain)
# --------------------------------------------------------------------------- #
@njit(cache=True, nogil=True)
def internal_clash(xyz, pairs, r0):
    """Repulsion-only intramolecular strain, in A^2 of summed overlap.

    Attraction is deliberately absent: an intramolecular well would reward the
    ligand for folding onto itself, which is not what a docking search should
    be optimising.  `r0[p]` is the contact distance of pair p.
    """
    total = 0.0
    for p in range(pairs.shape[0]):
        i = pairs[p, 0]
        j = pairs[p, 1]
        dx = xyz[i, 0] - xyz[j, 0]
        dy = xyz[i, 1] - xyz[j, 1]
        dz = xyz[i, 2] - xyz[j, 2]
        d2 = dx * dx + dy * dy + dz * dz
        if d2 >= r0[p] * r0[p]:
            continue
        d = r0[p] - np.sqrt(d2)
        total += d * d
    return total


@njit(cache=True, nogil=True)
def internal_ff(xyz, pairs, rii, epsii, charge):
    total = 0.0
    for p in range(pairs.shape[0]):
        i = pairs[p, 0]
        j = pairs[p, 1]
        dx = xyz[i, 0] - xyz[j, 0]
        dy = xyz[i, 1] - xyz[j, 1]
        dz = xyz[i, 2] - xyz[j, 2]
        r = np.sqrt(dx * dx + dy * dy + dz * dz)
        if r > 8.0:
            continue
        if r < 0.5:
            r = 0.5
        rij = 0.5 * (rii[i] + rii[j])
        eps = np.sqrt(epsii[i] * epsii[j])
        total += params.AD_W_VDW * gridmod._vdw_12_6(r, rij, eps, params.AD_SMOOTH)
        total += params.AD_W_ELEC * 332.06 * charge[i] * charge[j] / (
            r * gridmod._dielectric(r))
    return total


# --------------------------------------------------------------------------- #
# Receptor hydroxyl rotation
# --------------------------------------------------------------------------- #
@njit(cache=True, nogil=True)
def apply_rot_h(ref_xyz, out_xyz, genes, axis_a, axis_b, h_atoms, touched):
    """Rotate each flagged receptor OH/NH group to the angle in `genes`.

    Only the handful of hydrogens listed in `touched` are restored from the
    reference and rotated; the rest of `out_xyz` is already correct and must
    not be copied, or this becomes the most expensive part of a score.
    """
    for k in range(touched.shape[0]):
        i = touched[k]
        out_xyz[i, 0] = ref_xyz[i, 0]
        out_xyz[i, 1] = ref_xyz[i, 1]
        out_xyz[i, 2] = ref_xyz[i, 2]
    origin = np.empty(3)
    axis = np.empty(3)
    buf = np.empty(3, dtype=np.int64)
    for g in range(genes.shape[0]):
        a = axis_a[g]
        b = axis_b[g]
        origin[0] = out_xyz[b, 0]
        origin[1] = out_xyz[b, 1]
        origin[2] = out_xyz[b, 2]
        axis[0] = out_xyz[b, 0] - out_xyz[a, 0]
        axis[1] = out_xyz[b, 1] - out_xyz[a, 1]
        axis[2] = out_xyz[b, 2] - out_xyz[a, 2]
        m = 0
        for k in range(h_atoms.shape[1]):
            if h_atoms[g, k] >= 0:
                buf[m] = h_atoms[g, k]
                m += 1
        if m == 0:
            continue
        rotate_about_axis(out_xyz, buf[:m], origin, axis, genes[g])


# --------------------------------------------------------------------------- #
# Scorer
# --------------------------------------------------------------------------- #
@dataclass
class Scorer:
    kind: str                       # "empirical" | "forcefield"
    grid: gridmod.Grid
    lig: object
    rec: object

    # cached selections
    lig_map_id: np.ndarray = None
    lig_scored: np.ndarray = None   # atoms read from the grid
    lig_lipo: np.ndarray = None
    lig_dh: np.ndarray = None
    lig_acc: np.ndarray = None
    rec_acc: np.ndarray = None
    rec_dh: np.ndarray = None
    rec_metal: np.ndarray = None
    out_of_box: float = 100.0
    n_rot_h: int = 0
    _rec_buf: np.ndarray = None
    _last: dict = field(default_factory=dict)

    # weights, exposed so they can be tuned without touching code
    w_clash: float = 1.0
    w_internal: float = 1.0
    w_plp: float = 0.0        # set > 0 for a ChemPLP-like shape term

    def __post_init__(self):
        lig, rec, g = self.lig, self.rec, self.grid
        n = lig.n_atoms

        if self.kind == "forcefield":
            types = [t for t in g.map_names if t.startswith("vdw:")]
            tmap = {t.split(":", 1)[1]: i for i, t in enumerate(types)}
            self.lig_map_id = np.asarray(
                [tmap.get(lig.ad_type[i], -1) for i in range(n)], dtype=np.int64)
            self.lig_scored = np.where(self.lig_map_id >= 0)[0].astype(np.int64)
            # One selection per vdW map, computed once rather than per call.
            self._vdw_sel = [
                np.ascontiguousarray(
                    self.lig_scored[self.lig_map_id[self.lig_scored] == k])
                for k in range(len(types))
            ]
            self._n_vdw = len(types)
            self._q = np.ascontiguousarray(lig.charge[self.lig_scored])
            self._s_i = np.ascontiguousarray(
                lig.solpar[self.lig_scored] + gridmod.QSOLPAR * np.abs(self._q))
            self._vol = np.ascontiguousarray(lig.vol[self.lig_scored])
        else:
            self.lig_map_id = np.full(n, -1, dtype=np.int64)
            self.lig_scored = np.asarray(
                [i for i in range(n) if lig.element[i] != "H"], dtype=np.int64)
            polar = np.asarray([lig.element[i] in ("N", "O") for i in range(n)])
            self._lig_p = np.ascontiguousarray(
                np.asarray([i for i in self.lig_scored if polar[i]], dtype=np.int64))
            self._lig_np = np.ascontiguousarray(
                np.asarray([i for i in self.lig_scored if not polar[i]], dtype=np.int64))

        self.lig_lipo = np.ascontiguousarray(
            np.where(lig.is_lipophilic)[0].astype(np.int64))
        self.lig_dh = np.where(lig.is_donor_h)[0].astype(np.int64)
        self.lig_acc = np.where(lig.is_acceptor)[0].astype(np.int64)
        self.lig_scored = np.ascontiguousarray(self.lig_scored)

        # Intramolecular contact distances from element radii.
        vdw = np.asarray([params.ELEMENT_VDW.get(lig.element[i], 1.7)
                          for i in range(n)])
        p = lig.intra_pairs
        self._intra_r0 = np.ascontiguousarray(
            params.INTERNAL_SCALE * (vdw[p[:, 0]] + vdw[p[:, 1]])
        ) if p.shape[0] else np.zeros(0)

        # Only receptor polar atoms close to the box can ever contribute.
        lo = g.origin - (HB_CUTOFF + 1.0)
        hi = g.origin + (np.asarray(g.shape) - 1) * g.spacing + HB_CUTOFF + 1.0
        near = np.all((rec.coords >= lo) & (rec.coords <= hi), axis=1)
        self.rec_acc = np.where(near & rec.is_acceptor)[0].astype(np.int64)
        self.rec_dh = np.where(near & rec.is_donor_h)[0].astype(np.int64)
        self.rec_metal = np.where(near & rec.is_metal)[0].astype(np.int64)
        self._rec_buf = rec.coords.copy()

        # Receptor hydroxyls close enough to matter become extra genes.
        if rec.rot_h_axis_b.size:
            keep = np.asarray([bool(near[b]) for b in rec.rot_h_axis_b])
            self.rec = rec
            self._rot_a = rec.rot_h_axis_a[keep]
            self._rot_b = rec.rot_h_axis_b[keep]
            self._rot_h = np.ascontiguousarray(rec.rot_h_atoms[keep])
        else:
            self._rot_a = np.zeros(0, dtype=np.int64)
            self._rot_b = np.zeros(0, dtype=np.int64)
            self._rot_h = np.zeros((0, 3), dtype=np.int64)
        self.n_rot_h = int(self._rot_a.size)
        self._touched = np.ascontiguousarray(
            self._rot_h[self._rot_h >= 0].ravel().astype(np.int64))

    # ------------------------------------------------------------------ #
    def restrict_rot_h(self, points: np.ndarray, cutoff: float = 5.5) -> int:
        """Keep only the receptor hydroxyls that line the cavity.

        Every rotatable hydroxyl is a search dimension, and one on the far side
        of the box cannot affect the score no matter where it points.  Dropping
        those is the difference between a 50-gene and a 20-gene problem.
        """
        if self.n_rot_h == 0 or points is None or len(points) == 0:
            return self.n_rot_h
        pts = np.ascontiguousarray(points)
        keep = np.zeros(self.n_rot_h, dtype=bool)
        for k, b in enumerate(self._rot_b):
            d2 = np.min(np.sum((pts - self.rec.coords[b]) ** 2, axis=1))
            keep[k] = d2 <= cutoff * cutoff
        self._rot_a = np.ascontiguousarray(self._rot_a[keep])
        self._rot_b = np.ascontiguousarray(self._rot_b[keep])
        self._rot_h = np.ascontiguousarray(self._rot_h[keep])
        self.n_rot_h = int(self._rot_a.size)
        self._touched = np.ascontiguousarray(
            self._rot_h[self._rot_h >= 0].ravel().astype(np.int64))
        self._rec_buf = self.rec.coords.copy()
        return self.n_rot_h

    def receptor_coords(self, rot_h_genes: np.ndarray | None):
        if self.n_rot_h == 0 or rot_h_genes is None or rot_h_genes.size == 0:
            return self.rec.coords
        apply_rot_h(self.rec.coords, self._rec_buf, rot_h_genes,
                    self._rot_a, self._rot_b, self._rot_h, self._touched)
        return self._rec_buf

    def score(self, lig_xyz: np.ndarray, rot_h_genes=None,
              breakdown: bool = False):
        rec_xyz = self.receptor_coords(rot_h_genes)
        if self.kind == "empirical":
            return self._score_empirical(lig_xyz, rec_xyz, breakdown)
        return self._score_ff(lig_xyz, rec_xyz, breakdown)

    # ------------------------------------------------------------------ #
    def _out_of_box(self, lig_xyz):
        g = self.grid
        return gridmod.out_of_box(g.origin, 1.0 / g.spacing,
                                  g.shape[0], g.shape[1], g.shape[2],
                                  lig_xyz, self.lig_scored)

    def _score_empirical(self, lig_xyz, rec_xyz, breakdown):
        g = self.grid
        inv = 1.0 / g.spacing
        T = gridmod.trilinear_sel
        lipo = T(g.maps, 0, g.origin, inv, lig_xyz, self.lig_lipo) \
            if self.lig_lipo.size else 0.0
        clash = (T(g.maps, 1, g.origin, inv, lig_xyz, self._lig_np)
                 + T(g.maps, 3, g.origin, inv, lig_xyz, self._lig_np)
                 + T(g.maps, 2, g.origin, inv, lig_xyz, self._lig_p)
                 + T(g.maps, 4, g.origin, inv, lig_xyz, self._lig_p))
        plp = T(g.maps, 5, g.origin, inv, lig_xyz, self.lig_scored) \
            if self.w_plp else 0.0
        hb = hbond_chemscore(lig_xyz, rec_xyz,
                             self.lig_dh, self.lig.donor_parent,
                             self.lig_acc, self.lig.acc_antecedent,
                             self.rec_acc, self.rec.acc_antecedent,
                             self.rec_dh, self.rec.donor_parent)
        met = metal_term(lig_xyz, rec_xyz, self.lig_acc, self.rec_metal) \
            if self.rec_metal.size else 0.0
        intra = internal_clash(lig_xyz, self.lig.intra_pairs, self._intra_r0)
        oob = self._out_of_box(lig_xyz)

        total = (params.CS_DG0
                 + params.CS_DG_HBOND * hb
                 + params.CS_DG_METAL * met
                 + params.CS_DG_LIPO * lipo
                 + params.CS_DG_ROT * self.lig.n_frozen_rot
                 + params.CS_CLASH_WEIGHT * self.w_clash * clash
                 + params.CS_INTERNAL_WEIGHT * self.w_internal * intra
                 + self.w_plp * plp
                 + params.CS_OUT_OF_BOX * oob)
        if breakdown:
            self._last = dict(total=total, hbond=hb, metal=met, lipo=lipo,
                              clash=clash, plp_steric=plp, internal=intra,
                              rotors=self.lig.n_frozen_rot, outside_box=oob)
        return total

    def _score_ff(self, lig_xyz, rec_xyz, breakdown):
        g = self.grid
        inv = 1.0 / g.spacing
        T = gridmod.trilinear_sel
        n_vdw = self._n_vdw
        i_elec, i_dv, i_ds = n_vdw, n_vdw + 1, n_vdw + 2

        vdw = 0.0
        for k in range(n_vdw):
            sel = self._vdw_sel[k]
            if sel.size:
                vdw += T(g.maps, k, g.origin, inv, lig_xyz, sel)
        W = gridmod.trilinear_weighted
        elec = W(g.maps, i_elec, g.origin, inv, lig_xyz, self.lig_scored,
                 self._q, 0.0)
        desolv = (W(g.maps, i_dv, g.origin, inv, lig_xyz, self.lig_scored,
                    self._s_i, 0.0)
                  + W(g.maps, i_ds, g.origin, inv, lig_xyz, self.lig_scored,
                      self._vol, 0.0))
        hb = hbond_12_10(lig_xyz, rec_xyz,
                         self.lig_dh, self.lig.donor_parent,
                         self.lig_acc, self.lig.acc_antecedent,
                         self.rec_acc, self.rec.acc_antecedent,
                         self.rec_dh, self.rec.donor_parent)
        intra = internal_ff(lig_xyz, self.lig.intra_pairs,
                            self.lig.rii, self.lig.epsii, self.lig.charge)
        oob = self._out_of_box(lig_xyz)

        total = (params.AD_W_VDW * vdw
                 + params.AD_W_ELEC * elec
                 + params.AD_W_DESOLV * desolv
                 + params.AD_W_HBOND * hb
                 + params.AD_W_TORS * self.lig.n_frozen_rot
                 + self.w_internal * intra
                 + self.out_of_box * oob)
        if breakdown:
            self._last = dict(total=total, vdw=vdw, elec=elec, desolv=desolv,
                              hbond=hb, internal=intra,
                              rotors=self.lig.n_frozen_rot, outside_box=oob)
        return total

    def terms(self) -> dict:
        return dict(self._last)


def build_scorer(kind: str, rec, lig, centre, half_widths,
                 spacing: float = 0.4) -> Scorer:
    kind = kind.lower()
    if kind in ("empirical", "chemscore"):
        g = gridmod.build_empirical_grid(rec, centre, half_widths, spacing)
        return Scorer(kind="empirical", grid=g, lig=lig, rec=rec)
    if kind in ("chemplp", "plp"):
        g = gridmod.build_empirical_grid(rec, centre, half_widths, spacing)
        return Scorer(kind="empirical", grid=g, lig=lig, rec=rec, w_plp=1.0)
    if kind in ("forcefield", "ff", "autodock", "goldscore"):
        g = gridmod.build_ff_grid(rec, lig, centre, half_widths, spacing)
        return Scorer(kind="forcefield", grid=g, lig=lig, rec=rec)
    raise ValueError(f"unknown scoring function: {kind}")
