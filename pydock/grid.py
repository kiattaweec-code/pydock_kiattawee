"""Precomputed interaction grids.

Terms that depend only on (ligand atom type, position) are tabulated on a
regular lattice once per run and read back by trilinear interpolation.  That is
what makes a GA with ~10^6 score evaluations affordable in Python.

Hydrogen bonds are deliberately *not* gridded: their value depends on the
ligand atom's own direction, and in GOLD it also depends on rotatable receptor
hydroxyls.  Those are evaluated pairwise in :mod:`pydock.scoring`.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from numba import njit

from . import params
from .geometry import FASTMATH

QSOLPAR = 0.01097
DESOLV_SIGMA = 3.6
CUTOFF = 8.0


@dataclass
class Grid:
    origin: np.ndarray          # (3,) corner
    spacing: float
    shape: tuple[int, int, int]
    maps: np.ndarray            # (K, nx, ny, nz)
    map_names: list[str]

    @property
    def centre(self):
        n = np.asarray(self.shape, dtype=float) - 1.0
        return self.origin + 0.5 * n * self.spacing

    def index_of(self, name: str) -> int:
        return self.map_names.index(name)


def make_lattice(centre, half_widths, spacing: float = 0.4):
    half = np.asarray(half_widths, dtype=float)
    npts = (2 * np.ceil(half / spacing) + 1).astype(np.int64)
    origin = np.asarray(centre, dtype=float) - (npts - 1) / 2.0 * spacing
    return origin, tuple(int(x) for x in npts)


# --------------------------------------------------------------------------- #
# Interpolation
# --------------------------------------------------------------------------- #
@njit(**FASTMATH)
def trilinear(maps, k, origin, inv_sp, coords, out_of_box_penalty):
    """Sum map `k` over all `coords`.  Points outside the box are penalised."""
    nx, ny, nz = maps.shape[1], maps.shape[2], maps.shape[3]
    total = 0.0
    for a in range(coords.shape[0]):
        fx = (coords[a, 0] - origin[0]) * inv_sp
        fy = (coords[a, 1] - origin[1]) * inv_sp
        fz = (coords[a, 2] - origin[2]) * inv_sp
        if fx < 0.0 or fy < 0.0 or fz < 0.0 or fx > nx - 2 or fy > ny - 2 or fz > nz - 2:
            total += out_of_box_penalty
            continue
        i = int(fx)
        j = int(fy)
        k2 = int(fz)
        u = fx - i
        v = fy - j
        w = fz - k2
        c000 = maps[k, i, j, k2]
        c100 = maps[k, i + 1, j, k2]
        c010 = maps[k, i, j + 1, k2]
        c110 = maps[k, i + 1, j + 1, k2]
        c001 = maps[k, i, j, k2 + 1]
        c101 = maps[k, i + 1, j, k2 + 1]
        c011 = maps[k, i, j + 1, k2 + 1]
        c111 = maps[k, i + 1, j + 1, k2 + 1]
        total += ((1 - u) * ((1 - v) * ((1 - w) * c000 + w * c001)
                             + v * ((1 - w) * c010 + w * c011))
                  + u * ((1 - v) * ((1 - w) * c100 + w * c101)
                         + v * ((1 - w) * c110 + w * c111)))
    return total


@njit(**FASTMATH)
def trilinear_sel(maps, k, origin, inv_sp, coords, sel):
    """Sum map `k` over the subset `sel`, without copying coordinates."""
    nx, ny, nz = maps.shape[1], maps.shape[2], maps.shape[3]
    total = 0.0
    for s in range(sel.shape[0]):
        a = sel[s]
        fx = (coords[a, 0] - origin[0]) * inv_sp
        fy = (coords[a, 1] - origin[1]) * inv_sp
        fz = (coords[a, 2] - origin[2]) * inv_sp
        if fx < 0.0 or fy < 0.0 or fz < 0.0 or fx > nx - 2 or fy > ny - 2 or fz > nz - 2:
            continue
        i = int(fx)
        j = int(fy)
        k2 = int(fz)
        u = fx - i
        v = fy - j
        w = fz - k2
        c000 = maps[k, i, j, k2]
        c100 = maps[k, i + 1, j, k2]
        c010 = maps[k, i, j + 1, k2]
        c110 = maps[k, i + 1, j + 1, k2]
        c001 = maps[k, i, j, k2 + 1]
        c101 = maps[k, i + 1, j, k2 + 1]
        c011 = maps[k, i, j + 1, k2 + 1]
        c111 = maps[k, i + 1, j + 1, k2 + 1]
        total += ((1 - u) * ((1 - v) * ((1 - w) * c000 + w * c001)
                             + v * ((1 - w) * c010 + w * c011))
                  + u * ((1 - v) * ((1 - w) * c100 + w * c101)
                         + v * ((1 - w) * c110 + w * c111)))
    return total


@njit(**FASTMATH)
def out_of_box(origin, inv_sp, shape_x, shape_y, shape_z, coords, sel):
    """How many of the selected atoms fall outside the interpolable region."""
    n = 0
    for s in range(sel.shape[0]):
        a = sel[s]
        fx = (coords[a, 0] - origin[0]) * inv_sp
        fy = (coords[a, 1] - origin[1]) * inv_sp
        fz = (coords[a, 2] - origin[2]) * inv_sp
        if (fx < 0.0 or fy < 0.0 or fz < 0.0
                or fx > shape_x - 2 or fy > shape_y - 2 or fz > shape_z - 2):
            n += 1
    return n


@njit(**FASTMATH)
def trilinear_weighted(maps, k, origin, inv_sp, coords, sel, weights,
                       out_of_box_penalty):
    """Weighted sum of map `k` over the subset `sel` of `coords`."""
    nx, ny, nz = maps.shape[1], maps.shape[2], maps.shape[3]
    total = 0.0
    for s in range(sel.shape[0]):
        a = sel[s]
        fx = (coords[a, 0] - origin[0]) * inv_sp
        fy = (coords[a, 1] - origin[1]) * inv_sp
        fz = (coords[a, 2] - origin[2]) * inv_sp
        if fx < 0.0 or fy < 0.0 or fz < 0.0 or fx > nx - 2 or fy > ny - 2 or fz > nz - 2:
            total += out_of_box_penalty
            continue
        i = int(fx)
        j = int(fy)
        k2 = int(fz)
        u = fx - i
        v = fy - j
        w = fz - k2
        c000 = maps[k, i, j, k2]
        c100 = maps[k, i + 1, j, k2]
        c010 = maps[k, i, j + 1, k2]
        c110 = maps[k, i + 1, j + 1, k2]
        c001 = maps[k, i, j, k2 + 1]
        c101 = maps[k, i + 1, j, k2 + 1]
        c011 = maps[k, i, j + 1, k2 + 1]
        c111 = maps[k, i + 1, j + 1, k2 + 1]
        val = ((1 - u) * ((1 - v) * ((1 - w) * c000 + w * c001)
                          + v * ((1 - w) * c010 + w * c011))
               + u * ((1 - v) * ((1 - w) * c100 + w * c101)
                      + v * ((1 - w) * c110 + w * c111)))
        total += weights[s] * val
    return total


# --------------------------------------------------------------------------- #
# Pair potentials used while filling the maps
# --------------------------------------------------------------------------- #
@njit(**FASTMATH)
def _vdw_12_6(r, rij, eps, smooth):
    """Lennard-Jones 12-6 with the minimum flattened over +/- `smooth`."""
    rr = r
    if rr < rij - smooth:
        rr = rr + smooth
    elif rr > rij + smooth:
        rr = rr - smooth
    else:
        rr = rij
    if rr < 0.5:
        rr = 0.5
    ratio = rij / rr
    r6 = ratio ** 6
    return eps * (r6 * r6 - 2.0 * r6)


@njit(**FASTMATH)
def _dielectric(r):
    """Mehler-Solmajer sigmoidal distance-dependent dielectric."""
    lam_b = params.MS_LAMBDA * params.MS_B
    e = params.MS_A + params.MS_B / (1.0 + params.MS_K
                                     * np.exp(-lam_b * r))
    return e


@njit(**FASTMATH)
def _plp(r, a, b, c, d, well, wall):
    """Piecewise linear potential (Gehlhaar form).

    r < a          : linear wall rising to `wall`
    a <= r < b     : linear descent from 0 to `well`
    b <= r < c     : flat at `well`
    c <= r < d     : linear return to 0
    r >= d         : 0
    """
    if r >= d:
        return 0.0
    if r >= c:
        return well * (d - r) / (d - c)
    if r >= b:
        return well
    if r >= a:
        return well * (r - a) / (b - a)
    if r <= 0.1:
        return wall
    return wall * (a - r) / a


@njit(**FASTMATH)
def _lipo_block(r, r1, r2):
    if r <= r1:
        return 1.0
    if r >= r2:
        return 0.0
    return (r2 - r) / (r2 - r1)


# --------------------------------------------------------------------------- #
# Map construction (scatter from each receptor atom into its local sub-box)
# --------------------------------------------------------------------------- #
@njit(cache=True, nogil=True)
def _fill_vdw(maps, k, origin, spacing, rec_xyz, rec_rii, rec_eps,
              lig_rii, lig_eps, cutoff, smooth):
    nx, ny, nz = maps.shape[1], maps.shape[2], maps.shape[3]
    span = int(cutoff / spacing) + 1
    for j in range(rec_xyz.shape[0]):
        rij = 0.5 * (rec_rii[j] + lig_rii)
        eps = np.sqrt(rec_eps[j] * lig_eps)
        if eps <= 0.0:
            continue
        gx = (rec_xyz[j, 0] - origin[0]) / spacing
        gy = (rec_xyz[j, 1] - origin[1]) / spacing
        gz = (rec_xyz[j, 2] - origin[2]) / spacing
        i0 = max(0, int(gx) - span)
        i1 = min(nx, int(gx) + span + 2)
        j0 = max(0, int(gy) - span)
        j1 = min(ny, int(gy) + span + 2)
        k0 = max(0, int(gz) - span)
        k1 = min(nz, int(gz) + span + 2)
        for ii in range(i0, i1):
            dx = origin[0] + ii * spacing - rec_xyz[j, 0]
            for jj in range(j0, j1):
                dy = origin[1] + jj * spacing - rec_xyz[j, 1]
                d2xy = dx * dx + dy * dy
                if d2xy > cutoff * cutoff:
                    continue
                for kk in range(k0, k1):
                    dz = origin[2] + kk * spacing - rec_xyz[j, 2]
                    d2 = d2xy + dz * dz
                    if d2 > cutoff * cutoff:
                        continue
                    maps[k, ii, jj, kk] += _vdw_12_6(np.sqrt(d2), rij, eps, smooth)


@njit(cache=True, nogil=True)
def _fill_elec(maps, k, origin, spacing, rec_xyz, rec_q, cutoff):
    nx, ny, nz = maps.shape[1], maps.shape[2], maps.shape[3]
    span = int(cutoff / spacing) + 1
    for j in range(rec_xyz.shape[0]):
        if rec_q[j] == 0.0:
            continue
        gx = (rec_xyz[j, 0] - origin[0]) / spacing
        gy = (rec_xyz[j, 1] - origin[1]) / spacing
        gz = (rec_xyz[j, 2] - origin[2]) / spacing
        i0 = max(0, int(gx) - span)
        i1 = min(nx, int(gx) + span + 2)
        j0 = max(0, int(gy) - span)
        j1 = min(ny, int(gy) + span + 2)
        k0 = max(0, int(gz) - span)
        k1 = min(nz, int(gz) + span + 2)
        for ii in range(i0, i1):
            dx = origin[0] + ii * spacing - rec_xyz[j, 0]
            for jj in range(j0, j1):
                dy = origin[1] + jj * spacing - rec_xyz[j, 1]
                d2xy = dx * dx + dy * dy
                if d2xy > cutoff * cutoff:
                    continue
                for kk in range(k0, k1):
                    dz = origin[2] + kk * spacing - rec_xyz[j, 2]
                    d2 = d2xy + dz * dz
                    if d2 > cutoff * cutoff:
                        continue
                    r = np.sqrt(d2)
                    if r < 0.5:
                        r = 0.5
                    maps[k, ii, jj, kk] += 332.06 * rec_q[j] / (r * _dielectric(r))


@njit(cache=True, nogil=True)
def _fill_gaussian(maps, k, origin, spacing, rec_xyz, weight, sigma, cutoff):
    nx, ny, nz = maps.shape[1], maps.shape[2], maps.shape[3]
    span = int(cutoff / spacing) + 1
    two_s2 = 2.0 * sigma * sigma
    for j in range(rec_xyz.shape[0]):
        if weight[j] == 0.0:
            continue
        gx = (rec_xyz[j, 0] - origin[0]) / spacing
        gy = (rec_xyz[j, 1] - origin[1]) / spacing
        gz = (rec_xyz[j, 2] - origin[2]) / spacing
        i0 = max(0, int(gx) - span)
        i1 = min(nx, int(gx) + span + 2)
        j0 = max(0, int(gy) - span)
        j1 = min(ny, int(gy) + span + 2)
        k0 = max(0, int(gz) - span)
        k1 = min(nz, int(gz) + span + 2)
        for ii in range(i0, i1):
            dx = origin[0] + ii * spacing - rec_xyz[j, 0]
            for jj in range(j0, j1):
                dy = origin[1] + jj * spacing - rec_xyz[j, 1]
                d2xy = dx * dx + dy * dy
                if d2xy > cutoff * cutoff:
                    continue
                for kk in range(k0, k1):
                    dz = origin[2] + kk * spacing - rec_xyz[j, 2]
                    d2 = d2xy + dz * dz
                    if d2 > cutoff * cutoff:
                        continue
                    maps[k, ii, jj, kk] += weight[j] * np.exp(-d2 / two_s2)


@njit(cache=True, nogil=True)
def _fill_plp(maps, k, origin, spacing, rec_xyz, a, b, c, d, well, wall):
    nx, ny, nz = maps.shape[1], maps.shape[2], maps.shape[3]
    cutoff = d
    span = int(cutoff / spacing) + 1
    for j in range(rec_xyz.shape[0]):
        gx = (rec_xyz[j, 0] - origin[0]) / spacing
        gy = (rec_xyz[j, 1] - origin[1]) / spacing
        gz = (rec_xyz[j, 2] - origin[2]) / spacing
        i0 = max(0, int(gx) - span)
        i1 = min(nx, int(gx) + span + 2)
        j0 = max(0, int(gy) - span)
        j1 = min(ny, int(gy) + span + 2)
        k0 = max(0, int(gz) - span)
        k1 = min(nz, int(gz) + span + 2)
        for ii in range(i0, i1):
            dx = origin[0] + ii * spacing - rec_xyz[j, 0]
            for jj in range(j0, j1):
                dy = origin[1] + jj * spacing - rec_xyz[j, 1]
                d2xy = dx * dx + dy * dy
                if d2xy > cutoff * cutoff:
                    continue
                for kk in range(k0, k1):
                    dz = origin[2] + kk * spacing - rec_xyz[j, 2]
                    d2 = d2xy + dz * dz
                    if d2 > cutoff * cutoff:
                        continue
                    maps[k, ii, jj, kk] += _plp(np.sqrt(d2), a, b, c, d, well, wall)


@njit(cache=True, nogil=True)
def _fill_clash(maps, k, origin, spacing, rec_xyz, r0):
    """Quadratic steric repulsion: (r0 - r)^2 for r < r0, else zero."""
    nx, ny, nz = maps.shape[1], maps.shape[2], maps.shape[3]
    span = int(r0 / spacing) + 1
    for j in range(rec_xyz.shape[0]):
        gx = (rec_xyz[j, 0] - origin[0]) / spacing
        gy = (rec_xyz[j, 1] - origin[1]) / spacing
        gz = (rec_xyz[j, 2] - origin[2]) / spacing
        i0 = max(0, int(gx) - span)
        i1 = min(nx, int(gx) + span + 2)
        j0 = max(0, int(gy) - span)
        j1 = min(ny, int(gy) + span + 2)
        k0 = max(0, int(gz) - span)
        k1 = min(nz, int(gz) + span + 2)
        for ii in range(i0, i1):
            dx = origin[0] + ii * spacing - rec_xyz[j, 0]
            for jj in range(j0, j1):
                dy = origin[1] + jj * spacing - rec_xyz[j, 1]
                d2xy = dx * dx + dy * dy
                if d2xy > r0 * r0:
                    continue
                for kk in range(k0, k1):
                    dz = origin[2] + kk * spacing - rec_xyz[j, 2]
                    d2 = d2xy + dz * dz
                    if d2 > r0 * r0:
                        continue
                    d = r0 - np.sqrt(d2)
                    maps[k, ii, jj, kk] += d * d


@njit(cache=True, nogil=True)
def _fill_lipo(maps, k, origin, spacing, rec_xyz, r1, r2):
    nx, ny, nz = maps.shape[1], maps.shape[2], maps.shape[3]
    span = int(r2 / spacing) + 1
    for j in range(rec_xyz.shape[0]):
        gx = (rec_xyz[j, 0] - origin[0]) / spacing
        gy = (rec_xyz[j, 1] - origin[1]) / spacing
        gz = (rec_xyz[j, 2] - origin[2]) / spacing
        i0 = max(0, int(gx) - span)
        i1 = min(nx, int(gx) + span + 2)
        j0 = max(0, int(gy) - span)
        j1 = min(ny, int(gy) + span + 2)
        k0 = max(0, int(gz) - span)
        k1 = min(nz, int(gz) + span + 2)
        for ii in range(i0, i1):
            dx = origin[0] + ii * spacing - rec_xyz[j, 0]
            for jj in range(j0, j1):
                dy = origin[1] + jj * spacing - rec_xyz[j, 1]
                d2xy = dx * dx + dy * dy
                if d2xy > r2 * r2:
                    continue
                for kk in range(k0, k1):
                    dz = origin[2] + kk * spacing - rec_xyz[j, 2]
                    d2 = d2xy + dz * dz
                    if d2 > r2 * r2:
                        continue
                    maps[k, ii, jj, kk] += _lipo_block(np.sqrt(d2), r1, r2)


# --------------------------------------------------------------------------- #
# Public builders
# --------------------------------------------------------------------------- #
def _near_box(rec, origin, shape, spacing, margin=CUTOFF + 1.0):
    hi = origin + (np.asarray(shape) - 1) * spacing
    lo = origin - margin
    hi = hi + margin
    m = np.all((rec.coords >= lo) & (rec.coords <= hi), axis=1)
    return np.where(m)[0]


def build_ff_grid(rec, lig, centre, half_widths, spacing=0.4) -> Grid:
    """van der Waals maps per ligand type, plus electrostatic and desolvation."""
    origin, shape = make_lattice(centre, half_widths, spacing)
    near = _near_box(rec, origin, shape, spacing)
    xyz = np.ascontiguousarray(rec.coords[near])

    lig_types = sorted(set(t for t in lig.ad_type if t != "H"))
    names = [f"vdw:{t}" for t in lig_types] + ["elec", "desolv_v", "desolv_s"]
    maps = np.zeros((len(names),) + shape, dtype=np.float64)
    origin = np.ascontiguousarray(origin)

    for k, t in enumerate(lig_types):
        rii, eps, _, _ = params.AD_PARAMS.get(t, params.AD_PARAMS["C"])
        _fill_vdw(maps, k, origin, spacing, xyz,
                  np.ascontiguousarray(rec.rii[near]),
                  np.ascontiguousarray(rec.epsii[near]),
                  rii, eps, CUTOFF, params.AD_SMOOTH)

    _fill_elec(maps, len(lig_types), origin, spacing, xyz,
               np.ascontiguousarray(rec.charge[near]), CUTOFF)
    _fill_gaussian(maps, len(lig_types) + 1, origin, spacing, xyz,
                   np.ascontiguousarray(rec.vol[near]), DESOLV_SIGMA, CUTOFF)
    sj = rec.solpar[near] + QSOLPAR * np.abs(rec.charge[near])
    _fill_gaussian(maps, len(lig_types) + 2, origin, spacing, xyz,
                   np.ascontiguousarray(sj), DESOLV_SIGMA, CUTOFF)

    return Grid(origin=origin, spacing=spacing, shape=shape, maps=maps,
                map_names=names)


EMPIRICAL_MAPS = ["lipo",
                  "clash:recNP_ligNP", "clash:recNP_ligP",
                  "clash:recP_ligNP", "clash:recP_ligP",
                  "plp_steric"]


def build_empirical_grid(rec, centre, half_widths, spacing=0.4) -> Grid:
    """Lipophilic contact map, four polarity-resolved clash maps, PLP steric.

    Steric repulsion is split by the polarity of *both* partners because a
    polar contact pair legitimately approaches to ~2.6 A (that is a hydrogen
    bond) while an apolar pair at the same separation is a clash.  A single
    type-independent clash map would punish every hydrogen bond it found.
    """
    origin, shape = make_lattice(centre, half_widths, spacing)
    near = _near_box(rec, origin, shape, spacing)
    is_heavy = np.asarray([rec.element[i] != "H" for i in near])
    is_polar = np.asarray([rec.element[i] in ("N", "O") for i in near])
    heavy = near[is_heavy]
    rec_np = near[is_heavy & ~is_polar]
    rec_p = near[is_heavy & is_polar]
    lipo = near[np.asarray([bool(rec.is_lipophilic[i]) for i in near])]

    maps = np.zeros((len(EMPIRICAL_MAPS),) + shape, dtype=np.float64)
    origin = np.ascontiguousarray(origin)

    if lipo.size:
        _fill_lipo(maps, 0, origin, spacing,
                   np.ascontiguousarray(rec.coords[lipo]),
                   params.CS_LIPO_R1, params.CS_LIPO_R2)
    for k, (sel, r0) in enumerate((
            (rec_np, params.CLASH_R_NP_NP),
            (rec_np, params.CLASH_R_NP_P),
            (rec_p, params.CLASH_R_NP_P),
            (rec_p, params.CLASH_R_P_P)), start=1):
        if sel.size:
            _fill_clash(maps, k, origin, spacing,
                        np.ascontiguousarray(rec.coords[sel]), r0)

    a, b, c, d, well, wall = params.PLP_STERIC
    _fill_plp(maps, 5, origin, spacing,
              np.ascontiguousarray(rec.coords[heavy]), a, b, c, d, well, 0.0)

    return Grid(origin=origin, spacing=spacing, shape=shape, maps=maps,
                map_names=list(EMPIRICAL_MAPS))
