"""Cavity detection, LIGSITE-style.

Docking into a bounding box wastes almost all of the search effort on empty
solvent and on the protein interior.  GOLD gets around this by detecting the
cavity first and only ever placing the ligand inside it.  This module does the
same: it marks grid points that are both *free* (no protein atom within a probe
radius) and *buried* (enclosed by protein along most scan directions), keeps
the largest connected blob of them, and hands back those points.

Reference: Hendlich, Rippmann & Barnickel, J. Mol. Graph. Model. 1997, 15, 359.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from numba import njit

# 7 scan axes: 3 cartesian + 4 body diagonals.  Each is scanned both ways.
_SCAN = np.array([
    [1, 0, 0], [0, 1, 0], [0, 0, 1],
    [1, 1, 1], [1, 1, -1], [1, -1, 1], [-1, 1, 1],
], dtype=np.int64)


@dataclass
class Cavity:
    points: np.ndarray        # (K, 3) cavity point coordinates
    free: np.ndarray          # (nx, ny, nz) bool, no protein within the probe
    mask: np.ndarray          # (nx, ny, nz) bool, in the selected cavity
    origin: np.ndarray
    spacing: float

    @property
    def volume(self) -> float:
        return float(self.mask.sum()) * self.spacing ** 3

    def centre(self):
        return self.points.mean(axis=0)


@njit(cache=True, nogil=True)
def _mark_occupied(free, origin, spacing, xyz, radii):
    nx, ny, nz = free.shape
    for j in range(xyz.shape[0]):
        r = radii[j]
        span = int(r / spacing) + 1
        gx = int((xyz[j, 0] - origin[0]) / spacing)
        gy = int((xyz[j, 1] - origin[1]) / spacing)
        gz = int((xyz[j, 2] - origin[2]) / spacing)
        for ii in range(max(0, gx - span), min(nx, gx + span + 2)):
            dx = origin[0] + ii * spacing - xyz[j, 0]
            for jj in range(max(0, gy - span), min(ny, gy + span + 2)):
                dy = origin[1] + jj * spacing - xyz[j, 1]
                if dx * dx + dy * dy > r * r:
                    continue
                for kk in range(max(0, gz - span), min(nz, gz + span + 2)):
                    dz = origin[2] + kk * spacing - xyz[j, 2]
                    if dx * dx + dy * dy + dz * dz <= r * r:
                        free[ii, jj, kk] = False


@njit(cache=True, nogil=True)
def _burial(free, scan, max_steps, min_psp):
    """Count protein-solvent-protein events per free point."""
    nx, ny, nz = free.shape
    out = np.zeros(free.shape, dtype=np.bool_)
    for i in range(nx):
        for j in range(ny):
            for k in range(nz):
                if not free[i, j, k]:
                    continue
                psp = 0
                for s in range(scan.shape[0]):
                    dx = scan[s, 0]
                    dy = scan[s, 1]
                    dz = scan[s, 2]
                    hit_pos = False
                    x, y, z = i + dx, j + dy, k + dz
                    steps = 0
                    while 0 <= x < nx and 0 <= y < ny and 0 <= z < nz and steps < max_steps:
                        if not free[x, y, z]:
                            hit_pos = True
                            break
                        x += dx
                        y += dy
                        z += dz
                        steps += 1
                    if not hit_pos:
                        continue
                    hit_neg = False
                    x, y, z = i - dx, j - dy, k - dz
                    steps = 0
                    while 0 <= x < nx and 0 <= y < ny and 0 <= z < nz and steps < max_steps:
                        if not free[x, y, z]:
                            hit_neg = True
                            break
                        x -= dx
                        y -= dy
                        z -= dz
                        steps += 1
                    if hit_neg:
                        psp += 1
                if psp >= min_psp:
                    out[i, j, k] = True
    return out


@njit(cache=True, nogil=True)
def _label_components(mask):
    """Flood-fill every blob; returns (labels, size_per_label)."""
    nx, ny, nz = mask.shape
    label = np.zeros(mask.shape, dtype=np.int32)
    stack = np.empty((nx * ny * nz, 3), dtype=np.int64)
    sizes = np.zeros(nx * ny * nz // 4 + 2, dtype=np.int64)
    current = 0
    for si in range(nx):
        for sj in range(ny):
            for sk in range(nz):
                if not mask[si, sj, sk] or label[si, sj, sk] != 0:
                    continue
                current += 1
                if current >= sizes.shape[0]:
                    return label, sizes[:current]
                top = 0
                stack[top, 0] = si
                stack[top, 1] = sj
                stack[top, 2] = sk
                top = 1
                label[si, sj, sk] = current
                size = 0
                while top > 0:
                    top -= 1
                    x = stack[top, 0]
                    y = stack[top, 1]
                    z = stack[top, 2]
                    size += 1
                    for d in range(6):
                        if d == 0:
                            a, b, c = x + 1, y, z
                        elif d == 1:
                            a, b, c = x - 1, y, z
                        elif d == 2:
                            a, b, c = x, y + 1, z
                        elif d == 3:
                            a, b, c = x, y - 1, z
                        elif d == 4:
                            a, b, c = x, y, z + 1
                        else:
                            a, b, c = x, y, z - 1
                        if a < 0 or b < 0 or c < 0 or a >= nx or b >= ny or c >= nz:
                            continue
                        if mask[a, b, c] and label[a, b, c] == 0:
                            label[a, b, c] = current
                            stack[top, 0] = a
                            stack[top, 1] = b
                            stack[top, 2] = c
                            top += 1
                sizes[current] = size
    return label, sizes[:current + 1]


@njit(cache=True, nogil=True)
def _largest_component(mask):
    """Flood-fill; returns a mask holding only the biggest connected blob."""
    nx, ny, nz = mask.shape
    label = np.zeros(mask.shape, dtype=np.int32)
    stack = np.empty((nx * ny * nz, 3), dtype=np.int64)
    best_label = 0
    best_size = 0
    current = 0
    for si in range(nx):
        for sj in range(ny):
            for sk in range(nz):
                if not mask[si, sj, sk] or label[si, sj, sk] != 0:
                    continue
                current += 1
                top = 0
                stack[top, 0] = si
                stack[top, 1] = sj
                stack[top, 2] = sk
                top = 1
                label[si, sj, sk] = current
                size = 0
                while top > 0:
                    top -= 1
                    x = stack[top, 0]
                    y = stack[top, 1]
                    z = stack[top, 2]
                    size += 1
                    for d in range(6):
                        if d == 0:
                            a, b, c = x + 1, y, z
                        elif d == 1:
                            a, b, c = x - 1, y, z
                        elif d == 2:
                            a, b, c = x, y + 1, z
                        elif d == 3:
                            a, b, c = x, y - 1, z
                        elif d == 4:
                            a, b, c = x, y, z + 1
                        else:
                            a, b, c = x, y, z - 1
                        if a < 0 or b < 0 or c < 0 or a >= nx or b >= ny or c >= nz:
                            continue
                        if mask[a, b, c] and label[a, b, c] == 0:
                            label[a, b, c] = current
                            stack[top, 0] = a
                            stack[top, 1] = b
                            stack[top, 2] = c
                            top += 1
                if size > best_size:
                    best_size = size
                    best_label = current
    out = np.zeros(mask.shape, dtype=np.bool_)
    for i in range(nx):
        for j in range(ny):
            for k in range(nz):
                if label[i, j, k] == best_label:
                    out[i, j, k] = True
    return out


def detect_cavity(rec, origin, shape, spacing: float = 0.6,
                  probe: float = 1.4, min_psp: int = 4,
                  scan_range: float = 11.0,
                  seed_point=None, seed_radius: float = 11.0,
                  largest_only: bool = True) -> Cavity:
    """Find the cavity inside the given lattice.

    `probe` is added to each atom's van der Waals radius when marking occupied
    space, so `free` means "a ligand heavy atom could sit here".
    """
    from . import params

    origin = np.asarray(origin, dtype=np.float64)
    shape = tuple(int(s) for s in shape)
    hi = origin + (np.asarray(shape) - 1) * spacing

    lo_m = origin - (scan_range + 3.0)
    hi_m = hi + (scan_range + 3.0)
    sel = np.where(np.all((rec.coords >= lo_m) & (rec.coords <= hi_m), axis=1)
                   & np.asarray([e != "H" for e in rec.element]))[0]
    xyz = np.ascontiguousarray(rec.coords[sel])
    radii = np.ascontiguousarray(
        np.asarray([params.ELEMENT_VDW.get(rec.element[i], 1.7) for i in sel])
        + probe)

    free = np.ones(shape, dtype=np.bool_)
    _mark_occupied(free, origin, spacing, xyz, radii)
    buried = _burial(free, _SCAN, int(scan_range / spacing), min_psp)

    if seed_point is not None:
        gx, gy, gz = np.meshgrid(*[np.arange(s) for s in shape], indexing="ij")
        pts = np.stack([gx, gy, gz], axis=-1) * spacing + origin
        d = np.linalg.norm(pts - np.asarray(seed_point), axis=-1)
        buried &= d <= seed_radius

    mask = _largest_component(buried) if largest_only else buried
    idx = np.argwhere(mask)
    points = idx * spacing + origin
    return Cavity(points=np.ascontiguousarray(points), free=free, mask=mask,
                  origin=origin, spacing=spacing)


@dataclass
class Pocket:
    """One detected cavity, ranked by volume."""
    rank: int
    centre: np.ndarray
    volume: float
    n_points: int
    lower: np.ndarray          # bounding box of the pocket points
    upper: np.ndarray

    @property
    def size(self) -> np.ndarray:
        """A docking box around the pocket: enclosing it, but capped.

        Open channels and surface grooves can span 35 A or more, and a box
        that large makes the search far slower than it needs to be.  The cap
        gives a sane starting box that the user can widen if they mean to.
        """
        return np.clip(self.upper - self.lower + 6.0, 12.0, 30.0)


def find_pockets(rec, spacing: float = 0.8, probe: float = 1.4,
                 min_psp: int = 5, scan_range: float = 11.0,
                 max_pockets: int = 6, min_volume: float = 80.0,
                 margin: float = 3.0) -> list[Pocket]:
    """Scan the whole receptor and return its biggest pockets, best first.

    This is the search-free version of :func:`detect_cavity`: instead of
    confining the search to a box the caller already knows about, it looks over
    the entire structure so a user with no idea where the site is can be shown
    candidates to pick from.  A coarser default lattice keeps that affordable.
    """
    heavy = np.asarray([e != "H" for e in rec.element])
    xyz_all = rec.coords[heavy]
    lo = xyz_all.min(axis=0) - margin
    hi = xyz_all.max(axis=0) + margin
    shape = tuple(int(x) for x in np.ceil((hi - lo) / spacing) + 1)

    cav = detect_cavity(rec, lo, shape, spacing, probe=probe,
                        min_psp=min_psp, scan_range=scan_range,
                        largest_only=False)
    labels, sizes = _label_components(cav.mask)
    if sizes.shape[0] <= 1:
        return []

    order = np.argsort(sizes)[::-1]
    vox = spacing ** 3
    pockets: list[Pocket] = []
    for lab in order:
        if lab == 0 or sizes[lab] == 0:
            continue
        volume = float(sizes[lab]) * vox
        if volume < min_volume or len(pockets) >= max_pockets:
            if len(pockets) >= max_pockets:
                break
            continue
        pts = np.argwhere(labels == lab) * spacing + lo
        pockets.append(Pocket(
            rank=len(pockets) + 1,
            centre=pts.mean(axis=0),
            volume=volume,
            n_points=int(sizes[lab]),
            lower=pts.min(axis=0),
            upper=pts.max(axis=0),
        ))
    return pockets
