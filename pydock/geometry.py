"""Fast geometric primitives: quaternions, axis rotations, torsion-tree transforms.

All hot functions are Numba-compiled and operate on plain float64 arrays so the
GA inner loop never touches Python objects.
"""
from __future__ import annotations

import numpy as np
from numba import njit

FASTMATH = dict(cache=True, fastmath=True, nogil=True)


# --------------------------------------------------------------------------- #
# Rotations
# --------------------------------------------------------------------------- #
@njit(**FASTMATH)
def quat_to_matrix(q):
    """Unit-normalised quaternion (w, x, y, z) -> 3x3 rotation matrix."""
    n = np.sqrt(q[0] * q[0] + q[1] * q[1] + q[2] * q[2] + q[3] * q[3])
    if n < 1e-12:
        w, x, y, z = 1.0, 0.0, 0.0, 0.0
    else:
        w, x, y, z = q[0] / n, q[1] / n, q[2] / n, q[3] / n
    m = np.empty((3, 3))
    m[0, 0] = 1.0 - 2.0 * (y * y + z * z)
    m[0, 1] = 2.0 * (x * y - w * z)
    m[0, 2] = 2.0 * (x * z + w * y)
    m[1, 0] = 2.0 * (x * y + w * z)
    m[1, 1] = 1.0 - 2.0 * (x * x + z * z)
    m[1, 2] = 2.0 * (y * z - w * x)
    m[2, 0] = 2.0 * (x * z - w * y)
    m[2, 1] = 2.0 * (y * z + w * x)
    m[2, 2] = 1.0 - 2.0 * (x * x + y * y)
    return m


@njit(**FASTMATH)
def rotate_about_axis(coords, idx, origin, axis, angle):
    """Rotate coords[idx] about the line (origin, axis) by `angle` radians.

    Rodrigues' formula, applied in place.
    """
    n = np.sqrt(axis[0] ** 2 + axis[1] ** 2 + axis[2] ** 2)
    if n < 1e-9:
        return
    ux, uy, uz = axis[0] / n, axis[1] / n, axis[2] / n
    c = np.cos(angle)
    s = np.sin(angle)
    t = 1.0 - c

    # Pre-build the matrix once rather than per atom.
    m00 = t * ux * ux + c
    m01 = t * ux * uy - s * uz
    m02 = t * ux * uz + s * uy
    m10 = t * ux * uy + s * uz
    m11 = t * uy * uy + c
    m12 = t * uy * uz - s * ux
    m20 = t * ux * uz - s * uy
    m21 = t * uy * uz + s * ux
    m22 = t * uz * uz + c

    for k in range(idx.shape[0]):
        i = idx[k]
        px = coords[i, 0] - origin[0]
        py = coords[i, 1] - origin[1]
        pz = coords[i, 2] - origin[2]
        coords[i, 0] = origin[0] + m00 * px + m01 * py + m02 * pz
        coords[i, 1] = origin[1] + m10 * px + m11 * py + m12 * pz
        coords[i, 2] = origin[2] + m20 * px + m21 * py + m22 * pz


# --------------------------------------------------------------------------- #
# Torsion tree
# --------------------------------------------------------------------------- #
@njit(**FASTMATH)
def apply_pose(coords0, torsions, bond_a, bond_b, move_ptr, move_idx,
               quat, centre, out):
    """Build pose coordinates from a chromosome.

    Parameters
    ----------
    coords0 : (N, 3) reference internal geometry (never modified)
    torsions : (T,) torsion *deltas* in radians, applied to the reference
    bond_a, bond_b : (T,) atom indices defining each rotation axis a->b
    move_ptr : (T+1,) CSR-style offsets into move_idx
    move_idx : flat array of atom indices rotated by each bond
    quat : (4,) global orientation quaternion
    centre : (3,) target position of the ligand centroid
    out : (N, 3) output buffer

    Torsions must be ordered parent-before-child (BFS from the root fragment)
    so that a child's rotation axis is already in its final place.
    """
    n = coords0.shape[0]
    for i in range(n):
        out[i, 0] = coords0[i, 0]
        out[i, 1] = coords0[i, 1]
        out[i, 2] = coords0[i, 2]

    axis = np.empty(3)
    origin = np.empty(3)
    for t in range(torsions.shape[0]):
        a = bond_a[t]
        b = bond_b[t]
        origin[0] = out[a, 0]
        origin[1] = out[a, 1]
        origin[2] = out[a, 2]
        axis[0] = out[b, 0] - origin[0]
        axis[1] = out[b, 1] - origin[1]
        axis[2] = out[b, 2] - origin[2]
        rotate_about_axis(out, move_idx[move_ptr[t]:move_ptr[t + 1]],
                          origin, axis, torsions[t])

    # Global rigid-body placement: rotate about the current centroid, then
    # translate that centroid onto `centre`.
    cx = 0.0
    cy = 0.0
    cz = 0.0
    for i in range(n):
        cx += out[i, 0]
        cy += out[i, 1]
        cz += out[i, 2]
    cx /= n
    cy /= n
    cz /= n

    m = quat_to_matrix(quat)
    for i in range(n):
        px = out[i, 0] - cx
        py = out[i, 1] - cy
        pz = out[i, 2] - cz
        out[i, 0] = centre[0] + m[0, 0] * px + m[0, 1] * py + m[0, 2] * pz
        out[i, 1] = centre[1] + m[1, 0] * px + m[1, 1] * py + m[1, 2] * pz
        out[i, 2] = centre[2] + m[2, 0] * px + m[2, 1] * py + m[2, 2] * pz


@njit(**FASTMATH)
def rmsd(a, b):
    """In-place (no superposition) RMSD between two coordinate sets."""
    s = 0.0
    n = a.shape[0]
    for i in range(n):
        dx = a[i, 0] - b[i, 0]
        dy = a[i, 1] - b[i, 1]
        dz = a[i, 2] - b[i, 2]
        s += dx * dx + dy * dy + dz * dz
    return np.sqrt(s / n)


@njit(**FASTMATH)
def symmetry_rmsd(a, b, maps):
    """Best RMSD over a set of atom-index permutations (graph automorphisms)."""
    best = 1e30
    for k in range(maps.shape[0]):
        s = 0.0
        for i in range(maps.shape[1]):
            j = maps[k, i]
            dx = a[i, 0] - b[j, 0]
            dy = a[i, 1] - b[j, 1]
            dz = a[i, 2] - b[j, 2]
            s += dx * dx + dy * dy + dz * dz
        s = np.sqrt(s / maps.shape[1])
        if s < best:
            best = s
    return best


def random_quaternion(rng):
    """Uniformly distributed random rotation (Shoemake's method)."""
    u1, u2, u3 = rng.random(3)
    s1 = np.sqrt(1.0 - u1)
    s2 = np.sqrt(u1)
    return np.array([
        s1 * np.sin(2 * np.pi * u2),
        s1 * np.cos(2 * np.pi * u2),
        s2 * np.sin(2 * np.pi * u3),
        s2 * np.cos(2 * np.pi * u3),
    ])
