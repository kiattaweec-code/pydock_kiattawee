"""How close can the torsion tree get to the crystal pose at all?

The start conformer's rigid fragments (ring systems, bond angles) come from
ETKDG, not from the crystal, and docking never changes them.  So there is a
floor on the achievable RMSD that has nothing to do with the scoring function
or the search.  This script finds it by optimising the chromosome against RMSD
directly, which is the best any search could possibly do.
"""
import numpy as np
from pydock.optimize import nelder_mead

from pydock import ligand as L
from pydock.geometry import apply_pose, symmetry_rmsd

lig = L.load_ligand("test_data/ligand_start.sdf")
ref = L.load_ligand("test_data/reference.sdf")
assert lig.n_atoms == ref.n_atoms

heavy = lig.heavy_idx
target = np.ascontiguousarray(ref.coords0[ref.heavy_idx])
maps = np.ascontiguousarray(lig.sym_maps)
buf = lig.new_buffer()
n_t = lig.n_rot


def pose_rmsd(x):
    apply_pose(lig.coords0, x[7:7 + n_t], lig.bond_a, lig.bond_b,
               lig.move_ptr, lig.move_idx, x[3:7], x[0:3], buf)
    return symmetry_rmsd(np.ascontiguousarray(buf[heavy]), target, maps)


rng = np.random.default_rng(0)
best = np.inf
best_x = None
centre = target.mean(axis=0)
for trial in range(60):
    x0 = np.empty(7 + n_t)
    x0[0:3] = centre + rng.normal(0, 0.5, 3)
    q = rng.normal(0, 1, 4)
    x0[3:7] = q / np.linalg.norm(q)
    x0[7:] = rng.uniform(-np.pi, np.pi, n_t)
    r = nelder_mead(pose_rmsd, x0, maxiter=6000, maxfev=12000,
                    fatol=1e-4, xatol=1e-4)
    if r.fun < best:
        best, best_x = r.fun, r.x

print(f"best achievable heavy-atom RMSD with this start conformer: {best:.2f} A")

# How much of that is the rigid core alone?
core = lig.move_idx[lig.move_ptr[0]:lig.move_ptr[1]]
rigid = np.setdiff1d(heavy, core)
apply_pose(lig.coords0, best_x[7:7 + n_t], lig.bond_a, lig.bond_b,
           lig.move_ptr, lig.move_idx, best_x[3:7], best_x[0:3], buf)
d = np.linalg.norm(buf[heavy] - ref.coords0[ref.heavy_idx], axis=1)
print(f"per-atom deviation at the optimum: median {np.median(d):.2f} A, "
      f"90th pct {np.percentile(d, 90):.2f} A, max {d.max():.2f} A")
print(f"atoms within 1.0 A: {int((d < 1.0).sum())}/{len(d)}")
