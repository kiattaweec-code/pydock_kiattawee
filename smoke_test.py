"""Stage-by-stage smoke test of the pydock pipeline."""
import time

import numpy as np

import pydock
from pydock import ligand as ligmod
from pydock import receptor as recmod
from pydock import scoring as scoremod

t = time.time()
lig = ligmod.load_ligand("test_data/ligand_start.sdf")
ref = ligmod.load_ligand("test_data/reference.sdf")
print(f"[ligand] {lig.n_atoms} atoms, {lig.n_rot} rotatable genes "
      f"({lig.n_frozen_rot} heavy), {lig.intra_pairs.shape[0]} intra pairs, "
      f"{lig.sym_maps.shape[0]} symmetry maps  [{time.time()-t:.1f}s]")
print(f"         donors={int(lig.is_donor_h.sum())} "
      f"acceptors={int(lig.is_acceptor.sum())} "
      f"lipophilic={int(lig.is_lipophilic.sum())}")
print(f"         torsion tree bonds: {list(zip(lig.bond_a, lig.bond_b))}")

# Torsion tree sanity: changing a torsion must move atoms but preserve bonds.
from pydock.geometry import apply_pose  # noqa: E402

buf = lig.new_buffer()
zero = np.zeros(lig.n_rot)
apply_pose(lig.coords0, zero, lig.bond_a, lig.bond_b, lig.move_ptr,
           lig.move_idx, np.array([1.0, 0, 0, 0]),
           lig.coords0.mean(axis=0), buf)
print(f"[tree]   identity transform max deviation: "
      f"{np.abs(buf - lig.coords0).max():.2e} A")

bonds0 = [(b.GetBeginAtomIdx(), b.GetEndAtomIdx()) for b in lig.mol.GetBonds()]
d0 = np.array([np.linalg.norm(lig.coords0[i] - lig.coords0[j]) for i, j in bonds0])
rng = np.random.default_rng(0)
apply_pose(lig.coords0, rng.uniform(-np.pi, np.pi, lig.n_rot),
           lig.bond_a, lig.bond_b, lig.move_ptr, lig.move_idx,
           np.array([0.3, 0.5, 0.2, 0.1]), np.array([1.0, 2.0, 3.0]), buf)
d1 = np.array([np.linalg.norm(buf[i] - buf[j]) for i, j in bonds0])
print(f"[tree]   after random torsions: bond-length drift "
      f"max {np.abs(d1 - d0).max():.2e} A, "
      f"heavy-atom RMSD from start "
      f"{np.sqrt(((buf[lig.heavy_idx]-lig.coords0[lig.heavy_idx])**2).sum(1).mean()):.2f} A")

for rec_file in ("test_data/receptor.pdb", "test_data/receptor_noH.pdb"):
    t = time.time()
    rec = recmod.prepare_receptor(rec_file)
    print(f"[rec]    {rec_file}: {rec.n_atoms} atoms, "
          f"donorH={int(rec.is_donor_h.sum())} "
          f"acc={int(rec.is_acceptor.sum())} "
          f"lipo={int(rec.is_lipophilic.sum())} "
          f"rotOH={rec.rot_h_axis_b.size}  [{time.time()-t:.1f}s]")

centre, half = recmod.site_from_ligand(ref.coords0[ref.heavy_idx], padding=4.0)
print(f"[site]   centre {np.round(centre,2)} half {np.round(half,2)}")

rec = recmod.prepare_receptor("test_data/receptor_noH.pdb")
for kind in ("empirical", "chemplp", "forcefield"):
    t = time.time()
    sc = scoremod.build_scorer(kind, rec, lig, centre, half, spacing=0.4)
    tb = time.time() - t
    s = sc.score(ref.coords0, np.zeros(sc.n_rot_h), breakdown=True)
    t = time.time()
    n = 300
    for _ in range(n):
        sc.score(ref.coords0, np.zeros(sc.n_rot_h))
    per = (time.time() - t) / n * 1e6
    print(f"[score]  {kind:11s} grids {tb:5.1f}s  maps={len(sc.grid.map_names)} "
          f"{sc.grid.shape}  rotOH genes={sc.n_rot_h}")
    print(f"         crystal pose score = {s:9.3f}   "
          f"{per:6.0f} us/eval")
    print("         terms:", {k: round(float(v), 3) for k, v in sc.terms().items()})
