"""Redocking validation: both scoring functions, RMSD against the crystal pose."""
import sys
import time

import numpy as np

from pydock import DockingConfig, GAConfig, run_docking
from pydock import ligand as L
from pydock.geometry import symmetry_rmsd

REC = "test_data/receptor_noH.pdb"
LIG = "test_data/ligand_start.sdf"
REF = "test_data/reference.sdf"

n_evals = int(sys.argv[1]) if len(sys.argv) > 1 else 250_000
n_runs = int(sys.argv[2]) if len(sys.argv) > 2 else 8
kinds = sys.argv[3].split(",") if len(sys.argv) > 3 else ["empirical", "forcefield"]

ref = L.load_ligand(REF)
ref_heavy = np.ascontiguousarray(ref.coords0[ref.heavy_idx])

for kind in kinds:
    print("=" * 72)
    print(f"scoring function: {kind}   ({n_runs} runs x {n_evals} evals)")
    print("=" * 72)
    cfg = DockingConfig(
        receptor=REC, ligands=[LIG], reference_ligand=REF,
        scoring=kind, n_runs=n_runs, n_poses=n_runs,
        output_dir=f"out_{kind}",
        ga=GAConfig(pop_size=100, n_islands=5, seed=11,
                    autoscale=False, n_evals=n_evals,
                    early_stop=n_evals, ls_frequency=0.10),
    )
    t = time.time()
    res = run_docking(cfg)[0]
    wall = time.time() - t

    lig = L.load_ligand(LIG)
    rows = []
    for s in res.solutions:
        r = float(symmetry_rmsd(np.ascontiguousarray(s.coords[lig.heavy_idx]),
                                ref_heavy, np.ascontiguousarray(lig.sym_maps)))
        rows.append((s.run, s.score, r, s.terms))

    print(f"\n  {'run':>4} {'score':>10} {'rmsd':>7}   terms")
    for run, score, r, terms in rows:
        extra = " ".join(f"{k}={float(v):.2f}" for k, v in terms.items()
                         if k in ("hbond", "lipo", "clash", "vdw", "elec", "desolv"))
        print(f"  {run:>4} {score:10.3f} {r:7.2f}   {extra}")

    top_rmsd = rows[0][2]
    best_rmsd = min(r for _, _, r, _ in rows)
    print(f"\n  top-ranked pose RMSD : {top_rmsd:.2f} A"
          f"   {'SUCCESS' if top_rmsd < 2.0 else 'miss'} (< 2.0 A criterion)")
    print(f"  best pose found RMSD : {best_rmsd:.2f} A")
    print(f"  distinct binding modes: {len(res.clusters)}")
    print(f"  wall time            : {wall:.0f}s "
          f"({wall / n_runs:.0f}s per run)\n")
