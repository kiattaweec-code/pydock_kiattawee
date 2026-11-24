"""Split a protein-ligand complex into receptor.pdb + reference ligand SDF."""
import os
import sys

import numpy as np
from rdkit import Chem
from rdkit.Chem import rdDetermineBonds, rdMolDescriptors
from rdkit.Geometry import Point3D

SRC = sys.argv[1] if len(sys.argv) > 1 else "/tmp/extract/prolif/data/top.pdb"
OUT = sys.argv[2] if len(sys.argv) > 2 else "/home/claude/test_data"
LIGRES = sys.argv[3] if len(sys.argv) > 3 else "LIG"
SOLVENT = {"HOH", "WAT", "TIP3", "SOL", "NA", "CL", "SOD", "CLA", "POT", "K"}

os.makedirs(OUT, exist_ok=True)

lig_lines, rec_lines = [], []
with open(SRC) as fh:
    for line in fh:
        if line[:6] not in ("ATOM  ", "HETATM"):
            continue
        rn = line[17:20].strip().upper()
        if rn == LIGRES:
            lig_lines.append(line)
        elif rn not in SOLVENT:
            rec_lines.append(line)

print(f"receptor atoms {len(rec_lines)}, ligand atoms {len(lig_lines)}")

with open(f"{OUT}/receptor.pdb", "w") as fh:
    fh.writelines(rec_lines)
    fh.write("END\n")

# Heavy-atom-only copy, to exercise the polar-hydrogen placement path.
with open(f"{OUT}/receptor_noH.pdb", "w") as fh:
    for line in rec_lines:
        el = line[76:78].strip() or line[12:16].strip()[:1]
        if el.upper() != "H":
            fh.write(line)
    fh.write("END\n")

# --- ligand: perceive bond orders from the 3D geometry --------------------
els, xyz = [], []
for line in lig_lines:
    el = line[76:78].strip()
    if not el:
        el = "".join(c for c in line[12:16].strip() if c.isalpha())[:1]
    els.append(el.capitalize() if len(el) == 2 else el.upper())
    xyz.append([float(line[30:38]), float(line[38:46]), float(line[46:54])])

xyz_block = f"{len(els)}\nligand\n" + "\n".join(
    f"{e} {p[0]:.4f} {p[1]:.4f} {p[2]:.4f}" for e, p in zip(els, xyz)) + "\n"
mol = Chem.MolFromXYZBlock(xyz_block)
assert mol is not None

ok = False
for charge in (0, 1, -1, 2, -2):
    try:
        probe = Chem.Mol(mol)
        rdDetermineBonds.DetermineBonds(probe, charge=charge)
        Chem.SanitizeMol(probe)
        mol = probe
        ok = True
        print(f"bond perception succeeded with net charge {charge:+d}")
        break
    except Exception as exc:  # noqa: BLE001
        last = exc
if not ok:
    raise SystemExit(f"bond perception failed: {last}")

mol.SetProp("_Name", "reference")
conf = mol.GetConformer()
for i, p in enumerate(xyz):
    conf.SetAtomPosition(i, Point3D(*p))

w = Chem.SDWriter(f"{OUT}/reference.sdf")
w.write(mol)
w.close()

smi = Chem.MolToSmiles(Chem.RemoveHs(mol))
print("ligand SMILES:", smi)
print("formula:", rdMolDescriptors.CalcMolFormula(mol))
print("heavy atoms:", Chem.RemoveHs(mol).GetNumAtoms(),
      " rotatable:", rdMolDescriptors.CalcNumRotatableBonds(mol))

# A start structure with the same topology but a scrambled conformation, so
# the docking run cannot simply keep the answer it was handed.
start = Chem.Mol(mol)
from rdkit.Chem import AllChem  # noqa: E402

AllChem.EmbedMolecule(start, randomSeed=0xBADC0DE)
AllChem.MMFFOptimizeMolecule(start, maxIters=500)
start.SetProp("_Name", "ligand")
w = Chem.SDWriter(f"{OUT}/ligand_start.sdf")
w.write(start)
w.close()

ref_xyz = np.asarray(mol.GetConformer().GetPositions())
print("site centre:", np.round(ref_xyz.mean(axis=0), 2))
print("wrote receptor.pdb, receptor_noH.pdb, reference.sdf, ligand_start.sdf")
