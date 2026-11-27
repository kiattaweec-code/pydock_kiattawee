"""The MOL2 writer must produce files that read back as the same molecule.

RDKit has no MOL2 writer, so pydock.mol2 is ours, and a wrong Sybyl atom or
bond type is the kind of error that is invisible until someone opens the file
in another program.  Every case below is written, read back with RDKit's MOL2
reader, and compared.

    python packaging/test_mol2.py
"""
from __future__ import annotations

import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

from rdkit import Chem, RDLogger                                  # noqa: E402
from rdkit.Chem import AllChem                                    # noqa: E402

from pydock.mol2 import round_trips, sybyl_bond, sybyl_type, to_mol2  # noqa: E402

RDLogger.DisableLog("rdApp.*")

# Each is here because it exercises a type rule that is easy to get wrong.
CASES = [
    ("benzene",      "c1ccccc1"),                 # aromatic carbon
    ("pyridine",     "c1ccncc1"),                 # aromatic N, 6-membered
    ("imidazole",    "c1cnc[nH]1"),               # aromatic N, 5-membered
    ("indole",       "c1ccc2[nH]ccc2c1"),         # fused 5-membered
    ("thiophene",    "c1ccsc1"),                  # aromatic S
    ("aspirin",      "CC(=O)Oc1ccccc1C(=O)O"),    # ester + carboxylic acid
    ("benzoate",     "[O-]C(=O)c1ccccc1"),        # O.co2
    ("caffeine",     "CN1C=NC2=C1C(=O)N(C)C(=O)N2C"),
    ("acetanilide",  "CC(=O)Nc1ccc(N)cc1"),       # N.am and N.pl3 together
    ("arginine",     "NC(CCCNC(N)=[NH2+])C(=O)[O-]"),   # C.cat
    ("guanidine",    "NC(=[NH2+])N"),
    ("nitrobenzene", "O=[N+]([O-])c1ccccc1"),     # must not be O.co2
    ("sulfonamide",  "CS(=O)(=O)Nc1ccccc1"),      # S.o2
    ("sulfoxide",    "CS(=O)c1ccccc1"),           # S.o
    ("phosphate",    "OP(=O)(O)OCC"),
    ("alkyne",       "CC#CC(=O)O"),               # C.1
    ("imatinib",     "Cc1ccc(NC(=O)c2ccc(CN3CCN(C)CC3)cc2)cc1"
                     "Nc1nccc(-c2cccnc2)n1"),
]

ok, bad = [], []


def embed(smiles: str) -> Chem.Mol | None:
    mol = Chem.AddHs(Chem.MolFromSmiles(smiles))
    if AllChem.EmbedMolecule(mol, randomSeed=7) != 0:
        return None
    return mol


def main() -> int:
    for name, smiles in CASES:
        mol = embed(smiles)
        if mol is None:
            bad.append(f"{name}: could not be embedded in 3D")
            continue
        good, detail = round_trips(mol)
        note = f"  ({detail})" if good and detail.startswith("same skeleton") else ""
        (ok if good else bad).append(
            f"{name}{note}" if good else f"{name}: {detail}")

    # Specific types, so a silent regression to all-C.3 is caught.
    checks = [
        ("acetanilide amide N is N.am", "CC(=O)Nc1ccccc1",
         lambda m: "N.am" in {sybyl_type(a) for a in m.GetAtoms()}),
        ("benzene carbon is C.ar", "c1ccccc1",
         lambda m: {sybyl_type(a) for a in m.GetAtoms() if a.GetSymbol() == "C"} == {"C.ar"}),
        ("benzoate oxygens are O.co2", "[O-]C(=O)c1ccccc1",
         lambda m: sum(sybyl_type(a) == "O.co2" for a in m.GetAtoms()) == 2),
        ("nitro oxygens are not O.co2", "O=[N+]([O-])c1ccccc1",
         lambda m: not any(sybyl_type(a) == "O.co2" for a in m.GetAtoms())),
        ("arginine has C.cat", "NC(CCCNC(N)=[NH2+])C(=O)[O-]",
         lambda m: any(sybyl_type(a) == "C.cat" for a in m.GetAtoms())),
        ("sulfone S is S.o2", "CS(=O)(=O)N",
         lambda m: any(sybyl_type(a) == "S.o2" for a in m.GetAtoms())),
        ("amide bond is am", "CC(=O)Nc1ccccc1",
         lambda m: "am" in {sybyl_bond(b) for b in m.GetBonds()}),
        ("aromatic rings are written kekulised", "c1ccccc1",
         lambda m: "ar" not in {l.split()[3] for l in to_mol2(m).split("@<TRIPOS>BOND")[1].strip().splitlines()}),
    ]
    for label, smiles, test in checks:
        mol = embed(smiles)
        try:
            if mol is not None and test(mol):
                ok.append(label)
            else:
                bad.append(f"{label}: not satisfied")
        except Exception as exc:                                  # noqa: BLE001
            bad.append(f"{label}: {type(exc).__name__}: {exc}")

    # Negative control: the comparison must notice a corrupted file.
    mol = embed("CCO")
    corrupted = Chem.MolFromMol2Block(to_mol2(mol).replace("O.3", "N.3", 1),
                                      removeHs=False, sanitize=True)
    changed = corrupted is None or (
        Chem.MolToSmiles(Chem.RemoveHs(mol))
        != Chem.MolToSmiles(Chem.RemoveHs(corrupted)))
    (ok if changed else bad).append(
        "a corrupted atom type is detected" if changed
        else "a corrupted atom type went unnoticed")

    for line in ok:
        print("  ok   " + line)
    for line in bad:
        print("  FAIL " + line)
    print(f"\n{len(ok)} passed, {len(bad)} failed\n")
    return 1 if bad else 0


if __name__ == "__main__":
    raise SystemExit(main())
