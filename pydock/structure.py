"""Pull a downloaded complex apart: protein, ligand, waters, metals, ions.

A PDB entry straight from the database is not a docking receptor.  It carries
the bound ligand, crystallisation waters, buffer molecules and often several
copies of the whole assembly, and `pydock.receptor` keeps every HETATM it is
given -- so a ligand left in the file is treated as part of the protein and
fills the pocket it is meant to reveal.

This module does the splitting: :func:`inspect_structure` lists what is in a
file, :func:`extract_residue` lifts one group out as a molecule with bond
orders worked out from its geometry, and :func:`build_receptor` writes the
protein back without the parts you chose to drop.
"""
from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem, rdDetermineBonds, rdMolDescriptors

from .params import METALS
from .receptor import STANDARD_AA, canonical_resname

WATERS = {"HOH", "WAT", "DOD", "TIP3", "TIP4", "SOL", "H2O", "OH2", "T3P"}

NUCLEIC = {"A", "C", "G", "U", "T", "DA", "DC", "DG", "DT", "DU",
           "RA", "RC", "RG", "RU"}

# Common crystallisation additives: present in the file, almost never the
# thing being studied.  Flagged so they can be dropped without inspection.
BUFFERS = {
    "GOL", "EDO", "PEG", "PG4", "PGE", "1PE", "2PE", "P6G", "DMS", "DMSO",
    "MPD", "SO4", "PO4", "ACT", "ACY", "FMT", "EPE", "TRS", "MES", "CIT",
    "FLC", "TAR", "MLI", "OXL", "BME", "DTT", "IMD", "AZI", "NO3", "CO3",
    "BCT", "IOD", "BR", "CL", "FLU", "NH4", "ACE", "NAG", "BMA", "MAN",
    "XYP", "SIA", "FUC", "GAL", "GLC", "SUC", "BOG", "LDA", "C8E", "OLC",
}

# An ion is a single atom that is not one of the metals we keep for scoring.
HALIDES = {"CL", "BR", "I", "F", "IOD"}


@dataclass
class Component:
    """One residue-sized piece of the file."""

    key: str                      # "A/999/LIG" -- stable id for the UI
    resname: str
    chain: str
    resseq: int
    icode: str
    kind: str                     # polymer | ligand | water | metal | ion | buffer
    n_atoms: int
    n_heavy: int
    formula: str
    centre: tuple[float, float, float]
    elements: dict[str, int] = field(default_factory=dict)

    @property
    def label(self) -> str:
        seq = f"{self.resseq}{self.icode}".strip()
        return f"{self.resname} {self.chain}:{seq}"


@dataclass
class Structure:
    path: str
    records: list[dict]
    components: list[Component]
    chains: dict[str, dict]
    n_atoms: int
    n_hydrogens: int
    has_hydrogens: bool

    def ligand_candidates(self) -> list[Component]:
        """Plausible bound ligands, biggest first."""
        return sorted((c for c in self.components if c.kind == "ligand"),
                      key=lambda c: -c.n_heavy)

    def by_key(self, key: str) -> Component | None:
        return next((c for c in self.components if c.key == key), None)


# --------------------------------------------------------------------------- #
# Reading
# --------------------------------------------------------------------------- #
def _element_of(line: str) -> str:
    el = line[76:78].strip()
    if not el:
        name = line[12:16].strip()
        el = "".join(ch for ch in name if ch.isalpha())[:1]
    return el.capitalize() if len(el) == 2 else el.upper()


def read_records(path: str, model: int = 1) -> list[dict]:
    """Every ATOM/HETATM line of one model, keeping the raw text."""
    out: list[dict] = []
    current_model = 1
    with open(path, errors="replace") as fh:
        for line in fh:
            tag = line[:6]
            if tag == "MODEL ":
                try:
                    current_model = int(line[10:14])
                except ValueError:
                    current_model += 1
                continue
            if tag == "ENDMDL":
                if current_model >= model:
                    break
                continue
            if tag not in ("ATOM  ", "HETATM"):
                continue
            if current_model != model:
                continue
            altloc = line[16]
            if altloc not in (" ", "A", "1"):
                continue                      # keep one alternate conformation
            try:
                resseq = int(line[22:26])
                xyz = (float(line[30:38]), float(line[38:46]), float(line[46:54]))
            except ValueError:
                continue
            out.append({
                "line": line.rstrip("\n"),
                "hetatm": tag == "HETATM",
                "name": line[12:16].strip(),
                "resname": line[17:20].strip().upper(),
                "chain": (line[21].strip() or "A"),
                "resseq": resseq,
                "icode": line[26].strip(),
                "element": _element_of(line),
                "xyz": xyz,
            })
    if not out:
        raise ValueError(f"no atom records found in {path}")
    return out


def _classify(resname: str, records: list[dict]) -> str:
    if resname in WATERS:
        return "water"
    # Through canonical_resname first: force fields rename residues (HSD and
    # HSE for histidine, CYX for a disulfide cysteine), and without this a
    # CHARMM-prepared protein offers every histidine as a ligand candidate.
    if (resname in STANDARD_AA or canonical_resname(resname) in STANDARD_AA
            or resname in NUCLEIC):
        return "polymer"
    heavy = [r for r in records if r["element"] != "H"]
    if len(heavy) == 1:
        el = heavy[0]["element"].upper()
        if el in METALS or resname.upper() in METALS:
            return "metal"
        if el in HALIDES or resname.upper() in HALIDES:
            return "ion"
        return "ion"
    if resname in BUFFERS:
        return "buffer"
    return "ligand"


def inspect_structure(path: str, model: int = 1) -> Structure:
    """What is in this file, grouped by residue."""
    records = read_records(path, model=model)

    groups: dict[tuple, list[dict]] = {}
    for r in records:
        groups.setdefault((r["chain"], r["resseq"], r["icode"], r["resname"]),
                          []).append(r)

    components: list[Component] = []
    chains: dict[str, dict] = {}
    for (chain, resseq, icode, resname), rs in groups.items():
        kind = _classify(resname, rs)
        heavy = [r for r in rs if r["element"] != "H"]
        counts: dict[str, int] = {}
        for r in heavy:
            counts[r["element"]] = counts.get(r["element"], 0) + 1
        xyz = np.array([r["xyz"] for r in rs], dtype=float)
        components.append(Component(
            key=f"{chain}/{resseq}{icode}/{resname}",
            resname=resname, chain=chain, resseq=resseq, icode=icode,
            kind=kind, n_atoms=len(rs), n_heavy=len(heavy),
            formula="".join(f"{e}{n if n > 1 else ''}"
                            for e, n in sorted(counts.items())),
            centre=tuple(round(float(v), 2) for v in xyz.mean(axis=0)),
            elements=counts,
        ))
        info = chains.setdefault(chain, {"residues": 0, "atoms": 0,
                                         "polymer_residues": 0})
        info["residues"] += 1
        info["atoms"] += len(rs)
        if kind == "polymer":
            info["polymer_residues"] += 1

    n_h = sum(1 for r in records if r["element"] == "H")
    components.sort(key=lambda c: (c.kind != "ligand", -c.n_heavy, c.key))
    return Structure(path=path, records=records, components=components,
                     chains=chains, n_atoms=len(records), n_hydrogens=n_h,
                     has_hydrogens=n_h > 0)


# --------------------------------------------------------------------------- #
# Extracting one group as a molecule
# --------------------------------------------------------------------------- #
def _name_new_hydrogens(mol: Chem.Mol, comp: Component) -> None:
    """Give added hydrogens the residue of the atom they hang off.

    Without this they carry no residue information at all, and the MOL2
    writer files them under a second substructure -- one ligand arriving in
    another program as two.
    """
    counter: dict[int, int] = {}
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() != 1 or atom.GetPDBResidueInfo() is not None:
            continue
        heavy = next((n for n in atom.GetNeighbors()
                      if n.GetAtomicNum() != 1), None)
        parent = heavy.GetPDBResidueInfo() if heavy is not None else None
        idx = heavy.GetIdx() if heavy is not None else -1
        counter[idx] = counter.get(idx, 0) + 1
        base = (parent.GetName().strip() if parent is not None
                else comp.resname)[:2]
        info = Chem.AtomPDBResidueInfo()
        info.SetName(f"H{base}{counter[idx]}".ljust(4)[:4])
        info.SetResidueName(parent.GetResidueName() if parent is not None
                            else comp.resname)
        info.SetResidueNumber(parent.GetResidueNumber() if parent is not None
                              else comp.resseq)
        info.SetChainId(parent.GetChainId() if parent is not None
                        else comp.chain)
        info.SetIsHeteroAtom(True)
        atom.SetMonomerInfo(info)


def extract_residue(structure: Structure, key: str, smiles: str | None = None,
                    add_hydrogens: bool = True) -> Chem.Mol:
    """Lift one component out as a molecule with bond orders.

    A PDB file records no bond orders and no protonation, so they have to
    come from somewhere:

    * ``smiles`` -- the reliable route.  Bond orders are transferred onto the
      crystal coordinates from a template, so the result is exactly the
      compound you say it is.  Get the SMILES from the ligand's page on RCSB
      or from PubChem.
    * otherwise perception from geometry, which works for small, neutral,
      mostly-saturated ligands and fails on drug-sized ones: for a 43-heavy-
      atom ergoline no net charge gives a valid Kekule structure.

    The result carries ``pydock_bond_orders`` saying which route was used:
    ``template``, ``perceived`` or ``connectivity`` -- the last meaning every
    bond is single and the file is good for geometry and nothing else.
    """
    comp = structure.by_key(key)
    if comp is None:
        raise ValueError(f"no such component: {key}")

    rs = [r for r in structure.records
          if r["chain"] == comp.chain and r["resseq"] == comp.resseq
          and r["icode"] == comp.icode and r["resname"] == comp.resname]
    if not rs:
        raise ValueError(f"{key} has no atoms")

    # Through a PDB block, which keeps atom names and residue information and
    # perceives connectivity by distance.
    block = "\n".join(r["line"] for r in rs) + "\nEND\n"
    raw = Chem.MolFromPDBBlock(block, removeHs=False, sanitize=False,
                               proximityBonding=True)
    if raw is None or raw.GetNumAtoms() == 0:
        raise ValueError(f"could not read the coordinates of {comp.label}")

    mol, provenance, note = None, "connectivity", ""

    if smiles:
        template = Chem.MolFromSmiles(smiles)
        if template is None:
            raise ValueError("that SMILES could not be parsed")
        # Checked before use, because AssignBondOrdersFromTemplate is happy
        # with a template that matches only a fragment: ethanol against a
        # 43-atom ergoline "succeeds", assigns three bonds, and leaves the
        # rest single -- a confidently wrong molecule.
        n_template = Chem.RemoveHs(template).GetNumAtoms()
        if n_template != comp.n_heavy:
            raise ValueError(
                f"that SMILES has {n_template} heavy atoms but "
                f"{comp.label} has {comp.n_heavy} ({comp.formula}); "
                f"it is not the same compound.")
        try:
            mol = AllChem.AssignBondOrdersFromTemplate(template, raw)
            Chem.SanitizeMol(mol)
            provenance = "template"
            # And confirmed afterwards: the same atom counts can still be a
            # different compound, and the formula is cheap to compare.
            got = rdMolDescriptors.CalcMolFormula(Chem.RemoveHs(mol))
            want = rdMolDescriptors.CalcMolFormula(Chem.RemoveHs(template))
            if got != want:
                raise ValueError(
                    f"the molecule built from that SMILES is {got}, but the "
                    f"SMILES itself is {want}")
        except Exception as exc:                                  # noqa: BLE001
            raise ValueError(
                f"the SMILES does not match {comp.label}: {exc}. "
                f"It has {comp.n_heavy} heavy atoms "
                f"({comp.formula}); check you have the right compound.")

    if mol is None:
        last: Exception | None = None
        for charge in (0, -1, 1, -2, 2, -3, 3):
            try:
                probe = Chem.Mol(raw)
                rdDetermineBonds.DetermineBondOrders(probe, charge=charge)
                Chem.SanitizeMol(probe)
                mol = probe
                provenance = "perceived"
                break
            except Exception as exc:                              # noqa: BLE001
                last = exc
        if mol is None:
            # Keep the connectivity rather than failing: the geometry is
            # still right, and it is enough to place a docking box or to
            # measure an RMSD.  The caller is told not to trust the orders.
            mol = Chem.Mol(raw)
            Chem.SanitizeMol(mol, Chem.SanitizeFlags.SANITIZE_ALL
                             ^ Chem.SanitizeFlags.SANITIZE_PROPERTIES
                             ^ Chem.SanitizeFlags.SANITIZE_KEKULIZE,
                             catchErrors=True)
            note = (f"bond orders could not be worked out from the geometry "
                    f"({last}); every bond is written as single. Paste the "
                    f"ligand's SMILES to get this right.")

    # Not when the orders are unknown: saturating every single bond invents
    # hydrogens that are not there (57 instead of 36 on an ergoline), which
    # looks like a real molecule and is not one.
    if add_hydrogens and provenance != "connectivity":
        try:
            mol = Chem.AddHs(mol, addCoords=True)
            _name_new_hydrogens(mol, comp)
        except Exception:                                         # noqa: BLE001
            note = (note + " Hydrogens could not be added.").strip()

    # From the coordinates, not the template: the crystal structure is the
    # evidence for which stereoisomer is actually bound, and a template
    # SMILES is often drawn without stereochemistry or with the wrong one.
    try:
        Chem.AssignStereochemistryFrom3D(mol)
    except Exception:                                             # noqa: BLE001
        pass

    mol.SetProp("_Name", comp.label)
    mol.SetProp("pydock_source", structure.path.rsplit("/", 1)[-1])
    mol.SetProp("pydock_bond_orders", provenance)
    try:
        if provenance == "connectivity":
            # Counting implicit hydrogens here would report the valences of
            # an all-single-bond skeleton, not the atoms in the file.
            counts: dict[str, int] = {}
            for atom in mol.GetAtoms():
                counts[atom.GetSymbol()] = counts.get(atom.GetSymbol(), 0) + 1
            mol.SetProp("pydock_formula", "".join(
                f"{e}{n if n > 1 else ''}" for e, n in sorted(counts.items())))
        else:
            mol.SetProp("pydock_formula", rdMolDescriptors.CalcMolFormula(mol))
    except Exception:                                             # noqa: BLE001
        pass
    if note:
        mol.SetProp("pydock_warning", note)
    return mol


# --------------------------------------------------------------------------- #
# Writing the receptor back out
# --------------------------------------------------------------------------- #
def build_receptor(structure: Structure, drop_keys: set[str] | None = None,
                   remove_waters: bool = True, remove_hydrogens: bool = True,
                   keep_metals: bool = True, remove_buffers: bool = True,
                   remove_ions: bool = False,
                   chains: set[str] | None = None) -> tuple[str, dict]:
    """Return (PDB text, summary of what was dropped).

    Hydrogens are removed by default on purpose: pydock places polar
    hydrogens itself, and it only treats hydroxyls as rotatable when the
    file arrives without any -- so a protonated receptor silently loses that
    part of the search.
    """
    drop_keys = set(drop_keys or ())
    kept: list[str] = []
    dropped: dict[str, int] = {}
    kept_kinds: dict[str, int] = {}

    def note(bucket: dict, kind: str):
        bucket[kind] = bucket.get(kind, 0) + 1

    by_key = {c.key: c for c in structure.components}
    for r in structure.records:
        key = f"{r['chain']}/{r['resseq']}{r['icode']}/{r['resname']}"
        comp = by_key.get(key)
        kind = comp.kind if comp else "polymer"

        if r["element"] == "H" and remove_hydrogens:
            note(dropped, "hydrogen")
            continue
        if key in drop_keys:
            note(dropped, kind if kind != "polymer" else "selected residue")
            continue
        if kind == "water" and remove_waters:
            note(dropped, "water")
            continue
        if kind == "buffer" and remove_buffers:
            note(dropped, "buffer or cryoprotectant")
            continue
        if kind == "metal" and not keep_metals:
            note(dropped, "metal")
            continue
        if kind == "ion" and remove_ions:
            note(dropped, "ion")
            continue
        if chains and r["chain"] not in chains:
            note(dropped, "other chain")
            continue

        kept.append(r["line"])
        note(kept_kinds, kind)

    if not kept:
        raise ValueError("every atom was removed; loosen the options")

    text = "\n".join(kept) + "\nEND\n"
    return text, {"atoms_kept": len(kept), "dropped": dropped,
                  "kept": kept_kinds}


# --------------------------------------------------------------------------- #
# Putting a pose back with its protein
# --------------------------------------------------------------------------- #
def merge_complex(receptor_pdb: str, ligand: Chem.Mol, chain: str = "Z",
                  resname: str = "LIG", resseq: int = 900,
                  title: str = "") -> str:
    """One PDB holding the receptor and a docked pose, ready for PyMOL.

    The pose comes out of docking as an SDF and the protein stays a PDB, so
    looking at the two together means loading both and hoping they share a
    frame.  They do -- docking never moves the receptor -- but the pair is
    awkward to keep together, mail, or attach to a paper.

    The ligand is given its own chain and residue so it can be selected on
    its own, and its CONECT records are renumbered onto the merged serials;
    without them a viewer guesses the bonds by distance and draws aromatic
    rings as whatever the distances suggest.
    """
    lines = [ln.rstrip("\n") for ln in receptor_pdb.splitlines()
             if ln[:6] in ("ATOM  ", "HETATM")]
    if not lines:
        raise ValueError("the receptor has no atom records")

    last_serial = 0
    for ln in lines:
        try:
            last_serial = max(last_serial, int(ln[6:11]))
        except ValueError:
            pass

    # Residue information decides what the ligand is called on screen.
    probe = Chem.Mol(ligand)
    for i, atom in enumerate(probe.GetAtoms(), start=1):
        info = Chem.AtomPDBResidueInfo()
        info.SetName(f"{atom.GetSymbol()}{i}".ljust(4)[:4])
        info.SetResidueName(resname[:3].ljust(3))
        info.SetResidueNumber(resseq)
        info.SetChainId(chain[:1])
        info.SetIsHeteroAtom(True)
        info.SetOccupancy(1.0)
        info.SetTempFactor(0.0)
        atom.SetMonomerInfo(info)

    block = Chem.MolToPDBBlock(probe, flavor=4)   # 4: write CONECT records

    old_to_new: dict[int, int] = {}
    hetatm: list[str] = []
    conect: list[str] = []
    for ln in block.splitlines():
        if ln[:6] in ("ATOM  ", "HETATM"):
            old = int(ln[6:11])
            new = last_serial + len(hetatm) + 1
            old_to_new[old] = new
            hetatm.append("HETATM" + f"{new:5d}" + ln[11:].rstrip())
        elif ln.startswith("CONECT"):
            conect.append(ln.rstrip())

    remapped = []
    for ln in conect:
        serials = [ln[i:i + 5].strip() for i in range(6, len(ln.rstrip()), 5)]
        new = [old_to_new.get(int(s)) for s in serials if s]
        new = [n for n in new if n is not None]
        if len(new) > 1:
            remapped.append("CONECT" + "".join(f"{n:5d}" for n in new))

    out = []
    if title:
        out.append(f"REMARK   1 {title[:60]}")
    out.append("REMARK   1 receptor and docked pose, merged by pydock")
    out += lines
    out.append("TER")
    out += hetatm
    out += remapped
    out.append("END")
    return "\n".join(out) + "\n"
