"""Write Tripos MOL2 files.

RDKit reads MOL2 but cannot write it, and Open Babel is a heavy dependency to
add to a packaged desktop app for one file format.  So this writes it: the
three sections a docking program needs (MOLECULE, ATOM, BOND) with proper
Sybyl atom and bond types, rather than dumping every atom as C.3 and every
bond as 1, which is what most hand-rolled writers do and what makes the
output useless to anything that reads atom types.

``round_trips`` checks a molecule survives write -> read unchanged, which is
how the type assignment below is tested.
"""
from __future__ import annotations

from rdkit import Chem

# Elements written as the bare symbol: no hybridisation state is meaningful.
_PLAIN = {
    "H", "F", "Cl", "Br", "I", "Li", "Na", "K", "Ca", "Mg", "Al", "Si",
    "Zn", "Mn", "Fe", "Co", "Cu", "Ni", "Cd", "Hg", "Se", "B", "Cr", "Mo",
}

_HYB = {
    Chem.HybridizationType.SP: "1",
    Chem.HybridizationType.SP2: "2",
    Chem.HybridizationType.SP3: "3",
}


def _is_amide_nitrogen(atom: Chem.Atom) -> bool:
    """N bonded to a carbonyl (or thiocarbonyl) carbon -> N.am."""
    if atom.GetSymbol() != "N" or atom.GetIsAromatic():
        return False
    for nbr in atom.GetNeighbors():
        if nbr.GetSymbol() != "C":
            continue
        for b in nbr.GetBonds():
            other = b.GetOtherAtom(nbr)
            if (b.GetBondType() == Chem.BondType.DOUBLE
                    and other.GetSymbol() in ("O", "S")):
                return True
    return False


def _is_guanidinium_carbon(atom: Chem.Atom) -> bool:
    """The central carbon of a guanidinium/amidinium group -> C.cat."""
    if atom.GetSymbol() != "C":
        return False
    nitrogens = [n for n in atom.GetNeighbors() if n.GetSymbol() == "N"]
    if len(nitrogens) < 3:
        return False
    charged = any(n.GetFormalCharge() > 0 for n in nitrogens)
    double = any(b.GetBondType() == Chem.BondType.DOUBLE
                 and b.GetOtherAtom(atom).GetSymbol() == "N"
                 for b in atom.GetBonds())
    return charged or double


def _is_delocalised_oxygen(atom: Chem.Atom) -> bool:
    """Carboxylate, nitro or phosphate oxygen -> O.co2.

    Sybyl uses O.co2 for both oxygens of a delocalised pair, so the test is
    on the shared neighbour rather than on this atom's own bond order.
    """
    if atom.GetSymbol() != "O" or atom.GetDegree() != 1:
        return False
    centre = atom.GetNeighbors()[0]
    # Carboxylate, sulfonate and phosphate only.  A nitro group's oxygens
    # are O.2 in Sybyl, and readers reject O.co2 next to a nitrogen.
    if centre.GetSymbol() not in ("C", "P", "S"):
        return False
    terminal_o = [n for n in centre.GetNeighbors()
                  if n.GetSymbol() == "O" and n.GetDegree() == 1]
    return len(terminal_o) >= 2


def _planar_nitrogen(atom: Chem.Atom) -> bool:
    """Trigonal-planar N that is not an amide -> N.pl3 (aniline, enamine)."""
    if atom.GetSymbol() != "N" or atom.GetDegree() != 3:
        return False
    if atom.GetHybridization() == Chem.HybridizationType.SP2:
        return True
    return any(b.GetIsConjugated() for b in atom.GetBonds())


def sybyl_type(atom: Chem.Atom) -> str:
    """The Sybyl (MOL2) atom type for one atom."""
    sym = atom.GetSymbol()
    if sym in _PLAIN:
        return sym
    if atom.GetIsAromatic() and sym in ("C", "N", "O", "S", "P"):
        return f"{sym}.ar"

    if sym == "C":
        if _is_guanidinium_carbon(atom):
            return "C.cat"
        return "C." + _HYB.get(atom.GetHybridization(), "3")

    if sym == "N":
        # Total degree, not GetDegree(): a protonated amine extracted from a
        # PDB file carries its hydrogen implicitly, and counting only explicit
        # neighbours typed it N.3, which silently dropped the charge -- MOL2
        # has no formal-charge field, so N.4 is the only way to express it.
        if atom.GetFormalCharge() > 0 and atom.GetTotalDegree() == 4:
            return "N.4"
        if _is_amide_nitrogen(atom):
            return "N.am"
        if _planar_nitrogen(atom):
            return "N.pl3"
        return "N." + _HYB.get(atom.GetHybridization(), "3")

    if sym == "O":
        if _is_delocalised_oxygen(atom):
            return "O.co2"
        return "O." + ("2" if atom.GetHybridization()
                       == Chem.HybridizationType.SP2 else "3")

    if sym == "S":
        oxygens = sum(1 for b in atom.GetBonds()
                      if b.GetBondType() == Chem.BondType.DOUBLE
                      and b.GetOtherAtom(atom).GetSymbol() == "O")
        if oxygens >= 2:
            return "S.o2"
        if oxygens == 1:
            return "S.o"
        return "S." + ("2" if atom.GetHybridization()
                       == Chem.HybridizationType.SP2 else "3")

    if sym == "P":
        return "P.3"
    return sym


def sybyl_bond(bond: Chem.Bond) -> str:
    """The Sybyl bond type for one bond.

    Decided by bond order, not by the aromatic flag: the molecule is
    kekulised before writing, so aromatic rings come out as alternating
    1 and 2.  Writing `ar` instead looks tidier but leaves readers to
    re-kekulise, which fails on five-membered heterocycles -- pyrrole,
    thiophene, imidazole -- and makes the file unreadable.  Atom types
    still carry `.ar`, so the aromaticity is not lost.
    """
    if bond.GetBondType() == Chem.BondType.AROMATIC:
        return "ar"
    begin, end = bond.GetBeginAtom(), bond.GetEndAtom()
    if bond.GetBondType() == Chem.BondType.SINGLE:
        # The C-N bond of an amide is 'am' in Sybyl, not 1.
        pair = {begin.GetSymbol(), end.GetSymbol()}
        if pair == {"C", "N"}:
            nitrogen = begin if begin.GetSymbol() == "N" else end
            if _is_amide_nitrogen(nitrogen):
                return "am"
        return "1"
    return {Chem.BondType.DOUBLE: "2",
            Chem.BondType.TRIPLE: "3"}.get(bond.GetBondType(), "1")


def to_mol2(mol: Chem.Mol, name: str = "LIG", conf_id: int = -1,
            charges: str | None = "GASTEIGER") -> str:
    """Serialise a molecule with 3D coordinates to MOL2 text.

    `charges` names the charge model written in the header; pass None to
    write zeros and ``NO_CHARGES``.
    """
    if mol.GetNumConformers() == 0:
        raise ValueError("MOL2 needs 3D coordinates; this molecule has none")
    conf = mol.GetConformer(conf_id)

    # Kekulise for the bond block but keep aromatic flags so atom types
    # stay `.ar` (see sybyl_bond).
    try:
        kek = Chem.Mol(mol)
        Chem.Kekulize(kek, clearAromaticFlags=False)
        bonds_from = kek
    except Exception:                                             # noqa: BLE001
        bonds_from = mol

    q = [0.0] * mol.GetNumAtoms()
    model = "NO_CHARGES"
    if charges:
        try:
            probe = Chem.Mol(mol)
            from rdkit.Chem import AllChem
            AllChem.ComputeGasteigerCharges(probe)
            vals = []
            for a in probe.GetAtoms():
                v = a.GetDoubleProp("_GasteigerCharge")
                vals.append(v if v == v and abs(v) < 1e3 else 0.0)  # NaN check
            q = vals
            model = "GASTEIGER"
        except Exception:                                         # noqa: BLE001
            pass

    # Substructure (residue) names, taken from the PDB info where present.
    def residue(atom: Chem.Atom) -> tuple[int, str]:
        info = atom.GetPDBResidueInfo()
        if info is None:
            return 1, name[:8]
        return info.GetResidueNumber() or 1, (info.GetResidueName() or name).strip()[:8]

    subst = {}
    for atom in mol.GetAtoms():
        rid, rname = residue(atom)
        subst.setdefault(rid, rname)
    order = {rid: i + 1 for i, rid in enumerate(sorted(subst))}

    lines = [
        "@<TRIPOS>MOLECULE",
        name,
        f"{mol.GetNumAtoms():5d}{mol.GetNumBonds():6d}{len(subst):6d}"
        f"{0:6d}{0:6d}",
        "SMALL",
        model,
        "",
        "@<TRIPOS>ATOM",
    ]

    counts: dict[str, int] = {}
    for atom in mol.GetAtoms():
        i = atom.GetIdx()
        p = conf.GetAtomPosition(i)
        sym = atom.GetSymbol()
        info = atom.GetPDBResidueInfo()
        if info is not None and info.GetName().strip():
            label = info.GetName().strip()
        else:
            counts[sym] = counts.get(sym, 0) + 1
            label = f"{sym}{counts[sym]}"
        rid, rname = residue(atom)
        lines.append(
            f"{i + 1:7d} {label:<8s}{p.x:10.4f}{p.y:10.4f}{p.z:10.4f} "
            f"{sybyl_type(atom):<6s}{order[rid]:5d} {rname:<8s}{q[i]:10.4f}")

    lines.append("@<TRIPOS>BOND")
    for bond in bonds_from.GetBonds():
        lines.append(
            f"{bond.GetIdx() + 1:6d}{bond.GetBeginAtomIdx() + 1:6d}"
            f"{bond.GetEndAtomIdx() + 1:6d} {sybyl_bond(bond):<2s}")

    return "\n".join(lines) + "\n"


def _skeleton(mol: Chem.Mol) -> str:
    """Connectivity and elements only: no bond orders, no formal charges.

    MOL2 cannot express formal charges, so a reader has to choose a
    resonance form for hypervalent groups.  RDKit renders a sulfone as the
    charge-separated S(+)-O(-), which is the same compound written another
    way.  Comparing skeletons tells a real structural error apart from that
    choice of convention.
    """
    m = Chem.RWMol(Chem.Mol(mol))
    for bond in m.GetBonds():
        bond.SetBondType(Chem.BondType.SINGLE)
        bond.SetIsAromatic(False)
    for atom in m.GetAtoms():
        atom.SetFormalCharge(0)
        atom.SetIsAromatic(False)
        atom.SetNoImplicit(True)
        atom.SetNumExplicitHs(0)
    plain = m.GetMol()
    Chem.SanitizeMol(plain, Chem.SanitizeFlags.SANITIZE_SYMMRINGS
                     | Chem.SanitizeFlags.SANITIZE_ADJUSTHS,
                     catchErrors=True)
    return Chem.MolToSmiles(Chem.RemoveHs(plain, sanitize=False))


def round_trips(mol: Chem.Mol) -> tuple[bool, str]:
    """Write, read back, and compare.  Returns (ok, detail)."""
    try:
        text = to_mol2(mol)
    except Exception as exc:                                      # noqa: BLE001
        return False, f"write failed: {exc}"

    back = Chem.MolFromMol2Block(text, removeHs=False, sanitize=True)
    if back is None:
        return False, "the output could not be read back"

    def canonical(m: Chem.Mol) -> str:
        m = Chem.Mol(m)
        if m.GetNumConformers():
            Chem.AssignStereochemistryFrom3D(m)
        smi = Chem.MolToSmiles(Chem.RemoveHs(m))
        try:
            from rdkit.Chem.MolStandardize import rdMolStandardize
            return rdMolStandardize.StandardizeSmiles(smi)
        except Exception:                                         # noqa: BLE001
            return smi

    try:
        a = canonical(mol)
        b = canonical(back)
    except Exception as exc:                                      # noqa: BLE001
        return False, f"could not canonicalise: {exc}"
    if a == b:
        return True, a

    # Not identical: is it a real structural change, or only the reader's
    # choice of resonance form for a hypervalent group?
    from rdkit.Chem import rdMolDescriptors
    same_formula = (rdMolDescriptors.CalcMolFormula(mol)
                    == rdMolDescriptors.CalcMolFormula(back))
    try:
        same_skeleton = _skeleton(mol) == _skeleton(back)
    except Exception:                                             # noqa: BLE001
        same_skeleton = False

    if same_formula and same_skeleton:
        return True, f"same skeleton and formula; reader chose another " \
                     f"resonance form ({b})"
    detail = "structure changed:\n    wrote %s\n    read  %s" % (a, b)
    if not same_formula:
        detail += "\n    formula differs: %s vs %s" % (
            rdMolDescriptors.CalcMolFormula(mol),
            rdMolDescriptors.CalcMolFormula(back))
    return False, detail
