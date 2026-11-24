"""Ligand preparation: atom typing, rotatable-bond detection, torsion tree.

The torsion tree is the data structure that lets a short chromosome
(translation, orientation, torsion angles) be expanded into full Cartesian
coordinates.  Bonds are emitted in breadth-first order from the root rigid
fragment so that a child's rotation axis is already in its final position by
the time the child rotation is applied.
"""
from __future__ import annotations

from collections import deque
from dataclasses import dataclass, field

import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem, rdMolTransforms  # noqa: F401

from . import params

# Single, acyclic, not a triple/double bond, and not an amide C-N.
_ROT_SMARTS = Chem.MolFromSmarts(
    "[!$(*#*)&!D1]-&!@[!$(*#*)&!D1]"
)
_AMIDE = Chem.MolFromSmarts("[NX3][CX3](=[OX1])")
# Terminal polar rotors worth keeping even though one side is a lone hydrogen.
_POLAR_ROTOR = Chem.MolFromSmarts("[$([OX2H1]),$([SX2H1]),$([NX3;H1,H2]),$([NX4;H1,H2,H3])]-[!#1]")


@dataclass
class Ligand:
    """A docking-ready ligand."""

    mol: Chem.Mol
    name: str

    # Geometry
    coords0: np.ndarray                      # (N, 3) reference conformer
    heavy_idx: np.ndarray                    # indices of non-hydrogen atoms

    # Torsion tree
    bond_a: np.ndarray                       # (T,) axis tail atom
    bond_b: np.ndarray                       # (T,) axis head atom
    move_ptr: np.ndarray                     # (T+1,) CSR offsets
    move_idx: np.ndarray                     # flat moving-atom indices

    # Typing
    element: list[str] = field(default_factory=list)
    ad_type: list[str] = field(default_factory=list)
    charge: np.ndarray = None                # (N,) Gasteiger
    is_donor_h: np.ndarray = None            # (N,) bool, polar hydrogen
    donor_parent: np.ndarray = None          # (N,) int, heavy atom of a polar H
    is_acceptor: np.ndarray = None           # (N,) bool
    acc_antecedent: np.ndarray = None        # (N, 3) up to 3 neighbours, -1 pad
    is_lipophilic: np.ndarray = None         # (N,) bool

    # Force-field arrays, aligned with atom order
    rii: np.ndarray = None
    epsii: np.ndarray = None
    vol: np.ndarray = None
    solpar: np.ndarray = None

    # Intramolecular pairs whose distance can change under torsion changes
    intra_pairs: np.ndarray = None           # (P, 2)

    n_rot: int = 0
    n_frozen_rot: int = 0                    # for the ChemScore rotor penalty
    sym_maps: np.ndarray = None              # heavy-atom automorphisms

    @property
    def n_atoms(self) -> int:
        return self.mol.GetNumAtoms()

    def new_buffer(self) -> np.ndarray:
        return np.empty((self.n_atoms, 3), dtype=np.float64)


# --------------------------------------------------------------------------- #
# Typing helpers
# --------------------------------------------------------------------------- #
def _ad_type(atom: Chem.Atom) -> str:
    sym = atom.GetSymbol()
    if sym == "H":
        nbrs = atom.GetNeighbors()
        if nbrs and nbrs[0].GetSymbol() in ("N", "O", "S"):
            return "HD"
        return "H"
    if sym == "C":
        return "A" if atom.GetIsAromatic() else "C"
    if sym == "N":
        return "NA" if _is_acceptor(atom) else "N"
    if sym == "O":
        return "OA"
    if sym == "S":
        return "SA" if _is_acceptor(atom) else "S"
    if sym in params.AD_PARAMS:
        return sym
    return "C"


def _is_acceptor(atom: Chem.Atom) -> bool:
    """N/O/S carrying at least one lone pair available for an H-bond."""
    sym = atom.GetSymbol()
    if sym == "O":
        # Ethers, carbonyls, hydroxyls, carboxylates all accept.  A positively
        # charged oxygen does not.
        return atom.GetFormalCharge() <= 0
    if sym == "S":
        return atom.GetFormalCharge() <= 0 and atom.GetTotalDegree() <= 2
    if sym == "N":
        if atom.GetFormalCharge() > 0:
            return False
        # An amide / aniline / pyrrole nitrogen is planar and donates rather
        # than accepts: 3 sigma bonds plus conjugation uses up the lone pair.
        n_sigma = atom.GetTotalDegree()
        if n_sigma >= 3 and (atom.GetIsAromatic() or _conjugated_n(atom)):
            return False
        return n_sigma <= 3
    return False


def _conjugated_n(atom: Chem.Atom) -> bool:
    for b in atom.GetBonds():
        if b.GetIsConjugated():
            return True
    return False


def _is_lipophilic(atom: Chem.Atom) -> bool:
    """Carbon/sulfur/halogen not bonded to a heteroatom (ChemScore definition)."""
    sym = atom.GetSymbol()
    if sym in ("F", "Cl", "Br", "I"):
        return True
    if sym not in ("C", "S"):
        return False
    for nbr in atom.GetNeighbors():
        if nbr.GetSymbol() in ("N", "O", "P") or nbr.GetFormalCharge() != 0:
            return False
    return True


# --------------------------------------------------------------------------- #
# Rotatable bonds and the tree
# --------------------------------------------------------------------------- #
def find_rotatable_bonds(mol: Chem.Mol, rotate_polar_h: bool = True) -> list[tuple[int, int]]:
    amide_n = set()
    for match in mol.GetSubstructMatches(_AMIDE):
        amide_n.add((match[0], match[1]))
        amide_n.add((match[1], match[0]))

    bonds: list[tuple[int, int]] = []
    seen: set[tuple[int, int]] = set()
    for i, j in mol.GetSubstructMatches(_ROT_SMARTS):
        key = (min(i, j), max(i, j))
        if key in seen or (i, j) in amide_n:
            continue
        b = mol.GetBondBetweenAtoms(i, j)
        if b is None or b.IsInRing() or b.GetBondType() != Chem.BondType.SINGLE:
            continue
        # Skip rotations that only spin a methyl / tert-butyl style group:
        # all distal heavy neighbours are equivalent hydrogens.
        if _spins_only_hydrogens(mol, i, j) and _spins_only_hydrogens(mol, j, i):
            continue
        seen.add(key)
        bonds.append((i, j))

    if rotate_polar_h:
        for match in mol.GetSubstructMatches(_POLAR_ROTOR):
            i, j = match[0], match[1]
            key = (min(i, j), max(i, j))
            if key in seen:
                continue
            b = mol.GetBondBetweenAtoms(i, j)
            if b is None or b.IsInRing() or b.GetBondType() != Chem.BondType.SINGLE:
                continue
            if not any(n.GetSymbol() == "H" for n in mol.GetAtomWithIdx(i).GetNeighbors()):
                continue
            seen.add(key)
            bonds.append((j, i))  # rotate the OH/NH hydrogens
    return bonds


def _spins_only_hydrogens(mol: Chem.Mol, pivot: int, other: int) -> bool:
    """True if rotating about pivot->other moves only hydrogens on `pivot`."""
    a = mol.GetAtomWithIdx(pivot)
    heavy = [n.GetIdx() for n in a.GetNeighbors()
             if n.GetIdx() != other and n.GetSymbol() != "H"]
    return len(heavy) == 0


def build_torsion_tree(mol: Chem.Mol, rot_bonds: list[tuple[int, int]]):
    """Return (bond_a, bond_b, move_ptr, move_idx, frag_id, root_frag)."""
    n = mol.GetNumAtoms()
    rot_set = {(min(i, j), max(i, j)) for i, j in rot_bonds}

    # Rigid fragments = connected components after cutting rotatable bonds.
    adj: list[list[int]] = [[] for _ in range(n)]
    for b in mol.GetBonds():
        i, j = b.GetBeginAtomIdx(), b.GetEndAtomIdx()
        if (min(i, j), max(i, j)) in rot_set:
            continue
        adj[i].append(j)
        adj[j].append(i)

    frag_id = np.full(n, -1, dtype=np.int32)
    n_frag = 0
    for start in range(n):
        if frag_id[start] >= 0:
            continue
        q = deque([start])
        frag_id[start] = n_frag
        while q:
            u = q.popleft()
            for v in adj[u]:
                if frag_id[v] < 0:
                    frag_id[v] = n_frag
                    q.append(v)
        n_frag += 1

    # Root = fragment with the most heavy atoms (keeps the core still).
    heavy_count = np.zeros(n_frag, dtype=np.int32)
    for a in mol.GetAtoms():
        if a.GetSymbol() != "H":
            heavy_count[frag_id[a.GetIdx()]] += 1
    root_frag = int(np.argmax(heavy_count))

    # Fragment tree, edges carry the atom pair of the rotatable bond.
    fadj: list[list[tuple[int, int, int]]] = [[] for _ in range(n_frag)]
    for i, j in rot_bonds:
        fi, fj = int(frag_id[i]), int(frag_id[j])
        fadj[fi].append((fj, i, j))
        fadj[fj].append((fi, j, i))

    atoms_of_frag: list[list[int]] = [[] for _ in range(n_frag)]
    for a in range(n):
        atoms_of_frag[frag_id[a]].append(a)

    bond_a: list[int] = []
    bond_b: list[int] = []
    move_sets: list[list[int]] = []

    visited = np.zeros(n_frag, dtype=bool)
    visited[root_frag] = True
    queue = deque([root_frag])
    parent = {root_frag: -1}
    order: list[int] = []
    tree_children: dict[int, list[tuple[int, int, int]]] = {}
    while queue:
        f = queue.popleft()
        order.append(f)
        kids = []
        for g, pa, pb in fadj[f]:
            if visited[g]:
                continue
            visited[g] = True
            parent[g] = f
            kids.append((g, pa, pb))
            queue.append(g)
        tree_children[f] = kids

    # Subtree membership: process fragments in reverse BFS order.
    subtree: dict[int, list[int]] = {}
    for f in reversed(order):
        atoms = list(atoms_of_frag[f])
        for g, _, _ in tree_children.get(f, []):
            atoms.extend(subtree[g])
        subtree[f] = atoms

    for f in order:
        for g, pa, pb in tree_children.get(f, []):
            bond_a.append(pa)
            bond_b.append(pb)
            move_sets.append(subtree[g])

    move_ptr = np.zeros(len(move_sets) + 1, dtype=np.int64)
    for k, s in enumerate(move_sets):
        move_ptr[k + 1] = move_ptr[k] + len(s)
    move_idx = np.empty(int(move_ptr[-1]), dtype=np.int64)
    for k, s in enumerate(move_sets):
        move_idx[move_ptr[k]:move_ptr[k + 1]] = np.asarray(s, dtype=np.int64)

    return (np.asarray(bond_a, dtype=np.int64),
            np.asarray(bond_b, dtype=np.int64),
            move_ptr, move_idx, frag_id, root_frag)


def _intra_pairs(mol: Chem.Mol, frag_id: np.ndarray) -> np.ndarray:
    """Atom pairs in different rigid fragments and >3 bonds apart.

    Only these contribute a variable internal energy; everything else is
    constant for a fixed reference conformer.
    """
    n = mol.GetNumAtoms()
    dm = Chem.GetDistanceMatrix(mol)
    pairs = []
    for i in range(n):
        for j in range(i + 1, n):
            if frag_id[i] == frag_id[j]:
                continue
            if dm[i, j] < 4:
                continue
            pairs.append((i, j))
    if not pairs:
        return np.zeros((0, 2), dtype=np.int64)
    return np.asarray(pairs, dtype=np.int64)


def _symmetry_maps(mol: Chem.Mol, heavy_idx: np.ndarray, limit: int = 64) -> np.ndarray:
    """Heavy-atom graph automorphisms, for symmetry-corrected RMSD."""
    heavy = Chem.RWMol(mol)
    try:
        heavy_mol = Chem.RemoveHs(Chem.Mol(mol))
    except Exception:
        heavy_mol = Chem.Mol(mol)
    try:
        matches = heavy_mol.GetSubstructMatches(
            heavy_mol, uniquify=False, useChirality=False, maxMatches=limit)
    except Exception:
        matches = ()
    if not matches:
        return np.arange(len(heavy_idx), dtype=np.int64).reshape(1, -1)
    return np.asarray(matches, dtype=np.int64)


# --------------------------------------------------------------------------- #
# Entry point
# --------------------------------------------------------------------------- #
def prepare_ligand(mol: Chem.Mol, name: str = "ligand",
                   rotate_polar_h: bool = True,
                   embed_if_flat: bool = True,
                   seed: int = 0xC0FFEE) -> Ligand:
    """Type a molecule and build its torsion tree."""
    mol = Chem.Mol(mol)
    if mol.GetNumConformers() == 0:
        mol = Chem.AddHs(mol)
        p = AllChem.ETKDGv3()
        p.randomSeed = seed
        if AllChem.EmbedMolecule(mol, p) != 0:
            raise ValueError(f"could not generate a 3D conformer for {name}")
        AllChem.MMFFOptimizeMolecule(mol, maxIters=500)
    else:
        mol = Chem.AddHs(mol, addCoords=True)
        if embed_if_flat:
            c = mol.GetConformer().GetPositions()
            if np.ptp(c, axis=0).min() < 1e-3:
                p = AllChem.ETKDGv3()
                p.randomSeed = seed
                AllChem.EmbedMolecule(mol, p)
                AllChem.MMFFOptimizeMolecule(mol, maxIters=500)

    Chem.SanitizeMol(mol, Chem.SanitizeFlags.SANITIZE_ALL
                     ^ Chem.SanitizeFlags.SANITIZE_PROPERTIES,
                     catchErrors=True)
    AllChem.ComputeGasteigerCharges(mol)

    n = mol.GetNumAtoms()
    coords0 = np.asarray(mol.GetConformer().GetPositions(), dtype=np.float64)

    rot_bonds = find_rotatable_bonds(mol, rotate_polar_h=rotate_polar_h)
    bond_a, bond_b, move_ptr, move_idx, frag_id, _root = build_torsion_tree(mol, rot_bonds)

    element: list[str] = []
    adt: list[str] = []
    charge = np.zeros(n)
    is_donor_h = np.zeros(n, dtype=bool)
    donor_parent = np.full(n, -1, dtype=np.int64)
    is_acceptor = np.zeros(n, dtype=bool)
    acc_ante = np.full((n, 3), -1, dtype=np.int64)
    is_lipo = np.zeros(n, dtype=bool)
    rii = np.zeros(n)
    epsii = np.zeros(n)
    vol = np.zeros(n)
    solpar = np.zeros(n)

    for a in mol.GetAtoms():
        i = a.GetIdx()
        sym = a.GetSymbol()
        element.append(sym)
        t = _ad_type(a)
        adt.append(t)
        q = a.GetDoubleProp("_GasteigerCharge")
        if not np.isfinite(q):
            q = 0.0
        hq = a.GetDoubleProp("_GasteigerHCharge") if a.HasProp("_GasteigerHCharge") else 0.0
        charge[i] = q + (hq if np.isfinite(hq) else 0.0)

        if t == "HD":
            is_donor_h[i] = True
            donor_parent[i] = a.GetNeighbors()[0].GetIdx()
        if _is_acceptor(a):
            is_acceptor[i] = True
            nb = [x.GetIdx() for x in a.GetNeighbors() if x.GetSymbol() != "H"][:3]
            for k, v in enumerate(nb):
                acc_ante[i, k] = v
        is_lipo[i] = _is_lipophilic(a)

        p = params.AD_PARAMS.get(t, params.AD_PARAMS["C"])
        rii[i], epsii[i], vol[i], solpar[i] = p

    heavy_idx = np.asarray([a.GetIdx() for a in mol.GetAtoms()
                            if a.GetSymbol() != "H"], dtype=np.int64)

    lig = Ligand(
        mol=mol, name=name, coords0=coords0, heavy_idx=heavy_idx,
        bond_a=bond_a, bond_b=bond_b, move_ptr=move_ptr, move_idx=move_idx,
        element=element, ad_type=adt, charge=charge,
        is_donor_h=is_donor_h, donor_parent=donor_parent,
        is_acceptor=is_acceptor, acc_antecedent=acc_ante, is_lipophilic=is_lipo,
        rii=rii, epsii=epsii, vol=vol, solpar=solpar,
        intra_pairs=_intra_pairs(mol, frag_id),
        n_rot=len(bond_a),
        n_frozen_rot=len([1 for i, j in rot_bonds
                          if not _spins_only_hydrogens(mol, i, j)]),
        sym_maps=_symmetry_maps(mol, heavy_idx),
    )
    return lig


def read_molecules(path: str) -> list[tuple[Chem.Mol, str]]:
    """Read every molecule in a ligand file.

    SDF, SMILES and MOL2 files routinely hold whole libraries, so all records
    are returned, not just the first.  Each comes back with a name: the
    record's own title where it has one, otherwise the file name with an index.
    """
    low = path.lower()
    stem = path.rsplit("/", 1)[-1].rsplit(".", 1)[0] or "ligand"
    out: list[tuple[Chem.Mol, str]] = []

    if low.endswith((".sdf", ".mol")):
        supp = Chem.SDMolSupplier(path, removeHs=False, sanitize=True)
        for i, m in enumerate(supp):
            if m is None:
                continue                       # skip the bad record, keep going
            out.append((m, _record_name(m, stem, i, len(supp))))
    elif low.endswith(".mol2"):
        # RDKit has no multi-MOL2 supplier, so split on the record marker.
        text = open(path, errors="replace").read()
        chunks = [c for c in text.split("@<TRIPOS>MOLECULE") if c.strip()]
        for i, chunk in enumerate(chunks):
            m = Chem.MolFromMol2Block("@<TRIPOS>MOLECULE" + chunk,
                                      removeHs=False, sanitize=True)
            if m is not None:
                out.append((m, _record_name(m, stem, i, len(chunks))))
    elif low.endswith(".pdb"):
        m = Chem.MolFromPDBFile(path, removeHs=False, sanitize=True)
        if m is not None:
            out.append((m, _record_name(m, stem, 0, 1)))
    elif low.endswith((".smi", ".smiles")):
        with open(path, errors="replace") as fh:
            for i, line in enumerate(fh):
                line = line.strip()
                if not line or line.startswith("#"):
                    continue
                parts = line.split(None, 1)
                m = Chem.MolFromSmiles(parts[0])
                if m is None:
                    continue
                name = parts[1].strip() if len(parts) > 1 else ""
                out.append((m, name or f"{stem}_{i + 1}"))
    else:
        raise ValueError(f"unrecognised ligand format: {path}")

    if not out:
        raise ValueError(f"no readable molecules in {path}")
    return out


def _record_name(mol: Chem.Mol, stem: str, index: int, total: int) -> str:
    if mol.HasProp("_Name"):
        name = mol.GetProp("_Name").strip()
        if name:
            return name
    return stem if total == 1 else f"{stem}_{index + 1}"


def load_ligands(path: str, failures: list | None = None, **kw) -> list[Ligand]:
    """Prepare every molecule in a ligand file for docking.

    A record that cannot be prepared (no 3D embedding, broken valence) is
    skipped rather than killing the whole library.  Pass a list as `failures`
    to collect `(name, reason)` for the ones that were dropped.
    """
    prepared: list[Ligand] = []
    dropped: list[tuple[str, str]] = []
    for mol, name in read_molecules(path):
        try:
            prepared.append(prepare_ligand(mol, name=name, **kw))
        except Exception as exc:                                  # noqa: BLE001
            dropped.append((name, f"{type(exc).__name__}: {exc}"))
    if failures is not None:
        failures.extend(dropped)
    if not prepared:
        detail = "; ".join(f"{n}: {e}" for n, e in dropped[:3])
        raise ValueError(f"no ligand in {path} could be prepared ({detail})")
    return prepared


def load_ligand(path: str, **kw) -> Ligand:
    """Read a single ligand -- the first readable record in the file."""
    return load_ligands(path, **kw)[0]
