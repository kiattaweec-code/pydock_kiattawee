"""Receptor preparation from PDB: typing, charges, polar-hydrogen placement.

Crystal structures almost never ship hydrogens, and hydrogen-bond geometry is
the single most informative term in an empirical score, so this module places
polar hydrogens itself from idealised geometry.  Hydrogens whose position is
not determined by the heavy-atom frame -- hydroxyls, thiols, lysine ammonium --
are flagged as *rotatable* and handed to the GA as extra genes, which is how
GOLD treats them.
"""
from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np

from . import params

DEG = np.pi / 180.0

# --------------------------------------------------------------------------- #
# Polar-hydrogen placement rules
# --------------------------------------------------------------------------- #
# Hydrogens fixed by the heavy-atom frame: direction = -(sum of unit vectors
# towards the listed neighbours).
_BISECTOR_H = {
    ("*", "N"): ("CA", "C_prev"),
    ("TRP", "NE1"): ("CD1", "CE2"),
    ("HIS", "ND1"): ("CG", "CE1"),
}

# Planar NH2 groups: two hydrogens in the plane of (donor, parent, reference).
_SP2_NH2 = {
    ("ASN", "ND2"): ("CG", "OD1"),
    ("GLN", "NE2"): ("CD", "OE1"),
    ("ARG", "NH1"): ("CZ", "NE"),
    ("ARG", "NH2"): ("CZ", "NE"),
}
_SP2_NH = {
    ("ARG", "NE"): ("CD", "CZ"),
}

# Rotatable hydrogens: (parent, n_hydrogens, cone half-angle in degrees)
_ROTATABLE_H = {
    ("SER", "OG"): ("CB", 1, 109.5),
    ("THR", "OG1"): ("CB", 1, 109.5),
    ("TYR", "OH"): ("CZ", 1, 109.5),
    ("CYS", "SG"): ("CB", 1, 96.0),
    ("LYS", "NZ"): ("CE", 3, 109.5),
}

# Acceptor atoms by name.  Every oxygen accepts; nitrogen and sulfur only when
# listed here.
_ACCEPTOR_N = {("HIS", "NE2")}
_ACCEPTOR_S = {("MET", "SD"), ("CYS", "SG")}

# Carbons bonded to N/O, which disqualifies them as lipophilic (ChemScore).
_NON_LIPO_C = {
    "C", "CA",
    "ASP:CG", "GLU:CD", "ASN:CG", "GLN:CD",
    "ARG:CZ", "ARG:CD", "LYS:CE", "PRO:CD",
    "HIS:CG", "HIS:CD2", "HIS:CE1",
    "SER:CB", "THR:CB", "TYR:CZ",
    "TRP:CD1", "TRP:CE2",
}

STANDARD_AA = {
    "ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS", "ILE",
    "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL",
}

# Force-field and tautomer-specific spellings all collapse onto the standard
# three-letter names, so the typing rules below need only one entry each.
RESNAME_ALIASES = {
    "HID": "HIS", "HIE": "HIS", "HIP": "HIS",
    "HSD": "HIS", "HSE": "HIS", "HSP": "HIS",
    "CYX": "CYS", "CYM": "CYS", "CYSS": "CYS",
    "MSE": "MET", "LYN": "LYS", "ASH": "ASP", "GLH": "GLU", "TYM": "TYR",
    "ARN": "ARG", "HOH": "HOH", "CSO": "CYS",
}


def canonical_resname(rn: str) -> str:
    return RESNAME_ALIASES.get(rn, rn)


@dataclass
class Receptor:
    coords: np.ndarray                   # (M, 3) heavy atoms + placed polar H
    element: list[str]
    name: list[str]
    resname: list[str]
    resid: list[tuple[str, int]]
    ad_type: list[str] = field(default_factory=list)
    charge: np.ndarray = None
    is_donor_h: np.ndarray = None
    donor_parent: np.ndarray = None
    is_acceptor: np.ndarray = None
    acc_antecedent: np.ndarray = None
    is_lipophilic: np.ndarray = None
    is_metal: np.ndarray = None
    rii: np.ndarray = None
    epsii: np.ndarray = None
    vol: np.ndarray = None
    solpar: np.ndarray = None

    # Rotatable receptor hydrogens: one gene each.
    rot_h_axis_a: np.ndarray = None      # parent heavy atom index
    rot_h_axis_b: np.ndarray = None      # pivot (donor) heavy atom index
    rot_h_atoms: np.ndarray = None       # (R, 3) hydrogen indices, -1 padded
    rot_h_ref: np.ndarray = None         # reference coords of those hydrogens

    # What preparation did to the file, so callers can report it.
    supplied_hydrogens: int = 0
    replaced_hydrogens: int = 0

    @property
    def n_atoms(self) -> int:
        return self.coords.shape[0]


# --------------------------------------------------------------------------- #
# PDB parsing
# --------------------------------------------------------------------------- #
def parse_pdb(path: str, keep_waters: bool = False,
              keep_hetatm: bool = True,
              exclude_resnames: set[str] | None = None):
    """Minimal PDB reader.  Returns a list of atom dicts."""
    exclude = {r.upper() for r in (exclude_resnames or set())}
    atoms = []
    seen = set()
    with open(path) as fh:
        for line in fh:
            rec = line[:6]
            if rec not in ("ATOM  ", "HETATM"):
                if line.startswith("ENDMDL"):
                    break
                continue
            altloc = line[16]
            if altloc not in (" ", "A"):
                continue
            resname = line[17:20].strip().upper()
            if resname in exclude:
                continue
            resname = canonical_resname(resname)
            if resname in ("HOH", "WAT", "DOD") and not keep_waters:
                continue
            if rec == "HETATM" and not keep_hetatm and resname not in STANDARD_AA:
                continue
            name = line[12:16].strip()
            chain = line[21].strip() or "A"
            try:
                resseq = int(line[22:26])
            except ValueError:
                continue
            icode = line[26].strip()
            key = (chain, resseq, icode, name)
            if key in seen:
                continue
            seen.add(key)
            el = line[76:78].strip()
            if not el:
                el = "".join(c for c in name if c.isalpha())[:1]
            el = el.capitalize() if len(el) == 2 else el.upper()
            atoms.append(dict(
                name=name, resname=resname, chain=chain, resseq=resseq,
                icode=icode, element=el,
                xyz=np.array([float(line[30:38]), float(line[38:46]),
                              float(line[46:54])]),
                occ=float(line[54:60] or 1.0),
            ))
    if not atoms:
        raise ValueError(f"no atoms parsed from {path}")
    return atoms


def _unit(v):
    n = np.linalg.norm(v)
    return v / n if n > 1e-8 else np.array([0.0, 0.0, 1.0])


def _cone_position(pivot, parent_dir, angle_deg, azimuth, length):
    """A point at `angle_deg` from -parent_dir, rotated by `azimuth` about it."""
    axis = _unit(parent_dir)
    # Build an orthonormal frame around the axis.
    tmp = np.array([1.0, 0.0, 0.0])
    if abs(np.dot(tmp, axis)) > 0.9:
        tmp = np.array([0.0, 1.0, 0.0])
    e1 = _unit(np.cross(axis, tmp))
    e2 = np.cross(axis, e1)
    th = angle_deg * DEG
    d = (np.cos(np.pi - th) * axis
         + np.sin(np.pi - th) * (np.cos(azimuth) * e1 + np.sin(azimuth) * e2))
    return pivot + length * _unit(d)


_BOND_LEN = {"N": 1.01, "O": 0.97, "S": 1.34}


# --------------------------------------------------------------------------- #
# Main preparation
# --------------------------------------------------------------------------- #
def prepare_receptor(path: str, keep_waters: bool = False,
                     rotate_oh: bool = True,
                     exclude_resnames: set[str] | None = None,
                     strip_hydrogens: bool = True) -> Receptor:
    atoms = parse_pdb(path, keep_waters=keep_waters,
                      exclude_resnames=exclude_resnames)

    # Hydrogens that arrive in the file are replaced by default, because the
    # placement below is also what builds the rotatable OH/SH/NH3+ genes: a
    # protonated receptor used as-is scored with zero of them, losing 82
    # degrees of freedom on a protein this size without saying anything.
    # Pass strip_hydrogens=False to keep a protonation you chose deliberately,
    # and accept that those torsions stay fixed.
    supplied_h = sum(1 for a in atoms if a["element"] == "H")
    replaced_h = 0
    if supplied_h and strip_hydrogens:
        atoms = [a for a in atoms if a["element"] != "H"]
        replaced_h = supplied_h

    # Does the file already carry hydrogens?
    has_h = any(a["element"] == "H" for a in atoms)

    # Index atoms by residue for rule lookups.
    residues: dict[tuple, dict[str, int]] = {}
    for i, a in enumerate(atoms):
        key = (a["chain"], a["resseq"], a["icode"])
        residues.setdefault(key, {})[a["name"]] = i
    res_order = list(residues.keys())
    prev_c: dict[tuple, int] = {}
    for k in range(1, len(res_order)):
        a_key, b_key = res_order[k - 1], res_order[k]
        if a_key[0] != b_key[0]:
            continue
        if "C" in residues[a_key]:
            prev_c[b_key] = residues[a_key]["C"]

    added: list[dict] = []
    rot_groups: list[dict] = []

    if not has_h:
        for key, amap in residues.items():
            rn = atoms[next(iter(amap.values()))]["resname"]

            # --- fixed, frame-determined hydrogens -------------------------
            for (rule_res, dname), nbrs in _BISECTOR_H.items():
                if rule_res not in ("*", rn) or dname not in amap:
                    continue
                if rule_res == "*" and rn == "PRO":
                    continue
                d = amap[dname]
                pts = []
                for nb in nbrs:
                    if nb == "C_prev":
                        j = prev_c.get(key)
                    else:
                        j = amap.get(nb)
                    if j is None:
                        pts = []
                        break
                    pts.append(atoms[j]["xyz"])
                if len(pts) != 2:
                    continue
                dpos = atoms[d]["xyz"]
                direction = -_unit(_unit(pts[0] - dpos) + _unit(pts[1] - dpos))
                added.append(dict(
                    name="H" + dname, resname=rn, chain=key[0], resseq=key[1],
                    icode=key[2], element="H",
                    xyz=dpos + _BOND_LEN.get(atoms[d]["element"], 1.0) * direction,
                    parent=d, occ=1.0))

            for table, n_h in ((_SP2_NH2, 2), (_SP2_NH, 1)):
                for (rule_res, dname), (pname, rname) in table.items():
                    if rule_res != rn or dname not in amap:
                        continue
                    if pname not in amap or rname not in amap:
                        continue
                    d = amap[dname]
                    dpos = atoms[d]["xyz"]
                    ppos = atoms[amap[pname]]["xyz"]
                    rpos = atoms[amap[rname]]["xyz"]
                    normal = _unit(np.cross(ppos - dpos, rpos - ppos))
                    axis = _unit(dpos - ppos)
                    bl = _BOND_LEN.get(atoms[d]["element"], 1.0)
                    if n_h == 2:
                        for s in (+1, -1):
                            dirv = _unit(np.cos(60 * DEG) * axis
                                         + s * np.sin(60 * DEG)
                                         * _unit(np.cross(normal, axis)))
                            added.append(dict(
                                name=f"H{dname}{1 if s > 0 else 2}", resname=rn,
                                chain=key[0], resseq=key[1], icode=key[2],
                                element="H", xyz=dpos + bl * dirv,
                                parent=d, occ=1.0))
                    else:
                        # Single sp2 NH: anti to the reference substituent.
                        dirv = -_unit(_unit(ppos - dpos) + _unit(rpos - dpos))
                        added.append(dict(
                            name="H" + dname, resname=rn, chain=key[0],
                            resseq=key[1], icode=key[2], element="H",
                            xyz=dpos + bl * dirv, parent=d, occ=1.0))

            # --- rotatable hydrogens ---------------------------------------
            for (rule_res, dname), (pname, n_h, cone) in _ROTATABLE_H.items():
                if rule_res != rn or dname not in amap or pname not in amap:
                    continue
                d = amap[dname]
                p = amap[pname]
                dpos = atoms[d]["xyz"]
                bl = _BOND_LEN.get(atoms[d]["element"], 1.0)
                group = []
                for k in range(n_h):
                    az = 2 * np.pi * k / n_h
                    pos = _cone_position(dpos, dpos - atoms[p]["xyz"], cone, az, bl)
                    added.append(dict(
                        name=f"H{dname}{k + 1}" if n_h > 1 else "H" + dname,
                        resname=rn, chain=key[0], resseq=key[1], icode=key[2],
                        element="H", xyz=pos, parent=d, occ=1.0))
                    group.append(len(atoms) + len(added) - 1)
                if rotate_oh:
                    rot_groups.append(dict(parent=p, pivot=d, hydrogens=group))

    all_atoms = atoms + added
    m = len(all_atoms)
    coords = np.asarray([a["xyz"] for a in all_atoms], dtype=np.float64)

    element = [a["element"] for a in all_atoms]
    aname = [a["name"] for a in all_atoms]
    resname = [a["resname"] for a in all_atoms]
    resid = [(a["chain"], a["resseq"]) for a in all_atoms]

    adt: list[str] = []
    charge = np.zeros(m)
    is_donor_h = np.zeros(m, dtype=bool)
    donor_parent = np.full(m, -1, dtype=np.int64)
    is_acceptor = np.zeros(m, dtype=bool)
    acc_ante = np.full((m, 3), -1, dtype=np.int64)
    is_lipo = np.zeros(m, dtype=bool)
    is_metal = np.zeros(m, dtype=bool)
    rii = np.zeros(m)
    epsii = np.zeros(m)
    vol = np.zeros(m)
    solpar = np.zeros(m)

    # Neighbour lists for acceptor directionality, by distance (no CONECT).
    heavy = np.asarray([i for i in range(m) if element[i] != "H"], dtype=np.int64)
    heavy_xyz = coords[heavy]

    for i, a in enumerate(all_atoms):
        el = element[i]
        rn = resname[i]
        nm = aname[i]
        key2 = f"{rn}:{nm}"

        if el.upper() in params.METALS and rn not in STANDARD_AA:
            is_metal[i] = True
            t = el if el in params.AD_PARAMS else "Zn"
            adt.append(t)
            charge[i] = 2.0
        elif el == "H":
            parent = a.get("parent", -1)
            if parent < 0:
                # Hydrogens from the input file: polar if next to N/O/S.
                d2 = np.sum((heavy_xyz - coords[i]) ** 2, axis=1)
                j = int(heavy[np.argmin(d2)])
                parent = j if d2.min() < 1.6 ** 2 else -1
            polar = parent >= 0 and element[parent] in ("N", "O", "S")
            adt.append("HD" if polar else "H")
            if polar:
                is_donor_h[i] = True
                donor_parent[i] = parent
            charge[i] = _charge_of(rn, nm, el)
        else:
            acc = False
            if el == "O":
                acc = True
            elif el == "N" and ((rn, nm) in _ACCEPTOR_N or rn not in STANDARD_AA):
                acc = True
            elif el == "S" and (rn, nm) in _ACCEPTOR_S:
                acc = True
            is_acceptor[i] = acc
            if acc:
                d2 = np.sum((heavy_xyz - coords[i]) ** 2, axis=1)
                order = np.argsort(d2)
                k = 0
                for o in order:
                    j = int(heavy[o])
                    if j == i or d2[o] > 1.95 ** 2:
                        continue
                    acc_ante[i, k] = j
                    k += 1
                    if k == 3:
                        break
            if el in ("C", "S"):
                is_lipo[i] = (nm not in _NON_LIPO_C) and (key2 not in _NON_LIPO_C)
            elif el in ("F", "Cl", "Br", "I"):
                is_lipo[i] = True

            if el == "C":
                t = "A" if _aromatic_carbon(rn, nm) else "C"
            elif el == "N":
                t = "NA" if acc else "N"
            elif el == "O":
                t = "OA"
            elif el == "S":
                t = "SA" if acc else "S"
            else:
                t = el if el in params.AD_PARAMS else "C"
            adt.append(t)
            charge[i] = _charge_of(rn, nm, el)

        p = params.AD_PARAMS.get(adt[i], params.AD_PARAMS["C"])
        rii[i], epsii[i], vol[i], solpar[i] = p

    rec = Receptor(
        coords=coords, element=element, name=aname, resname=resname, resid=resid,
        ad_type=adt, charge=charge, is_donor_h=is_donor_h,
        donor_parent=donor_parent, is_acceptor=is_acceptor,
        acc_antecedent=acc_ante, is_lipophilic=is_lipo, is_metal=is_metal,
        rii=rii, epsii=epsii, vol=vol, solpar=solpar,
    )

    if rot_groups:
        rec.rot_h_axis_a = np.asarray([g["parent"] for g in rot_groups], dtype=np.int64)
        rec.rot_h_axis_b = np.asarray([g["pivot"] for g in rot_groups], dtype=np.int64)
        hidx = np.full((len(rot_groups), 3), -1, dtype=np.int64)
        for k, g in enumerate(rot_groups):
            for j, h in enumerate(g["hydrogens"][:3]):
                hidx[k, j] = h
        rec.rot_h_atoms = hidx
        rec.rot_h_ref = coords.copy()
    else:
        rec.rot_h_axis_a = np.zeros(0, dtype=np.int64)
        rec.rot_h_axis_b = np.zeros(0, dtype=np.int64)
        rec.rot_h_atoms = np.zeros((0, 3), dtype=np.int64)
        rec.rot_h_ref = coords.copy()

    rec.supplied_hydrogens = supplied_h
    rec.replaced_hydrogens = replaced_h
    return rec


_AROMATIC = {
    "PHE": {"CG", "CD1", "CD2", "CE1", "CE2", "CZ"},
    "TYR": {"CG", "CD1", "CD2", "CE1", "CE2", "CZ"},
    "TRP": {"CG", "CD1", "CD2", "CE2", "CE3", "CZ2", "CZ3", "CH2"},
    "HIS": {"CG", "CD2", "CE1"},
}


def _aromatic_carbon(resname: str, name: str) -> bool:
    return name in _AROMATIC.get(resname, ())


def _charge_of(resname: str, name: str, element: str) -> float:
    sc = params.SIDECHAIN_CHARGES.get(resname)
    if sc and name in sc:
        return sc[name]
    if name in params.BACKBONE_CHARGES and resname in STANDARD_AA:
        return params.BACKBONE_CHARGES[name]
    if element == "H":
        return 0.15 if name.startswith(("HN", "HO", "HG", "HH", "HE", "HD", "HZ")) else 0.06
    if element == "O":
        return -0.45
    if element == "N":
        return -0.40
    if element == "S":
        return -0.20
    return 0.0


# --------------------------------------------------------------------------- #
# Binding site
# --------------------------------------------------------------------------- #
def site_from_ligand(coords: np.ndarray, padding: float = 4.0):
    """Centre and half-widths of a box around a reference ligand."""
    centre = coords.mean(axis=0)
    half = (coords.max(axis=0) - coords.min(axis=0)) / 2.0 + padding
    return centre, np.maximum(half, 5.0)


def site_from_residues(rec: Receptor, residues: list[tuple[str, int]],
                       padding: float = 5.0):
    mask = np.asarray([r in set(residues) for r in rec.resid])
    if not mask.any():
        raise ValueError("none of the requested residues were found")
    sel = rec.coords[mask]
    return sel.mean(axis=0), (sel.max(axis=0) - sel.min(axis=0)) / 2.0 + padding
