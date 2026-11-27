"""Structure preparation: turn a downloaded complex into docking inputs.

A PDB entry is not a receptor.  It holds the bound ligand, waters, buffer
molecules and often several chains, and the docking engine keeps every
hetero atom it is given -- so a ligand left in the file fills the pocket it
was meant to reveal, and crystallographic hydrogens switch off the hydroxyl
rotation the search relies on.

These routes let someone do the splitting in a browser: upload, look at the
structure, click the ligand, and get back a clean receptor plus the ligand as
MOL2 -- then send both straight into docking without downloading anything.
"""
from __future__ import annotations

import json
import time
import uuid
from pathlib import Path

from fastapi import APIRouter, File, Form, HTTPException, UploadFile
from fastapi.responses import FileResponse, JSONResponse
from rdkit import Chem

from pydock import mol2 as mol2mod
from pydock import structure as structmod

router = APIRouter(prefix="/api/prep", tags=["prepare"])

MAX_UPLOAD = 40 * 1024 * 1024
KEEP_FOR = 24 * 3600          # a prepared structure is scratch, not storage
ROOT: Path | None = None      # set by app.py, which owns the data directory


def configure(data_root: Path) -> None:
    global ROOT
    ROOT = data_root / "prep"
    ROOT.mkdir(parents=True, exist_ok=True)


def _dir(prep_id: str) -> Path:
    if ROOT is None:
        raise HTTPException(500, "preparation storage is not configured")
    # uuid4().hex only, so a crafted id cannot walk out of the directory.
    if not prep_id or len(prep_id) != 32 or not all(
            c in "0123456789abcdef" for c in prep_id):
        raise HTTPException(400, "bad preparation id")
    path = ROOT / prep_id
    if not path.is_dir():
        raise HTTPException(404, "that prepared structure has expired; "
                                 "upload the file again")
    return path


def _sweep() -> None:
    """Drop preparations older than a day, so the directory cannot grow."""
    if ROOT is None or not ROOT.is_dir():
        return
    cutoff = time.time() - KEEP_FOR
    for path in ROOT.iterdir():
        try:
            if path.is_dir() and path.stat().st_mtime < cutoff:
                for child in path.iterdir():
                    child.unlink(missing_ok=True)
                path.rmdir()
        except OSError:
            pass


def _flag(value) -> bool:
    return str(value).strip().lower() in ("1", "true", "yes", "on")


def _viewer_copy(st: structmod.Structure, path: Path) -> int:
    """A light copy for the 3D viewer: no hydrogens, nothing else removed.

    The viewer has to show the structure exactly as uploaded -- the ligand
    still in place, since picking it is the whole point -- but hydrogens
    triple the atom count for nothing on screen.
    """
    lines = [r["line"] for r in st.records if r["element"] != "H"]
    path.write_text("\n".join(lines) + "\nEND\n")
    return len(lines)


# --------------------------------------------------------------------------- #
# Inspect
# --------------------------------------------------------------------------- #
@router.post("/inspect")
def inspect(structure: UploadFile = File(...)):
    """Upload a complex and get back what is in it."""
    _sweep()
    name = Path(structure.filename or "structure.pdb").name
    if Path(name).suffix.lower() not in (".pdb", ".ent", ".txt"):
        raise HTTPException(400, f"{name}: expected a .pdb file")
    data = structure.file.read(MAX_UPLOAD + 1)
    if not data:
        raise HTTPException(400, f"{name} is empty")
    if len(data) > MAX_UPLOAD:
        raise HTTPException(413, f"{name} is larger than "
                                 f"{MAX_UPLOAD // (1024 * 1024)} MB")

    prep_id = uuid.uuid4().hex
    wd = (ROOT or Path(".")) / prep_id
    wd.mkdir(parents=True, exist_ok=True)
    source = wd / "source.pdb"
    source.write_bytes(data)

    try:
        st = structmod.inspect_structure(str(source))
    except Exception as exc:                                      # noqa: BLE001
        raise HTTPException(400, f"could not read {name}: {exc}")

    _viewer_copy(st, wd / "viewer.pdb")

    def pack(c: structmod.Component) -> dict:
        return {"key": c.key, "label": c.label, "resname": c.resname,
                "chain": c.chain, "resseq": c.resseq, "kind": c.kind,
                "n_atoms": c.n_atoms, "n_heavy": c.n_heavy,
                "formula": c.formula, "centre": list(c.centre)}

    ligands = [pack(c) for c in st.ligand_candidates()]
    others = [pack(c) for c in st.components
              if c.kind in ("metal", "ion", "buffer")]
    waters = sum(1 for c in st.components if c.kind == "water")
    polymer = sum(1 for c in st.components if c.kind == "polymer")

    meta = {
        "prep_id": prep_id, "name": name, "created": time.time(),
        "n_atoms": st.n_atoms, "n_hydrogens": st.n_hydrogens,
        "has_hydrogens": st.has_hydrogens,
        "polymer_residues": polymer, "waters": waters,
        "chains": {k: v for k, v in sorted(st.chains.items())},
        "ligands": ligands, "others": others,
    }
    (wd / "meta.json").write_text(json.dumps(meta))

    notes = []
    if not ligands:
        notes.append("No bound ligand was found. If this protein came without "
                     "one, it is ready to use as a receptor as it is.")
    if st.has_hydrogens:
        notes.append(f"The file carries {st.n_hydrogens} hydrogens. Removing "
                     f"them is recommended: pydock places polar hydrogens "
                     f"itself and only treats hydroxyls as rotatable when the "
                     f"receptor arrives without any.")
    if waters:
        notes.append(f"{waters} waters will be removed unless you keep them.")

    return JSONResponse({**meta, "notes": notes})


# --------------------------------------------------------------------------- #
# Build
# --------------------------------------------------------------------------- #
@router.post("/build")
def build(prep_id: str = Form(...),
          ligand_key: str = Form(""),
          drop_keys: str = Form(""),
          smiles: str = Form(""),
          remove_hydrogens: str = Form("true"),
          remove_waters: str = Form("true"),
          remove_buffers: str = Form("true"),
          keep_metals: str = Form("true"),
          remove_ions: str = Form("false"),
          chains: str = Form("")):
    """Write the cleaned receptor, and the chosen ligand as MOL2 and SDF."""
    wd = _dir(prep_id)
    try:
        st = structmod.inspect_structure(str(wd / "source.pdb"))
    except Exception as exc:                                      # noqa: BLE001
        raise HTTPException(400, f"could not re-read the upload: {exc}")

    extra = {k.strip() for k in drop_keys.split(",") if k.strip()}
    keep_chains = {c.strip() for c in chains.split(",") if c.strip()}
    unknown = [k for k in extra | ({ligand_key} if ligand_key else set())
               if st.by_key(k) is None]
    if unknown:
        raise HTTPException(400, f"no such component: {', '.join(unknown)}")

    warnings: list[str] = []
    ligand_info = None

    if ligand_key:
        try:
            lig = structmod.extract_residue(st, ligand_key,
                                            smiles=smiles.strip() or None)
        except ValueError as exc:
            raise HTTPException(400, str(exc))
        if lig.HasProp("pydock_warning"):
            warnings.append(lig.GetProp("pydock_warning"))
        try:
            (wd / "ligand.mol2").write_text(mol2mod.to_mol2(lig))
        except Exception as exc:                                  # noqa: BLE001
            raise HTTPException(500, f"could not write the MOL2: {exc}")
        writer = Chem.SDWriter(str(wd / "ligand.sdf"))
        writer.write(lig)
        writer.close()
        orders = lig.GetProp("pydock_bond_orders")
        ligand_info = {
            "key": ligand_key,
            "label": lig.GetProp("_Name"),
            "n_atoms": lig.GetNumAtoms(),
            "formula": lig.GetProp("pydock_formula")
            if lig.HasProp("pydock_formula") else "",
            "bond_orders": orders,
            "trusted": orders == "template",
            "mol2": f"/api/prep/{prep_id}/ligand.mol2",
            "sdf": f"/api/prep/{prep_id}/ligand.sdf",
        }
        extra.add(ligand_key)

    try:
        text, summary = structmod.build_receptor(
            st, drop_keys=extra,
            remove_waters=_flag(remove_waters),
            remove_hydrogens=_flag(remove_hydrogens),
            keep_metals=_flag(keep_metals),
            remove_buffers=_flag(remove_buffers),
            remove_ions=_flag(remove_ions),
            chains=keep_chains or None)
    except ValueError as exc:
        raise HTTPException(400, str(exc))
    (wd / "receptor.pdb").write_text(text)

    if not ligand_key:
        warnings.append("No ligand was extracted, so nothing marks the "
                        "binding site; you will have to place the docking "
                        "box by hand.")

    return JSONResponse({
        "prep_id": prep_id,
        "receptor": {"atoms": summary["atoms_kept"],
                     "dropped": summary["dropped"],
                     "kept": summary["kept"],
                     "pdb": f"/api/prep/{prep_id}/receptor.pdb"},
        "ligand": ligand_info,
        "warnings": warnings,
    })


# --------------------------------------------------------------------------- #
# Files
# --------------------------------------------------------------------------- #
FILES = {
    "viewer.pdb": ("chemical/x-pdb", "structure.pdb"),
    "receptor.pdb": ("chemical/x-pdb", "receptor_clean.pdb"),
    "ligand.mol2": ("chemical/x-mol2", "ligand.mol2"),
    "ligand.sdf": ("chemical/x-mdl-sdfile", "ligand.sdf"),
}


@router.get("/{prep_id}/{filename}")
def download(prep_id: str, filename: str):
    if filename not in FILES:
        raise HTTPException(404, "no such file")
    wd = _dir(prep_id)
    path = wd / filename
    if not path.is_file():
        raise HTTPException(404, f"{filename} has not been built yet")
    media, nice = FILES[filename]
    return FileResponse(str(path), media_type=media, filename=nice)


def prepared_inputs(prep_id: str) -> tuple[Path, Path | None]:
    """The receptor and ligand a docking job should use. For app.py."""
    wd = _dir(prep_id)
    receptor = wd / "receptor.pdb"
    if not receptor.is_file():
        raise HTTPException(400, "that structure has not been prepared yet")
    ligand = wd / "ligand.sdf"
    return receptor, (ligand if ligand.is_file() else None)
