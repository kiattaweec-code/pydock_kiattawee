"""FastAPI service wrapping the pydock engine.

    uvicorn pydock_web.app:app --host 127.0.0.1 --port 8000

Design notes
------------
* Docking is CPU-bound, so jobs are serialised through one worker thread
  (``jobs.JobManager``) and the request handlers never block on a run.
* Uploads are written to a per-job directory and parsed only inside the worker,
  so a malformed file fails that job and nothing else.
* Results are returned as JSON metadata plus SDF/PDB text fetched separately,
  which keeps the status poll small while a 40-pose run is in flight.
"""
from __future__ import annotations

import os
import sys
import tempfile
from pathlib import Path

import numpy as np
from fastapi import FastAPI, File, Form, HTTPException, UploadFile
from fastapi.responses import FileResponse, JSONResponse, PlainTextResponse
from fastapi.staticfiles import StaticFiles

from rdkit import Chem

from pydock import DockingConfig, GAConfig, __version__, run_docking
from pydock import cavity as cavmod
from pydock import ligand as ligmod
from pydock import receptor as recmod
from pydock.dock import ligand_efficiency, resolve_site

from . import prep as prepmod
from .jobs import JobManager, Status

HERE = Path(__file__).resolve().parent


def _bundled(*parts: str) -> Path:
    """Locate a shipped file, whether running from source or frozen.

    PyInstaller unpacks data files to sys._MEIPASS, so a path built from
    __file__ points into the wrong place inside a packaged app.
    """
    meipass = getattr(sys, "_MEIPASS", None)
    if meipass:
        p = Path(meipass).joinpath(*parts)
        if p.exists():
            return p
    return HERE.parent.joinpath(*parts)


STATIC = _bundled("pydock_web", "static") if getattr(sys, "_MEIPASS", None) \
    else HERE / "static"
DATA_ROOT = Path(os.environ["PYDOCK_DATA"]) if os.environ.get("PYDOCK_DATA") \
    else Path(tempfile.gettempdir()) / "pydock_jobs"
EXAMPLE_DIR = Path(os.environ["PYDOCK_EXAMPLES"]) if os.environ.get("PYDOCK_EXAMPLES") \
    else _bundled("test_data")

MAX_UPLOAD = int(os.environ.get("PYDOCK_MAX_UPLOAD", 40 * 1024 * 1024))
MAX_LIGANDS = int(os.environ.get("PYDOCK_MAX_LIGANDS", 25))

RECEPTOR_EXT = {".pdb", ".ent"}
LIGAND_EXT = {".sdf", ".mol", ".mol2", ".pdb", ".smi", ".smiles"}

SCORING = {"empirical", "chemscore", "chemplp", "plp", "forcefield", "ff",
           "autodock"}

app = FastAPI(title="pydock", version=__version__,
              description="Flexible-ligand molecular docking in the browser")


# --------------------------------------------------------------------------- #
# The worker body
# --------------------------------------------------------------------------- #
def run_job(job):
    p = job.params
    wd = job.workdir
    poses_dir = wd / "poses"

    receptor = str(wd / job.receptor_name)
    ligands = [str(wd / n) for n in job.ligand_names]
    reference = str(wd / p["reference_name"]) if p.get("reference_name") else None

    cfg = DockingConfig(
        receptor=receptor,
        ligands=ligands,
        reference_ligand=reference,
        centre=tuple(p["centre"]) if p.get("centre") else None,
        radius=p.get("radius"),
        box_size=tuple(p["box_size"]) if p.get("box_size") else None,
        site_residues=[tuple(r) for r in p["site_residues"]]
        if p.get("site_residues") else None,
        padding=p.get("padding", 4.0),
        scoring=p["scoring"],
        spacing=p.get("spacing", 0.4),
        n_runs=p["n_runs"],
        n_poses=p["n_poses"],
        keep_waters=p.get("keep_waters", False),
        rotate_receptor_oh=p.get("rotate_receptor_oh", True),
        detect_cavity=p.get("detect_cavity", True),
        output_dir=str(poses_dir),
        verbose=False,
        progress=job.report,
        cancelled=job.cancelled,
        ga=GAConfig(
            pop_size=p.get("pop_size", 100),
            n_islands=p.get("n_islands", 5),
            seed=p.get("seed", 1),
            autoscale=not p.get("n_evals"),
            n_evals=p.get("n_evals") or 60_000,
            early_stop=p.get("n_evals") or 60_000,
            ls_frequency=p.get("ls_frequency", 0.10),
        ),
    )
    job.status = Status.RUNNING
    results = run_docking(cfg)

    # Remember where the site ended up, for the viewer.
    rec = recmod.prepare_receptor(receptor, keep_waters=cfg.keep_waters)
    out: list[dict] = []
    for i, r in enumerate(results):
        # Take the name the engine actually wrote: library records get an
        # index prefix, so recomputing it here would drift out of step.
        sdf_path = Path(r.pose_file) if r.pose_file else None
        sdf_name = sdf_path.name if sdf_path else None
        centre = r.solutions[0].coords[r.ligand.heavy_idx].mean(axis=0)
        out.append({
            "index": i,
            "name": r.name,
            "n_atoms": int(r.ligand.n_atoms),
            "n_heavy": int(len(r.ligand.heavy_idx)),
            "n_rotatable": int(r.ligand.n_rot),
            "score": round(float(r.solutions[0].score), 3),
            "efficiency": round(ligand_efficiency(r), 4),
            "seconds": round(r.seconds, 1),
            "n_modes": len(r.clusters),
            "rmsd_to_reference": (None if r.rmsd_to_reference is None
                                  else round(r.rmsd_to_reference, 2)),
            "centre": [round(float(x), 2) for x in centre],
            "sdf": sdf_name if (sdf_path and sdf_path.exists()) else None,
            "poses": [
                {
                    "rank": k + 1,
                    "run": int(s.run),
                    "score": round(float(s.score), 3),
                    "rmsd": None if s.rmsd is None else round(float(s.rmsd), 2),
                    "terms": {kk: round(float(vv), 3)
                              for kk, vv in s.terms.items() if kk != "total"},
                }
                for k, s in enumerate(r.solutions[:cfg.n_poses])
            ],
        })
    job.results = out
    job._payload["receptor_path"] = receptor
    job._payload["n_receptor_atoms"] = int(rec.n_atoms)
    try:
        (centre, half), _ref = resolve_site(cfg, rec)
        job._payload["site"] = {
            "centre": [round(float(v), 2) for v in centre],
            "size": [round(float(v), 2) for v in 2 * np.asarray(half)],
        }
        job._payload["site_centre"] = job._payload["site"]["centre"]
    except Exception:                                             # noqa: BLE001
        job._payload["site_centre"] = out[0]["centre"] if out else [0.0, 0.0, 0.0]
    if len(out) > 1 and (poses_dir / "ranking.csv").exists():
        job._payload["ranking"] = "ranking.csv"


def _safe(name: str) -> str:
    return "".join(c if c.isalnum() or c in "-_." else "_" for c in name)[:60] \
        or "ligand"


manager = JobManager(DATA_ROOT, run_job)

prepmod.configure(DATA_ROOT)
app.include_router(prepmod.router)


# --------------------------------------------------------------------------- #
# Upload validation
# --------------------------------------------------------------------------- #
def _store(upload: UploadFile, wd: Path, allowed: set[str]) -> str:
    name = Path(upload.filename or "file").name
    ext = Path(name).suffix.lower()
    if ext not in allowed:
        raise HTTPException(
            400, f"{name}: unsupported extension {ext or '(none)'}; "
                 f"expected one of {', '.join(sorted(allowed))}")
    data = upload.file.read(MAX_UPLOAD + 1)
    if not data:
        raise HTTPException(400, f"{name} is empty")
    if len(data) > MAX_UPLOAD:
        raise HTTPException(
            413, f"{name} is larger than {MAX_UPLOAD // (1024 * 1024)} MB")
    safe = _safe(name)
    (wd / safe).write_bytes(data)
    return safe


def _as_int(value, lo, hi, default, label):
    if value in (None, ""):
        return default
    try:
        n = int(float(value))
    except (TypeError, ValueError):
        raise HTTPException(400, f"{label} must be a number")
    if not lo <= n <= hi:
        raise HTTPException(400, f"{label} must be between {lo} and {hi}")
    return n


def _as_float(value, lo, hi, default, label):
    if value in (None, ""):
        return default
    try:
        x = float(value)
    except (TypeError, ValueError):
        raise HTTPException(400, f"{label} must be a number")
    if not lo <= x <= hi:
        raise HTTPException(400, f"{label} must be between {lo} and {hi}")
    return x


# --------------------------------------------------------------------------- #
# Routes
# --------------------------------------------------------------------------- #
@app.get("/api/health")
def health():
    import numba
    import rdkit

    return {
        "status": "ok",
        "pydock": __version__,
        "rdkit": rdkit.__version__,
        "numba": numba.__version__,
        "queued": manager.pending(),
        "scoring_functions": ["empirical", "chemplp", "forcefield"],
        "example_available": EXAMPLE_DIR.is_dir()
        and (EXAMPLE_DIR / "receptor_noH.pdb").exists(),
    }


@app.post("/api/prepare")
def prepare(receptor: UploadFile = File(...),
            reference: UploadFile | None = File(None),
            find_pockets: str = Form("true")):
    """Parse a receptor before docking so the browser can draw the site.

    Returns the structure to display, its bounds, the pockets found across the
    whole protein, and the box implied by a reference ligand if one is given.
    Nothing is kept on the server: the browser re-sends the files on submit,
    which keeps this endpoint stateless and leaves no orphaned uploads.
    """
    import shutil
    import tempfile as tf

    wd = Path(tf.mkdtemp(prefix="pydock_prep_"))
    try:
        rec_name = _store(receptor, wd, RECEPTOR_EXT)
        rec_path = wd / rec_name
        try:
            rec = recmod.prepare_receptor(str(rec_path))
        except Exception as exc:                                  # noqa: BLE001
            raise HTTPException(
                400, f"could not read {rec_name} as a PDB: {exc}")

        heavy = np.asarray([e != "H" for e in rec.element])
        if not heavy.any():
            raise HTTPException(400, f"{rec_name} has no heavy atoms")
        xyz = rec.coords[heavy]
        lower, upper = xyz.min(axis=0), xyz.max(axis=0)

        out: dict = {
            "receptor_name": rec_name,
            "n_atoms": int(rec.n_atoms),
            "n_residues": len(set(rec.resid)),
            "bounds": {
                "lower": [round(float(v), 2) for v in lower],
                "upper": [round(float(v), 2) for v in upper],
                "centre": [round(float(v), 2) for v in (lower + upper) / 2],
                "size": [round(float(v), 2) for v in (upper - lower)],
            },
            "pockets": [],
            "reference": None,
        }

        # Strip hydrogens from what the viewer receives: they treble the file
        # size and the cartoon and sticks do not need them.
        lines = [ln for ln in rec_path.read_text(errors="replace").splitlines(True)
                 if ln[:6] in ("ATOM  ", "HETATM")
                 and (ln[76:78].strip() or ln[12:16].strip()[:1]).upper() != "H"]
        out["pdb"] = "".join(lines) + "END\n"

        if reference is not None and reference.filename:
            ref_name = _store(reference, wd, LIGAND_EXT)
            try:
                ref = ligmod.load_ligand(str(wd / ref_name))
            except Exception as exc:                              # noqa: BLE001
                raise HTTPException(
                    400, f"could not read the reference ligand: {exc}")
            rc = ref.coords0[ref.heavy_idx]
            centre, half = recmod.site_from_ligand(rc, padding=4.0)
            out["reference"] = {
                "name": ref.name,
                "centre": [round(float(v), 2) for v in centre],
                "size": [round(float(v), 2) for v in 2 * half],
                "sdf": Chem.MolToMolBlock(ref.mol) + "$$$$\n",
            }

        if find_pockets.lower() in ("1", "true", "yes", "on"):
            try:
                for p in cavmod.find_pockets(rec):
                    out["pockets"].append({
                        "rank": p.rank,
                        "centre": [round(float(v), 2) for v in p.centre],
                        "size": [round(float(v), 2) for v in p.size],
                        "volume": round(p.volume),
                    })
            except Exception as exc:                              # noqa: BLE001
                out["pocket_error"] = f"{type(exc).__name__}: {exc}"
        return out
    finally:
        shutil.rmtree(wd, ignore_errors=True)


@app.post("/api/jobs")
def submit(
    receptor: UploadFile | None = File(None),
    ligands: list[UploadFile] = File([]),
    reference: UploadFile | None = File(None),
    prep_id: str = Form(""),
    dock_prepared_ligand: str = Form("false"),
    scoring: str = Form("empirical"),
    n_runs: str = Form("10"),
    n_poses: str = Form("5"),
    n_evals: str = Form(""),
    pop_size: str = Form("100"),
    n_islands: str = Form("5"),
    spacing: str = Form("0.4"),
    padding: str = Form("4.0"),
    seed: str = Form("1"),
    centre_x: str = Form(""),
    centre_y: str = Form(""),
    centre_z: str = Form(""),
    radius: str = Form(""),
    size_x: str = Form(""),
    size_y: str = Form(""),
    size_z: str = Form(""),
    site_residues: str = Form(""),
    keep_waters: str = Form("false"),
    detect_cavity: str = Form("true"),
    rotate_receptor_oh: str = Form("true"),
):
    if scoring not in SCORING:
        raise HTTPException(400, f"unknown scoring function: {scoring}")
    ligands = [f for f in ligands if f and f.filename]
    prep_id = prep_id.strip()

    # A structure prepared on the preparation page supplies the receptor, and
    # its extracted ligand marks the binding site -- so neither has to make
    # the round trip through the user's downloads folder.
    prep_receptor = prep_ligand = None
    if prep_id:
        prep_receptor, prep_ligand = prepmod.prepared_inputs(prep_id)
    elif receptor is None or not receptor.filename:
        raise HTTPException(400, "a receptor file is required")

    dock_prepared = (prep_ligand is not None
                     and dock_prepared_ligand.lower() in ("1", "true", "yes", "on"))
    if not ligands and not dock_prepared:
        raise HTTPException(
            400, "at least one ligand file is required"
            + (" (or tick the box to redock the prepared ligand)"
               if prep_ligand is not None else ""))
    if len(ligands) > MAX_LIGANDS:
        raise HTTPException(400, f"at most {MAX_LIGANDS} ligand files per job")

    params: dict = {
        "scoring": scoring,
        "n_runs": _as_int(n_runs, 1, 50, 10, "runs"),
        "n_poses": _as_int(n_poses, 1, 50, 5, "poses"),
        "n_evals": _as_int(n_evals, 0, 2_000_000, 0, "evaluations"),
        "pop_size": _as_int(pop_size, 10, 500, 100, "population"),
        "n_islands": _as_int(n_islands, 1, 20, 5, "islands"),
        "spacing": _as_float(spacing, 0.2, 1.0, 0.4, "grid spacing"),
        "padding": _as_float(padding, 0.0, 15.0, 4.0, "padding"),
        "seed": _as_int(seed, 0, 2 ** 31 - 1, 1, "seed"),
        "keep_waters": keep_waters.lower() in ("1", "true", "yes", "on"),
        "detect_cavity": detect_cavity.lower() in ("1", "true", "yes", "on"),
        "rotate_receptor_oh": rotate_receptor_oh.lower() in ("1", "true", "yes", "on"),
    }
    params["n_poses"] = min(params["n_poses"], params["n_runs"])

    if site_residues.strip():
        res = []
        for tok in site_residues.replace(",", " ").split():
            try:
                chain, num = (tok.split(":") if ":" in tok
                              else (tok[0], tok[1:]))
                res.append([chain, int(num)])
            except (ValueError, IndexError):
                raise HTTPException(
                    400, f"could not read residue '{tok}'; use CHAIN:NUMBER")
        params["site_residues"] = res

    if centre_x.strip() or centre_y.strip() or centre_z.strip():
        if not (centre_x.strip() and centre_y.strip() and centre_z.strip()):
            raise HTTPException(400, "give all three centre coordinates, or none")
        params["centre"] = [
            _as_float(centre_x, -1e4, 1e4, 0.0, "centre x"),
            _as_float(centre_y, -1e4, 1e4, 0.0, "centre y"),
            _as_float(centre_z, -1e4, 1e4, 0.0, "centre z"),
        ]
        if size_x.strip() or size_y.strip() or size_z.strip():
            params["box_size"] = [
                _as_float(size_x, 6.0, 80.0, 22.0, "box size x"),
                _as_float(size_y, 6.0, 80.0, 22.0, "box size y"),
                _as_float(size_z, 6.0, 80.0, 22.0, "box size z"),
            ]
        else:
            params["radius"] = _as_float(radius, 3.0, 40.0, 11.0, "radius")

    job = manager.create(params, "pending", [])
    try:
        if prep_receptor is not None:
            job.receptor_name = "receptor_prepared.pdb"
            (job.workdir / job.receptor_name).write_bytes(
                prep_receptor.read_bytes())
        else:
            job.receptor_name = _store(receptor, job.workdir, RECEPTOR_EXT)

        job.ligand_names = [_store(f, job.workdir, LIGAND_EXT) for f in ligands]
        if dock_prepared:
            (job.workdir / "prepared_ligand.sdf").write_bytes(
                prep_ligand.read_bytes())
            job.ligand_names.append("prepared_ligand.sdf")

        if reference is not None and reference.filename:
            params["reference_name"] = _store(reference, job.workdir, LIGAND_EXT)
        elif prep_ligand is not None:
            # The ligand that came out of the crystal structure is the best
            # possible definition of the site, and it gives an RMSD to
            # measure a redock against.
            params["reference_name"] = "reference_prepared.sdf"
            (job.workdir / "reference_prepared.sdf").write_bytes(
                prep_ligand.read_bytes())
    except HTTPException:
        manager.delete(job.id)
        raise

    if not (params.get("reference_name") or params.get("centre")
            or params.get("site_residues")):
        manager.delete(job.id)
        raise HTTPException(
            400, "define the binding site: upload a reference ligand, "
                 "or give a centre and radius, or list site residues")

    manager.enqueue(job)
    return JSONResponse({"id": job.id, "status": job.status.value}, status_code=202)


@app.post("/api/jobs/example")
def submit_example(scoring: str = Form("empirical"),
                   n_runs: str = Form("4"),
                   n_evals: str = Form("120000")):
    """Queue the bundled ergotamine redocking case."""
    rec = EXAMPLE_DIR / "receptor_noH.pdb"
    lig = EXAMPLE_DIR / "ligand_start.sdf"
    ref = EXAMPLE_DIR / "reference.sdf"
    if not (rec.exists() and lig.exists() and ref.exists()):
        raise HTTPException(404, f"example data not found in {EXAMPLE_DIR}")
    if scoring not in SCORING:
        raise HTTPException(400, f"unknown scoring function: {scoring}")

    runs = _as_int(n_runs, 1, 50, 4, "runs")
    params = {
        "scoring": scoring, "n_runs": runs, "n_poses": runs,
        "n_evals": _as_int(n_evals, 0, 2_000_000, 120_000, "evaluations"),
        "pop_size": 100, "n_islands": 5, "spacing": 0.4, "padding": 4.0,
        "seed": 11, "keep_waters": False, "detect_cavity": True,
        "rotate_receptor_oh": True, "reference_name": "reference.sdf",
    }
    job = manager.create(params, "receptor_noH.pdb", ["ligand_start.sdf"])
    for src in (rec, lig, ref):
        (job.workdir / src.name).write_bytes(src.read_bytes())
    manager.enqueue(job)
    return JSONResponse({"id": job.id, "status": job.status.value}, status_code=202)


@app.get("/api/jobs")
def list_jobs():
    return {"jobs": manager.list(), "queued": manager.pending()}


@app.get("/api/jobs/{jid}")
def job_status(jid: str):
    job = manager.get(jid)
    if job is None:
        raise HTTPException(404, "no such job")
    return job.public()


@app.post("/api/jobs/{jid}/cancel")
def cancel(jid: str):
    job = manager.get(jid)
    if job is None:
        raise HTTPException(404, "no such job")
    if job.status in (Status.DONE, Status.FAILED, Status.CANCELLED):
        return {"id": jid, "status": job.status.value, "cancelled": False}
    job.cancel()
    return {"id": jid, "status": job.status.value, "cancelled": True}


@app.delete("/api/jobs/{jid}")
def delete(jid: str):
    if not manager.delete(jid):
        raise HTTPException(404, "no such job")
    return {"deleted": jid}


@app.get("/api/jobs/{jid}/receptor.pdb")
def receptor_pdb(jid: str, around: float = 20.0):
    """Receptor for the viewer, trimmed to residues near the docked site."""
    job = manager.get(jid)
    if job is None:
        raise HTTPException(404, "no such job")
    path = job._payload.get("receptor_path") or str(job.workdir / job.receptor_name)
    if not Path(path).exists():
        raise HTTPException(404, "receptor file is gone")
    centre = np.asarray(job._payload.get("site_centre") or [0.0, 0.0, 0.0])
    around = max(6.0, min(float(around), 60.0))

    keep_res: set[tuple] = set()
    rows: list[str] = []
    with open(path) as fh:
        for line in fh:
            if line[:6] not in ("ATOM  ", "HETATM"):
                continue
            rows.append(line)
            try:
                xyz = np.array([float(line[30:38]), float(line[38:46]),
                                float(line[46:54])])
            except ValueError:
                continue
            if np.linalg.norm(xyz - centre) <= around:
                keep_res.add((line[21], line[22:27]))
    out = [r for r in rows if (r[21], r[22:27]) in keep_res]
    if not out:
        out = rows
    return PlainTextResponse("".join(out) + "END\n",
                             media_type="chemical/x-pdb")


@app.get("/api/jobs/{jid}/ranking.csv")
def ranking_csv(jid: str):
    """The screening table: one row per ligand, best score first."""
    job = manager.get(jid)
    if job is None:
        raise HTTPException(404, "no such job")
    path = job.workdir / "poses" / "ranking.csv"
    if not path.exists():
        raise HTTPException(404, "this job docked a single ligand")
    return FileResponse(path, media_type="text/csv",
                        filename=f"{job.id}_ranking.csv")


@app.get("/api/jobs/{jid}/reference.sdf")
def reference_sdf(jid: str):
    """The reference ligand, for overlaying on the docked pose."""
    job = manager.get(jid)
    if job is None:
        raise HTTPException(404, "no such job")
    name = job.params.get("reference_name")
    if not name:
        return PlainTextResponse("", media_type="chemical/x-mdl-sdfile")
    path = job.workdir / name
    if not path.exists():
        return PlainTextResponse("", media_type="chemical/x-mdl-sdfile")
    if path.suffix.lower() in (".sdf", ".mol"):
        return FileResponse(path, media_type="chemical/x-mdl-sdfile")
    # Convert mol2/pdb references so the viewer only ever handles SDF.
    from rdkit import Chem
    loader = (Chem.MolFromMol2File if path.suffix.lower() == ".mol2"
              else Chem.MolFromPDBFile)
    mol = loader(str(path), removeHs=False, sanitize=True)
    if mol is None:
        return PlainTextResponse("", media_type="chemical/x-mdl-sdfile")
    return PlainTextResponse(Chem.MolToMolBlock(mol) + "$$$$\n",
                             media_type="chemical/x-mdl-sdfile")


@app.get("/api/jobs/{jid}/poses/{index}.sdf")
def poses_sdf(jid: str, index: int):
    job = manager.get(jid)
    if job is None:
        raise HTTPException(404, "no such job")
    if index < 0 or index >= len(job.results):
        raise HTTPException(404, "no such ligand in this job")
    name = job.results[index].get("sdf")
    path = job.workdir / "poses" / (name or "")
    if not name or not path.exists():
        raise HTTPException(404, "poses have not been written")
    return FileResponse(path, media_type="chemical/x-mdl-sdfile",
                        filename=f"{job.id}_{name}")


app.mount("/", StaticFiles(directory=str(STATIC), html=True), name="static")
