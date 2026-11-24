# pydock web

A browser front end for the [pydock](../README.md) docking engine: upload a
receptor and ligands, define the site, watch the genetic algorithm work, and
inspect the poses in 3D.

![the result view](../shot_viewer.png)

## Run it

```bash
pip install -e .                       # installs pydock
pip install fastapi "uvicorn[standard]" python-multipart
uvicorn pydock_web.app:app --host 127.0.0.1 --port 8000
```

Open <http://127.0.0.1:8000>. The **Run bundled example** button redocks
ergotamine into the receptor site it came from, which is also the engine's
validation case — about 75 seconds for 4 runs on one core.

3Dmol.js is vendored in `static/vendor/`, so the viewer works with no internet
access.

### Configuration

| Variable | Default | Meaning |
|---|---|---|
| `PYDOCK_DATA` | system temp dir | where job directories are written |
| `PYDOCK_EXAMPLES` | `../test_data` | folder holding the bundled example |
| `PYDOCK_MAX_UPLOAD` | 40 MB | per-file upload ceiling |
| `PYDOCK_MAX_LIGANDS` | 25 | ligand files per job |

## API

Interactive docs are at `/docs`. The whole surface:

| Method | Path | Purpose |
|---|---|---|
| `GET` | `/api/health` | engine and library versions, queue depth |
| `POST` | `/api/prepare` | parse a receptor and return its structure, bounds and detected pockets |
| `POST` | `/api/jobs` | submit (multipart: receptor, ligands, reference, settings) |
| `GET` | `/api/jobs/{id}/ranking.csv` | screening table, for multi-ligand jobs |
| `POST` | `/api/jobs/example` | queue the bundled validation case |
| `GET` | `/api/jobs` | recent jobs, newest first |
| `GET` | `/api/jobs/{id}` | status, fraction, rolling log, results |
| `POST` | `/api/jobs/{id}/cancel` | stop a running job between GA runs |
| `DELETE` | `/api/jobs/{id}` | delete the job and its files |
| `GET` | `/api/jobs/{id}/receptor.pdb` | receptor, trimmed around the site |
| `GET` | `/api/jobs/{id}/reference.sdf` | the reference ligand, if one was given |
| `GET` | `/api/jobs/{id}/poses/{i}.sdf` | ranked poses for ligand *i* |

Scripted use:

```bash
JOB=$(curl -s -X POST localhost:8000/api/jobs \
        -F receptor=@protein.pdb \
        -F ligands=@ligand.sdf \
        -F reference=@crystal.sdf \
        -F scoring=empirical -F n_runs=10 -F n_evals=250000 \
      | python3 -c 'import sys,json;print(json.load(sys.stdin)["id"])')

while :; do
  curl -s localhost:8000/api/jobs/$JOB |
    python3 -c 'import sys,json;j=json.load(sys.stdin);print(j["status"],j["message"])'
  sleep 5
done

curl -s -o poses.sdf localhost:8000/api/jobs/$JOB/poses/0.sdf
```

Every pose carries its score, the GA run that produced it, the per-term
breakdown, and — when a reference was supplied — its RMSD.

## Placing the docking box

Upload a receptor and the **Binding site** panel appears before you run
anything. It shows the structure with the docking box drawn as a wireframe
cage, and four ways to place it:

- **Detected pockets.** The whole protein is scanned (LIGSITE burial, the same
  code the engine uses to confine its search) and the biggest enclosed pockets
  come back as chips, largest first. One click moves the box onto a pocket and
  sizes it to fit.
- **A reference ligand.** *Fit reference* snaps the box around a known binder,
  which is also what supplies the RMSD column in the results.
- **Click to centre.** Tick the box, click any atom, and the box recentres
  there — the quickest way to aim at a residue you already know.
- **Sliders and numbers.** Six sliders for centre and size, mirrored by the
  numeric fields in the sidebar; either drives the other.

The readout shows the box volume and warns when it passes ~27,000 Å<sup>3</sup>,
because search cost grows with the volume and an oversized box is the most
common reason a run is slow and misses.

Touching the box takes over the site: the site mode switches to *Coordinates*
and the exact centre and size you see are what gets submitted. A reference
ligand stays attached either way, so you keep the RMSD even when you place
the box by hand.

## Docking more than one ligand

Three ways, which combine freely:

- select several ligand files at once,
- give one SDF or SMILES file holding a whole library — every record is read
  and docked, not just the first,
- or both.

A record that cannot be prepared is skipped with a line in the log rather than
failing the job. Results come back as a **Screening** table ranked by best
score, with a toggle to rank by **ligand efficiency** (score per heavy atom)
instead — worth using, because raw docking scores grow with molecular size and
a ranking by score alone rewards whatever is biggest. Clicking a row loads
that ligand's poses into the viewer, and the whole table downloads as CSV.

For a library the empirical grids are built once and reused across every
ligand, since they depend only on the receptor and the box.

## How it is put together

```
pydock_web/
  app.py              FastAPI routes, upload validation, the worker body
  jobs.py             in-process job queue, progress, cancellation
  static/
    index.html        the page
    style.css         the look
    app.js            submission, polling, results table, 3D viewer
    vendor/3Dmol-min.js
```

Docking is CPU-bound and single-threaded, so **one worker thread drains a
queue** and the request handlers never block. Progress reaches the browser
because the engine takes a `progress(fraction, message)` callback, which the
job object turns into a percentage and a rolling log; the client polls once a
second. Cancellation is a flag the engine checks between GA runs.

Uploads are written to a per-job directory and parsed only inside the worker,
so a malformed file fails that one job and nothing else — verified by sending
a text file renamed `.sdf`, which produced a clean `failed` status with the
parse error while the server stayed up.

The receptor is trimmed to residues within 20 Å of the docked site before
being sent to the viewer (193 KB → 92 KB on the example), which keeps the
first paint quick.

## Verified behaviour

The backend was exercised over HTTP and the front end driven in a real
browser and in jsdom:

- bundled example submitted through the API finished in 75 s and reproduced
  the engine's validated result — top pose **−14.199, 1.19 Å RMSD**, the
  correct binding mode ranked first with a clear margin over the next pose
  (−9.7)
- input validation returns a usable message for every bad input tried: wrong
  extension, missing site definition, out-of-range numbers, a partial centre,
  a malformed residue spec, an unknown scoring function
- a corrupt ligand file fails only its own job; the server stays healthy
- cancel stops a 20-run job within one GA run
- 29 front-end assertions pass in jsdom with no page errors, covering
  site-mode switching, submission, live progress, the pose table, pose
  switching in the viewer, the display toggles, and the recent-jobs list
- 25 further assertions pass in headless Chrome against the live server,
  covering pocket detection and the chips, the box drawing, sliders, numeric
  fields, click-to-centre, *Fit reference*, the oversized-box warning, and the
  screening table including re-ranking and ligand switching

Two bugs that these caught and that are now fixed: box sizes with one decimal
failed the numeric inputs' `step` validation, which made Chrome refuse to
submit the form *silently*; and pose SDFs 404'd for library jobs because the
indexed filenames the engine writes no longer matched what the API looked for.

## Limitations

- **Jobs live in memory.** A restart loses the job list, though the files on
  disk survive. Fine for one box; a real service wants a database.
- **One job at a time.** Deliberate on a small machine, since each job
  saturates a core. For throughput, move the worker to a process pool or
  Celery and key jobs by id in Redis.
- **No authentication and no quotas.** Bind to `127.0.0.1`, or put it behind a
  reverse proxy with auth before exposing it. Anyone who can reach the port
  can queue CPU-heavy work and upload files.
- **Uploads are trusted to be chemistry.** They are size- and
  extension-checked and parsed in the worker, but a hostile file is only as
  contained as RDKit's parsers.
- **Old jobs are evicted after 60**, oldest finished first, taking their
  directories with them. There is no scheduled cleanup beyond that.
