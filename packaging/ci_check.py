"""Prove a freshly built app actually works, on whatever OS built it.

Starts the frozen binary, waits for the local server, docks the bundled
example, and checks the result is the pose we expect.  Exits non-zero if
anything fails, so a broken build never reaches a release.
"""
from __future__ import annotations

import json
import subprocess
import sys
import time
import urllib.error
import urllib.request
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
DIST = ROOT / "dist"
PORT = 8765
BASE = f"http://127.0.0.1:{PORT}"
HOME = ROOT / "build" / "ci_home"


def binary() -> Path:
    if sys.platform == "darwin":
        p = DIST / "pydock.app" / "Contents" / "MacOS" / "pydock"
        if p.exists():
            return p
    name = "pydock.exe" if sys.platform == "win32" else "pydock"
    p = DIST / "pydock" / name
    if not p.exists():
        raise SystemExit(f"no built app at {p}")
    return p


def get(path: str, timeout: float = 10):
    with urllib.request.urlopen(BASE + path, timeout=timeout) as r:
        return json.loads(r.read().decode())


def post(path: str, fields: dict[str, str]):
    boundary = "----pydockci"
    body = b""
    for k, v in fields.items():
        body += (f"--{boundary}\r\nContent-Disposition: form-data; "
                 f'name="{k}"\r\n\r\n{v}\r\n').encode()
    body += f"--{boundary}--\r\n".encode()
    req = urllib.request.Request(
        BASE + path, data=body, method="POST",
        headers={"Content-Type": f"multipart/form-data; boundary={boundary}"})
    with urllib.request.urlopen(req, timeout=30) as r:
        return json.loads(r.read().decode())


def main() -> int:
    exe = binary()
    print(f"launching {exe}", flush=True)
    import os

    env = dict(os.environ, PYDOCK_HOME=str(HOME))
    proc = subprocess.Popen([str(exe), "--no-open", "--port", str(PORT)],
                            env=env, stdout=subprocess.PIPE,
                            stderr=subprocess.STDOUT, text=True)
    try:
        deadline = time.time() + 240
        health = None
        while time.time() < deadline:
            if proc.poll() is not None:
                print(proc.stdout.read() if proc.stdout else "")
                raise SystemExit(f"the app exited early ({proc.returncode})")
            try:
                health = get("/api/health", timeout=3)
                break
            except (urllib.error.URLError, OSError):
                time.sleep(1)
        if health is None:
            raise SystemExit("the server never became reachable")
        print(f"  health: pydock {health['pydock']}, rdkit {health['rdkit']}, "
              f"numba {health['numba']}", flush=True)

        if not health.get("example_available"):
            raise SystemExit("the bundled example data is missing from the build")

        with urllib.request.urlopen(BASE + "/", timeout=10) as r:
            page = r.read()
        assert b"pydock" in page and len(page) > 4000, "the UI did not render"
        print(f"  web UI: {len(page)} bytes", flush=True)

        job = post("/api/jobs/example",
                   {"scoring": "empirical", "n_runs": "2", "n_evals": "60000"})
        print(f"  queued job {job['id']}", flush=True)

        deadline = time.time() + 1200
        state = {}
        while time.time() < deadline:
            state = get(f"/api/jobs/{job['id']}")
            if state["status"] in ("done", "failed", "cancelled"):
                break
            time.sleep(5)

        if state.get("status") != "done":
            print(json.dumps(state, indent=2)[:2000])
            raise SystemExit(f"the docking job ended as {state.get('status')}")

        result = state["results"][0]
        best = result["poses"][0]
        print(f"  docked {result['name']}: score {best['score']}, "
              f"rmsd {result['rmsd_to_reference']} A", flush=True)

        if best["score"] > 0:
            raise SystemExit(f"implausible score {best['score']}")
        if not result["poses"]:
            raise SystemExit("no poses came back")

        print("\nOK — the built app serves the UI and docks correctly.\n")
        return 0
    finally:
        proc.terminate()
        try:
            proc.wait(timeout=20)
        except subprocess.TimeoutExpired:
            proc.kill()


if __name__ == "__main__":
    raise SystemExit(main())
