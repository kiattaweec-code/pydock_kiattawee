"""Job manager: a serialised queue of docking runs with live progress.

Docking is CPU-bound and single-threaded, so jobs run one at a time in a
background worker thread while the HTTP layer stays responsive.  Each job keeps
its own working directory, a rolling log, and a cancellation flag the engine
polls between GA runs.

This is deliberately in-process: no Redis, no Celery, no database.  It is the
right shape for a single-box deployment and a prototype.  The notes in
README.md say what to change for a real multi-worker service.
"""
from __future__ import annotations

import shutil
import threading
import time
import traceback
import uuid
from collections import deque
from dataclasses import dataclass, field
from datetime import datetime, timezone
from enum import Enum
from pathlib import Path
from queue import Queue
from typing import Any

MAX_LOG_LINES = 400


class Status(str, Enum):
    QUEUED = "queued"
    PREPARING = "preparing"
    RUNNING = "running"
    DONE = "done"
    FAILED = "failed"
    CANCELLED = "cancelled"


@dataclass
class Job:
    id: str
    workdir: Path
    params: dict[str, Any]
    receptor_name: str
    ligand_names: list[str]
    status: Status = Status.QUEUED
    fraction: float = 0.0
    message: str = "queued"
    log: deque = field(default_factory=lambda: deque(maxlen=MAX_LOG_LINES))
    results: list[dict] = field(default_factory=list)
    error: str | None = None
    created: float = field(default_factory=time.time)
    started: float | None = None
    finished: float | None = None
    _cancel: threading.Event = field(default_factory=threading.Event)
    _lock: threading.Lock = field(default_factory=threading.Lock)
    # Kept so poses can be re-serialised on request without redocking.
    _payload: dict = field(default_factory=dict)

    # -- progress sink handed to the engine ----------------------------- #
    def report(self, fraction: float | None, message: str):
        with self._lock:
            if fraction is not None:
                self.fraction = max(0.0, min(1.0, float(fraction)))
            text = (message or "").strip()
            if text:
                self.message = text.split("\n")[0][:200]
                stamp = datetime.now(timezone.utc).strftime("%H:%M:%S")
                for line in text.split("\n"):
                    if line.strip():
                        self.log.append(f"{stamp}  {line.rstrip()}")

    def cancelled(self) -> bool:
        return self._cancel.is_set()

    def cancel(self):
        self._cancel.set()

    @property
    def elapsed(self) -> float:
        end = self.finished or time.time()
        return end - (self.started or self.created)

    def public(self, with_log: bool = True) -> dict:
        with self._lock:
            out = {
                "id": self.id,
                "status": self.status.value,
                "fraction": round(self.fraction, 4),
                "message": self.message,
                "receptor": self.receptor_name,
                "ligands": self.ligand_names,
                "params": self.params,
                "results": self.results,
                "error": self.error,
                "elapsed": round(self.elapsed, 1),
                "queued_at": self.created,
                "site": self._payload.get("site"),
                "ranking": self._payload.get("ranking"),
            }
            if with_log:
                out["log"] = list(self.log)
            return out


class JobManager:
    """Single-worker queue.  Thread-safe for the handful of operations used."""

    def __init__(self, root: Path, runner, max_jobs: int = 60):
        self.root = Path(root)
        self.root.mkdir(parents=True, exist_ok=True)
        self.runner = runner
        self.max_jobs = max_jobs
        self._jobs: dict[str, Job] = {}
        self._order: list[str] = []
        self._queue: Queue[str] = Queue()
        self._lock = threading.Lock()
        self._worker = threading.Thread(target=self._loop, daemon=True,
                                        name="pydock-worker")
        self._worker.start()

    # ------------------------------------------------------------------ #
    def create(self, params: dict, receptor_name: str,
               ligand_names: list[str]) -> Job:
        jid = uuid.uuid4().hex[:12]
        wd = self.root / jid
        wd.mkdir(parents=True, exist_ok=True)
        job = Job(id=jid, workdir=wd, params=params,
                  receptor_name=receptor_name, ligand_names=ligand_names)
        with self._lock:
            self._jobs[jid] = job
            self._order.append(jid)
            self._evict_locked()
        return job

    def enqueue(self, job: Job):
        job.report(0.0, f"queued behind {self.pending()} job(s)")
        self._queue.put(job.id)

    def get(self, jid: str) -> Job | None:
        with self._lock:
            return self._jobs.get(jid)

    def list(self) -> list[dict]:
        with self._lock:
            jobs = [self._jobs[j] for j in reversed(self._order)
                    if j in self._jobs]
        return [j.public(with_log=False) for j in jobs]

    def pending(self) -> int:
        return self._queue.qsize()

    def delete(self, jid: str) -> bool:
        with self._lock:
            job = self._jobs.pop(jid, None)
            if jid in self._order:
                self._order.remove(jid)
        if job is None:
            return False
        job.cancel()
        shutil.rmtree(job.workdir, ignore_errors=True)
        return True

    def _evict_locked(self):
        """Drop the oldest finished jobs once the cap is exceeded."""
        while len(self._order) > self.max_jobs:
            for jid in list(self._order):
                job = self._jobs.get(jid)
                if job and job.status in (Status.DONE, Status.FAILED,
                                          Status.CANCELLED):
                    self._order.remove(jid)
                    self._jobs.pop(jid, None)
                    shutil.rmtree(job.workdir, ignore_errors=True)
                    break
            else:
                break

    # ------------------------------------------------------------------ #
    def _loop(self):
        while True:
            jid = self._queue.get()
            job = self.get(jid)
            if job is None or job.cancelled():
                if job is not None:
                    job.status = Status.CANCELLED
                    job.finished = time.time()
                self._queue.task_done()
                continue
            job.started = time.time()
            job.status = Status.PREPARING
            try:
                self.runner(job)
                if job.cancelled():
                    job.status = Status.CANCELLED
                    job.report(job.fraction, "cancelled")
                else:
                    job.status = Status.DONE
                    job.fraction = 1.0
                    job.report(1.0, "finished")
            except Exception as exc:  # noqa: BLE001
                if job.cancelled():
                    job.status = Status.CANCELLED
                    job.report(job.fraction, "cancelled")
                else:
                    job.status = Status.FAILED
                    job.error = f"{type(exc).__name__}: {exc}"
                    job.report(job.fraction, f"failed: {job.error}")
                    job.log.append(traceback.format_exc()[-1500:])
            finally:
                job.finished = time.time()
                self._queue.task_done()
