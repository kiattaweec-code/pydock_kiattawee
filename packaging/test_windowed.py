"""Regression test: the app must start with no stdout and no stderr.

A windowed PyInstaller build has no console, so Python sets ``sys.stdout``
and ``sys.stderr`` to ``None``.  Anything that inspects them then fails on
attribute access, and 0.1.0 shipped exactly that bug: uvicorn's log
formatter calls ``sys.stdout.isatty()`` while installing its logging
config, so the app died with "Unable to configure formatter 'default'"
before the server started.  It only happened when the app was launched by
double-clicking -- running it from a terminal or a pipe gives it real
streams, which is why it passed CI.

This reproduces the condition in-process and is cheap enough to run on
every platform before the build.

    python packaging/test_windowed.py
"""
from __future__ import annotations

import sys
import urllib.request
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

failures: list[str] = []
notes: list[str] = []
real_stdout = sys.stdout


class Discard:
    """A stream with no destination, as the runtime hook falls back to."""

    def write(self, text):
        return len(text)

    def flush(self):
        pass

    def isatty(self):
        return False

    def writable(self):
        return True


def check(label, fn):
    try:
        fn()
        notes.append(f"ok   {label}")
    except BaseException as exc:                                  # noqa: BLE001
        failures.append(f"{label}: {type(exc).__name__}: {exc}")


def serves_with_stdout_none():
    """The real windowed condition: both streams are None, nothing replaced.

    This has to go through the app's own Server, not a hand-written uvicorn
    Config, or it tests the test instead of the program.
    """
    from pydock_app.launcher import Server, free_port, setup_logging

    sys.stdout = None
    sys.stderr = None
    try:
        setup_logging(False)
        port = free_port(0)
        server = Server(port)
        server.start()                   # 0.1.0 raised ValueError here
        if not server.wait_until_ready(timeout=300):
            raise RuntimeError("the server never became reachable")
        try:
            base = f"http://127.0.0.1:{port}"
            with urllib.request.urlopen(base + "/api/health", timeout=10) as r:
                if r.status != 200:
                    raise RuntimeError(f"health returned {r.status}")
            with urllib.request.urlopen(base, timeout=10) as r:
                if len(r.read()) < 4000:
                    raise RuntimeError("the UI did not render")
        finally:
            server.stop()
    finally:
        sys.stdout = real_stdout
        sys.stderr = sys.__stderr__


def serves_with_replacement_streams():
    """And with the runtime hook's fallback stream in place of a console."""
    from pydock_app.launcher import Server, free_port, setup_logging

    sys.stdout = Discard()
    sys.stderr = Discard()
    try:
        setup_logging(False)
        port = free_port(0)
        server = Server(port)
        server.start()
        if not server.wait_until_ready(timeout=300):
            raise RuntimeError("the server never became reachable")
        server.stop()
    finally:
        sys.stdout = real_stdout
        sys.stderr = sys.__stderr__


def printing_is_harmless():
    """print() must not raise when there is nowhere to print."""
    sys.stdout = None
    try:
        print("this goes nowhere")           # a no-op when stdout is None
    finally:
        sys.stdout = real_stdout


if __name__ == "__main__":
    check("print() with sys.stdout = None", printing_is_harmless)
    check("server starts and serves with sys.stdout = None", serves_with_stdout_none)
    check("server starts with the fallback streams", serves_with_replacement_streams)

    for line in notes:
        print("  " + line)
    for line in failures:
        print("  FAIL " + line)
    print(f"\n{len(notes)} passed, {len(failures)} failed\n")
    raise SystemExit(1 if failures else 0)
