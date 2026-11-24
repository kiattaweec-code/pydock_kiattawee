"""Where things live, whether running from source or from a frozen bundle.

PyInstaller unpacks data files to a temporary directory and points
``sys._MEIPASS`` at it, so any path derived from ``__file__`` has to be
resolved through here instead of being assumed.  Everything the app *writes*
goes to a per-user directory, because the bundle itself is read-only (and on
macOS it lives inside a quarantined .app).
"""
from __future__ import annotations

import os
import sys
from pathlib import Path

APP_NAME = "pydock"


def frozen() -> bool:
    return getattr(sys, "frozen", False) and hasattr(sys, "_MEIPASS")


def resource(*parts: str) -> Path:
    """A read-only file shipped with the app."""
    if frozen():
        return Path(sys._MEIPASS).joinpath(*parts)          # type: ignore[attr-defined]
    return Path(__file__).resolve().parent.parent.joinpath(*parts)


def data_dir() -> Path:
    """Per-user writable directory, following each platform's convention."""
    override = os.environ.get("PYDOCK_HOME")
    if override:
        base = Path(override).expanduser()
    elif sys.platform == "win32":
        base = Path(os.environ.get("LOCALAPPDATA")
                    or Path.home() / "AppData" / "Local") / APP_NAME
    elif sys.platform == "darwin":
        base = Path.home() / "Library" / "Application Support" / APP_NAME
    else:
        base = Path(os.environ.get("XDG_DATA_HOME")
                    or Path.home() / ".local" / "share") / APP_NAME
    base.mkdir(parents=True, exist_ok=True)
    return base


def jobs_dir() -> Path:
    p = data_dir() / "jobs"
    p.mkdir(parents=True, exist_ok=True)
    return p


def cache_dir() -> Path:
    p = data_dir() / "cache"
    p.mkdir(parents=True, exist_ok=True)
    return p


def log_file() -> Path:
    return data_dir() / "pydock.log"


def configure_environment() -> dict[str, str]:
    """Point every library that wants to write somewhere at a writable place.

    Numba is the one that matters: ``@njit(cache=True)`` tries to write its
    cache next to the source file, which inside a bundle is either missing or
    read-only.  Without this the app either warns on every function or fails
    outright, and recompiles ~40 kernels on every launch.
    """
    env: dict[str, str] = {}

    numba_cache = cache_dir() / "numba"
    numba_cache.mkdir(parents=True, exist_ok=True)
    env["NUMBA_CACHE_DIR"] = str(numba_cache)

    # Matplotlib and fontconfig also write caches; harmless to redirect.
    env["MPLCONFIGDIR"] = str(cache_dir() / "mpl")
    Path(env["MPLCONFIGDIR"]).mkdir(parents=True, exist_ok=True)

    env["PYDOCK_DATA"] = str(jobs_dir())

    examples = resource("test_data")
    if examples.is_dir():
        env["PYDOCK_EXAMPLES"] = str(examples)

    # setdefault, so anything the user set on the command line still wins.
    for k, v in env.items():
        os.environ.setdefault(k, v)
    return env
