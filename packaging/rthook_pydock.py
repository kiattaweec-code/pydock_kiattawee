"""Runtime hook: runs inside the frozen app before anything else imports.

Two things have to be true before Numba or RDKit are touched:

* the Numba cache has to point somewhere writable, because the bundle is
  read-only and on macOS it sits inside a quarantined .app;
* the bundled ``pydock`` sources have to be importable from the unpacked
  directory, so Numba can find the source files its cache keys on.
"""
import os
import sys
from pathlib import Path

if getattr(sys, "frozen", False) and hasattr(sys, "_MEIPASS"):
    base = Path(sys._MEIPASS)

    if str(base) not in sys.path:
        sys.path.insert(0, str(base))

    # Mirrors pydock_app.paths.data_dir(), which cannot be imported this
    # early.  PYDOCK_HOME has to be honoured here too, or the cache lands
    # somewhere other than the directory the rest of the app uses.
    override = os.environ.get("PYDOCK_HOME")
    if override:
        home = Path(override).expanduser()
    elif sys.platform == "win32":
        home = Path(os.environ.get("LOCALAPPDATA")
                    or Path.home() / "AppData" / "Local") / "pydock"
    elif sys.platform == "darwin":
        home = Path.home() / "Library" / "Application Support" / "pydock"
    else:
        home = Path(os.environ.get("XDG_DATA_HOME")
                    or Path.home() / ".local" / "share") / "pydock"

    try:
        cache = home / "cache" / "numba"
        cache.mkdir(parents=True, exist_ok=True)
        os.environ.setdefault("NUMBA_CACHE_DIR", str(cache))
    except OSError:
        # A read-only home only costs us the cache, not the run.
        os.environ.setdefault("NUMBA_CACHE_DIR", "")
