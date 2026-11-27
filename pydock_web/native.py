"""Writing files where the person actually wants them.

The browser can only drop a download into whatever folder it was told to
use, one file at a time. But this server runs on the person's own machine,
bound to 127.0.0.1, so it can write straight into a folder they choose --
and when the app is in its own window, it can even ask with a real folder
dialog. Everything here is about that: picking a destination, and putting
files in it without overwriting their work.
"""
from __future__ import annotations

import os
import shutil
from pathlib import Path

# Set by the launcher once a native window exists; a plain browser session
# leaves it None and the caller falls back to a typed or default folder.
_folder_picker = None


def set_folder_picker(fn) -> None:
    global _folder_picker
    _folder_picker = fn


def has_folder_picker() -> bool:
    return _folder_picker is not None


def pick_folder() -> str | None:
    """Ask with a native dialog. None if cancelled or unavailable."""
    if _folder_picker is None:
        return None
    try:
        return _folder_picker()
    except Exception:                                             # noqa: BLE001
        return None


def default_folder() -> Path:
    """Somewhere obvious, created on demand.

    Documents is where people look for their own files; a Downloads folder
    is a staging area they clear out, which is the wrong place for work
    they prepared deliberately.
    """
    home = Path.home()
    for candidate in (home / "Documents", home / "เอกสาร", home):
        if candidate.is_dir():
            return candidate / "pydock"
    return home / "pydock"


def resolve_folder(folder: str | None) -> Path:
    """Where to write, as an absolute path, created if need be."""
    path = Path(os.path.expanduser(folder.strip())) if folder and folder.strip() \
        else default_folder()
    try:
        path = path.resolve()
    except OSError as exc:
        raise ValueError(f"that folder cannot be used: {exc}")
    if path.exists() and not path.is_dir():
        raise ValueError(f"{path} is a file, not a folder")
    try:
        path.mkdir(parents=True, exist_ok=True)
    except OSError as exc:
        raise ValueError(f"could not create {path}: {exc}")
    if not os.access(path, os.W_OK):
        raise ValueError(f"{path} cannot be written to")
    return path


def free_name(folder: Path, name: str) -> Path:
    """A path in `folder` that is not already taken.

    Preparing a second structure should not quietly replace the first, and
    a file the person has open in PyMOL should not change under them.
    """
    target = folder / name
    if not target.exists():
        return target
    stem, suffix = Path(name).stem, Path(name).suffix
    for n in range(2, 1000):
        candidate = folder / f"{stem}_{n}{suffix}"
        if not candidate.exists():
            return candidate
    raise ValueError(f"too many files named like {name} in {folder}")


def copy_into(folder: Path, pairs: list[tuple[Path, str]]) -> list[dict]:
    """Copy (source, wanted name) pairs in, reporting where each landed."""
    written: list[dict] = []
    for source, name in pairs:
        if not source.is_file():
            continue
        target = free_name(folder, name)
        shutil.copyfile(source, target)
        written.append({"name": target.name, "path": str(target),
                        "bytes": target.stat().st_size,
                        "renamed": target.name != name})
    if not written:
        raise ValueError("there was nothing to save")
    return written
