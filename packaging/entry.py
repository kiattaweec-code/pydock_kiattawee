"""Frozen entry point.

PyInstaller runs its entry script as ``__main__`` with no parent package, so
``pydock_app/__main__.py`` and its relative imports cannot be used directly.
This module imports absolutely instead, and is the only thing the spec points
at.  ``python -m pydock_app`` keeps working through the other file.
"""
import multiprocessing
import sys


def run() -> int:
    from pydock_app.launcher import main
    return main()


if __name__ == "__main__":
    # On Windows every spawned process re-executes this file; without
    # freeze_support that means a new copy of the app each time.
    multiprocessing.freeze_support()
    sys.exit(run())
