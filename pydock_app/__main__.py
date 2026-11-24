"""`python -m pydock_app`, and the entry point PyInstaller freezes."""
import multiprocessing
import sys

from .launcher import main

if __name__ == "__main__":
    # Without this a frozen app on Windows re-runs the launcher in every
    # worker process it spawns, which looks like the app opening itself
    # forever.  Harmless elsewhere.
    multiprocessing.freeze_support()
    sys.exit(main())
