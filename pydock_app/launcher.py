"""Desktop entry point: run the server locally and show it in a window.

The app is the same FastAPI service the browser version uses, started on a
free loopback port inside this process.  What the user sees is a native window
pointed at it -- pywebview where that is available (WKWebView on macOS, the
Edge WebView2 runtime on Windows), and otherwise their default browser.

Nothing listens on an external interface and nothing is uploaded anywhere:
every structure stays in the per-user data directory.
"""
from __future__ import annotations

import argparse
import logging
import socket
import sys
import threading
import time
import urllib.error
import urllib.request
import webbrowser
from pathlib import Path

from . import paths

WINDOW_TITLE = "pydock"
STARTUP_TIMEOUT = 180.0          # first launch compiles Numba kernels
log = logging.getLogger("pydock.app")


# --------------------------------------------------------------------------- #
def free_port(preferred: int = 0) -> int:
    """Ask the OS for a free loopback port (or check the one requested)."""
    with socket.socket(socket.AF_INET, socket.SOCK_STREAM) as s:
        s.setsockopt(socket.SOL_SOCKET, socket.SO_REUSEADDR, 1)
        try:
            s.bind(("127.0.0.1", preferred))
        except OSError:
            if preferred:
                log.warning("port %d is in use; picking another", preferred)
                s.bind(("127.0.0.1", 0))
            else:
                raise
        return s.getsockname()[1]


def setup_logging(verbose: bool):
    handlers: list[logging.Handler] = [logging.StreamHandler(sys.stderr)]
    try:
        handlers.append(logging.FileHandler(paths.log_file(), encoding="utf-8"))
    except OSError:
        pass                                   # a read-only home is survivable
    logging.basicConfig(
        level=logging.DEBUG if verbose else logging.INFO,
        format="%(asctime)s  %(levelname)-7s %(name)s  %(message)s",
        handlers=handlers, force=True)


class Server:
    """uvicorn on a background thread, with a clean shutdown."""

    def __init__(self, port: int):
        self.port = port
        self.url = f"http://127.0.0.1:{port}"
        self._server = None
        self._thread: threading.Thread | None = None
        self._error: BaseException | None = None

    def start(self):
        import uvicorn

        from pydock_web.app import app

        # Pin the plain stack: the optional uvloop/httptools/websockets
        # extras buy nothing on a loopback socket and add ~30 MB to the
        # packaged app, so they are excluded from the bundle entirely.
        config = uvicorn.Config(app, host="127.0.0.1", port=self.port,
                                log_level="warning", access_log=False,
                                loop="asyncio", http="h11", ws="none")
        self._server = uvicorn.Server(config)
        # The app owns the process lifetime, not uvicorn: without this a
        # Ctrl-C handler from uvicorn fights the window event loop.
        self._server.install_signal_handlers = lambda: None

        def run():
            try:
                self._server.run()
            except BaseException as exc:                      # noqa: BLE001
                self._error = exc
                log.exception("the local server stopped")

        self._thread = threading.Thread(target=run, daemon=True,
                                        name="pydock-server")
        self._thread.start()

    def wait_until_ready(self, timeout: float = STARTUP_TIMEOUT) -> bool:
        deadline = time.time() + timeout
        while time.time() < deadline:
            if self._error is not None:
                return False
            try:
                with urllib.request.urlopen(self.url + "/api/health", timeout=2) as r:
                    if r.status == 200:
                        return True
            except (urllib.error.URLError, OSError, ConnectionError):
                time.sleep(0.25)
        return False

    def stop(self):
        if self._server is not None:
            self._server.should_exit = True
        if self._thread is not None:
            self._thread.join(timeout=8)


# --------------------------------------------------------------------------- #
def show_window(url: str) -> bool:
    """Native window via pywebview.  False if it is unavailable."""
    try:
        import webview                                        # noqa: PLC0415
    except Exception:                                         # noqa: BLE001
        return False
    try:
        webview.create_window(WINDOW_TITLE, url, width=1480, height=980,
                              min_size=(900, 640), text_select=True)
        webview.start()                 # blocks until the window is closed
        return True
    except Exception:                                         # noqa: BLE001
        log.warning("could not open a native window; falling back to the browser",
                    exc_info=True)
        return False


def wait_for_exit(url: str):
    """Browser fallback: keep the process alive until interrupted."""
    print(f"\n  pydock is running at {url}")
    print("  Your browser should have opened. If not, paste that address in.")
    print("  Leave this window open while you work; press Ctrl-C to quit.\n")
    try:
        while True:
            time.sleep(0.5)
    except KeyboardInterrupt:
        print("\n  shutting down…")


def selftest() -> int:
    """Print what is installed, where it writes, and whether caching works.

    The kernel cache is the one thing that silently degrades: when it misses,
    the app still works but spends a minute recompiling at every launch, so
    it is worth being able to check directly.
    """
    import inspect
    import os
    import time

    print(f"\npydock self-test")
    print(f"  frozen            {paths.frozen()}")
    print(f"  executable        {sys.executable}")
    if paths.frozen():
        print(f"  bundle            {sys._MEIPASS}")   # type: ignore[attr-defined]
    print(f"  data directory    {paths.data_dir()}")
    print(f"  numba cache dir   {os.environ.get('NUMBA_CACHE_DIR') or '(unset)'}")

    import numba
    import rdkit

    from pydock import __version__, geometry
    print(f"\n  pydock {__version__} · rdkit {rdkit.__version__} "
          f"· numba {numba.__version__}")

    src = inspect.getfile(geometry)
    print(f"\n  geometry module   {src}")
    print(f"  source present    {os.path.exists(src)}")
    py = os.path.splitext(src)[0] + ".py"
    print(f"  .py beside it     {os.path.exists(py)}  ({py})")

    from pydock_web.app import EXAMPLE_DIR, STATIC
    print(f"\n  static files      {STATIC}  exists={STATIC.exists()}")
    print(f"  example data      {EXAMPLE_DIR}  exists={EXAMPLE_DIR.exists()}")

    print("\n  compiling a kernel…")
    import numpy as np
    a = np.zeros((4, 3))
    b = np.ones((4, 3))
    t = time.time()
    geometry.rmsd(a, b)
    dt = time.time() - t

    cached = list(Path(os.environ.get("NUMBA_CACHE_DIR") or ".").rglob("*.nbi")) \
        if os.environ.get("NUMBA_CACHE_DIR") else []
    print(f"  first call        {dt:.2f}s")
    print(f"  cache files       {len(cached)}")
    if cached:
        print("  kernel cache      WORKING — later launches start quickly")
        ok = True
    else:
        print("  kernel cache      NOT WORKING — every launch will recompile")
        ok = False
    print()
    return 0 if ok else 2


# --------------------------------------------------------------------------- #
def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(
        prog="pydock", description="Flexible-ligand molecular docking.")
    ap.add_argument("--port", type=int, default=0,
                    help="port to serve on (default: pick a free one)")
    ap.add_argument("--browser", action="store_true",
                    help="open in your default browser instead of a window")
    ap.add_argument("--no-open", action="store_true",
                    help="start the server but open nothing")
    ap.add_argument("--data-dir", help="where jobs and caches are stored")
    ap.add_argument("-v", "--verbose", action="store_true")
    ap.add_argument("--version", action="store_true")
    ap.add_argument("--selftest", action="store_true",
                    help="check the installation and print a diagnosis")
    args = ap.parse_args(argv)

    if args.data_dir:
        import os
        os.environ["PYDOCK_HOME"] = args.data_dir

    env = paths.configure_environment()
    setup_logging(args.verbose)

    if args.version:
        from pydock import __version__
        print(f"pydock {__version__}")
        return 0

    log.info("data directory: %s", paths.data_dir())
    log.debug("environment: %s", env)

    if args.selftest:
        return selftest()

    port = free_port(args.port)
    server = Server(port)
    print(f"  starting pydock…  (first launch compiles the search kernels, "
          f"which takes a minute)")
    server.start()

    if not server.wait_until_ready():
        log.error("the local server did not come up in time")
        print("\n  pydock could not start. The log is at "
              f"{paths.log_file()}\n", file=sys.stderr)
        server.stop()
        return 1

    log.info("serving on %s", server.url)
    try:
        if args.no_open:
            wait_for_exit(server.url)
        elif args.browser or not show_window(server.url):
            webbrowser.open(server.url)
            wait_for_exit(server.url)
    finally:
        server.stop()
        log.info("stopped")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
