# PyInstaller spec, shared by Windows, macOS and Linux.
#
#     pyinstaller packaging/pydock.spec --noconfirm
#
# Two choices here are deliberate and worth keeping.
#
# *onedir, not onefile.*  Numba caches compiled kernels under a key derived
# from the source file's path.  A onefile build unpacks to a fresh temporary
# directory on every launch, so that key changes every time and the cache
# never hits -- the app would spend a minute recompiling at each start.  A
# onedir build has a stable path, so the cache written on first run is reused
# for good.  It also avoids unpacking ~300 MB on every launch.
#
# *The .py sources ship too.*  Numba will only cache a function whose source
# file it can find on disk.  Bundling the handful of modules that carry
# @njit kernels is a few hundred kilobytes and is what makes the above work.

import sys
from pathlib import Path

from PyInstaller.utils.hooks import collect_all, collect_submodules

ROOT = Path(SPECPATH).resolve().parent
IS_WIN = sys.platform == "win32"
IS_MAC = sys.platform == "darwin"

datas = [
    (str(ROOT / "pydock_web" / "static"), "pydock_web/static"),
    (str(ROOT / "test_data" / "receptor_noH.pdb"), "test_data"),
    (str(ROOT / "test_data" / "ligand_start.sdf"), "test_data"),
    (str(ROOT / "test_data" / "reference.sdf"), "test_data"),
    (str(ROOT / "test_data" / "library.sdf"), "test_data"),
]

# Source files for the Numba-compiled modules, so caching works (see above).
for mod in ("geometry", "grid", "scoring", "cavity", "ga", "params",
            "ligand", "receptor", "dock", "optimize", "__init__"):
    src = ROOT / "pydock" / f"{mod}.py"
    if src.exists():
        datas.append((str(src), "pydock"))

binaries = []
hiddenimports = [
    "pydock", "pydock_web", "pydock_web.app", "pydock_web.jobs",
    "pydock_web.prep", "pydock.structure", "pydock.mol2",
    "pydock_app", "pydock_app.launcher", "pydock_app.paths",
    "uvicorn.logging", "uvicorn.loops.asyncio",
    "uvicorn.protocols.http.h11_impl", "uvicorn.lifespan.on",
    "multipart", "python_multipart", "email.mime.multipart",
]
hiddenimports += collect_submodules("pydock")

# RDKit, Numba and llvmlite all load things dynamically; let the hooks find them.
for pkg in ("rdkit", "numba", "llvmlite"):
    d, b, h = collect_all(pkg)
    datas += d
    binaries += b
    # collect_all also reports every module that `excludes` drops -- the
    # depiction stack, and numba's whole CUDA tree.  Left in, each one logs a
    # "hidden import not found" error, and 200 of those bury a real problem.
    hiddenimports += [
        m for m in h
        if not m.startswith(("rdkit.Chem.Draw", "numba.cuda", "numba.roc"))
        and m != "rdkit.Chem.inchi"
    ]

# pywebview gives the app a real window instead of opening a browser tab.
# It is optional on purpose: it pulls in platform GUI bindings (pyobjc on
# macOS, pythonnet on Windows) that can be absent or refuse to bundle, and
# the launcher falls back to the default browser by itself.  So include it
# when the build machine has it, and carry on when it does not.
import importlib.util

if importlib.util.find_spec("webview") is not None:
    try:
        d, b, h = collect_all("webview")
        datas += d
        binaries += b
        hiddenimports += h
        hiddenimports += ["webview.platforms.cocoa"] if IS_MAC else (
            ["webview.platforms.winforms", "webview.platforms.edgechromium",
             "clr", "clr_loader"] if IS_WIN else ["webview.platforms.gtk"])
        print("[pydock] bundling pywebview: the app opens in its own window")
    except Exception as exc:                                      # noqa: BLE001
        print(f"[pydock] pywebview not bundled ({exc}); "
              f"the app will open in the default browser")
else:
    print("[pydock] pywebview is not installed; "
          "the app will open in the default browser")

# Trim what a docking app never touches.  SciPy is gone from the codebase
# entirely (pydock.optimize replaced its one use), and these would otherwise
# ride along through transitive imports.
excludes = [
    "scipy", "matplotlib", "tkinter", "PyQt5", "PyQt6", "PySide2", "PySide6",
    "IPython", "jupyter", "notebook", "pytest", "sphinx", "setuptools",
    "pip", "wheel", "pandas", "sqlalchemy", "test", "tests",
    "numba.cuda", "numba.roc",
    # RDKit's 2D depiction stack pulls in cairo, FreeType and the whole X11
    # client chain.  pydock never draws a molecule -- the browser does that --
    # so none of it is reachable.
    "rdkit.Chem.Draw", "rdkit.Chem.Draw.IPythonConsole",
    "rdkit.Chem.Draw.rdMolDraw2D", "rdkit.Chem.inchi",
    "PIL", "Pillow",
    # Only the icon generator uses Pillow, and that runs at build time.
    # Nothing here speaks TLS or WebSocket: the server is loopback-only and
    # the uvicorn Config below pins the plain asyncio/h11 stack.
    "cryptography", "websockets", "uvloop", "httptools", "watchfiles",
    "uvicorn.protocols.websockets", "uvicorn.loops.uvloop",
]

block_cipher = None

a = Analysis(
    [str(ROOT / "packaging" / "entry.py")],
    pathex=[str(ROOT)],
    binaries=binaries,
    datas=datas,
    hiddenimports=sorted(set(hiddenimports)),
    hookspath=[str(ROOT / "packaging" / "hooks")],
    runtime_hooks=[str(ROOT / "packaging" / "rthook_pydock.py")],
    excludes=excludes,
    cipher=block_cipher,
    noarchive=False,
)

# Drop the GUI shared libraries RDKit declares but pydock never loads.  They
# arrive as dependencies of the depiction extensions, which are excluded
# above; without the extension, nothing ever resolves against them.
_GUI_LIBS = ("libx11", "libxcb", "libxau", "libxdmcp", "libxext", "libxrender",
             "libcairo", "libpango", "libgdk", "libgtk", "libatk",
             "libfreetype", "libfontconfig", "libharfbuzz", "libpixman",
             "libepoxy", "libwayland", "libgio-", "libgobject", "libgthread")


def _keep(entry) -> bool:
    name = Path(entry[0]).name.lower()
    return not name.startswith(_GUI_LIBS)


_before = len(a.binaries)
a.binaries = TOC([e for e in a.binaries if _keep(e)])
if _before != len(a.binaries):
    print(f"[pydock] dropped {_before - len(a.binaries)} GUI-only shared libraries")

pyz = PYZ(a.pure, a.zipped_data, cipher=block_cipher)

exe = EXE(
    pyz,
    a.scripts,
    [],
    exclude_binaries=True,
    name="pydock",
    debug=False,
    bootloader_ignore_signals=False,
    # Stripping symbols takes a large bite out of libllvmlite (170 MB
    # unstripped).  Windows binaries are not stripped this way.
    strip=not IS_WIN,
    upx=False,                       # UPX corrupts some RDKit/LLVM libraries
    console=not (IS_WIN or IS_MAC),  # no console window on desktop platforms
    disable_windowed_traceback=False,
    icon=str(ROOT / "packaging" / "icon.ico") if IS_WIN else
         (str(ROOT / "packaging" / "icon.icns") if IS_MAC else None),
    # Gives the .exe a real publisher and description in Properties and in
    # the SmartScreen prompt, instead of blank fields.
    version=str(ROOT / "packaging" / "version_info.txt") if IS_WIN else None,
)

coll = COLLECT(
    exe,
    a.binaries,
    a.datas,
    strip=not IS_WIN,
    upx=False,
    name="pydock",
)

if IS_MAC:
    app = BUNDLE(
        coll,
        name="pydock.app",
        icon=str(ROOT / "packaging" / "icon.icns"),
        bundle_identifier="org.pydock.app",
        info_plist={
            "CFBundleName": "pydock",
            "CFBundleDisplayName": "pydock",
            "CFBundleShortVersionString": "0.1.1",
            "CFBundleVersion": "0.1.1",
            "NSHighResolutionCapable": True,
            # The app talks only to 127.0.0.1, which still counts as a local
            # network connection on recent macOS.
            "NSAppTransportSecurity": {"NSAllowsLocalNetworking": True},
            "LSMinimumSystemVersion": "11.0",
            "NSHumanReadableCopyright": "pydock",
        },
    )
