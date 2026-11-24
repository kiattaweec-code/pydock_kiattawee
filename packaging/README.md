# Building the desktop app

pydock ships as a normal desktop application: double-click it and a window
opens with the full interface. Inside, it is the same FastAPI service the
browser version uses, started on a loopback port in the same process and
displayed in a native webview. Nothing listens on an external interface and
no structure is uploaded anywhere.

## The one hard constraint

**PyInstaller cannot cross-compile.** A Windows `.exe` must be built on
Windows, a macOS `.app` on macOS. There is no way around this: the bundle
embeds the platform's Python, its C extensions, and its shared libraries.

So there are two routes to installers for both platforms:

1. **GitHub Actions** — `.github/workflows/build-apps.yml` builds all four
   targets on their own runners. Push a tag (`git tag v0.1.0 && git push
   --tags`) or run it from the Actions tab; it attaches the installers to a
   draft release. This is the route to use.
2. **Build on each machine** — if you have a Windows PC and a Mac to hand,
   run the same one command on each.

## Building locally

```bash
pip install -e . && pip install fastapi uvicorn python-multipart
pip install pyinstaller pillow pywebview

python packaging/build.py --package
```

`build.py` generates the icons, runs PyInstaller, runs the frozen binary's
own self-test, and then produces the right artefact for the OS it is on:

| Built on | Produces |
|---|---|
| Windows | `pydock-0.1.0-windows-x86_64-setup.exe` (Inno Setup) and a portable `.zip` |
| macOS | `pydock-0.1.0-macos-<arch>.dmg` |
| Linux | `pydock-0.1.0-linux-x86_64.tar.gz` |

The Windows installer needs Inno Setup (`winget install JRSoftware.InnoSetup`);
without it the portable zip is still produced.

## Verifying a build

```bash
dist/pydock/pydock --selftest     # paths, versions, and kernel caching
python packaging/ci_check.py      # launch it, serve the UI, dock the example
```

`ci_check.py` is what CI runs before publishing: it starts the built binary,
waits for the server, docks the bundled example end to end and checks the
result. A build that cannot dock never reaches a release.

## Two decisions worth knowing about

**onedir, not onefile.** Numba caches compiled kernels under a key derived
from the source file's path. A onefile build unpacks to a new temporary
directory on each launch, so the key changes every time and the cache never
hits. onedir has a stable path, so the kernels compiled on first run are
reused for good — and the app does not unpack 400 MB at every start. Measured
on the example job: 27.3 s cold, 24.2 s warm, with the server reachable about
a second after launch either way.

**The `.py` sources ship alongside the bytecode.** Numba only caches a
function whose source file it can find on disk. Bundling the dozen modules
that carry `@njit` kernels costs a few hundred kilobytes and is what makes
the caching above work at all. `--selftest` reports whether it is working, so
this does not silently regress.

## Size

| | |
|---|---|
| Download | ~125 MB compressed |
| Installed | ~360 MB |

Most of that is `libllvmlite.so` at 153 MB — the LLVM compiler Numba uses to
turn the scoring kernels into machine code, and not something that can be
trimmed without building llvmlite from source. RDKit is a further ~100 MB.

The build already strips what a docking app never touches: SciPy went when
`pydock.optimize` replaced its one use (−109 MB), RDKit's 2D depiction stack
and the X11/cairo/FreeType chain it drags in (−47 MB), Pillow, `cryptography`,
and uvicorn's optional `uvloop`/`httptools`/`websockets` extras, which buy
nothing on a loopback socket.

## Code signing

The CI builds are **unsigned**, so:

- **Windows** shows a SmartScreen warning. *More info* → *Run anyway*.
- **macOS** refuses a double-click. Right-click the app → *Open*, or
  `xattr -dr com.apple.quarantine /Applications/pydock.app`.

To sign properly:

- **macOS**: set the `MACOS_SIGN_IDENTITY` secret to your Developer ID
  Application identity and `build.py` signs the bundle with a hardened
  runtime. Notarisation needs a further `xcrun notarytool submit` step with
  an App Store Connect key — worth adding if this is distributed widely,
  since without it every user sees the Gatekeeper prompt.
- **Windows**: sign `dist/pydock/pydock.exe` and the finished setup
  executable with `signtool` before upload. An OV certificate still accrues
  SmartScreen reputation slowly; an EV certificate does not.

## Where the app keeps its data

| OS | Path |
|---|---|
| Windows | `%LOCALAPPDATA%\pydock` |
| macOS | `~/Library/Application Support/pydock` |
| Linux | `~/.local/share/pydock` |

Holding `jobs/` (uploaded structures and poses), `cache/numba/` (compiled
kernels) and `pydock.log`. Override the lot with `PYDOCK_HOME=/some/path`.
The Windows uninstaller removes the cache and asks before deleting jobs.

## Command line

The app is also a normal CLI:

```
pydock                      # open the window
pydock --browser            # use the default browser instead
pydock --no-open --port 8000  # headless, e.g. on a shared workstation
pydock --selftest           # diagnose an installation
pydock --data-dir D:\work   # keep jobs somewhere else
```
