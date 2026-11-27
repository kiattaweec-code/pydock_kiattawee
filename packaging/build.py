"""Cross-platform build driver.

    python packaging/build.py            # build the app for this OS
    python packaging/build.py --package  # ...and wrap it in an installer

PyInstaller cannot cross-compile: a Windows .exe has to be built on Windows
and a macOS .app on macOS.  This script is the same on all three so the CI
matrix and a local build run identical steps.
"""
from __future__ import annotations

import argparse
import os
import platform
import shutil
import subprocess
import sys
import time
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
DIST = ROOT / "dist"
VERSION = "0.2.0"

IS_WIN = sys.platform == "win32"
IS_MAC = sys.platform == "darwin"
ARCH = platform.machine().lower().replace("amd64", "x86_64")


def run(cmd: list[str], **kw) -> int:
    print("  $", " ".join(str(c) for c in cmd), flush=True)
    return subprocess.run(cmd, check=True, **kw).returncode


def human(path: Path) -> str:
    if path.is_file():
        n = path.stat().st_size
    else:
        n = sum(f.stat().st_size for f in path.rglob("*") if f.is_file())
    return f"{n / 1048576:.0f} MB"


# --------------------------------------------------------------------------- #
def build_app() -> Path:
    print(f"\n== building pydock {VERSION} for {sys.platform}/{ARCH} ==\n")
    try:
        run([sys.executable, str(HERE / "make_icons.py")])
    except Exception as exc:                                      # noqa: BLE001
        print(f"  (icons not regenerated: {exc})")

    for d in (DIST, ROOT / "build"):
        shutil.rmtree(d, ignore_errors=True)

    t0 = time.time()
    run([sys.executable, "-m", "PyInstaller", str(HERE / "pydock.spec"),
         "--noconfirm", "--distpath", str(DIST), "--workpath", str(ROOT / "build")])

    out = DIST / ("pydock.app" if IS_MAC else "pydock")
    if not out.exists():
        raise SystemExit(f"build produced nothing at {out}")
    print(f"\n  built {out.name}  {human(out)}  in {time.time() - t0:.0f}s")
    return out


def smoke_test(app: Path):
    """Run the frozen binary's own self-test before packaging it."""
    exe = (app / "Contents" / "MacOS" / "pydock") if IS_MAC else \
          (app / ("pydock.exe" if IS_WIN else "pydock"))
    print("\n== self-test ==\n")
    r = subprocess.run([str(exe), "--selftest"], capture_output=True, text=True,
                       timeout=600)
    print(r.stdout or r.stderr)
    if r.returncode == 2:
        print("  WARNING: the kernel cache is not working; the app will be "
              "slow to start every time")
    elif r.returncode != 0:
        raise SystemExit(f"the built app failed its self-test "
                         f"(exit {r.returncode})")


# --------------------------------------------------------------------------- #
FIRST_RUN_WIN = """\
Starting pydock the first time
==============================

This build is not code-signed, so Windows shows a blue
"Windows protected your PC" box the first time.

  1. Click "More info".
  2. Click "Run anyway".

That appears once. After installing, start pydock from the Start menu.

The first launch takes about a minute: pydock compiles its search kernels and
caches them. Later launches start in a second or two.

If you downloaded the ZIP instead of the installer, right-click the ZIP ->
Properties -> tick "Unblock" BEFORE extracting it. Then run pydock.exe from
the extracted folder.

On a Windows-on-ARM PC (Surface with a Snapdragon chip) this x64 build runs
through the built-in emulation. It works; docking is somewhat slower.

Everything runs on your PC. No structure is uploaded anywhere.
"""


def package_windows(app: Path) -> list[Path]:
    made = []
    print("\n== portable zip ==")
    (app / "READ ME FIRST.txt").write_text(FIRST_RUN_WIN, encoding="utf-8")
    portable = DIST / f"pydock-{VERSION}-windows-{ARCH}-portable.zip"
    shutil.make_archive(str(portable.with_suffix("")), "zip", DIST, "pydock")
    made.append(portable)

    iscc = shutil.which("iscc") or shutil.which("ISCC")
    if not iscc:
        for guess in (r"C:\Program Files (x86)\Inno Setup 6\ISCC.exe",
                      r"C:\Program Files\Inno Setup 6\ISCC.exe"):
            if Path(guess).exists():
                iscc = guess
                break
    if iscc:
        print("\n== installer (Inno Setup) ==")
        run([iscc, f"/DMyAppVersion={VERSION}", str(HERE / "windows_installer.iss")])
        for f in DIST.glob("pydock-*-setup.exe"):
            made.append(f)
    else:
        print("\n  Inno Setup not found; skipping the installer.")
        print("  Install it with:  winget install JRSoftware.InnoSetup")
    return made


FIRST_RUN_MAC = """\
Opening pydock the first time
=============================

Step 1 — drag pydock onto the Applications folder in this window.

Step 2 — this build is not notarised by Apple, so the first launch is
blocked. Pick whichever you prefer; both are one-time.

  FASTEST (one command)
    Open Terminal and paste:

        xattr -dr com.apple.quarantine /Applications/pydock.app

    Then open pydock normally.

  NO TERMINAL
    1. Double-click pydock. macOS refuses it. Click Done.
    2. Open System Settings > Privacy & Security.
    3. Scroll to Security. There is a line about pydock being blocked.
       Click "Open Anyway", then confirm with your password.
    4. Open pydock again and click Open.

    Note: Control-clicking the app and choosing Open no longer works on
    macOS Sequoia (15) and later. Apple removed that shortcut, so use one
    of the two routes above.

After this, pydock opens with a normal double-click.

The first launch takes about a minute: pydock compiles its search kernels
and caches them. Later launches start in a second or two.

Everything runs on your Mac. No structure is uploaded anywhere.
"""


def sign_macos(app: Path) -> str:
    """Sign the bundle, and never leave it unsigned.

    An ad-hoc signature is not optional on Apple silicon: the kernel refuses
    to run an arm64 bundle that has no signature at all, and the user is told
    the app is "damaged" -- which right-click > Open does not fix.  So when
    there is no Developer ID we still sign with `-` (ad-hoc), which makes the
    app launchable after the usual Gatekeeper prompt.
    """
    identity = os.environ.get("MACOS_SIGN_IDENTITY")
    if identity:
        print(f"  signing with: {identity}")
        run(["codesign", "--deep", "--force", "--options", "runtime",
             "--timestamp", "--sign", identity, str(app)])
        kind = "Developer ID"
    else:
        print("  MACOS_SIGN_IDENTITY not set: ad-hoc signing")
        run(["codesign", "--deep", "--force", "--sign", "-", str(app)])
        kind = "ad-hoc"

    # Verify, because a broken signature fails at the user's end, not here.
    check = subprocess.run(["codesign", "--verify", "--deep", "--strict",
                            str(app)], capture_output=True, text=True)
    if check.returncode != 0:
        raise SystemExit("the signature did not verify:\n"
                         + (check.stderr or check.stdout))
    print(f"  signature verified ({kind})")
    return kind


def _layout_dmg(volume: str):
    """Size the window and place the icons, so the drag target is obvious."""
    background = "" if not (HERE / "dmg_background.png").exists() else f'''
            set background picture of theViewOptions to file ".background:bg.png"'''
    script = f'''
    tell application "Finder"
      tell disk "{volume}"
        open
        set theXY to {{200, 120}}
        set current view of container window to icon view
        set toolbar visible of container window to false
        set statusbar visible of container window to false
        set bounds of container window to {{200, 120, 940, 580}}
        set theViewOptions to the icon view options of container window
        set arrangement of theViewOptions to not arranged
        set icon size of theViewOptions to 128{background}
        set position of item "pydock.app" of container window to {{180, 210}}
        set position of item "Applications" of container window to {{560, 210}}
        set position of item "READ ME FIRST.txt" of container window to {{370, 380}}
        close
        open
        update without registering applications
        delay 2
      end tell
    end tell
    '''
    r = subprocess.run(["osascript", "-e", script], capture_output=True, text=True)
    if r.returncode != 0:
        # Finder is not always available on a CI runner; the image still works,
        # it just opens with the default window.
        print(f"  (could not lay out the window: "
              f"{(r.stderr or '').strip()[:120]})")


def package_macos(app: Path) -> list[Path]:
    print("\n== disk image ==")
    sign_macos(app)

    dmg = DIST / f"pydock-{VERSION}-macos-{ARCH}.dmg"
    tmp_dmg = DIST / "pydock-tmp.dmg"
    for f in (dmg, tmp_dmg):
        f.unlink(missing_ok=True)

    staging = DIST / "dmg"
    shutil.rmtree(staging, ignore_errors=True)
    staging.mkdir(parents=True)
    # ditto, not copytree: it is the only copy that reliably preserves the
    # signature and the bundle's extended attributes.
    run(["ditto", str(app), str(staging / app.name)])
    os.symlink("/Applications", staging / "Applications")
    (staging / "READ ME FIRST.txt").write_text(FIRST_RUN_MAC)

    bg = HERE / "dmg_background.png"
    if bg.exists():
        (staging / ".background").mkdir(exist_ok=True)
        shutil.copy(bg, staging / ".background" / "bg.png")

    volume = f"pydock {VERSION}"
    if not _pretty_dmg(staging, tmp_dmg, dmg, volume):
        # The laid-out route failed; build the image in one shot instead.
        # Plainer window, same contents, and it cannot get stuck on a mount.
        print("  falling back to a plain disk image")
        dmg.unlink(missing_ok=True)
        run(["hdiutil", "create", "-volname", volume, "-srcfolder",
             str(staging), "-ov", "-format", "UDZO", str(dmg)])

    tmp_dmg.unlink(missing_ok=True)
    shutil.rmtree(staging, ignore_errors=True)
    return [dmg]


def _detach(device: str | None, volume: str) -> bool:
    """Detach the image and wait for the kernel to let go of it.

    `hdiutil convert` fails with "Resource temporarily unavailable" if the
    read-write image is still attached, and detach does not always take
    effect by the time it returns -- so this retries and then confirms the
    mount point is really gone.
    """
    target = device or f"/Volumes/{volume}"
    for attempt in range(6):
        subprocess.run(["hdiutil", "detach", target, "-force"],
                       capture_output=True, text=True)
        time.sleep(1.5)
        if not Path(f"/Volumes/{volume}").exists():
            return True
    return not Path(f"/Volumes/{volume}").exists()


def _pretty_dmg(staging: Path, tmp_dmg: Path, dmg: Path, volume: str) -> bool:
    """Build a disk image with the icons arranged.  False if that cannot be done.

    Needs a usable Finder, so it works on a desktop Mac and generally not on
    a CI runner, where the caller falls back to a plain image.
    """
    device = None
    try:
        run(["hdiutil", "create", "-volname", volume, "-srcfolder",
             str(staging), "-ov", "-fs", "HFS+", "-format", "UDRW",
             str(tmp_dmg)])
        mount = subprocess.run(
            ["hdiutil", "attach", str(tmp_dmg), "-readwrite", "-noverify",
             "-noautoopen"], capture_output=True, text=True)
        if mount.returncode != 0:
            return False
        for line in mount.stdout.splitlines():
            if line.strip().startswith("/dev/disk"):
                device = line.split()[0]
                break
        _layout_dmg(volume)
    except Exception as exc:                                      # noqa: BLE001
        print(f"  ({exc})")
        _detach(device, volume)
        return False

    if not _detach(device, volume):
        print("  the disk image stayed mounted")
        return False

    for attempt in range(4):
        r = subprocess.run(["hdiutil", "convert", str(tmp_dmg), "-format",
                            "UDZO", "-imagekey", "zlib-level=9", "-ov",
                            "-o", str(dmg)], capture_output=True, text=True)
        if r.returncode == 0:
            return True
        print(f"  convert attempt {attempt + 1} failed: "
              f"{(r.stderr or r.stdout).strip().splitlines()[-1][:90]}")
        time.sleep(4)
    return False


def package_linux(app: Path) -> list[Path]:
    print("\n== tarball ==")
    tar = DIST / f"pydock-{VERSION}-linux-{ARCH}.tar.gz"
    shutil.make_archive(str(tar).removesuffix(".tar.gz"), "gztar", DIST, "pydock")
    return [tar]


# --------------------------------------------------------------------------- #
def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--package", action="store_true",
                    help="also produce the installer/disk image for this OS")
    ap.add_argument("--skip-test", action="store_true")
    args = ap.parse_args()

    app = build_app()
    if not args.skip_test:
        smoke_test(app)

    if args.package:
        made = (package_windows(app) if IS_WIN else
                package_macos(app) if IS_MAC else
                package_linux(app))
        print("\n== artefacts ==")
        for f in made:
            print(f"  {f.name:<48} {human(f)}")
    print("\ndone.\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
