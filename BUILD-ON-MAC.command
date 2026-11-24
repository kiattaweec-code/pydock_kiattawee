#!/bin/bash
#
#  Double-click this file in Finder to build the macOS app.
#
#  It sets up an isolated Python environment, installs what pydock needs,
#  builds pydock.app, checks that the built app really docks, and wraps it
#  in a .dmg you can send to people.
#
#  Nothing is installed system-wide: everything lands in .venv next to this
#  file, and deleting the folder removes all of it.

set -u
cd "$(dirname "$0")" || exit 1

BOLD=$'\033[1m'; DIM=$'\033[2m'; RED=$'\033[31m'; GRN=$'\033[32m'; OFF=$'\033[0m'
say()  { printf '%s\n' "$*"; }
step() { printf '\n%s==> %s%s\n' "$BOLD" "$*" "$OFF"; }
die()  { printf '\n%sX  %s%s\n\n' "$RED" "$*" "$OFF"; say "Press return to close."; read -r _; exit 1; }

clear
say "${BOLD}Building pydock for macOS${OFF}"
say "${DIM}$(sw_vers -productName 2>/dev/null) $(sw_vers -productVersion 2>/dev/null) · $(uname -m)${OFF}"

# --- 1. find a usable Python -------------------------------------------------
step "Looking for Python 3.10 or newer"
PY=""
for c in python3.13 python3.12 python3.11 python3.10 python3; do
    p=$(command -v "$c" 2>/dev/null) || continue
    "$p" -c 'import sys; raise SystemExit(0 if (3,10) <= sys.version_info < (3,14) else 1)' 2>/dev/null \
        && { PY="$p"; break; }
done

if [ -z "$PY" ]; then
    say ""
    say "${RED}No suitable Python found.${OFF}"
    say ""
    say "Install one, then run this again:"
    say ""
    say "  ${BOLD}brew install python@3.12${OFF}"
    say ""
    say "or download the installer from"
    say "  https://www.python.org/downloads/macos/"
    say ""
    say "(The Python that ships with macOS cannot install packages, so it"
    say " will not work for this.)"
    say ""
    say "Press return to close."; read -r _; exit 1
fi
say "  using $PY  ($($PY -V 2>&1))"

# Xcode command line tools provide codesign and hdiutil.
if ! command -v codesign >/dev/null 2>&1; then
    step "Installing the Xcode command line tools"
    say "  A system dialog will appear. Finish it, then run this script again."
    xcode-select --install 2>/dev/null
    say ""; say "Press return to close."; read -r _; exit 1
fi

# --- 2. environment ----------------------------------------------------------
step "Setting up the build environment  (a few minutes the first time)"
if [ ! -d .venv ]; then
    "$PY" -m venv .venv || die "could not create the virtual environment"
fi
# shellcheck disable=SC1091
source .venv/bin/activate || die "could not activate the virtual environment"

python -m pip install --upgrade pip --quiet || die "could not update pip"
say "  installing pydock and its build tools…"
python -m pip install -e ".[app,build]" --quiet \
    || die "dependency installation failed — scroll up for the reason"

# --- 3. build ----------------------------------------------------------------
step "Building the app  (10-15 minutes)"
python packaging/build.py --package || die "the build failed — scroll up for the reason"

# --- 4. done -----------------------------------------------------------------
DMG=$(ls -t dist/*.dmg 2>/dev/null | head -1)
if [ -z "$DMG" ]; then
    die "the build finished but produced no .dmg"
fi

step "${GRN}Done${OFF}"
say ""
say "  ${BOLD}$(basename "$DMG")${OFF}  ($(du -h "$DMG" | cut -f1))"
say ""
say "  Send that file to anyone with a Mac of this kind ($(uname -m))."
say "  Opening the folder for you now."
say ""
open dist 2>/dev/null

say "Press return to close."
read -r _
