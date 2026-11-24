@echo off
REM
REM  Double-click this file to build the Windows app.
REM
REM  It sets up an isolated Python environment, installs what pydock needs,
REM  builds pydock.exe, checks that the built app really docks, and wraps it
REM  in an installer you can send to people.
REM
REM  Nothing is installed system-wide: everything lands in .venv next to this
REM  file, and deleting the folder removes all of it.

setlocal EnableDelayedExpansion
cd /d "%~dp0"

REM Captured before any parenthesised block: the ")" inside the variable
REM name would otherwise close the block early and break the script.
set "PF86=%ProgramFiles(x86)%"
set "PF64=%ProgramFiles%"
title Building pydock for Windows
cls

echo.
echo  Building pydock for Windows
echo  ---------------------------
echo.

REM --- 1. find a usable Python ----------------------------------------------
echo  ==^> Looking for Python 3.10 or newer
set "PY="
for %%V in (3.13 3.12 3.11 3.10) do (
    if not defined PY (
        py -%%V -c "import sys" >nul 2>&1 && set "PY=py -%%V"
    )
)
if not defined PY (
    python -c "import sys; raise SystemExit(0 if (3,10) <= sys.version_info < (3,14) else 1)" >nul 2>&1 && set "PY=python"
)

if not defined PY (
    echo.
    echo  No suitable Python found.
    echo.
    echo  Install it, then run this again:
    echo.
    echo      winget install Python.Python.3.12
    echo.
    echo  Or download it from https://www.python.org/downloads/windows/
    echo  and tick "Add python.exe to PATH" during setup.
    echo.
    echo  After installing, CLOSE this window and open a new one, so that
    echo  Windows picks up the new PATH.
    echo.
    pause
    exit /b 1
)
for /f "delims=" %%O in ('%PY% -V 2^>^&1') do echo     using %%O

REM 64-bit matters: RDKit and Numba publish no 32-bit wheels.
%PY% -c "import sys; raise SystemExit(0 if sys.maxsize > 2**32 else 1)" >nul 2>&1
if errorlevel 1 (
    echo.
    echo  That Python is 32-bit. RDKit and Numba need a 64-bit Python.
    echo  Install the 64-bit build and run this again.
    echo.
    pause
    exit /b 1
)

REM --- 2. environment --------------------------------------------------------
echo.
echo  ==^> Setting up the build environment  ^(a few minutes the first time^)
if not exist ".venv\Scripts\python.exe" (
    %PY% -m venv .venv
    if errorlevel 1 goto :failed_env
)
set "VPY=.venv\Scripts\python.exe"

"%VPY%" -m pip install --upgrade pip --quiet
if errorlevel 1 goto :failed_env

echo     installing pydock and its build tools...
"%VPY%" -m pip install -e ".[app,build]" --quiet
if errorlevel 1 goto :failed_deps

REM --- 3. Inno Setup (optional, makes the .exe installer) --------------------
echo.
echo  ==^> Checking for Inno Setup
set "HAVE_ISCC="
where iscc >nul 2>&1 && set "HAVE_ISCC=1"
if not defined HAVE_ISCC if exist "%PF86%\Inno Setup 6\ISCC.exe" set "HAVE_ISCC=1"
if not defined HAVE_ISCC if exist "%PF64%\Inno Setup 6\ISCC.exe" set "HAVE_ISCC=1"

if not defined HAVE_ISCC (
    echo     not installed - trying to install it with winget
    winget install --id JRSoftware.InnoSetup --silent --accept-package-agreements --accept-source-agreements >nul 2>&1
    if exist "!PF86!\Inno Setup 6\ISCC.exe" (
        set "HAVE_ISCC=1"
    ) else (
        echo     could not install it - you will still get the portable ZIP
    )
) else (
    echo     found
)

REM --- 4. build --------------------------------------------------------------
echo.
echo  ==^> Building the app  ^(10-15 minutes^)
echo.
"%VPY%" packaging\build.py --package
if errorlevel 1 goto :failed_build

REM --- 5. done ---------------------------------------------------------------
echo.
echo  ==^> Done
echo.
for %%F in (dist\pydock-*-setup.exe) do echo     %%~nxF   ^(%%~zF bytes^)
for %%F in (dist\pydock-*-portable.zip) do echo     %%~nxF   ^(%%~zF bytes^)
echo.
echo     Send the setup .exe to anyone with a 64-bit Windows PC.
echo     ^(It also runs on Windows-on-ARM through the built-in emulation.^)
echo.
start "" "dist"
pause
exit /b 0

:failed_env
echo.
echo  Could not create the Python environment. Scroll up for the reason.
pause
exit /b 1

:failed_deps
echo.
echo  Installing the dependencies failed. Scroll up for the reason.
echo  The usual cause is a Python version outside 3.10-3.13, or no network.
pause
exit /b 1

:failed_build
echo.
echo  The build failed. Scroll up for the reason.
pause
exit /b 1
