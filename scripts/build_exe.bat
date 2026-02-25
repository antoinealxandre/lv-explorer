@echo off
REM =========================================================================
REM  LV Explorer — Windows build script
REM  Produces: dist\LVExplorer\LVExplorer.exe
REM
REM  Prerequisites:
REM    pip install pyinstaller
REM    (All other dependencies must also be installed in the active venv)
REM =========================================================================

SETLOCAL

SET SCRIPT_DIR=%~dp0
SET ROOT=%SCRIPT_DIR%..
SET DIST_DIR=%ROOT%\dist
SET BUILD_DIR=%ROOT%\build

echo.
echo ===================================================
echo   LV Explorer — Building standalone .exe
echo ===================================================
echo.

cd /d "%ROOT%"

REM --- Clean previous build ---
if exist "%BUILD_DIR%" (
    echo [1/4] Cleaning previous build...
    rmdir /s /q "%BUILD_DIR%"
)
if exist "%DIST_DIR%\LVExplorer" (
    rmdir /s /q "%DIST_DIR%\LVExplorer"
)

REM --- Run PyInstaller ---
echo [2/4] Running PyInstaller...
pyinstaller lv_explorer.spec --noconfirm
if ERRORLEVEL 1 (
    echo.
    echo ERROR: PyInstaller failed. See output above.
    exit /b 1
)

REM --- Create ZIP archive ---
echo [3/4] Creating ZIP archive...
SET ZIP_NAME=LVExplorer_v0.1.0_alpha_win64.zip
cd "%DIST_DIR%"
powershell -Command "Compress-Archive -Path 'LVExplorer' -DestinationPath '%ZIP_NAME%' -Force"
cd /d "%ROOT%"

REM --- Done ---
echo [4/4] Done!
echo.
echo Output:
echo   Folder : dist\LVExplorer\LVExplorer.exe
echo   Archive: dist\%ZIP_NAME%
echo.
ENDLOCAL
