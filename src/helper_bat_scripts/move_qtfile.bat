@echo off
rem ------------------------------------------------------------
rem  Move latest QtPlasKin output to Qt_1 or Qt_2 automatically
rem ------------------------------------------------------------
setlocal enabledelayedexpansion

REM ----------------------------------------
REM Define base project path relative to this script
REM ----------------------------------------

REM Get the directory where this script is located
set "SCRIPT_DIR=%~dp0"
REM Go up two levels to reach Zukini/
pushd "%SCRIPT_DIR%\..\.."
set "PROJECT_ROOT=%CD%"
popd

rem -----------------------------------------------------------------
rem Paths (edit if your layout is different)
rem -----------------------------------------------------------------
set "SOURCE=%PROJECT_ROOT%\resources\ZDPlaskin\QT_files_simulation"
set "DEST1=%PROJECT_ROOT%\data\QtFiles\Qt_1"
set "DEST2=%PROJECT_ROOT%\data\QtFiles\Qt_2"

rem Default destination is Qt_1 …
set "DESTINATION=%DEST1%"

rem … unless Qt_1 already has *any* files:
dir /a-d /b "%DEST1%" >nul 2>&1 && (
    echo Qt_1 is not empty - using Qt_2 instead.
    set "DESTINATION=%DEST2%"
)

rem -----------------------------------------------------------------
rem Make sure the destination folder exists
rem -----------------------------------------------------------------
if not exist "%DESTINATION%" (
    mkdir "%DESTINATION%"
)

rem -----------------------------------------------------------------
rem Move qt_*.txt (only if at least one file exists)
rem -----------------------------------------------------------------
echo Moving files from "%SOURCE%" to "%DESTINATION%" …
set "MOVED_ANY=0"
for %%F in ("%SOURCE%\qt_*.txt") do (
    if exist "%%~F" (
        move /Y "%%~F" "%DESTINATION%" >nul
        set "MOVED_ANY=1"
    )
)

if "!MOVED_ANY!"=="0" (
    echo [WARNING] No qt_*.txt files found - nothing moved.
) else (
    echo Done.
)

rem -----------------------------------------------------------------
rem Copy kinet.inp beside the moved Qt files (optional)
rem -----------------------------------------------------------------
rem copy /Y "ZDPlaskin\kinet.inp" "%DESTINATION%\kinet.inp" nul
rem echo kinet.inp copied to "%DESTINATION%".

exit /b