@echo off
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

REM Set full paths for source and destination
set "SOURCE=%PROJECT_ROOT%\output\kinet.inp"
set "DESTINATION=%PROJECT_ROOT%\resources\ZDPlaskin"

REM === Check if source file exists ===
if not exist "%SOURCE%" (
    echo Error: Source file not found: %SOURCE%
    exit /b 1
)

REM === Move the file ===
move /Y "%SOURCE%" "%DESTINATION%\kinet.inp" >nul
if errorlevel 1 (
    echo Error: Failed to move file to %DESTINATION%
    exit /b 2
)

echo kinet.inp successfully moved to %DESTINATION%
