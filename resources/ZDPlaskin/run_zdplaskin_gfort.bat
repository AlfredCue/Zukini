@echo off
setlocal enabledelayedexpansion

REM ============================
REM CONFIGURACIÓN GENERAL
REM ============================
set PREPROCESSOR=preprocessor.exe
set KINET_FILE=kinet.inp
set SOURCE_FILES=dvode_f90_m.F90 zdplaskin_m.F90 dbd0D.F90
set LIBRARY_FILE=bolsig_g.dll
set RESULTS_DIR=QT_files_simulation
set EXECUTABLE=dvode_f90_m.exe

set TOTAL_STEPS=4

echo ================================================================
echo     ZDPlasKin Simulation Workflow [GFORTRAN]
echo ================================================================

REM ============================
REM STEP 1: Run Preprocessor
REM ============================
echo [1/%TOTAL_STEPS%] Running preprocessor on %KINET_FILE%...
if not exist %PREPROCESSOR% (
    echo ERROR: Preprocessor not found: %PREPROCESSOR%
    exit /b 1
)
%PREPROCESSOR% %KINET_FILE% zdplaskin_m.F90 > preprocessor_output.txt 2>&1
if errorlevel 1 (
    echo ERROR: Preprocessor failed. See preprocessor_output.txt
    exit /b 1
)
echo [1/%TOTAL_STEPS%] Preprocessor completed successfully.

REM ============================
REM STEP 2: Compile with GFortran
REM ============================
echo [2/%TOTAL_STEPS%] Compiling Fortran sources with GFortran...
gfortran -O3 -march=native -funroll-loops -fopenmp -fno-second-underscore -ffree-line-length-none -Wall dvode_f90_m.F90 zdplaskin_m.F90 dbd0D.F90 bolsig_x86_64_g.dll -o %EXECUTABLE% > compile_output.txt 2>&1
if errorlevel 1 (
    echo ERROR: Compilation failed. See compile_output.txt
    exit /b 2
)
echo [2/%TOTAL_STEPS%] Compilation successful.

REM ============================
REM STEP 3: Run Simulation
REM ============================
echo [3/%TOTAL_STEPS%] Running simulation...
if not exist %EXECUTABLE% (
    echo ERROR: Executable not found: %EXECUTABLE%
    exit /b 3
)

REM --- Run via PowerShell so we can “tee” stdout+stderr --------------
REM     ▸ &  : calls the EXE
REM     ▸ 2>&1: merges stderr into stdout
REM     ▸ Tee-Object writes to the file *and* passes through to console
REM     ▸ exit $LASTEXITCODE: bubbles the EXE’s real exit-code back to CMD
powershell -NoLogo -NoProfile -Command ^
  "& { & '%CD%\%EXECUTABLE%' 2>&1 | Tee-Object -Variable out | Out-File -FilePath 'simulation_output.txt' -Encoding UTF8; $out }"
if errorlevel 1 (
    echo ERROR: Simulation failed. See simulation_output.txt
    exit /b 3
)
echo [3/%TOTAL_STEPS%] Simulation completed successfully.

REM ============================
REM STEP 4: Organize results
REM ============================
echo [4/%TOTAL_STEPS%] Collecting results into %RESULTS_DIR%...
mkdir "%RESULTS_DIR%" >nul 2>&1
move qt_*.txt "%RESULTS_DIR%" >nul 2>&1
echo [4/%TOTAL_STEPS%] Results moved successfully.

echo ================================================================
echo       ZDPlasKin workflow finished without errors.
echo       Results directory: %RESULTS_DIR%
echo ================================================================
