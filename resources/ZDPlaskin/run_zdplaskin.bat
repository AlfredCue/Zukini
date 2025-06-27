@echo off
setlocal enabledelayedexpansion

REM ============================
REM CONFIGURACIÓN GENERAL
REM ============================
set PREPROCESSOR=preprocessor.exe
set KINET_FILE=kinet.inp
set SOURCE_FILES=dvode_f90_m.F90 zdplaskin_m.F90 dbd0D.F90
set LIBRARY_FILE=bolsig_x86_64.lib
set RESULTS_DIR=QT_files_simulation
set EXECUTABLE=zdplaskin_sim.exe

set TOTAL_STEPS=4

echo ================================================================
echo     ZDPlasKin Simulation Workflow
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
REM STEP 2: Compile Fortran code
REM ============================
echo [2/%TOTAL_STEPS%] Compiling Fortran sources with Intel Fortran...
ifort /O3 /arch:CORE_AVX2 /Qipo /Qfma %SOURCE_FILES% %LIBRARY_FILE% /Fe:%EXECUTABLE% > compile_output.txt 2>&1
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
%EXECUTABLE% 2>&1 | tee simulation_output.txt
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
