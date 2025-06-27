@echo off
setlocal

rem -- absolute paths -------------
set "ROOT=%~dp0"                           rem ...\pumpkin-1.4\
set "PUMPKIN_EXE=%ROOT%src\pumpkin.exe"
set "CASE_DIR=%ROOT%src\Examples\Sim"      rem contains *.qtf + input file

rem -- run PumpKin -----------------

echo -1 | "%PUMPKIN_EXE%" "%CASE_DIR%"  > "%ROOT%pumpkin_output.txt"

if errorlevel 1 (
    echo [ERROR] PumpKin execution failed.
    exit /b 1
)

echo [INFO] PumpKin completed
