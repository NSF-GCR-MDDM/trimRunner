@echo off
setlocal EnableDelayedExpansion

call "C:\root_v6.34.08\bin\thisroot.bat"

set MAX_JOBS=6
set JOBS_RUNNING=0

call :run python trimRunner.py 1H 500000
call :run python trimRunner.py 2H 300000
call :run python trimRunner.py 3H 250000
call :run python trimRunner.py 4He 200000
call :run python trimRunner.py 6He 100000
call :run python trimRunner.py 6Li 100000
call :run python trimRunner.py 7Li 100000
call :run python trimRunner.py 8Li 100000
call :run python trimRunner.py 16N 75000
call :run python trimRunner.py 18O 75000
call :run python trimRunner.py 19F 75000
call :run python trimRunner.py 19O 75000
call :run python trimRunner.py 20F 75000

REM Wait for all remaining jobs
:wait_all
if %JOBS_RUNNING% GTR 0 (
    timeout /t 1 >nul
    goto wait_all
)

echo All ion jobs finished.
pause
exit /b

:run
:wait_slot
if %JOBS_RUNNING% GEQ %MAX_JOBS% (
    timeout /t 1 >nul
    goto wait_slot
)

set /a JOBS_RUNNING+=1
start "" /b cmd /c "%* & set /a JOBS_RUNNING-=1"
exit /b
