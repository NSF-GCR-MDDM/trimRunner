@echo off
REM === ROOT + MSVC 64-bit Build Script ===

REM Step 0: Set up environment
call "C:\root_v6.34.08\bin\thisroot.bat"

REM Step 1: Clean previous build
echo Cleaning previous build...
rmdir /s /q build 2>nul
mkdir build
cd build

REM Step 2: Configure with CMake using NMake
echo Running CMake configuration...
cmake .. -G "NMake Makefiles" -DCMAKE_BUILD_TYPE=Release

IF %ERRORLEVEL% NEQ 0 (
    echo CMake configuration failed. Exiting.
    exit /b %ERRORLEVEL%
)

REM Step 3: Build using NMake
echo Building project...
nmake

IF %ERRORLEVEL% NEQ 0 (
    echo Build failed. Exiting.
    exit /b %ERRORLEVEL%
)

REM Step 4: Done
echo Build complete.
cd ..
pause
