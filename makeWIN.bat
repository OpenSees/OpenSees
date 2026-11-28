@echo off
REM
REM Build OpenSees, OpenSeesPy (and OpenSeesSP/OpenSeesMP if MUMPS given) on Windows.
REM
REM Usage: makeWIN.bat [--mumps <mumps build dir>] [--prefix <install dir>]
REM                    [--jobs <n>] [--no-clean]
REM
REM Requires: Visual Studio (Desktop development with C++), Intel oneAPI Base &
REM   HPC toolkits (MKL, MPI, ifx), Python, and conan >= 2 (pip install conan).
REM Run from a Command Prompt (not PowerShell).

setlocal EnableDelayedExpansion

REM always run from the OpenSees source folder (where this script lives)
cd /d "%~dp0"

set "BUILD_DIR=build"
set "OUTPUT_DIR=%USERPROFILE%\bin\OpenSeesLatest"
set "MUMPS_DIR="
set "JOBS=%NUMBER_OF_PROCESSORS%"
set "CLEAN=1"

:parse
if "%~1"=="" goto endparse
if /i "%~1"=="--mumps" (
    if not exist "%~2\" ( echo FAIL: MUMPS folder "%~2" does not exist & exit /b 1 )
    set "MUMPS_DIR=%~f2"
    shift & shift & goto parse
)
if /i "%~1"=="--prefix" (
    set "OUTPUT_DIR=%~f2"
    shift & shift & goto parse
)
if /i "%~1"=="--jobs" (
    set "JOBS=%~2"
    shift & shift & goto parse
)
if /i "%~1"=="--no-clean" (
    set "CLEAN=0"
    shift & goto parse
)
if /i "%~1"=="-h" goto usage
if /i "%~1"=="--help" goto usage
echo Unknown argument: %~1
goto usage
:endparse

REM set up the oneAPI (and Visual Studio) environment, once per Command Prompt
if not defined SETVARS_COMPLETED (
    if not defined ONEAPI_ROOT set "ONEAPI_ROOT=%ProgramFiles(x86)%\Intel\oneAPI"
    call "!ONEAPI_ROOT!\setvars.bat" intel64
    if errorlevel 1 ( echo FAIL: oneAPI setvars.bat failed & exit /b 1 )
)

REM check tools
where conan >nul 2>&1 || ( echo FAIL: conan not found - pip install conan & exit /b 1 )
where cmake >nul 2>&1 || ( echo FAIL: cmake not found & exit /b 1 )
where ninja >nul 2>&1 || ( echo FAIL: ninja not found - install VS "C++ CMake tools for Windows" & exit /b 1 )
where ifx   >nul 2>&1 || ( echo FAIL: Intel Fortran compiler ifx not found & exit /b 1 )

REM some conan recipes (tcl, hdf5) have cmake_minimum_required < 3.5, which CMake 4 rejects
set "CMAKE_ARGS=-DCMAKE_POLICY_VERSION_MINIMUM=3.5"

call conan profile detect --exist-ok
if errorlevel 1 ( echo FAIL: conan profile detect failed & exit /b 1 )

echo OpenSees build: build=%BUILD_DIR% install=%OUTPUT_DIR% jobs=%JOBS% mumps=%MUMPS_DIR%

if "%CLEAN%"=="1" (
    echo Cleaning up old build directory...
    if exist "%BUILD_DIR%" rd /s /q "%BUILD_DIR%"
)

REM
REM run conan to install dependencies (into build\Release via cmake_layout)
REM

call conan install . -s arch=x86_64 -s build_type=Release -s compiler.runtime=static --build=missing -c tools.cmake.cmaketoolchain:generator=Ninja
if errorlevel 1 ( echo FAIL: Conan install failed. & exit /b 1 )

REM
REM run cmake to generate build files, then build & finally install
REM

set "CMAKE_OPTS=-DCMAKE_BUILD_TYPE=Release -DCMAKE_Fortran_COMPILER=ifx"
set "CMAKE_OPTS=%CMAKE_OPTS% -DBLA_STATIC=ON -DMKL_LINK=static -DMKL_INTERFACE_FULL=intel_lp64"
set "CMAKE_OPTS=%CMAKE_OPTS% -DCMAKE_NINJA_FORCE_RESPONSE_FILE=ON"
set "CMAKE_OPTS=%CMAKE_OPTS% -DCMAKE_C_USE_RESPONSE_FILE_FOR_OBJECTS=ON -DCMAKE_CXX_USE_RESPONSE_FILE_FOR_OBJECTS=ON"
if defined MUMPS_DIR set "CMAKE_OPTS=%CMAKE_OPTS% -DMUMPS_DIR="%MUMPS_DIR%""

cmake -S . -B %BUILD_DIR%\Release -G Ninja -DCMAKE_TOOLCHAIN_FILE=%BUILD_DIR%\Release\generators\conan_toolchain.cmake %CMAKE_OPTS%
if errorlevel 1 ( echo FAIL: CMake configure failed. & exit /b 1 )

cmake --build %BUILD_DIR%\Release --parallel %JOBS%
if errorlevel 1 ( echo FAIL: CMake build failed. & exit /b 1 )

cmake --build %BUILD_DIR%\Release --target OpenSeesPy --parallel %JOBS%
if errorlevel 1 ( echo FAIL: OpenSeesPy build failed. & exit /b 1 )

cmake --install %BUILD_DIR%\Release --prefix "%OUTPUT_DIR%"
if errorlevel 1 ( echo FAIL: CMake install failed. & exit /b 1 )

echo Done. Installed to %OUTPUT_DIR%\bin
exit /b 0

:usage
echo Usage: makeWIN.bat [--mumps ^<mumps build dir^>] [--prefix ^<install dir^>] [--jobs ^<n^>] [--no-clean]
echo.
echo   --mumps    MUMPS build folder; needed (with MPI) for OpenSeesSP ^& OpenSeesMP.
echo   --prefix   install folder. Default: %%USERPROFILE%%\bin\OpenSeesLatest
echo   --jobs     parallel compile jobs. Default: number of cores.
echo   --no-clean keep the existing build folder (faster rebuilds).
exit /b 1
