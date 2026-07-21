@echo off
setlocal EnableDelayedExpansion

REM ---------------------------------------------------------------------------
REM OpenSeesMP Windows build (Ninja + Conan + Intel oneAPI).
REM Edit these paths for your machine.
REM ---------------------------------------------------------------------------
set "BUILD_DIR=build-mp"
set "MUMPS_DIR=%~dp0..\mumps\build"
set "METIS5_DIR=%~dp0..\metis-5.1.0\install"
set "ONEAPI_SETVARS=C:\Program Files (x86)\Intel\oneAPI\setvars.bat"
set "JOBS=10"

cd /d "%~dp0"

if not exist "!ONEAPI_SETVARS!" (
  echo ERROR: Intel oneAPI setvars.bat not found:
  echo   !ONEAPI_SETVARS!
  exit /b 1
)
if not exist "!MUMPS_DIR!\dmumps.lib" (
  echo ERROR: MUMPS libraries not found under:
  echo   !MUMPS_DIR!
  exit /b 1
)
if not exist "!METIS5_DIR!\include\metis.h" (
  echo ERROR: METIS 5 headers not found under:
  echo   !METIS5_DIR!
  exit /b 1
)

call "!ONEAPI_SETVARS!" intel64 mod
if errorlevel 1 exit /b 1

REM Conan provides Tcl, HDF5, Eigen and zlib.
conan install . -of "!BUILD_DIR!" ^
  -s build_type=Release ^
  -s arch=x86_64 ^
  -s compiler.runtime=static ^
  --build=missing ^
  -c tools.cmake.cmaketoolchain:generator=Ninja
if errorlevel 1 exit /b 1

REM Conan 2 layouts can place the toolchain in either location.
set "TOOLCHAIN=!BUILD_DIR!\build\Release\generators\conan_toolchain.cmake"
if not exist "!TOOLCHAIN!" set "TOOLCHAIN=!BUILD_DIR!\Release\generators\conan_toolchain.cmake"
if not exist "!TOOLCHAIN!" (
  echo ERROR: Conan CMake toolchain was not generated.
  exit /b 1
)

cmake -S . -B "!BUILD_DIR!\Release" -G Ninja ^
  -DCMAKE_TOOLCHAIN_FILE="!TOOLCHAIN!" ^
  -DCMAKE_BUILD_TYPE=Release ^
  -DCMAKE_MSVC_RUNTIME_LIBRARY=MultiThreaded ^
  -DCMAKE_EXE_LINKER_FLAGS="/FORCE:MULTIPLE" ^
  -DCMAKE_SHARED_LINKER_FLAGS="/FORCE:MULTIPLE" ^
  -DCMAKE_Fortran_COMPILER=ifx ^
  -DBLA_STATIC=ON ^
  -DMKL_LINK=static ^
  -DMKL_INTERFACE_FULL=intel_lp64 ^
  -DMUMPS_DIR="!MUMPS_DIR!" ^
  -DMETIS5_DIR="!METIS5_DIR!" ^
  -UOPENSEES_METIS5_LIBRARY ^
  -UOPENMPI ^
  -DPARALLEL_PROCESSING=OFF ^
  -DCMAKE_NINJA_FORCE_RESPONSE_FILE=ON ^
  -DCMAKE_C_USE_RESPONSE_FILE_FOR_OBJECTS=ON ^
  -DCMAKE_CXX_USE_RESPONSE_FILE_FOR_OBJECTS=ON
if errorlevel 1 exit /b 1

cmake --build "!BUILD_DIR!\Release" --target OpenSeesMP --parallel !JOBS!
if errorlevel 1 exit /b 1

REM OpenSeesMP.exe looks for Tcl's init.tcl under build-mp\lib at runtime.
REM Must match conanfile.py tcl/8.6.11 (not an older cached tcl/8.6.10 package).
set "TCL_VERSION=8.6.11"
set "TCL_PKG_LIB="
for /d %%d in ("%USERPROFILE%\.conan2\p\b\tcl*") do (
  if exist "%%d\p\lib\tcl8.6\init.tcl" (
    findstr /C:"package require -exact Tcl !TCL_VERSION!" "%%d\p\lib\tcl8.6\init.tcl" >nul 2>&1 && (
      set "TCL_PKG_LIB=%%d\p\lib"
      goto tcl_found
    )
  )
)
echo ERROR: Conan Tcl !TCL_VERSION! not found under %USERPROFILE%\.conan2\p\b.
exit /b 1

:tcl_found
if not exist "!BUILD_DIR!\lib" mkdir "!BUILD_DIR!\lib"
robocopy "!TCL_PKG_LIB!\tcl8.6" "!BUILD_DIR!\lib\tcl8.6" /E /NFL /NDL /NJH /NJS /NC /NS >nul
if errorlevel 8 exit /b 1
robocopy "!TCL_PKG_LIB!\tcl8" "!BUILD_DIR!\lib\tcl8" /E /NFL /NDL /NJH /NJS /NC /NS >nul
if errorlevel 8 exit /b 1
echo Tcl !TCL_VERSION! copied to !BUILD_DIR!\lib\

echo.
echo OpenSeesMP built successfully:
echo   %CD%\!BUILD_DIR!\Release\OpenSeesMP.exe
endlocal
