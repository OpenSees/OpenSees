#!/usr/bin/env bash
set -euo pipefail

# ---------------------------------------------------------------------------
# OpenSeesMP Ubuntu build (Ninja + Conan + Intel MPI/MKL).
# Edit these paths for your machine or override them in the environment:
#   IMPI_ROOT=/opt/intel/oneapi/mpi/latest \
#   MUMPS_DIR=/path/to/mumps/build ./makeOpenSeesMP_Ubuntu.sh
# ---------------------------------------------------------------------------
ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
BUILD_DIR="${BUILD_DIR:-${ROOT_DIR}/build-mp}"
MUMPS_DIR="${MUMPS_DIR:-${ROOT_DIR}/../mumps/build}"
METIS5_DIR="${METIS5_DIR:-/usr}"
IMPI_ROOT="${IMPI_ROOT:-/opt/intel/oneapi/mpi/2021.16}"
MKL_LIB="${MKL_LIB:-/opt/intel/oneapi/mkl/2025.2/lib}"
JOBS="${JOBS:-$(nproc)}"

cd "${ROOT_DIR}"

if [[ -f /opt/intel/oneapi/setvars.sh ]]; then
  set +u
  # shellcheck source=/dev/null
  source /opt/intel/oneapi/setvars.sh --force
  set -u
fi

if [[ ! -x "${IMPI_ROOT}/bin/mpigcc" ]]; then
  echo "ERROR: Intel MPI not found under: ${IMPI_ROOT}" >&2
  exit 1
fi
if [[ ! -f "${MUMPS_DIR}/libdmumps.a" ]]; then
  echo "ERROR: MUMPS libraries not found under: ${MUMPS_DIR}" >&2
  exit 1
fi
if [[ ! -f "${METIS5_DIR}/include/metis.h" ]]; then
  echo "ERROR: METIS 5 headers not found under: ${METIS5_DIR}" >&2
  echo "Install libmetis-dev or set METIS5_DIR." >&2
  exit 1
fi

# MUMPS was built against Intel MPI and LP64 MKL ScaLAPACK.
SCALAPACK_LIBRARIES="${SCALAPACK_LIBRARIES:-\
${MKL_LIB}/libmkl_scalapack_lp64.so;\
${MKL_LIB}/libmkl_gf_lp64.so;\
${MKL_LIB}/libmkl_gnu_thread.so;\
${MKL_LIB}/libmkl_core.so;\
${MKL_LIB}/libmkl_blacs_intelmpi_lp64.so}"

# Conan provides Tcl, HDF5, Eigen and zlib.
conan install . -of "${BUILD_DIR}" \
  -s build_type=Release \
  -s arch=x86_64 \
  --build=missing \
  -c tools.cmake.cmaketoolchain:generator=Ninja

# Conan 2 layouts can place the toolchain in either location.
TOOLCHAIN="${BUILD_DIR}/build/Release/generators/conan_toolchain.cmake"
if [[ ! -f "${TOOLCHAIN}" ]]; then
  TOOLCHAIN="${BUILD_DIR}/Release/generators/conan_toolchain.cmake"
fi
if [[ ! -f "${TOOLCHAIN}" ]]; then
  echo "ERROR: Conan CMake toolchain was not generated." >&2
  exit 1
fi

cmake -S . -B "${BUILD_DIR}/Release" -G Ninja \
  -DCMAKE_TOOLCHAIN_FILE="${TOOLCHAIN}" \
  -DCMAKE_BUILD_TYPE=Release \
  -DMPI_C_COMPILER="${IMPI_ROOT}/bin/mpigcc" \
  -DMPI_CXX_COMPILER="${IMPI_ROOT}/bin/mpigxx" \
  -DMPI_Fortran_COMPILER="${IMPI_ROOT}/bin/mpif90" \
  -DMUMPS_DIR="${MUMPS_DIR}" \
  -DMETIS5_DIR="${METIS5_DIR}" \
  -DSCALAPACK_LIBRARIES="${SCALAPACK_LIBRARIES}" \
  -UOPENMPI \
  -DPARALLEL_PROCESSING=OFF

cmake --build "${BUILD_DIR}/Release" \
  --target OpenSeesMP \
  --parallel "${JOBS}"

# OpenSeesMP looks for Tcl's init.tcl under build-mp/lib at runtime.
# Must match conanfile.py tcl/8.6.11 (not an older cached tcl/8.6.10 package).
TCL_VERSION="8.6.11"
TCL_STAGED=0
for TCL_PKG_LIB in "${HOME}/.conan2/p/b"/tcl*/p/lib; do
  INIT_TCL="${TCL_PKG_LIB}/tcl8.6/init.tcl"
  if [[ -f "${INIT_TCL}" ]] && grep -q "package require -exact Tcl ${TCL_VERSION}" "${INIT_TCL}"; then
    mkdir -p "${BUILD_DIR}/lib"
    cp -a "${TCL_PKG_LIB}/tcl8.6" "${TCL_PKG_LIB}/tcl8" "${BUILD_DIR}/lib/"
    TCL_STAGED=1
    break
  fi
done
if [[ "${TCL_STAGED}" -ne 1 ]]; then
  echo "ERROR: Conan Tcl ${TCL_VERSION} runtime was not found under ${HOME}/.conan2/p/b." >&2
  exit 1
fi

echo
echo "OpenSeesMP built successfully:"
echo "  ${BUILD_DIR}/Release/OpenSeesMP"
