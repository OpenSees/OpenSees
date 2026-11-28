#!/bin/bash
#
# Build OpenSees, OpenSeesPy (and OpenSeesSP/OpenSeesMP if MUMPS given) on Ubuntu.
#
# Usage: makeUbuntu.sh [--mumps <mumps build dir>] [--prefix <install dir>]
#                      [--jobs <n>] [--no-clean]
#
#   --mumps    MUMPS build folder; needed (with MPI) for OpenSeesSP & OpenSeesMP.
#   --prefix   install folder. Default: $HOME/bin/OpenSees
#   --jobs     parallel compile jobs. Default: number of cores.
#   --no-clean keep the existing build folder (faster rebuilds).
#
# Requires: conan >= 2 (pip install conan) and the apt packages:
#   cmake gcc g++ gfortran python3-pip liblapack-dev
#   (for OpenSeesSP/MP also: libopenmpi-dev libscalapack-openmpi-dev)
# see OpenSeesAWS-Ubuntu22.04.sh for a full machine setup.

usage() {
    sed -n '5,11p' "$0" | sed 's/^# \{0,1\}//'
    exit 1
}

die() {
    echo "FAIL: $1"
    exit 1
}

# always run from the OpenSees source folder (where this script lives)
cd "$(dirname "$0")" || die "cannot cd to source folder"

export CC=gcc
export CXX=g++
export FC=gfortran

# some conan recipes (tcl, hdf5) have cmake_minimum_required < 3.5, which CMake 4 rejects
export CMAKE_ARGS=-DCMAKE_POLICY_VERSION_MINIMUM=3.5

BUILD_DIR="build"
OUTPUT_DIR="$HOME/bin/OpenSees"
MUMPS_DIR=""
JOBS="$(nproc 2>/dev/null || echo 4)"
CLEAN=1

while [[ $# -gt 0 ]]; do
    case "$1" in
        --mumps|-m)
            [[ -d "$2" ]] || die "MUMPS folder '$2' does not exist"
            MUMPS_DIR="$(cd "$2" && pwd)"   # must be absolute
            shift 2
            ;;
        --prefix|-p)
            OUTPUT_DIR="$2"
            shift 2
            ;;
        --jobs|-j)
            JOBS="$2"
            shift 2
            ;;
        --no-clean)
            CLEAN=0
            shift
            ;;
        -h|--help)
            usage
            ;;
        *)
            echo "Unknown argument: $1"
            usage
            ;;
    esac
done

# check tools & conan profile
command -v conan >/dev/null || die "conan not found (pip install conan; may need \$HOME/.local/bin in PATH)"
command -v cmake >/dev/null || die "cmake not found (sudo apt install cmake)"
command -v "$FC"  >/dev/null || die "fortran compiler '$FC' not found (sudo apt install gfortran)"

conan profile detect --exist-ok || die "conan profile detect failed"

echo "OpenSees build: build=$BUILD_DIR install=$OUTPUT_DIR jobs=$JOBS mumps=${MUMPS_DIR:-none}"

if [[ $CLEAN -eq 1 ]]; then
    echo "Cleaning up old build directory..."
    rm -rf "${BUILD_DIR}"
fi

#
# run conan to install dependencies
#

conan install . --output-folder="${BUILD_DIR}" --build=missing -s build_type=Release \
    || die "Conan install failed."

#
# run cmake to generate build files, then build & finally install
#

CMAKE_OPTS=(-DCMAKE_BUILD_TYPE=Release)
CMAKE_OPTS+=(-DCMAKE_TOOLCHAIN_FILE="${BUILD_DIR}/build/Release/generators/conan_toolchain.cmake")
[[ -n "$MUMPS_DIR" ]] && CMAKE_OPTS+=(-DMUMPS_DIR="${MUMPS_DIR}")

cmake -S . -B "${BUILD_DIR}" "${CMAKE_OPTS[@]}" || die "CMake configure failed."

cmake --build "${BUILD_DIR}" --parallel "${JOBS}" || die "CMake build failed."
cmake --build "${BUILD_DIR}" --target OpenSeesPy --parallel "${JOBS}" || die "OpenSeesPy build failed."

cmake --install "${BUILD_DIR}" --prefix "${OUTPUT_DIR}" || die "CMake install failed."

echo "Done. Installed to ${OUTPUT_DIR}/bin"
