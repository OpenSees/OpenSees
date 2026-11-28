#!/bin/bash
#
# Build OpenSees, OpenSeesPy (and OpenSeesSP/OpenSeesMP if MUMPS given) on MacOS.
#
# Usage: makeMac.sh [--arch arm64|x86_64] [--mumps <mumps build dir>]
#                   [--prefix <install dir>] [--jobs <n>] [--no-clean]
#
#   --arch     build for a specific architecture using conan profile
#              macos-arm64 or macos-x86_64-gcc (must exist in ~/.conan2/profiles).
#              Default: the native arch, using the conan "default" profile.
#   --mumps    MUMPS build folder; needed (with MPI) for OpenSeesSP & OpenSeesMP.
#   --prefix   install folder. Default: $HOME/bin/OpenSeesLatest[_<arch>]
#   --jobs     parallel compile jobs. Default: number of cores.
#   --no-clean keep the existing build folder (faster rebuilds).
#
# Requires: conan >= 2, cmake, Xcode command line tools, gfortran (brew install gcc).
#
# NOTE if fails in nataf_gsa on Mac it is because clang still not working with omp,
# SOLN: you need a brew install of libomp, issue:
#           brew install libomp

usage() {
    sed -n '5,14p' "$0" | sed 's/^# \{0,1\}//'
    exit 1
}

die() {
    echo "FAIL: $1"
    exit 1
}

# always run from the OpenSees source folder (where this script lives)
cd "$(dirname "$0")" || die "cannot cd to source folder"

# default compilers
export CC=clang
export CXX=clang++
export FC=gfortran

# some conan recipes (tcl, hdf5) have cmake_minimum_required < 3.5, which CMake 4 rejects
export CMAKE_ARGS=-DCMAKE_POLICY_VERSION_MINIMUM=3.5

CONAN_PROFILE="default"
BUILD_DIR="build"
OUTPUT_DIR=""
ARCH=""
MUMPS_DIR=""
JOBS="$(sysctl -n hw.ncpu 2>/dev/null || echo 4)"
CLEAN=1

while [[ $# -gt 0 ]]; do
    case "$1" in
        --arch|-a)
            ARCH="$2"
            case "$ARCH" in
                arm64)
                    CONAN_PROFILE="macos-arm64"
                    ;;
                x86_64)
                    # x86_64 build uses Intel homebrew gcc (under Rosetta) in /usr/local
                    CONAN_PROFILE="macos-x86_64-gcc"
                    export FC=/usr/local/bin/gfortran
                    export CC=/usr/local/bin/gcc-15
                    export CXX=/usr/local/bin/g++-15
                    ;;
                *)
                    die "invalid architecture '$ARCH' (use arm64 or x86_64)"
                    ;;
            esac
            BUILD_DIR="build_$ARCH"
            shift 2
            ;;
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

if [[ -z "$OUTPUT_DIR" ]]; then
    OUTPUT_DIR="$HOME/bin/OpenSeesLatest${ARCH:+_$ARCH}"
fi

# check tools & conan profile
command -v conan >/dev/null || die "conan not found (pip install conan)"
command -v cmake >/dev/null || die "cmake not found (brew install cmake)"
command -v "$FC"  >/dev/null || die "fortran compiler '$FC' not found (brew install gcc)"

if [[ "$CONAN_PROFILE" == "default" ]]; then
    conan profile detect --exist-ok || die "conan profile detect failed"
else
    conan profile path "$CONAN_PROFILE" >/dev/null 2>&1 \
        || die "conan profile '$CONAN_PROFILE' not found in ~/.conan2/profiles"
fi

echo "OpenSees build: arch=${ARCH:-native} profile=$CONAN_PROFILE build=$BUILD_DIR"
echo "                install=$OUTPUT_DIR jobs=$JOBS mumps=${MUMPS_DIR:-none}"

if [[ $CLEAN -eq 1 ]]; then
    echo "Cleaning up old build directory..."
    rm -rf "${BUILD_DIR}"
fi

#
# run conan to install dependencies
#

conan install . --output-folder="${BUILD_DIR}" --build=missing -s build_type=Release \
    -pr "${CONAN_PROFILE}" || die "Conan install failed."

#
# run cmake to generate build files, then build & finally install
#

CMAKE_OPTS=(-DCMAKE_BUILD_TYPE=Release)
CMAKE_OPTS+=(-DCMAKE_TOOLCHAIN_FILE="${BUILD_DIR}/build/Release/generators/conan_toolchain.cmake")
[[ -n "$ARCH" ]]      && CMAKE_OPTS+=(-DCMAKE_OSX_ARCHITECTURES="${ARCH}")
[[ -n "$MUMPS_DIR" ]] && CMAKE_OPTS+=(-DMUMPS_DIR="${MUMPS_DIR}")

cmake -S . -B "${BUILD_DIR}" "${CMAKE_OPTS[@]}" || die "CMake configure failed."

cmake --build "${BUILD_DIR}" --parallel "${JOBS}" || die "CMake build failed."
cmake --build "${BUILD_DIR}" --target OpenSeesPy --parallel "${JOBS}" || die "OpenSeesPy build failed."

cmake --install "${BUILD_DIR}" --prefix "${OUTPUT_DIR}" || die "CMake install failed."

echo "Done. Installed to ${OUTPUT_DIR}/bin"
