#!/bin/bash

# Usage: makeMac.sh [--arch <x86_64 or arm64>]
# Defaults: arch = ""

# NOTE if fails in nataf_gsa on Mac it is because clang still not working with omp,
# SOLN: you need a brew install of libomp, issue:
#           brew install libomp

# define the compilers
export CC=clang
export CXX=clang++
export FC=gfortran

export CMAKE_ARGS=-DCMAKE_POLICY_VERSION_MINIMUM=3.5

CONAN_PROFILE="default"
BUILD_DIR="build"
OUTPUT_DIR="$HOME/bin/OpenSeesLatest"
ARCH=""

while [[ $# -gt 0 ]]; do
    case "$1" in
        --arch|-a)
            # 1. Capture the value
            ARCH="$2"
            
            # 2. VALIDATION: Check if it's one of the allowed types
            if [[ "$ARCH" != "x86_64" && "$ARCH" != "arm64" ]]; then
                echo "Error: Invalid architecture '$ARCH'."
                echo "Supported architectures: x86_64, arm64"
                exit 1
            fi

            # 3. If valid, set the variables
            CONAN_PROFILE="macos-$ARCH"
            BUILD_DIR="build_$ARCH"
            OUTPUT_DIR="${OUTPUT_DIR}_$ARCH"	    

	    #. set fortran compiler, 
	    if [ "${ARCH}" = "x86_64" ]; then
		CONAN_PROFILE="macos-$ARCH-gcc"		
		export FC=/usr/local/bin/gfortran
		export CC=/usr/local/bin/gcc-15
		export CXX=/usr/local/bin/g++-15		
	    fi

	    
            shift 2
            ;;
        *)
            echo "Unknown argument: $1"
            echo "Usage: $0 --arch <arm64|x86_64>"
            exit 1
            ;;
    esac
done

# Define a helper to print error and stop without closing the terminal
die() {
    echo "$1"
    exit 1
}

echo "Cleaning up old build directory..."
rm -rf "${BUILD_DIR}"

#
# run conan to install dependencies
#

conan install . --output-folder="${BUILD_DIR}" --build missing -s build_type=Release -pr "${CONAN_PROFILE}" || die "FAIL: Conan install failed."

#
# run cmake to generate CMakefiles, then build & finally install

CMAKE_ARCH_FLAG=""
CMAKE_FORTRAN_LIBDIR=""
if [ -n "${ARCH}" ]; then
    CMAKE_ARCH_FLAG="-DCMAKE_OSX_ARCHITECTURES=${ARCH}"
fi

cmake -B "${BUILD_DIR}" -S . \
    -DCMAKE_TOOLCHAIN_FILE="${BUILD_DIR}/build/Release/generators/conan_toolchain.cmake" \
    ${CMAKE_ARCH_FLAG} \
    -DCMAKE_BUILD_TYPE=Release || die "FAIL: CMake configure failed."

cmake --build "${BUILD_DIR}" --parallel 8 || die "FAIL: CMake build failed."
cmake --build "${BUILD_DIR}" --target OpenSeesPy || die "FAIL: CMake build failed."

cmake --install "${BUILD_DIR}" --prefix "${OUTPUT_DIR}" || die "FAIL: CMake install failed."

#rm -fr build
#conan install . --build missing
#cmake -S . -B build/Release -DCMAKE_TOOLCHAIN_FILE=build/Release/generators/conan_toolchain.cmake -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX=$HOME/bin/OpenSees3.8.0
#cd build/Release
#cmake --build . --parallel 10
#cmake --build . --target OpenSeesPy
#cmake --install .
