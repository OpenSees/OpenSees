#!/bin/bash

#
# NOTE: assumes conan>2, cmake, gcc, g++ and gfortran, python3-pip and liblapack-dev are apt installed
#  and that conan profile detect has been run and available in current env .. see OpenSeesAWS-Ubuntu22.04.sh if lost

export CC=gcc
export CXX=g++
export FC=gfortran

export CMAKE_ARGS=-DCMAKE_POLICY_VERSION_MINIMUM=3.5

BUILD_DIR="build"
OUTPUT_DIR="$HOME/bin/OpenSees"

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

conan install . --output-folder="${BUILD_DIR}" --build missing -s build_type=Release || die "FAIL: Conan install failed."

#
# run cmake to generate CMakefiles, then build & finally install
#

cmake -B "${BUILD_DIR}" -S . \
    -DCMAKE_TOOLCHAIN_FILE="${BUILD_DIR}/build/Release/generators/conan_toolchain.cmake" \
    -DCMAKE_BUILD_TYPE=Release || die "FAIL: CMake configure failed."

cmake --build "${BUILD_DIR}" --parallel 4 || die "FAIL: CMake build failed."
cmake --build "${BUILD_DIR}" --target OpenSeesPy || die "FAIL: CMake build failed."

cmake --install "${BUILD_DIR}" --prefix "${OUTPUT_DIR}" || die "FAIL: CMake install failed."

