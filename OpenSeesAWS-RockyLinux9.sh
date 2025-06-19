#!/bin/bash
#
# Script to build OpenSees on AWS Rocky Linux 9 instance
#  - Build in t2.medium with about 124GB RAM
#  - Minimal Linux serial build (no MPI, no parallel features)
#
# written: fmk (adapted for Rocky Linux 9)

#
# get some packages
#

sudo dnf update -y
sudo dnf groupinstall -y "Development Tools"
sudo dnf install -y cmake
sudo dnf install -y gcc gcc-c++ gcc-gfortran
sudo dnf install -y libstdc++-devel # Ensure libstdc++ development package is installed
sudo dnf install -y python3-pip python3-virtualenv
sudo dnf install -y lapack-devel blas-devel

# Enable PowerTools/CRB repository for additional packages
sudo dnf config-manager --set-enabled crb

# Install EPEL repository for additional packages
sudo dnf install -y epel-release

# Install Intel MKL (if available) or use alternative
# Note: Intel MKL may not be directly available in Rocky Linux 9 repos
# We'll use the standard BLAS/LAPACK libraries instead
sudo dnf install -y flexiblas-devel

#
# hdf5
#

git clone --depth 1 --branch hdf5-1_12_2 https://github.com/HDFGroup/hdf5.git
cd hdf5
./configure --prefix=/usr/local/hdf5
make -j 4
sudo make install
cd ..

#
# install conan & perform initial setup
#

python3 -m venv conan
source ./conan/bin/activate
python3 -m pip install conan
conan profile detect --force

#
# OpenSees Executables and Python lib
#

git clone https://github.com/OpenSees/OpenSees.git
cd OpenSees/

# Create build directory and install dependencies with conan
mkdir -p build
cd build
conan install .. --build missing

# Get Python version and paths for Rocky Linux 9 (typically Python 3.9)
PYTHON_VERSION=$(python3 -c "import sys; print(f'{sys.version_info.major}.{sys.version_info.minor}')")
PYTHON_LIB_PATH="/usr/lib64/libpython${PYTHON_VERSION}.so"
PYTHON_INCLUDE_PATH="/usr/include/python${PYTHON_VERSION}"

# Find the correct library paths for Rocky Linux 9
LAPACK_LIBS="/lib64/liblapack.so"

# Use the Conan toolchain file directly instead of CMake preset
cmake .. \
    -DCMAKE_BUILD_TYPE=Release \
    -DPython_LIBRARIES="$PYTHON_LIB_PATH" \
    -DPython_INCLUDE_DIRS="$PYTHON_INCLUDE_PATH" \
    -DLAPACK_LIBRARIES="$LAPACK_LIBS"

# Build only the serial targets (removed OpenSeesSP and OpenSeesMP)
make OpenSees -j32
make OpenSeesPy -j32
