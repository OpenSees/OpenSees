#
# Script to build OpenSees on AWS Ubuntu 22.04 instance
#  - Build in t2.medium with about 124GB RAM
#
# written: fmk

#
# get some packages
#

sudo apt-get update
sudo apt install -y cmake
sudo apt install -y gcc-11 g++-11 gfortran-11
sudo apt install -y python3-pip python3-venv
sudo apt install -y liblapack-dev
sudo apt install -y libopenmpi-dev
sudo apt install -y libmkl-rt
sudo apt install -y libmkl-blacs-openmpi-lp64
sudo apt install -y libscalapack-openmpi-dev

#
# build mumps
#

git clone https://github.com/OpenSees/mumps.git
cd mumps
mkdir build
cd build
cmake .. -Darith=d
cmake --build . --config Release --parallel 4
cd ../..

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
mkdir build
conan install conan.txt --output-folder=build --build=missing
cd build
cmake .. -DCMAKE_TOOLCHAIN_FILE=conan_toolchain.cmake -DCMAKE_BUILD_TYPE=Release  -DMUMPS_DIR=$PWD/../../mumps/build -DOPENMPI=TRUE -DSCALAPACK_LIBRARIES="/usr/lib/x86_64-linux-gnu/libmkl_blacs_openmpi_lp64.so;/usr/lib/x86_64-linux-gnu/libscalapack-openmpi.so"  -DPython_LIBRARIES=/usr/lib/x86_64-linux-gnu/libpython3.10.so -DPython_INCLUDE_DIRS=/usr/include/python3.10 -DLAPACK_LIBRARIES="-m64  -Wl,--start-group /usr/lib/x86_64-linux-gnu/libmkl_intel_ilp64.so /usr/lib/x86_64-linux-gnu/libmkl_sequential.so /usr/lib/x86_64-linux-gnu/libmkl_core.so -Wl,--end-group -lpthread"

cmake --build . --config Release --target OpenSees   --parallel 4
cmake --build . --config Release --target OpenSeesPy --parallel 4
cmake --build . --config Release --target OpenSeesSP --parallel 4
cmake --build . --config Release --target OpenSeesMP --parallel 4

