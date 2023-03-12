#! /bin/bash
module load gcc
module load hdf5/1.12.0
mkdir -p build
cd build
export CC=gcc; export CXX=g++; export FC=gfortran;
cmake .. -DLAPACK_LIBRARIES="$HOME/lib/libLapack.a;$HOME/lib/libBlas.a;-lz;-lpthread"
cmake --build . --config Release --target OpenSeesPy 
cp OpenSeesPy.so opensees
cd ..


