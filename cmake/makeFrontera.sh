#! /bin/bash
cd ..
module load intel
module load hdf5/1.12.0
mkdir -p build
cd build
MKL_LIB=${TACC_MKL_LIB}
export CC=icc; export CXX=icpc; export FC=ifort;
cmake .. -DBLA_STATIC=ON -DBLA_VENDOR=Intel10_64lp -DSCALAPACK_LIBRARIES="-Wl,--start-group $TACC_MKL_LIB/libmkl_scalapack_lp64.a $TACC_MKL_LIB/libmkl_core.a $TACC_MKL_LIB/libmkl_blacs_intelmpi_lp64.a $TACC_MKL_LIB/libmkl_intel_lp64.a $TACC_MKL_LIB/libmkl_sequential.a $TACC_MKL_LIB/libmkl_core.a $TACC_MKL_LIB/libmkl_sequential.a $TACC_MKL_LIB/libmkl_core.a -Wl,--end-group" -DLAPACK_LIBRARIES="-Wl,--start-group  $TACC_MKL_LIB/libmkl_intel_lp64.a $TACC_MKL_LIB/libmkl_sequential.a $TACC_MKL_LIB/libmkl_core.a  -Wl,--end-group -lpthread" -DMUMPS_DIR=/work2/00477/tg457427/frontera/mumps/build
cmake --build . --config Release --target OpenSees --parallel 4
cmake --build . --config Release --target OpenSeesPy 
cmake --build . --config Release --target OpenSeesMP
cmake --build . --config Release --target OpenSeesSP
cd ..


