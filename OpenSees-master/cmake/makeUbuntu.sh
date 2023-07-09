#! /bin/bash
cd ..
mkdir build                                                            
cd build                                                               
/home/ubuntu/.local/bin/conan install .. --build missing               
cmake .. -DMUMPS_DIR=$PWD/../../mumps/build -DOPENMPI=TRUE -DSCALAPACK_LIBRARIES="/usr/lib/x86_64-linux-gnu/libmkl_blacs_openmpi_lp64.so;/usr/lib/x86_64-linux-gnu/libscalapack-openmpi.so.2.1"                                      
cmake --build . --config Release --target OpenSees --parallel 4        
cmake --build . --config Release --target OpenSeesPy                   
cmake --build . --config Release --target OpenSeesMP                   
cmake --build . --config Release --target OpenSeesSP                   
mv ./OpenSeesPy.so ./opensees.so
