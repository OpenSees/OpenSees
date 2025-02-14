# Command to run to install all deps:

### Deps that must be installed:


- MUMPS
- Scalapack
- UMFPACK contained in Suite-Sparse 
- SUPERLU
- SUPERLUMT
- SUPERLUDIST
- Openblas
- Parmetis
- ARPACK
- Libevent
- GCC
- TCL (8.x)
- Python3

Update file:  /opt/local/etc/macports/sources.conf with the content on this repos.

Run the follow commands:

```
mkdir -p ~/ports

sudo port self-update
sudo port update
sudo port install mumps
sudo port install superlu 
sudo port install superlu_dist
sudo port install arpack
sudo port install SuiteSparse_UMFPACK
sudo port install libevent
sudo port install tcl@8.6.16
sudo port install sqlite3-tcl
sudo port install eigen3
sudo port install python310
sudo port select --set python python310
sudo port select --set python3 python310
```
**Install hdf5@1.12**

```
mkdir -p ~/GIT/
cd ~/GIT/
git clone --single-branch https://github.com/macports/macports-ports.git
cd macports-ports
git checkout 7bd523e
cd science/hdf5
sudo port install
```

**Note on Superlu_MT**

`mkdir -p ~/ports/math/superlu_mt`

copy Portfile into the folder and run portindex from folder /Users/devops/ports
```
cp Portfiles_superlu_mt ~/ports/math/superlu_mt/Portfile
cd ~/ports
sudo portindex
Creating port index in /Users/devops/ports
Adding port math/superlu_mt

Total number of ports parsed:	1
Ports successfully parsed:	1
Ports failed:			0
Up-to-date ports skipped:	0
```

After you can install: `sudo port install superlu_mt`


After that you have installed all deps follow the build instruction for MacOS
and replace the CMakeLists.txt on the root of the repos with the file that you find in the
Portfile folder.

To select the right version of Python use the follow cmake command:

```
#! /bin/bash
cd ..
mkdir build                                                            
cd build                                                               
cmake .. -DMUMPS_DIR=/opt/local/lib -DOPENMPI=TRUE -DSCALAPACK_LIBRARIES=/opt/local/lib/libscalapack.dylib -DPython_LIBRARIES=/opt/local/Library/Frameworks/Python.framework/Versions/3.10/lib/libpython3.10.dylib -DPython_INCLUDE_DIRS=/opt/local/Library/Frameworks/Python.framework/Versions/3.10/include/python3.10/ -DPython_VERSION=3.10
cmake --build . --target OpenSees -j4
cmake --build . --target OpenSeesPy -j4
mv ./OpenSeesPy.dylib ./opensees.so
```
