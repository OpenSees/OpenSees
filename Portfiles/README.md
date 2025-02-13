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
sudo port install qt5-qtcreator
sudo port select --set python python311
sudo port select --set python3 python311
```
**Note on Superlu_MT**

`mkdir -p ~/ports/math/superlu_mt`

copy Portfile into the folder and run portindex from folder /Users/devops/ports

After you can install: `sudo port install superlu_mt`

After that you have installed all deps follow the build instruction for MacOS
and replace the CMakeLists.txt on the root of the repos with the file that you find in the
Portfile folder
