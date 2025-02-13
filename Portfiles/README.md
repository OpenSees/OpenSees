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
mkdir -p /Users/devops/ports

port self-update
port update
port install mumps
port install superlu 
port install superlu_dist
port install arpack
port install SuiteSparse_UMFPACK
port install libevent
port install tcl@8.6.16
port install sqlite3-tcl
port install qt5-qtcreator
port select --set python python311
port select --set python3 python311
```
**Note on Superlu_MT**

`mkdir -p ~/ports/math/superlu_mt`

copy Portfile into the folder and run portindex from folder /Users/devops/ports

After you can install: `port install superlu_mt`

After that you have installed all deps follow the build instruction for MacOS
and replace the CMakeLists.txt on the root of the repos with the file that you find in the
Portfile folder
