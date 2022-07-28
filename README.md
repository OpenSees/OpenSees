# OpenSees Source Code Repository

This git repository contains all revisions to OpenSees source code since Version 2.3.2.

Older revisions of the code are available upon request.

If you plan on collaborating or even using OpenSees as your base code it is highly recommended that
you FORK this repo to your own account and work on it there. We will not allow anybody to write to
this repo. Only PULL requests will be considered. To fork the repo click on the FORK at the top of this page.

For a brief outline on forking we suggest:
https://www.atlassian.com/git/tutorials/comparing-workflows/forking-workflow

For a brief introduction to using your new FORK we suggest:
https://www.atlassian.com/git/tutorials/saving-changes

## Documentation
The documentation for OpenSees is being moved to a parellel github repo:
https://github.com/OpenSees/OpenSeesDocumentation

The documentation (in its present form) can be viewed in the browser using the following url:
https://OpenSees.github.io/OpenSeesDocumentation


Linux: (from a terminal)

cmake -DCMAKE_C_COMPILER=icc -DCMAKE_CXX_COMPILER=icpc -DCMAKE_Fortran_COMPILER=ifort ..

Windows: (From an Intel Compiler console:)

cmake -G "NMake Makefiles" -DCMAKE_C_COMPILER=icl -DCMAKE_CXX_COMPILER=icl -

## How to Build OpenSeesPy for Mac

Kazuki Ichinohe 2022/07/27

1. `xcode-select install`

    Install Command Line Tools, which make AppleClang and Git available.

   if `xcode-select: note: install requested for command line developer tools` appears, continue installation with GUI.

   if `xcode-select: error: command line tools are already installed, use "Software Update" to install updates` appears, skip because it's already installed.

2. `/bin/bash -c "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/master/install.sh)"`

    Install HomeBrew a package manager.

    Enter your password when it requires sudo authorization and press enter if needed.

3. `brew install cmake`
4. `brew install gfortran`
5. `brew install open-mpi`
6. `brew install python`
7. `pip3 install conan`

    cmake: A tool controling the software compilation process.

    gfortran: A compiler for Fortran.

    open-mpi: A library for multi processing.

    python: Python3.9

    >Since /usr/bin/python3 is Python3.8, install Python3.9 via HomeBrew. It requires reflesh the shell to reflect a new path setting. python@3.10 may be OK, but it sometimes fails to set the path automatically and need to set following a suggestion by HomeBrew.

    conan: A package manager for C and C++.

    >I think installing libraries via HomeBrew will work, but they enable a management with conan, so it might be preferable.

8. `git clone https://github.com/OpenSees/OpenSees.git`

    Command Line Tools include Git, so git clone the repository of OpenSees. "OpenSees" directory will be made at the working directory.

9. `cd OpenSees`
10. `mkdir build`

    Make a directory for build files. I heard it's common practice when using cmake. It doesn't contaminate the source directory.

11. `cd build`
12. `conan install ..`

    Load conanfile.py in .., install packages properly and make files needed after at the current directory.

13. `cmake -S .. -B ./`

    Prepare for build to ./ from the source ..

14. `mv ../OTHER/METIS/version ../OTHER/METIS/metis_version`

    Rename the file in order not to be included wrongly when AppleClang includes \<version\>.

15. `cmake --build ./ --target OpenSeesPy`

    Execute building with the target of OpenSeesPy.

16. `mv -f ./lib/OpenSeesPy.dylib /usr/local/lib/python3.9/site-packages/openseespymac/opensees.so`

    Replace the OpenSeesPy's library with OpenSeesPy.dylib generated at the ./lib. Its extension is also altered from .dylib to .so but it has no affect. The directory where OpenSeesPy's library is can be different.

