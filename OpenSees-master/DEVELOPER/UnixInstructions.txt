
for Mac users you need to edit the first part of the Makefile.def, uncomment
the first line and comment or delete the second line!


there are a number of example classes and example scripts to test these
classes. In addition for elements and materials there are example c and 
fortran routines.

typing make in the DEVELOPER directory will build them all. the resulting
shared object libraries (dynamic link libraries) will be found in the same
directory as the source code.

to test them you need to ensure that the LD_LIBRARY_PATH is set to include ./

to test them run 'OpenSees example1.tcl', note OpenSees must be on your path.

