/*
  This is the include file for source code that needs to know if the code is being run under valgrind
*/
#if !defined(__PETSCVALGRIND_H)
#define __PETSCVALGRIND_H

#if defined(PETSC_HAVE_VALGRIND)
#  include <valgrind/valgrind.h>
#  define PETSC_RUNNING_ON_VALGRIND RUNNING_ON_VALGRIND
#else
#  define PETSC_RUNNING_ON_VALGRIND PETSC_FALSE
#endif

#endif
