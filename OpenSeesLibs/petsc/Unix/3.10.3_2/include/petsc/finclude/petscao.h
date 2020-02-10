!
!
!  Include file for Fortran use of the AO (application ordering) package in PETSc
!
#if !defined (__PETSCAODEF_H)
#define __PETSCAODEF_H

#include "petsc/finclude/petscis.h"

#define AO PetscFortranAddr
#define AOType character*(80)

#define AOBASIC           'basic'
#define AOADVANCED        'advanced'
#define AOMAPPING         'mapping'
#define AOMEMORYSCALABLE  'memoryscalable'

#endif
