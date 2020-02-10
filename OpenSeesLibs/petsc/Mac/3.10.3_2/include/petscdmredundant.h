/* DM for redundant globally coupled degrees of freedom */
#if !defined(__PETSCDMREDUNDANT_H)
#define __PETSCDMREDUNDANT_H

#include <petscdm.h>

PETSC_EXTERN PetscErrorCode DMRedundantCreate(MPI_Comm,PetscMPIInt,PetscInt,DM*);
PETSC_EXTERN PetscErrorCode DMRedundantSetSize(DM,PetscMPIInt,PetscInt);
PETSC_EXTERN PetscErrorCode DMRedundantGetSize(DM,PetscMPIInt*,PetscInt*);

#endif
