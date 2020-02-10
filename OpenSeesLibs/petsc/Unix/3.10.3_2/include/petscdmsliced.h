/* Very minimal unstructured DM */
#if !defined(__PETSCDMSLICED_H)
#define __PETSCDMSLICED_H

#include <petscdm.h>

PETSC_EXTERN PetscErrorCode DMSlicedCreate(MPI_Comm,PetscInt,PetscInt,PetscInt,const PetscInt[],const PetscInt[],const PetscInt[],DM*);
PETSC_EXTERN PetscErrorCode DMSlicedSetPreallocation(DM,PetscInt,const PetscInt[],PetscInt,const PetscInt[]);
PETSC_EXTERN PetscErrorCode DMSlicedSetBlockFills(DM,const PetscInt*,const PetscInt*);
PETSC_EXTERN PetscErrorCode DMSlicedSetGhosts(DM,PetscInt,PetscInt,PetscInt,const PetscInt[]);

#endif
