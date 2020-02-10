#ifndef __PETSCMATHYPRE_H
#define __PETSCMATHYPRE_H

#include <petscmat.h>
#include <_hypre_parcsr_mv.h>

PETSC_EXTERN PetscErrorCode MatCreateFromParCSR(hypre_ParCSRMatrix*,MatType,PetscCopyMode,Mat*);
PETSC_EXTERN PetscErrorCode MatHYPREGetParCSR(Mat,hypre_ParCSRMatrix**);

#endif
