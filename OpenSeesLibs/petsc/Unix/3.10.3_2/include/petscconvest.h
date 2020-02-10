/*
      Objects which encapsulate discretizations+continuum residuals
*/
#if !defined(__PETSCCE_H)
#define __PETSCCE_H
#include <petscsnes.h>

/*S
  PetscConvEst - Provides an estimated convergence rate for a discretized problem

  Level: developer

.seealso:  PetscConvEstCreate(), PetscConvEstDestroy()
S*/
typedef struct _p_PetscConvEst *PetscConvEst;

PETSC_EXTERN PetscErrorCode PetscConvEstCreate(MPI_Comm, PetscConvEst *);
PETSC_EXTERN PetscErrorCode PetscConvEstDestroy(PetscConvEst *);
PETSC_EXTERN PetscErrorCode PetscConvEstView(PetscConvEst, PetscViewer);
PETSC_EXTERN PetscErrorCode PetscConvEstSetFromOptions(PetscConvEst);
PETSC_EXTERN PetscErrorCode PetscConvEstGetSolver(PetscConvEst, SNES *);
PETSC_EXTERN PetscErrorCode PetscConvEstSetSolver(PetscConvEst, SNES);
PETSC_EXTERN PetscErrorCode PetscConvEstSetUp(PetscConvEst);
PETSC_EXTERN PetscErrorCode PetscConvEstGetConvRate(PetscConvEst, PetscReal[]);
PETSC_EXTERN PetscErrorCode PetscConvEstRateView(PetscConvEst, PetscReal[], PetscViewer);

#endif
