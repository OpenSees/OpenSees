/*
      Objects which encapsulate mesh adaptation operation
*/
#if !defined(__PETSCDMADAPTOR_H)
#define __PETSCDMADAPTOR_H
#include <petscdm.h>
#include <petscconvest.h>

/*S
  DMAdaptor - The adaptor constructs a DMLabel or metric Vec that can be used to modify the DM.

  Level: developer

.seealso:  PetscConvEstCreate(), PetscConvEstDestroy()
S*/
typedef struct _p_DMAdaptor *DMAdaptor;

PETSC_EXTERN PetscErrorCode DMAdaptorCreate(MPI_Comm, DMAdaptor *);
PETSC_EXTERN PetscErrorCode DMAdaptorDestroy(DMAdaptor *);
PETSC_EXTERN PetscErrorCode DMAdaptorView(DMAdaptor, PetscViewer);
PETSC_EXTERN PetscErrorCode DMAdaptorSetFromOptions(DMAdaptor);
PETSC_EXTERN PetscErrorCode DMAdaptorGetSolver(DMAdaptor, SNES *);
PETSC_EXTERN PetscErrorCode DMAdaptorSetSolver(DMAdaptor, SNES);
PETSC_EXTERN PetscErrorCode DMAdaptorGetSequenceLength(DMAdaptor, PetscInt *);
PETSC_EXTERN PetscErrorCode DMAdaptorSetSequenceLength(DMAdaptor, PetscInt);
PETSC_EXTERN PetscErrorCode DMAdaptorSetUp(DMAdaptor);
PETSC_EXTERN PetscErrorCode DMAdaptorGetTransferFunction(DMAdaptor, PetscErrorCode (**)(DMAdaptor, DM, Vec, DM, Vec, void *));
PETSC_EXTERN PetscErrorCode DMAdaptorSetTransferFunction(DMAdaptor, PetscErrorCode (*)(DMAdaptor, DM, Vec, DM, Vec, void *));
PETSC_EXTERN PetscErrorCode DMAdaptorAdapt(DMAdaptor, Vec, DMAdaptationStrategy, DM *, Vec *);

#endif
