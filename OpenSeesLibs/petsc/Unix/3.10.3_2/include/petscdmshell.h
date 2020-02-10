#if !defined(__PETSCDMSHELL_H)
#define __PETSCDMSHELL_H

#include <petscdm.h>

PETSC_EXTERN PetscErrorCode DMShellCreate(MPI_Comm,DM*);
PETSC_EXTERN PetscErrorCode DMShellSetContext(DM,void *);
PETSC_EXTERN PetscErrorCode DMShellGetContext(DM,void **);
PETSC_EXTERN PetscErrorCode DMShellSetMatrix(DM,Mat);
PETSC_EXTERN PetscErrorCode DMShellSetGlobalVector(DM,Vec);
PETSC_EXTERN PetscErrorCode DMShellSetLocalVector(DM,Vec);
PETSC_EXTERN PetscErrorCode DMShellSetCreateGlobalVector(DM,PetscErrorCode (*)(DM,Vec*));
PETSC_EXTERN PetscErrorCode DMShellSetCreateLocalVector(DM,PetscErrorCode (*)(DM,Vec*));
PETSC_EXTERN PetscErrorCode DMShellSetGlobalToLocal(DM,PetscErrorCode (*)(DM,Vec,InsertMode,Vec),PetscErrorCode (*)(DM,Vec,InsertMode,Vec));
PETSC_EXTERN PetscErrorCode DMShellSetGlobalToLocalVecScatter(DM,VecScatter);
PETSC_EXTERN PetscErrorCode DMShellSetLocalToGlobal(DM,PetscErrorCode (*)(DM,Vec,InsertMode,Vec),PetscErrorCode (*)(DM,Vec,InsertMode,Vec));
PETSC_EXTERN PetscErrorCode DMShellSetLocalToGlobalVecScatter(DM,VecScatter);
PETSC_EXTERN PetscErrorCode DMShellSetLocalToLocal(DM,PetscErrorCode (*)(DM,Vec,InsertMode,Vec),PetscErrorCode (*)(DM,Vec,InsertMode,Vec));
PETSC_EXTERN PetscErrorCode DMShellSetLocalToLocalVecScatter(DM,VecScatter);
PETSC_EXTERN PetscErrorCode DMShellSetCreateMatrix(DM,PetscErrorCode (*)(DM,Mat*));
PETSC_EXTERN PetscErrorCode DMShellSetCoarsen(DM,PetscErrorCode (*)(DM,MPI_Comm,DM*));
PETSC_EXTERN PetscErrorCode DMShellSetRefine(DM,PetscErrorCode (*)(DM,MPI_Comm,DM*));
PETSC_EXTERN PetscErrorCode DMShellSetCreateInterpolation(DM,PetscErrorCode (*)(DM,DM,Mat*,Vec*));
PETSC_EXTERN PetscErrorCode DMShellSetCreateRestriction(DM, PetscErrorCode (*)(DM,DM,Mat*));
PETSC_EXTERN PetscErrorCode DMShellSetCreateInjection(DM,PetscErrorCode (*)(DM,DM,Mat*));
PETSC_EXTERN PetscErrorCode DMShellSetCreateFieldDecomposition(DM,PetscErrorCode (*)(DM,PetscInt*,char***,IS**,DM**));
PETSC_EXTERN PetscErrorCode DMShellSetCreateDomainDecomposition(DM,PetscErrorCode (*)(DM,PetscInt*,char***,IS**,IS**,DM**));
PETSC_EXTERN PetscErrorCode DMShellSetCreateDomainDecompositionScatters(DM,PetscErrorCode (*)(DM,PetscInt,DM*,VecScatter**,VecScatter**,VecScatter**));
PETSC_EXTERN PetscErrorCode DMShellSetCreateSubDM(DM,PetscErrorCode (*)(DM,PetscInt,const PetscInt[],IS*,DM*));
PETSC_EXTERN PetscErrorCode DMGlobalToLocalBeginDefaultShell(DM,Vec,InsertMode,Vec);
PETSC_EXTERN PetscErrorCode DMGlobalToLocalEndDefaultShell(DM,Vec,InsertMode,Vec);
PETSC_EXTERN PetscErrorCode DMLocalToGlobalBeginDefaultShell(DM,Vec,InsertMode,Vec);
PETSC_EXTERN PetscErrorCode DMLocalToGlobalEndDefaultShell(DM,Vec,InsertMode,Vec);
PETSC_EXTERN PetscErrorCode DMLocalToLocalBeginDefaultShell(DM,Vec,InsertMode,Vec);
PETSC_EXTERN PetscErrorCode DMLocalToLocalEndDefaultShell(DM,Vec,InsertMode,Vec);

#endif
