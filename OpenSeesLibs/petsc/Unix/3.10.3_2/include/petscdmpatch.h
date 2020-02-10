/*
  DMPatch, for domains covered by sets of patches.
*/
#if !defined(__PETSCDMPATCH_H)
#define __PETSCDMPATCH_H
#include <petscdm.h>

/*S
  DMPATCH - DM object that encapsulates a domain divided into many patches

  Level: intermediate

  Concepts: grids, grid refinement

.seealso:  DM, DMPatchCreate()
S*/
PETSC_EXTERN PetscErrorCode DMPatchCreate(MPI_Comm, DM*);

PETSC_EXTERN PetscErrorCode DMPatchZoom(DM,Vec,MatStencil,MatStencil,MPI_Comm,DM*,PetscSF*,PetscSF*);
PETSC_EXTERN PetscErrorCode DMPatchSolve(DM);
PETSC_EXTERN PetscErrorCode DMPatchGetPatchSize(DM,MatStencil*);
PETSC_EXTERN PetscErrorCode DMPatchSetPatchSize(DM,MatStencil);
PETSC_EXTERN PetscErrorCode DMPatchGetCommSize(DM,MatStencil*);
PETSC_EXTERN PetscErrorCode DMPatchSetCommSize(DM,MatStencil);
PETSC_EXTERN PetscErrorCode DMPatchGetCoarse(DM,DM*);
PETSC_EXTERN PetscErrorCode DMPatchCreateGrid(MPI_Comm,PetscInt,MatStencil,MatStencil,MatStencil,DM*);

/*
 * We want each patch to consist of an entire DM, DMDA at first
 - We cannot afford to store much more than the data from a single patch in memory
   - No global PetscSection, only PetscLayout
   - Optional scatters
   * There is a storable coarse level, which will also be a traditional DM (DMDA here)
   * The local and global vectors correspond to a ghosted patch
 * Need a way to activate a patch
   * Jack in sizes for l/g vectors
 - Need routine for viewing a full global vector
 - Jed handles solver
*/

#endif
