#if !defined(_PATCHIMPL_H)
#define _PATCHIMPL_H

#include <petscmat.h>     /*I      "petscmat.h"        I*/
#include <petscdmpatch.h> /*I      "petscdmpatch.h"    I*/
#include <petsc/private/dmimpl.h>

typedef struct {
  PetscInt   refct;
  DM         dmCoarse;
  MatStencil patchSize;
  MatStencil commSize;
} DM_Patch;

#endif /* _PATCHIMPL_H */
