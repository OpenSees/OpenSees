#if !defined(_PETSCCEIMPL_H)
#define _PETSCCEIMPL_H

#include <petscconvest.h>
#include <petsc/private/petscimpl.h>

typedef struct _PetscConvEstOps *PetscConvEstOps;
struct _PetscConvEstOps {
  PetscErrorCode (*setfromoptions)(PetscConvEst);
  PetscErrorCode (*setup)(PetscConvEst);
  PetscErrorCode (*view)(PetscConvEst,PetscViewer);
  PetscErrorCode (*destroy)(PetscConvEst);
};

struct _p_PetscConvEst
{
  PETSCHEADER(struct _PetscConvEstOps);
  /* Inputs */
  DM                idm;  /* Initial grid */
  SNES              snes; /* Solver */
  PetscInt          Nr;   /* The number of refinements */
  PetscInt          Nf;   /* The number of fields in the DM */
  PetscErrorCode (**initGuess)(PetscInt, PetscReal, const PetscReal[], PetscInt, PetscScalar[], void *);
  PetscErrorCode (**exactSol)(PetscInt, PetscReal, const PetscReal[], PetscInt, PetscScalar[], void *);
  /* Outputs */
  PetscBool  monitor;
  PetscReal *errors;
};

#endif
