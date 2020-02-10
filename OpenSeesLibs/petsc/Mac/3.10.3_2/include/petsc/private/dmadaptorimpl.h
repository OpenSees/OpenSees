#if !defined(_DMADAPTORIMPL_H)
#define _DMADAPTORIMPL_H

#include <petscdmadaptor.h>
#include <petsc/private/petscimpl.h>

typedef struct _DMAdaptorOps *DMAdaptorOps;
struct _DMAdaptorOps {
  PetscErrorCode (*setfromoptions)(DMAdaptor);
  PetscErrorCode (*setup)(DMAdaptor);
  PetscErrorCode (*view)(DMAdaptor,PetscViewer);
  PetscErrorCode (*destroy)(DMAdaptor);
  PetscErrorCode (*transfersolution)(DMAdaptor,DM,Vec,DM,Vec,void*);
  PetscErrorCode (*computeerrorindicator)(DMAdaptor,PetscInt,PetscInt,const PetscScalar*,const PetscScalar*,const PetscFVCellGeom*,PetscReal*,void*);
};

struct _p_DMAdaptor {
  PETSCHEADER(struct _DMAdaptorOps);
  /* Inputs */
  DM                 idm;  /* Initial grid */
  SNES               snes; /* Solver */
  VecTagger          refineTag, coarsenTag; /* Criteria for adaptivity */
  /*   control */
  DMAdaptationCriterion adaptCriterion;
  PetscBool          femType;
  PetscInt           numSeq;            /* Number of sequential adaptations */
  PetscInt           Nadapt;            /* Target number of vertices */
  PetscReal          refinementFactor;  /* N_adapt = r^dim N_orig */
  PetscReal          h_min, h_max;      /* Min and max Hessian eigenvalues */
  /*   FVM support */
  PetscBool          computeGradient;
  DM                 cellDM, gradDM;
  Vec                cellGeom, faceGeom, cellGrad; /* Local vectors */
  const PetscScalar *cellGeomArray, *cellGradArray;
  /* Outputs */
  PetscBool          monitor;
  /* Auxiliary objects */
  PetscLimiter       limiter;
  PetscErrorCode  (**exactSol)(PetscInt, PetscReal, const PetscReal[], PetscInt, PetscScalar[], void *);
};

#endif
