#if !defined(_FORESTIMPL_H)
#define _FORESTIMPL_H

#include <petscmat.h>       /*I      "petscmat.h"          I*/
#include <petscdmforest.h> /*I      "petscdmforest.h"    I*/
#include <petscbt.h>
#include <petsc/private/dmimpl.h>

typedef struct {
  PetscInt                   refct;
  void                       *data;
  PetscErrorCode             (*clearadaptivityforest)(DM);
  PetscErrorCode             (*getadaptivitysuccess)(DM,PetscBool*);
  PetscErrorCode             (*transfervec)(DM,Vec,DM,Vec,PetscBool,PetscReal);
  PetscErrorCode             (*createcellchart)(DM,PetscInt*,PetscInt*);
  PetscErrorCode             (*createcellsf)(DM,PetscSF*);
  PetscErrorCode             (*destroy)(DM);
  PetscErrorCode             (*ftemplate)(DM,DM);
  PetscBool                  setfromoptionscalled;
  PetscBool                  computeAdaptSF;
  PetscErrorCode             (*mapcoordinates)(DM,PetscInt,PetscInt,const PetscReal[],PetscReal[],void*);
  void                       *mapcoordinatesctx;
  DMForestTopology           topology;
  DM                         base;
  DM                         adapt;
  DMAdaptFlag                adaptPurpose;
  PetscInt                   adjDim;
  PetscInt                   overlap;
  PetscInt                   minRefinement;
  PetscInt                   maxRefinement;
  PetscInt                   initRefinement;
  PetscInt                   cStart;
  PetscInt                   cEnd;
  PetscSF                    cellSF;
  PetscSF                    preCoarseToFine;
  PetscSF                    coarseToPreFine;
  DMLabel                    adaptLabel;
  DMForestAdaptivityStrategy adaptStrategy;
  PetscInt                   gradeFactor;
  PetscReal                  *cellWeights;
  PetscCopyMode              cellWeightsCopyMode;
  PetscReal                  weightsFactor;
  PetscReal                  weightCapacity;
} DM_Forest;

PETSC_EXTERN PetscErrorCode DMCreate_Forest(DM);
PETSC_EXTERN PetscErrorCode DMClone_Forest(DM,DM*);
PETSC_EXTERN PetscErrorCode DMSetFromOptions_Forest(PetscOptionItems*,DM);

#endif /* _FORESTIMPL_H */
