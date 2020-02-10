/*
  DMFOREST, for parallel, hierarchically refined, distributed mesh problems.
*/
#if !defined(__PETSCDMFOREST_H)
#define __PETSCDMFOREST_H

#include <petscdm.h>

/*J
    DMForestTopology - String with the name of a PETSc DMFOREST base mesh topology. The topology is a string (e.g.
  "cube", "shell") and can be interpreted by subtypes of DMFOREST) to construct the base DM of a forest during
  DMSetUp().

   Level: beginner

.seealso: DMForestSetTopology(), DMForestGetTopology(), DMFOREST
J*/
typedef const char* DMForestTopology;

/* Just a name for the shape of the domain */
PETSC_EXTERN PetscErrorCode DMForestSetTopology(DM, DMForestTopology);
PETSC_EXTERN PetscErrorCode DMForestGetTopology(DM, DMForestTopology *);

/* this is the coarsest possible forest: can be any DM which we can
 * convert to a DMForest (right now: plex) */
PETSC_EXTERN PetscErrorCode DMForestSetBaseDM(DM, DM);
PETSC_EXTERN PetscErrorCode DMForestGetBaseDM(DM, DM *);
PETSC_EXTERN PetscErrorCode DMForestSetBaseCoordinateMapping(DM, PetscErrorCode(*)(DM,PetscInt,PetscInt,const PetscReal[],PetscReal[],void*),void*);
PETSC_EXTERN PetscErrorCode DMForestGetBaseCoordinateMapping(DM, PetscErrorCode(**)(DM,PetscInt,PetscInt,const PetscReal[],PetscReal[],void*),void*);

/* this is the forest from which we adapt */
PETSC_EXTERN PetscErrorCode DMForestSetAdaptivityForest(DM, DM);
PETSC_EXTERN PetscErrorCode DMForestGetAdaptivityForest(DM, DM *);

PETSC_EXTERN PetscErrorCode DMForestSetAdaptivityPurpose(DM, DMAdaptFlag);
PETSC_EXTERN PetscErrorCode DMForestGetAdaptivityPurpose(DM, DMAdaptFlag*);

/* what we consider adjacent, for the purposes of cell grading, overlap, etc. */
PETSC_EXTERN PetscErrorCode DMForestSetAdjacencyDimension(DM, PetscInt);
PETSC_EXTERN PetscErrorCode DMForestGetAdjacencyDimension(DM, PetscInt *);
PETSC_EXTERN PetscErrorCode DMForestSetAdjacencyCodimension(DM, PetscInt);
PETSC_EXTERN PetscErrorCode DMForestGetAdjacencyCodimension(DM, PetscInt *);

PETSC_EXTERN PetscErrorCode DMForestSetPartitionOverlap(DM, PetscInt);
PETSC_EXTERN PetscErrorCode DMForestGetPartitionOverlap(DM, PetscInt *);

PETSC_EXTERN PetscErrorCode DMForestSetMinimumRefinement(DM, PetscInt);
PETSC_EXTERN PetscErrorCode DMForestGetMinimumRefinement(DM, PetscInt *);

PETSC_EXTERN PetscErrorCode DMForestSetMaximumRefinement(DM, PetscInt);
PETSC_EXTERN PetscErrorCode DMForestGetMaximumRefinement(DM, PetscInt *);

PETSC_EXTERN PetscErrorCode DMForestSetInitialRefinement(DM, PetscInt);
PETSC_EXTERN PetscErrorCode DMForestGetInitialRefinement(DM, PetscInt *);

PETSC_EXTERN PetscErrorCode DMForestGetCellChart(DM, PetscInt *, PetscInt *);
PETSC_EXTERN PetscErrorCode DMForestGetCellSF(DM, PetscSF *);


/* flag each cell with an adaptivity count: should match the cell section */
PETSC_EXTERN PetscErrorCode DMForestSetAdaptivityLabel(DM, DMLabel);
PETSC_EXTERN PetscErrorCode DMForestGetAdaptivityLabel(DM, DMLabel*);

/*J
    DMForestAdaptivityStrategy - String with the name of a PETSc DMForest adaptivity strategy

   Level: intermediate

.seealso: DMForestSetType(), DMForest
J*/
typedef const char* DMForestAdaptivityStrategy;
#define DMFORESTADAPTALL "all"
#define DMFORESTADAPTANY "any"

/* how to combine: -flags         from multiple processes,
 *                 -coarsen flags from multiple children
 */
PETSC_EXTERN PetscErrorCode DMForestSetAdaptivityStrategy(DM, DMForestAdaptivityStrategy);
PETSC_EXTERN PetscErrorCode DMForestGetAdaptivityStrategy(DM, DMForestAdaptivityStrategy *);

PETSC_EXTERN PetscErrorCode DMForestSetComputeAdaptivitySF(DM, PetscBool);
PETSC_EXTERN PetscErrorCode DMForestGetComputeAdaptivitySF(DM, PetscBool *);

PETSC_EXTERN PetscErrorCode DMForestGetAdaptivitySF(DM, PetscSF *, PetscSF *);

PETSC_EXTERN PetscErrorCode DMForestGetAdaptivitySuccess(DM, PetscBool *);

PETSC_EXTERN PetscErrorCode DMForestTransferVec(DM, Vec, DM, Vec, PetscBool, PetscReal);

/* for a quadtree/octree mesh, this is the x:1 condition: 1 indicates a uniform mesh,
 *                                                        2 indicates typical 2:1,
 */
PETSC_EXTERN PetscErrorCode DMForestSetGradeFactor(DM, PetscInt);
PETSC_EXTERN PetscErrorCode DMForestGetGradeFactor(DM, PetscInt *);

/* weights for repartitioning */
PETSC_EXTERN PetscErrorCode DMForestSetCellWeights(DM, PetscReal[], PetscCopyMode);
PETSC_EXTERN PetscErrorCode DMForestGetCellWeights(DM, PetscReal *[]);

/* weight multiplier for refinement level: useful for sub-cycling time steps */
PETSC_EXTERN PetscErrorCode DMForestSetCellWeightFactor(DM, PetscReal);
PETSC_EXTERN PetscErrorCode DMForestGetCellWeightFactor(DM, PetscReal *);

/* this process's capacity when redistributing the cells */
PETSC_EXTERN PetscErrorCode DMForestSetWeightCapacity(DM, PetscReal);
PETSC_EXTERN PetscErrorCode DMForestGetWeightCapacity(DM, PetscReal *);

PETSC_EXTERN PetscErrorCode DMForestGetFineProjector(DM,Mat *);
PETSC_EXTERN PetscErrorCode DMForestGetCoarseRestrictor(DM,Mat *);

/* miscellaneous */
PETSC_EXTERN PetscErrorCode DMForestTemplate(DM,MPI_Comm,DM*);

/* type management */
PETSC_EXTERN PetscErrorCode DMForestRegisterType(DMType);
PETSC_EXTERN PetscErrorCode DMIsForest(DM,PetscBool*);

/* p4est */
PETSC_EXTERN PetscErrorCode DMP4estGetPartitionForCoarsening(DM,PetscBool *);
PETSC_EXTERN PetscErrorCode DMP4estSetPartitionForCoarsening(DM,PetscBool);
PETSC_EXTERN PetscErrorCode DMP8estGetPartitionForCoarsening(DM,PetscBool *);
PETSC_EXTERN PetscErrorCode DMP8estSetPartitionForCoarsening(DM,PetscBool);

#endif
