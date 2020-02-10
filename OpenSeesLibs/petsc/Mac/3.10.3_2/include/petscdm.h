/*
      Objects to manage the interactions between the mesh data structures and the algebraic objects
*/
#if !defined(__PETSCDM_H)
#define __PETSCDM_H
#include <petscmat.h>
#include <petscdmtypes.h>
#include <petscfetypes.h>
#include <petscdstypes.h>
#include <petscdmlabel.h>

PETSC_EXTERN PetscErrorCode DMInitializePackage(void);

PETSC_EXTERN PetscClassId DM_CLASSID;

#define DMLOCATEPOINT_POINT_NOT_FOUND -367

/*J
    DMType - String with the name of a PETSc DM

   Level: beginner

.seealso: DMSetType(), DM
J*/
typedef const char* DMType;
#define DMDA        "da"
#define DMCOMPOSITE "composite"
#define DMSLICED    "sliced"
#define DMSHELL     "shell"
#define DMPLEX      "plex"
#define DMREDUNDANT "redundant"
#define DMPATCH     "patch"
#define DMMOAB      "moab"
#define DMNETWORK   "network"
#define DMFOREST    "forest"
#define DMP4EST     "p4est"
#define DMP8EST     "p8est"
#define DMSWARM     "swarm"

PETSC_EXTERN const char *const DMBoundaryTypes[];
PETSC_EXTERN PetscFunctionList DMList;
PETSC_EXTERN PetscErrorCode DMCreate(MPI_Comm,DM*);
PETSC_EXTERN PetscErrorCode DMClone(DM,DM*);
PETSC_EXTERN PetscErrorCode DMSetType(DM, DMType);
PETSC_EXTERN PetscErrorCode DMGetType(DM, DMType *);
PETSC_EXTERN PetscErrorCode DMRegister(const char[],PetscErrorCode (*)(DM));
PETSC_EXTERN PetscErrorCode DMRegisterDestroy(void);

PETSC_EXTERN PetscErrorCode DMView(DM,PetscViewer);
PETSC_EXTERN PetscErrorCode DMLoad(DM,PetscViewer);
PETSC_EXTERN PetscErrorCode DMDestroy(DM*);
PETSC_EXTERN PetscErrorCode DMCreateGlobalVector(DM,Vec*);
PETSC_EXTERN PetscErrorCode DMCreateLocalVector(DM,Vec*);
PETSC_EXTERN PetscErrorCode DMGetLocalVector(DM,Vec *);
PETSC_EXTERN PetscErrorCode DMRestoreLocalVector(DM,Vec *);
PETSC_EXTERN PetscErrorCode DMGetGlobalVector(DM,Vec *);
PETSC_EXTERN PetscErrorCode DMRestoreGlobalVector(DM,Vec *);
PETSC_EXTERN PetscErrorCode DMClearGlobalVectors(DM);
PETSC_EXTERN PetscErrorCode DMClearLocalVectors(DM);
PETSC_EXTERN PetscErrorCode DMHasNamedGlobalVector(DM,const char*,PetscBool*);
PETSC_EXTERN PetscErrorCode DMGetNamedGlobalVector(DM,const char*,Vec*);
PETSC_EXTERN PetscErrorCode DMRestoreNamedGlobalVector(DM,const char*,Vec*);
PETSC_EXTERN PetscErrorCode DMHasNamedLocalVector(DM,const char*,PetscBool*);
PETSC_EXTERN PetscErrorCode DMGetNamedLocalVector(DM,const char*,Vec*);
PETSC_EXTERN PetscErrorCode DMRestoreNamedLocalVector(DM,const char*,Vec*);
PETSC_EXTERN PetscErrorCode DMGetLocalToGlobalMapping(DM,ISLocalToGlobalMapping*);
PETSC_EXTERN PetscErrorCode DMCreateFieldIS(DM,PetscInt*,char***,IS**);
PETSC_EXTERN PetscErrorCode DMGetBlockSize(DM,PetscInt*);
PETSC_EXTERN PetscErrorCode DMCreateColoring(DM,ISColoringType,ISColoring*);
PETSC_EXTERN PetscErrorCode DMCreateMatrix(DM,Mat*);
PETSC_EXTERN PetscErrorCode DMSetMatrixPreallocateOnly(DM,PetscBool);
PETSC_EXTERN PetscErrorCode DMSetMatrixStructureOnly(DM,PetscBool);
PETSC_EXTERN PetscErrorCode DMCreateInterpolation(DM,DM,Mat*,Vec*);
PETSC_EXTERN PetscErrorCode DMCreateRestriction(DM,DM,Mat*);
PETSC_EXTERN PetscErrorCode DMCreateInjection(DM,DM,Mat*);
PETSC_EXTERN PetscErrorCode DMCreateMassMatrix(DM,DM,Mat*);
PETSC_EXTERN PetscErrorCode DMGetWorkArray(DM,PetscInt,MPI_Datatype,void*);
PETSC_EXTERN PetscErrorCode DMRestoreWorkArray(DM,PetscInt,MPI_Datatype,void*);
PETSC_EXTERN PetscErrorCode DMRefine(DM,MPI_Comm,DM*);
PETSC_EXTERN PetscErrorCode DMCoarsen(DM,MPI_Comm,DM*);
PETSC_EXTERN PetscErrorCode DMGetCoarseDM(DM,DM*);
PETSC_EXTERN PetscErrorCode DMSetCoarseDM(DM,DM);
PETSC_EXTERN PetscErrorCode DMGetFineDM(DM,DM*);
PETSC_EXTERN PetscErrorCode DMSetFineDM(DM,DM);
PETSC_EXTERN PetscErrorCode DMRefineHierarchy(DM,PetscInt,DM[]);
PETSC_EXTERN PetscErrorCode DMCoarsenHierarchy(DM,PetscInt,DM[]);
PETSC_EXTERN PetscErrorCode DMCoarsenHookAdd(DM,PetscErrorCode (*)(DM,DM,void*),PetscErrorCode (*)(DM,Mat,Vec,Mat,DM,void*),void*);
PETSC_EXTERN PetscErrorCode DMCoarsenHookRemove(DM,PetscErrorCode (*)(DM,DM,void*),PetscErrorCode (*)(DM,Mat,Vec,Mat,DM,void*),void*);
PETSC_EXTERN PetscErrorCode DMRefineHookAdd(DM,PetscErrorCode (*)(DM,DM,void*),PetscErrorCode (*)(DM,Mat,DM,void*),void*);
PETSC_EXTERN PetscErrorCode DMRefineHookRemove(DM,PetscErrorCode (*)(DM,DM,void*),PetscErrorCode (*)(DM,Mat,DM,void*),void*);
PETSC_EXTERN PetscErrorCode DMRestrict(DM,Mat,Vec,Mat,DM);
PETSC_EXTERN PetscErrorCode DMInterpolate(DM,Mat,DM);
PETSC_EXTERN PetscErrorCode DMSetFromOptions(DM);
PETSC_STATIC_INLINE PetscErrorCode DMViewFromOptions(DM A,PetscObject obj,const char name[]) {return PetscObjectViewFromOptions((PetscObject)A,obj,name);}

PETSC_EXTERN PetscErrorCode DMAdaptLabel(DM,DMLabel,DM*);
PETSC_EXTERN PetscErrorCode DMAdaptMetric(DM, Vec, DMLabel, DM *);

PETSC_EXTERN PetscErrorCode DMSetUp(DM);
PETSC_EXTERN PetscErrorCode DMCreateInterpolationScale(DM,DM,Mat,Vec*);
PETSC_EXTERN PetscErrorCode DMCreateAggregates(DM,DM,Mat*);
PETSC_EXTERN PetscErrorCode DMGlobalToLocalHookAdd(DM,PetscErrorCode (*)(DM,Vec,InsertMode,Vec,void*),PetscErrorCode (*)(DM,Vec,InsertMode,Vec,void*),void*);
PETSC_EXTERN PetscErrorCode DMLocalToGlobalHookAdd(DM,PetscErrorCode (*)(DM,Vec,InsertMode,Vec,void*),PetscErrorCode (*)(DM,Vec,InsertMode,Vec,void*),void*);
PETSC_EXTERN PetscErrorCode DMGlobalToLocalBegin(DM,Vec,InsertMode,Vec);
PETSC_EXTERN PetscErrorCode DMGlobalToLocalEnd(DM,Vec,InsertMode,Vec);
PETSC_EXTERN PetscErrorCode DMLocalToGlobalBegin(DM,Vec,InsertMode,Vec);
PETSC_EXTERN PetscErrorCode DMLocalToGlobalEnd(DM,Vec,InsertMode,Vec);
PETSC_EXTERN PetscErrorCode DMLocalToLocalBegin(DM,Vec,InsertMode,Vec);
PETSC_EXTERN PetscErrorCode DMLocalToLocalEnd(DM,Vec,InsertMode,Vec);
PETSC_EXTERN PetscErrorCode DMConvert(DM,DMType,DM*);

/* Topology support */
PETSC_EXTERN PetscErrorCode DMGetDimension(DM,PetscInt*);
PETSC_EXTERN PetscErrorCode DMSetDimension(DM,PetscInt);
PETSC_EXTERN PetscErrorCode DMGetDimPoints(DM,PetscInt,PetscInt*,PetscInt*);
PETSC_EXTERN PetscErrorCode DMGetUseNatural(DM,PetscBool*);
PETSC_EXTERN PetscErrorCode DMSetUseNatural(DM,PetscBool);

/* Coordinate support */
PETSC_EXTERN PetscErrorCode DMGetCoordinateDM(DM,DM*);
PETSC_EXTERN PetscErrorCode DMSetCoordinateDM(DM,DM);
PETSC_EXTERN PetscErrorCode DMGetCoordinateDim(DM,PetscInt*);
PETSC_EXTERN PetscErrorCode DMSetCoordinateDim(DM,PetscInt);
PETSC_EXTERN PetscErrorCode DMGetCoordinateSection(DM,PetscSection*);
PETSC_EXTERN PetscErrorCode DMSetCoordinateSection(DM,PetscInt,PetscSection);
PETSC_EXTERN PetscErrorCode DMGetCoordinates(DM,Vec*);
PETSC_EXTERN PetscErrorCode DMSetCoordinates(DM,Vec);
PETSC_EXTERN PetscErrorCode DMGetCoordinatesLocal(DM,Vec*);
PETSC_EXTERN PetscErrorCode DMSetCoordinatesLocal(DM,Vec);
PETSC_EXTERN PetscErrorCode DMLocatePoints(DM,Vec,DMPointLocationType,PetscSF*);
PETSC_EXTERN PetscErrorCode DMGetPeriodicity(DM,PetscBool*,const PetscReal**,const PetscReal**,const DMBoundaryType**);
PETSC_EXTERN PetscErrorCode DMSetPeriodicity(DM,PetscBool,const PetscReal[],const PetscReal[],const DMBoundaryType[]);
PETSC_EXTERN PetscErrorCode DMLocalizeCoordinate(DM, const PetscScalar[], PetscBool, PetscScalar[]);
PETSC_EXTERN PetscErrorCode DMLocalizeCoordinates(DM);
PETSC_EXTERN PetscErrorCode DMGetCoordinatesLocalized(DM,PetscBool*);
PETSC_EXTERN PetscErrorCode DMGetNeighbors(DM,PetscInt*,const PetscMPIInt**);
PETSC_EXTERN PetscErrorCode DMGetCoordinateField(DM,DMField*);
PETSC_EXTERN PetscErrorCode DMSetCoordinateField(DM,DMField);

/* block hook interface */
PETSC_EXTERN PetscErrorCode DMSubDomainHookAdd(DM,PetscErrorCode (*)(DM,DM,void*),PetscErrorCode (*)(DM,VecScatter,VecScatter,DM,void*),void*);
PETSC_EXTERN PetscErrorCode DMSubDomainHookRemove(DM,PetscErrorCode (*)(DM,DM,void*),PetscErrorCode (*)(DM,VecScatter,VecScatter,DM,void*),void*);
PETSC_EXTERN PetscErrorCode DMSubDomainRestrict(DM,VecScatter,VecScatter,DM);

PETSC_EXTERN PetscErrorCode DMSetOptionsPrefix(DM,const char []);
PETSC_EXTERN PetscErrorCode DMAppendOptionsPrefix(DM,const char []);
PETSC_EXTERN PetscErrorCode DMGetOptionsPrefix(DM,const char*[]);
PETSC_EXTERN PetscErrorCode DMSetVecType(DM,VecType);
PETSC_EXTERN PetscErrorCode DMGetVecType(DM,VecType*);
PETSC_EXTERN PetscErrorCode DMSetMatType(DM,MatType);
PETSC_EXTERN PetscErrorCode DMGetMatType(DM,MatType*);
PETSC_EXTERN PetscErrorCode DMSetISColoringType(DM,ISColoringType);
PETSC_EXTERN PetscErrorCode DMGetISColoringType(DM,ISColoringType*);
PETSC_EXTERN PetscErrorCode DMSetApplicationContext(DM,void*);
PETSC_EXTERN PetscErrorCode DMSetApplicationContextDestroy(DM,PetscErrorCode (*)(void**));
PETSC_EXTERN PetscErrorCode DMGetApplicationContext(DM,void*);
PETSC_EXTERN PetscErrorCode DMSetVariableBounds(DM,PetscErrorCode (*)(DM,Vec,Vec));
PETSC_EXTERN PetscErrorCode DMHasVariableBounds(DM,PetscBool *);
PETSC_EXTERN PetscErrorCode DMHasColoring(DM,PetscBool *);
PETSC_EXTERN PetscErrorCode DMHasCreateRestriction(DM,PetscBool *);
PETSC_EXTERN PetscErrorCode DMHasCreateInjection(DM,PetscBool *);
PETSC_EXTERN PetscErrorCode DMComputeVariableBounds(DM,Vec,Vec);

PETSC_EXTERN PetscErrorCode DMCreateSubDM(DM, PetscInt, const PetscInt[], IS *, DM *);
PETSC_EXTERN PetscErrorCode DMCreateSuperDM(DM[], PetscInt, IS **, DM *);
PETSC_EXTERN PetscErrorCode DMCreateFieldDecomposition(DM,PetscInt*,char***,IS**,DM**);
PETSC_EXTERN PetscErrorCode DMCreateDomainDecomposition(DM,PetscInt*,char***,IS**,IS**,DM**);
PETSC_EXTERN PetscErrorCode DMCreateDomainDecompositionScatters(DM,PetscInt,DM*,VecScatter**,VecScatter**,VecScatter**);

PETSC_EXTERN PetscErrorCode DMGetRefineLevel(DM,PetscInt*);
PETSC_EXTERN PetscErrorCode DMSetRefineLevel(DM,PetscInt);
PETSC_EXTERN PetscErrorCode DMGetCoarsenLevel(DM,PetscInt*);
PETSC_EXTERN PetscErrorCode DMFinalizePackage(void);

PETSC_EXTERN PetscErrorCode VecGetDM(Vec, DM*);
PETSC_EXTERN PetscErrorCode VecSetDM(Vec, DM);
PETSC_EXTERN PetscErrorCode MatGetDM(Mat, DM*);
PETSC_EXTERN PetscErrorCode MatSetDM(Mat, DM);
PETSC_EXTERN PetscErrorCode MatFDColoringUseDM(Mat,MatFDColoring);

typedef struct NLF_DAAD* NLF;

#define DM_FILE_CLASSID 1211221

/* FEM support */
PETSC_EXTERN PetscErrorCode DMPrintCellVector(PetscInt, const char [], PetscInt, const PetscScalar []);
PETSC_EXTERN PetscErrorCode DMPrintCellMatrix(PetscInt, const char [], PetscInt, PetscInt, const PetscScalar []);
PETSC_EXTERN PetscErrorCode DMPrintLocalVec(DM, const char [], PetscReal, Vec);

PETSC_EXTERN PetscErrorCode DMGetSection(DM, PetscSection *);
PETSC_EXTERN PetscErrorCode DMSetSection(DM, PetscSection);
PETSC_EXTERN PetscErrorCode DMGetDefaultConstraints(DM, PetscSection *, Mat *);
PETSC_EXTERN PetscErrorCode DMSetDefaultConstraints(DM, PetscSection, Mat);
PETSC_EXTERN PetscErrorCode DMGetGlobalSection(DM, PetscSection *);
PETSC_EXTERN PetscErrorCode DMSetGlobalSection(DM, PetscSection);
PETSC_EXTERN PetscErrorCode DMGetDefaultSF(DM, PetscSF *);
PETSC_EXTERN PetscErrorCode DMSetDefaultSF(DM, PetscSF);
PETSC_EXTERN PetscErrorCode DMCreateDefaultSF(DM, PetscSection, PetscSection);
PETSC_EXTERN PetscErrorCode DMGetPointSF(DM, PetscSF *);
PETSC_EXTERN PetscErrorCode DMSetPointSF(DM, PetscSF);

PETSC_STATIC_INLINE PetscErrorCode DMGetDefaultSection(DM dm, PetscSection *s) {return DMGetSection(dm,s);}
PETSC_STATIC_INLINE PetscErrorCode DMSetDefaultSection(DM dm, PetscSection s) {return DMSetSection(dm,s);}
PETSC_STATIC_INLINE PetscErrorCode DMGetDefaultGlobalSection(DM dm, PetscSection *s) {return DMGetGlobalSection(dm,s);}
PETSC_STATIC_INLINE PetscErrorCode DMSetDefaultGlobalSection(DM dm, PetscSection s) {return DMSetGlobalSection(dm,s);}


PETSC_EXTERN PetscErrorCode DMGetOutputDM(DM, DM *);
PETSC_EXTERN PetscErrorCode DMGetOutputSequenceNumber(DM, PetscInt *, PetscReal *);
PETSC_EXTERN PetscErrorCode DMSetOutputSequenceNumber(DM, PetscInt, PetscReal);
PETSC_EXTERN PetscErrorCode DMOutputSequenceLoad(DM, PetscViewer, const char *, PetscInt, PetscReal *);

PETSC_EXTERN PetscErrorCode DMGetDS(DM, PetscDS *);
PETSC_EXTERN PetscErrorCode DMSetDS(DM, PetscDS);
PETSC_EXTERN PetscErrorCode DMGetNumFields(DM, PetscInt *);
PETSC_EXTERN PetscErrorCode DMSetNumFields(DM, PetscInt);
PETSC_EXTERN PetscErrorCode DMGetField(DM, PetscInt, PetscObject *);
PETSC_EXTERN PetscErrorCode DMSetField(DM, PetscInt, PetscObject);

struct _DMInterpolationInfo {
  MPI_Comm   comm;
  PetscInt   dim;    /*1 The spatial dimension of points */
  PetscInt   nInput; /* The number of input points */
  PetscReal *points; /* The input point coordinates */
  PetscInt  *cells;  /* The cell containing each point */
  PetscInt   n;      /* The number of local points */
  Vec        coords; /* The point coordinates */
  PetscInt   dof;    /* The number of components to interpolate */
};
typedef struct _DMInterpolationInfo *DMInterpolationInfo;

PETSC_EXTERN PetscErrorCode DMInterpolationCreate(MPI_Comm, DMInterpolationInfo *);
PETSC_EXTERN PetscErrorCode DMInterpolationSetDim(DMInterpolationInfo, PetscInt);
PETSC_EXTERN PetscErrorCode DMInterpolationGetDim(DMInterpolationInfo, PetscInt *);
PETSC_EXTERN PetscErrorCode DMInterpolationSetDof(DMInterpolationInfo, PetscInt);
PETSC_EXTERN PetscErrorCode DMInterpolationGetDof(DMInterpolationInfo, PetscInt *);
PETSC_EXTERN PetscErrorCode DMInterpolationAddPoints(DMInterpolationInfo, PetscInt, PetscReal[]);
PETSC_EXTERN PetscErrorCode DMInterpolationSetUp(DMInterpolationInfo, DM, PetscBool);
PETSC_EXTERN PetscErrorCode DMInterpolationGetCoordinates(DMInterpolationInfo, Vec *);
PETSC_EXTERN PetscErrorCode DMInterpolationGetVector(DMInterpolationInfo, Vec *);
PETSC_EXTERN PetscErrorCode DMInterpolationRestoreVector(DMInterpolationInfo, Vec *);
PETSC_EXTERN PetscErrorCode DMInterpolationEvaluate(DMInterpolationInfo, DM, Vec, Vec);
PETSC_EXTERN PetscErrorCode DMInterpolationDestroy(DMInterpolationInfo *);

PETSC_EXTERN PetscErrorCode DMCreateLabel(DM, const char []);
PETSC_EXTERN PetscErrorCode DMGetLabelValue(DM, const char[], PetscInt, PetscInt *);
PETSC_EXTERN PetscErrorCode DMSetLabelValue(DM, const char[], PetscInt, PetscInt);
PETSC_EXTERN PetscErrorCode DMClearLabelValue(DM, const char[], PetscInt, PetscInt);
PETSC_EXTERN PetscErrorCode DMGetLabelSize(DM, const char[], PetscInt *);
PETSC_EXTERN PetscErrorCode DMGetLabelIdIS(DM, const char[], IS *);
PETSC_EXTERN PetscErrorCode DMGetStratumSize(DM, const char [], PetscInt, PetscInt *);
PETSC_EXTERN PetscErrorCode DMGetStratumIS(DM, const char [], PetscInt, IS *);
PETSC_EXTERN PetscErrorCode DMSetStratumIS(DM, const char [], PetscInt, IS);
PETSC_EXTERN PetscErrorCode DMClearLabelStratum(DM, const char[], PetscInt);
PETSC_EXTERN PetscErrorCode DMGetLabelOutput(DM, const char[], PetscBool *);
PETSC_EXTERN PetscErrorCode DMSetLabelOutput(DM, const char[], PetscBool);

PETSC_EXTERN PetscErrorCode DMGetNumLabels(DM, PetscInt *);
PETSC_EXTERN PetscErrorCode DMGetLabelName(DM, PetscInt, const char **);
PETSC_EXTERN PetscErrorCode DMHasLabel(DM, const char [], PetscBool *);
PETSC_EXTERN PetscErrorCode DMGetLabel(DM, const char *, DMLabel *);
PETSC_EXTERN PetscErrorCode DMGetLabelByNum(DM, PetscInt, DMLabel *);
PETSC_EXTERN PetscErrorCode DMAddLabel(DM, DMLabel);
PETSC_EXTERN PetscErrorCode DMRemoveLabel(DM, const char [], DMLabel *);
PETSC_EXTERN PetscErrorCode DMCopyLabels(DM, DM);

PETSC_EXTERN PetscErrorCode DMAddBoundary(DM, DMBoundaryConditionType, const char[], const char[], PetscInt, PetscInt, const PetscInt *, void (*)(void), PetscInt, const PetscInt *, void *);
PETSC_EXTERN PetscErrorCode DMGetNumBoundary(DM, PetscInt *);
PETSC_EXTERN PetscErrorCode DMGetBoundary(DM, PetscInt, DMBoundaryConditionType *, const char **, const char **, PetscInt *, PetscInt *, const PetscInt **, void (**)(void), PetscInt *, const PetscInt **, void **);
PETSC_EXTERN PetscErrorCode DMIsBoundaryPoint(DM, PetscInt, PetscBool *);
PETSC_EXTERN PetscErrorCode DMCopyBoundary(DM, DM);

PETSC_EXTERN PetscErrorCode DMProjectFunction(DM,PetscReal,PetscErrorCode(**)(PetscInt,PetscReal,const PetscReal[],PetscInt,PetscScalar *,void *),void**,InsertMode,Vec);
PETSC_EXTERN PetscErrorCode DMProjectFunctionLocal(DM,PetscReal,PetscErrorCode(**)(PetscInt,PetscReal,const PetscReal[],PetscInt,PetscScalar *,void *),void**,InsertMode,Vec);
PETSC_EXTERN PetscErrorCode DMProjectFunctionLabel(DM, PetscReal, DMLabel, PetscInt, const PetscInt[], PetscInt, const PetscInt[], PetscErrorCode (**)(PetscInt, PetscReal, const PetscReal [], PetscInt, PetscScalar *, void *), void **, InsertMode, Vec);
PETSC_EXTERN PetscErrorCode DMProjectFunctionLabelLocal(DM,PetscReal,DMLabel,PetscInt,const PetscInt[],PetscInt,const PetscInt[],PetscErrorCode(**)(PetscInt,PetscReal,const PetscReal[],PetscInt,PetscScalar *,void *),void **,InsertMode,Vec);
PETSC_EXTERN PetscErrorCode DMProjectFieldLocal(DM,PetscReal,Vec,void (**)(PetscInt,PetscInt,PetscInt,const PetscInt[],const PetscInt[],const PetscScalar[],const PetscScalar[],const PetscScalar[],const PetscInt[],const PetscInt[],const PetscScalar[],const PetscScalar[],const PetscScalar[],PetscReal,const PetscReal[],PetscInt,const PetscScalar[],PetscScalar[]),InsertMode,Vec);
PETSC_EXTERN PetscErrorCode DMProjectFieldLabelLocal(DM,PetscReal,DMLabel,PetscInt,const PetscInt[],PetscInt,const PetscInt[],Vec,void (**)(PetscInt,PetscInt,PetscInt,const PetscInt[],const PetscInt[],const PetscScalar[],const PetscScalar[],const PetscScalar[],const PetscInt[],const PetscInt[],const PetscScalar[],const PetscScalar[],const PetscScalar[],PetscReal,const PetscReal[],PetscInt,const PetscScalar[],PetscScalar[]),InsertMode,Vec);
PETSC_EXTERN PetscErrorCode DMComputeL2Diff(DM,PetscReal,PetscErrorCode(**)(PetscInt,PetscReal,const PetscReal[],PetscInt,PetscScalar *,void *),void **,Vec,PetscReal *);
PETSC_EXTERN PetscErrorCode DMComputeL2GradientDiff(DM, PetscReal, PetscErrorCode (**)(PetscInt, PetscReal, const PetscReal [], const PetscReal [], PetscInt, PetscScalar *, void *), void **, Vec, const PetscReal [], PetscReal *);
PETSC_EXTERN PetscErrorCode DMComputeL2FieldDiff(DM,PetscReal,PetscErrorCode(**)(PetscInt,PetscReal,const PetscReal[],PetscInt,PetscScalar *,void *),void **,Vec,PetscReal *);

PETSC_EXTERN PetscErrorCode DMGetNullSpaceConstructor(DM, PetscInt, PetscErrorCode (**)(DM, PetscInt, MatNullSpace *));
PETSC_EXTERN PetscErrorCode DMSetNullSpaceConstructor(DM, PetscInt, PetscErrorCode (*)(DM, PetscInt, MatNullSpace *));

PETSC_EXTERN PetscErrorCode DMGetCompatibility(DM,DM,PetscBool*,PetscBool*);
#endif
