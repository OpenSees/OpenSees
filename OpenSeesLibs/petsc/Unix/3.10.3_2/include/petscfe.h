/*
      Objects which encapsulate finite element spaces and operations
*/
#if !defined(__PETSCFE_H)
#define __PETSCFE_H
#include <petscdm.h>
#include <petscdt.h>
#include <petscfetypes.h>
#include <petscdstypes.h>

typedef struct _n_PetscFEGeom {
  const PetscReal *xi;
  PetscReal *v;
  PetscReal *J;
  PetscReal *invJ;
  PetscReal *detJ;
  PetscReal *n;
  PetscInt  (*face)[2];
  PetscReal *suppInvJ[2];
  PetscInt  dim;
  PetscInt  dimEmbed;
  PetscInt  numCells;
  PetscInt  numPoints;
  PetscBool isAffine;
} PetscFEGeom;

PETSC_EXTERN PetscErrorCode PetscFEInitializePackage(void);

PETSC_EXTERN PetscClassId PETSCSPACE_CLASSID;

/*J
  PetscSpaceType - String with the name of a PETSc linear space

  Level: beginner

.seealso: PetscSpaceSetType(), PetscSpace
J*/
typedef const char *PetscSpaceType;
#define PETSCSPACEPOLYNOMIAL "poly"
#define PETSCSPACETENSOR     "tensor"
#define PETSCSPACEPOINT      "point"
#define PETSCSPACESUBSPACE   "subspace"

PETSC_EXTERN PetscFunctionList PetscSpaceList;
PETSC_EXTERN PetscErrorCode PetscSpaceCreate(MPI_Comm, PetscSpace *);
PETSC_EXTERN PetscErrorCode PetscSpaceDestroy(PetscSpace *);
PETSC_EXTERN PetscErrorCode PetscSpaceSetType(PetscSpace, PetscSpaceType);
PETSC_EXTERN PetscErrorCode PetscSpaceGetType(PetscSpace, PetscSpaceType *);
PETSC_EXTERN PetscErrorCode PetscSpaceSetUp(PetscSpace);
PETSC_EXTERN PetscErrorCode PetscSpaceSetFromOptions(PetscSpace);
PETSC_STATIC_INLINE PetscErrorCode PetscSpaceViewFromOptions(PetscSpace A,PetscObject B,const char name[]) {return PetscObjectViewFromOptions((PetscObject)A,B,name);}

PETSC_EXTERN PetscErrorCode PetscSpaceView(PetscSpace,PetscViewer);
PETSC_EXTERN PetscErrorCode PetscSpaceRegister(const char [], PetscErrorCode (*)(PetscSpace));
PETSC_EXTERN PetscErrorCode PetscSpaceRegisterDestroy(void);

PETSC_EXTERN PetscErrorCode PetscSpaceGetDimension(PetscSpace, PetscInt *);
PETSC_EXTERN PetscErrorCode PetscSpaceSetNumComponents(PetscSpace, PetscInt);
PETSC_EXTERN PetscErrorCode PetscSpaceGetNumComponents(PetscSpace, PetscInt *);
PETSC_EXTERN PetscErrorCode PetscSpaceSetNumVariables(PetscSpace, PetscInt);
PETSC_EXTERN PetscErrorCode PetscSpaceGetNumVariables(PetscSpace, PetscInt *);
PETSC_EXTERN PetscErrorCode PetscSpaceSetDegree(PetscSpace, PetscInt, PetscInt);
PETSC_EXTERN PetscErrorCode PetscSpaceGetDegree(PetscSpace, PetscInt *, PetscInt *);
PETSC_EXTERN PetscErrorCode PetscSpaceEvaluate(PetscSpace, PetscInt, const PetscReal[], PetscReal[], PetscReal[], PetscReal[]);
PETSC_EXTERN PetscErrorCode PetscSpaceGetHeightSubspace(PetscSpace, PetscInt, PetscSpace *);

PETSC_EXTERN PetscErrorCode PetscSpacePolynomialSetSymmetric(PetscSpace, PetscBool);
PETSC_EXTERN PetscErrorCode PetscSpacePolynomialGetSymmetric(PetscSpace, PetscBool *);
PETSC_EXTERN PetscErrorCode PetscSpacePolynomialSetTensor(PetscSpace, PetscBool);
PETSC_EXTERN PetscErrorCode PetscSpacePolynomialGetTensor(PetscSpace, PetscBool *);

PETSC_EXTERN PetscErrorCode PetscSpaceTensorSetNumSubspaces(PetscSpace, PetscInt);
PETSC_EXTERN PetscErrorCode PetscSpaceTensorGetNumSubspaces(PetscSpace, PetscInt *);
PETSC_EXTERN PetscErrorCode PetscSpaceTensorSetSubspace(PetscSpace, PetscInt, PetscSpace);
PETSC_EXTERN PetscErrorCode PetscSpaceTensorGetSubspace(PetscSpace, PetscInt, PetscSpace *);

PETSC_EXTERN PetscErrorCode PetscSpacePointGetPoints(PetscSpace, PetscQuadrature *);
PETSC_EXTERN PetscErrorCode PetscSpacePointSetPoints(PetscSpace, PetscQuadrature);

PETSC_EXTERN PetscErrorCode PetscSpaceCreateSubspace(PetscSpace, PetscDualSpace, PetscReal *, PetscReal *, PetscReal *, PetscReal *, PetscCopyMode, PetscSpace *);

PETSC_EXTERN PetscClassId PETSCDUALSPACE_CLASSID;

/*J
  PetscDualSpaceType - String with the name of a PETSc dual space

  Level: beginner

.seealso: PetscDualSpaceSetType(), PetscDualSpace
J*/
typedef const char *PetscDualSpaceType;
#define PETSCDUALSPACELAGRANGE "lagrange"
#define PETSCDUALSPACESIMPLE   "simple"

PETSC_EXTERN PetscFunctionList PetscDualSpaceList;
PETSC_EXTERN PetscErrorCode PetscDualSpaceCreate(MPI_Comm, PetscDualSpace *);
PETSC_EXTERN PetscErrorCode PetscDualSpaceDestroy(PetscDualSpace *);
PETSC_EXTERN PetscErrorCode PetscDualSpaceDuplicate(PetscDualSpace, PetscDualSpace *);
PETSC_EXTERN PetscErrorCode PetscDualSpaceSetType(PetscDualSpace, PetscDualSpaceType);
PETSC_EXTERN PetscErrorCode PetscDualSpaceGetType(PetscDualSpace, PetscDualSpaceType *);
PETSC_EXTERN PetscErrorCode PetscDualSpaceGetNumDof(PetscDualSpace, const PetscInt **);
PETSC_EXTERN PetscErrorCode PetscDualSpaceCreateSection(PetscDualSpace, PetscSection *);
PETSC_EXTERN PetscErrorCode PetscDualSpaceSetUp(PetscDualSpace);
PETSC_EXTERN PetscErrorCode PetscDualSpaceSetFromOptions(PetscDualSpace);
PETSC_STATIC_INLINE PetscErrorCode PetscDualSpaceViewFromOptions(PetscDualSpace A,PetscObject B,const char name[]) {return PetscObjectViewFromOptions((PetscObject)A,B,name);}

PETSC_EXTERN PetscErrorCode PetscDualSpaceView(PetscDualSpace,PetscViewer);
PETSC_EXTERN PetscErrorCode PetscDualSpaceRegister(const char [], PetscErrorCode (*)(PetscDualSpace));
PETSC_EXTERN PetscErrorCode PetscDualSpaceRegisterDestroy(void);

PETSC_EXTERN PetscErrorCode PetscDualSpaceGetDimension(PetscDualSpace, PetscInt *);
PETSC_EXTERN PetscErrorCode PetscDualSpaceSetNumComponents(PetscDualSpace, PetscInt);
PETSC_EXTERN PetscErrorCode PetscDualSpaceGetNumComponents(PetscDualSpace, PetscInt *);
PETSC_EXTERN PetscErrorCode PetscDualSpaceSetOrder(PetscDualSpace, PetscInt);
PETSC_EXTERN PetscErrorCode PetscDualSpaceGetOrder(PetscDualSpace, PetscInt *);
PETSC_EXTERN PetscErrorCode PetscDualSpaceSetDM(PetscDualSpace, DM);
PETSC_EXTERN PetscErrorCode PetscDualSpaceGetDM(PetscDualSpace, DM *);
PETSC_EXTERN PetscErrorCode PetscDualSpaceGetFunctional(PetscDualSpace, PetscInt, PetscQuadrature *);
PETSC_EXTERN PetscErrorCode PetscDualSpaceCreateReferenceCell(PetscDualSpace, PetscInt, PetscBool, DM *);
PETSC_EXTERN PetscErrorCode PetscDualSpaceGetSymmetries(PetscDualSpace, const PetscInt ****, const PetscScalar ****);

PETSC_EXTERN PetscErrorCode PetscFEGeomCreate(PetscQuadrature,PetscInt,PetscInt,PetscBool,PetscFEGeom**);
PETSC_EXTERN PetscErrorCode PetscFEGeomGetChunk(PetscFEGeom*,PetscInt,PetscInt,PetscFEGeom**);
PETSC_EXTERN PetscErrorCode PetscFEGeomRestoreChunk(PetscFEGeom*,PetscInt,PetscInt,PetscFEGeom**);
PETSC_EXTERN PetscErrorCode PetscFEGeomComplete(PetscFEGeom*);
PETSC_EXTERN PetscErrorCode PetscFEGeomDestroy(PetscFEGeom**);

PETSC_EXTERN PetscErrorCode PetscDualSpaceApply(PetscDualSpace, PetscInt, PetscReal, PetscFEGeom *, PetscInt, PetscErrorCode (*)(PetscInt, PetscReal, const PetscReal [], PetscInt, PetscScalar *, void *), void *, PetscScalar *);
PETSC_EXTERN PetscErrorCode PetscDualSpaceApplyDefault(PetscDualSpace, PetscInt, PetscReal, PetscFEGeom *, PetscInt, PetscErrorCode (*)(PetscInt, PetscReal, const PetscReal [], PetscInt, PetscScalar *, void *), void *, PetscScalar *);

PETSC_EXTERN PetscErrorCode PetscDualSpaceGetAllPoints(PetscDualSpace, PetscQuadrature *);
PETSC_EXTERN PetscErrorCode PetscDualSpaceCreateAllPointsDefault(PetscDualSpace, PetscQuadrature *);

PETSC_EXTERN PetscErrorCode PetscDualSpaceApplyAll(PetscDualSpace, const PetscScalar *, PetscScalar *);
PETSC_EXTERN PetscErrorCode PetscDualSpaceApplyAllDefault(PetscDualSpace, const PetscScalar *, PetscScalar *);

PETSC_EXTERN PetscErrorCode PetscDualSpaceLagrangeGetContinuity(PetscDualSpace, PetscBool *);
PETSC_EXTERN PetscErrorCode PetscDualSpaceLagrangeSetContinuity(PetscDualSpace, PetscBool);
PETSC_EXTERN PetscErrorCode PetscDualSpaceLagrangeGetTensor(PetscDualSpace, PetscBool *);
PETSC_EXTERN PetscErrorCode PetscDualSpaceLagrangeSetTensor(PetscDualSpace, PetscBool);

PETSC_EXTERN PetscErrorCode PetscDualSpaceGetHeightSubspace(PetscDualSpace,PetscInt,PetscDualSpace *);
PETSC_EXTERN PetscErrorCode PetscDualSpaceGetPointSubspace(PetscDualSpace,PetscInt,PetscDualSpace *);
PETSC_EXTERN PetscErrorCode PetscDualSpaceSimpleSetDimension(PetscDualSpace, PetscInt);
PETSC_EXTERN PetscErrorCode PetscDualSpaceSimpleSetFunctional(PetscDualSpace, PetscInt, PetscQuadrature);

PETSC_EXTERN PetscClassId PETSCFE_CLASSID;

/*J
  PetscFEType - String with the name of a PETSc finite element space

  Level: beginner

  Note: Currently, the classes are concerned with the implementation of element integration

.seealso: PetscFESetType(), PetscFE
J*/
typedef const char *PetscFEType;
#define PETSCFEBASIC     "basic"
#define PETSCFEOPENCL    "opencl"
#define PETSCFECOMPOSITE "composite"

PETSC_EXTERN PetscFunctionList PetscFEList;
PETSC_EXTERN PetscErrorCode PetscFECreate(MPI_Comm, PetscFE *);
PETSC_EXTERN PetscErrorCode PetscFEDestroy(PetscFE *);
PETSC_EXTERN PetscErrorCode PetscFESetType(PetscFE, PetscFEType);
PETSC_EXTERN PetscErrorCode PetscFEGetType(PetscFE, PetscFEType *);
PETSC_EXTERN PetscErrorCode PetscFESetUp(PetscFE);
PETSC_EXTERN PetscErrorCode PetscFESetFromOptions(PetscFE);
PETSC_STATIC_INLINE PetscErrorCode PetscFEViewFromOptions(PetscFE A,PetscObject B,const char name[]) {return PetscObjectViewFromOptions((PetscObject)A,B,name);}

PETSC_EXTERN PetscErrorCode PetscFEView(PetscFE,PetscViewer);
PETSC_EXTERN PetscErrorCode PetscFERegister(const char [], PetscErrorCode (*)(PetscFE));
PETSC_EXTERN PetscErrorCode PetscFERegisterDestroy(void);
PETSC_EXTERN PetscErrorCode PetscFECreateDefault(MPI_Comm, PetscInt, PetscInt, PetscBool, const char [], PetscInt, PetscFE *);

PETSC_EXTERN PetscErrorCode PetscFEGetDimension(PetscFE, PetscInt *);
PETSC_EXTERN PetscErrorCode PetscFEGetSpatialDimension(PetscFE, PetscInt *);
PETSC_EXTERN PetscErrorCode PetscFESetNumComponents(PetscFE, PetscInt);
PETSC_EXTERN PetscErrorCode PetscFEGetNumComponents(PetscFE, PetscInt *);
PETSC_EXTERN PetscErrorCode PetscFEGetTileSizes(PetscFE, PetscInt *, PetscInt *, PetscInt *, PetscInt *);
PETSC_EXTERN PetscErrorCode PetscFESetTileSizes(PetscFE, PetscInt, PetscInt, PetscInt, PetscInt);
PETSC_EXTERN PetscErrorCode PetscFESetBasisSpace(PetscFE, PetscSpace);
PETSC_EXTERN PetscErrorCode PetscFEGetBasisSpace(PetscFE, PetscSpace *);
PETSC_EXTERN PetscErrorCode PetscFESetDualSpace(PetscFE, PetscDualSpace);
PETSC_EXTERN PetscErrorCode PetscFEGetDualSpace(PetscFE, PetscDualSpace *);
PETSC_EXTERN PetscErrorCode PetscFESetQuadrature(PetscFE, PetscQuadrature);
PETSC_EXTERN PetscErrorCode PetscFEGetQuadrature(PetscFE, PetscQuadrature *);
PETSC_EXTERN PetscErrorCode PetscFESetFaceQuadrature(PetscFE, PetscQuadrature);
PETSC_EXTERN PetscErrorCode PetscFEGetFaceQuadrature(PetscFE, PetscQuadrature *);
PETSC_EXTERN PetscErrorCode PetscFEGetNumDof(PetscFE, const PetscInt **);
PETSC_EXTERN PetscErrorCode PetscFEGetDefaultTabulation(PetscFE, PetscReal **, PetscReal **, PetscReal **);
PETSC_EXTERN PetscErrorCode PetscFEGetFaceTabulation(PetscFE, PetscReal **, PetscReal **, PetscReal **);
PETSC_EXTERN PetscErrorCode PetscFEGetFaceCentroidTabulation(PetscFE, PetscReal **);
PETSC_EXTERN PetscErrorCode PetscFEGetTabulation(PetscFE, PetscInt, const PetscReal[], PetscReal **, PetscReal **, PetscReal **);
PETSC_EXTERN PetscErrorCode PetscFERestoreTabulation(PetscFE, PetscInt, const PetscReal[], PetscReal **, PetscReal **, PetscReal **);
PETSC_EXTERN PetscErrorCode PetscFERefine(PetscFE, PetscFE *);
PETSC_EXTERN PetscErrorCode PetscFEGetHeightSubspace(PetscFE, PetscInt, PetscFE *);

PETSC_EXTERN PetscErrorCode PetscFEIntegrate(PetscFE, PetscDS, PetscInt, PetscInt, PetscFEGeom *, const PetscScalar[], PetscDS, const PetscScalar[], PetscScalar[]);
PETSC_EXTERN PetscErrorCode PetscFEIntegrateBd(PetscFE, PetscDS, PetscInt,
                                               void (*)(PetscInt, PetscInt, PetscInt,
                                                        const PetscInt[], const PetscInt[], const PetscScalar[], const PetscScalar[], const PetscScalar[],
                                                        const PetscInt[], const PetscInt[], const PetscScalar[], const PetscScalar[], const PetscScalar[],
                                                        PetscReal, const PetscReal[], const PetscReal[], PetscInt, const
                                                        PetscScalar[], PetscScalar[]),
                                               PetscInt, PetscFEGeom *, const PetscScalar[], PetscDS, const PetscScalar[], PetscScalar[]);
PETSC_EXTERN PetscErrorCode PetscFEIntegrateResidual(PetscFE, PetscDS, PetscInt, PetscInt, PetscFEGeom *, const PetscScalar[], const PetscScalar[], PetscDS, const PetscScalar[], PetscReal, PetscScalar[]);
PETSC_EXTERN PetscErrorCode PetscFEIntegrateBdResidual(PetscFE, PetscDS, PetscInt, PetscInt, PetscFEGeom *, const PetscScalar[], const PetscScalar[], PetscDS, const PetscScalar[], PetscReal, PetscScalar[]);
PETSC_EXTERN PetscErrorCode PetscFEIntegrateJacobian(PetscFE, PetscDS, PetscFEJacobianType, PetscInt, PetscInt, PetscInt, PetscFEGeom *, const PetscScalar[], const PetscScalar[], PetscDS, const PetscScalar[], PetscReal, PetscReal, PetscScalar[]);
PETSC_EXTERN PetscErrorCode PetscFEIntegrateBdJacobian(PetscFE, PetscDS, PetscInt, PetscInt, PetscInt, PetscFEGeom *, const PetscScalar[], const PetscScalar[], PetscDS, const PetscScalar[], PetscReal, PetscReal, PetscScalar[]);

PETSC_EXTERN PetscErrorCode PetscFECompositeGetMapping(PetscFE, PetscInt *, const PetscReal *[], const PetscReal *[], const PetscReal *[]);

PETSC_EXTERN PetscErrorCode PetscFECreateHeightTrace(PetscFE, PetscInt, PetscFE *);
PETSC_EXTERN PetscErrorCode PetscFECreatePointTrace(PetscFE, PetscInt, PetscFE *);

PETSC_EXTERN PetscErrorCode PetscFEOpenCLSetRealType(PetscFE, PetscDataType);
PETSC_EXTERN PetscErrorCode PetscFEOpenCLGetRealType(PetscFE, PetscDataType *);

#endif
