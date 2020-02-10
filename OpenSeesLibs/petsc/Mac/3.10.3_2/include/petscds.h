/*
      Objects which encapsulate discretizations+continuum residuals
*/
#if !defined(__PETSCDS_H)
#define __PETSCDS_H
#include <petscfe.h>
#include <petscfv.h>
#include <petscdstypes.h>

PETSC_EXTERN PetscErrorCode PetscDSInitializePackage(void);

PETSC_EXTERN PetscClassId PETSCDS_CLASSID;

/*J
  PetscDSType - String with the name of a PETSc discrete system

  Level: beginner

.seealso: PetscDSSetType(), PetscDS
J*/
typedef const char *PetscDSType;
#define PETSCDSBASIC "basic"

typedef void (*PetscPointFunc)(PetscInt, PetscInt, PetscInt,
                               const PetscInt[], const PetscInt[], const PetscScalar[], const PetscScalar[], const PetscScalar[],
                               const PetscInt[], const PetscInt[], const PetscScalar[], const PetscScalar[], const PetscScalar[],
                               PetscReal, const PetscReal[], PetscInt, const PetscScalar[], PetscScalar[]);
typedef void (*PetscPointJac)(PetscInt, PetscInt, PetscInt,
                              const PetscInt[], const PetscInt[], const PetscScalar[], const PetscScalar[], const PetscScalar[],
                              const PetscInt[], const PetscInt[], const PetscScalar[], const PetscScalar[], const PetscScalar[],
                              PetscReal, PetscReal, const PetscReal[], PetscInt, const PetscScalar[], PetscScalar[]);
typedef void (*PetscBdPointFunc)(PetscInt, PetscInt, PetscInt,
                                 const PetscInt[], const PetscInt[], const PetscScalar[], const PetscScalar[], const PetscScalar[],
                                 const PetscInt[], const PetscInt[], const PetscScalar[], const PetscScalar[], const PetscScalar[],
                                 PetscReal, const PetscReal[], const PetscReal[], PetscInt, const PetscScalar[], PetscScalar[]);
typedef void (*PetscBdPointJac)(PetscInt, PetscInt, PetscInt,
                                const PetscInt[], const PetscInt[], const PetscScalar[], const PetscScalar[], const PetscScalar[],
                                const PetscInt[], const PetscInt[], const PetscScalar[], const PetscScalar[], const PetscScalar[],
                                PetscReal, PetscReal, const PetscReal[], const PetscReal[], PetscInt, const PetscScalar[], PetscScalar[]);
typedef void (*PetscRiemannFunc)(PetscInt, PetscInt, const PetscReal[], const PetscReal[], const PetscScalar[], const PetscScalar[], PetscInt, const PetscScalar[], PetscScalar[], void *);
typedef PetscErrorCode (*PetscSimplePointFunc)(PetscInt, PetscReal, const PetscReal[], PetscInt, PetscScalar[], void *);

PETSC_EXTERN PetscFunctionList PetscDSList;
PETSC_EXTERN PetscErrorCode PetscDSCreate(MPI_Comm, PetscDS *);
PETSC_EXTERN PetscErrorCode PetscDSDestroy(PetscDS *);
PETSC_EXTERN PetscErrorCode PetscDSSetType(PetscDS, PetscDSType);
PETSC_EXTERN PetscErrorCode PetscDSGetType(PetscDS, PetscDSType *);
PETSC_EXTERN PetscErrorCode PetscDSSetUp(PetscDS);
PETSC_EXTERN PetscErrorCode PetscDSSetFromOptions(PetscDS);
PETSC_STATIC_INLINE PetscErrorCode PetscDSViewFromOptions(PetscDS A,PetscObject B,const char name[]) {return PetscObjectViewFromOptions((PetscObject)A,B,name);}

PETSC_EXTERN PetscErrorCode PetscDSView(PetscDS,PetscViewer);
PETSC_EXTERN PetscErrorCode PetscDSRegister(const char [], PetscErrorCode (*)(PetscDS));
PETSC_EXTERN PetscErrorCode PetscDSRegisterDestroy(void);

PETSC_EXTERN PetscErrorCode PetscDSGetHeightSubspace(PetscDS, PetscInt, PetscDS *);
PETSC_EXTERN PetscErrorCode PetscDSGetSpatialDimension(PetscDS, PetscInt *);
PETSC_EXTERN PetscErrorCode PetscDSGetCoordinateDimension(PetscDS, PetscInt *);
PETSC_EXTERN PetscErrorCode PetscDSSetCoordinateDimension(PetscDS, PetscInt);
PETSC_EXTERN PetscErrorCode PetscDSGetNumFields(PetscDS, PetscInt *);
PETSC_EXTERN PetscErrorCode PetscDSGetTotalDimension(PetscDS, PetscInt *);
PETSC_EXTERN PetscErrorCode PetscDSGetTotalComponents(PetscDS, PetscInt *);
PETSC_EXTERN PetscErrorCode PetscDSGetFieldIndex(PetscDS, PetscObject, PetscInt *);
PETSC_EXTERN PetscErrorCode PetscDSGetFieldSize(PetscDS, PetscInt, PetscInt *);
PETSC_EXTERN PetscErrorCode PetscDSGetFieldOffset(PetscDS, PetscInt, PetscInt *);
PETSC_EXTERN PetscErrorCode PetscDSGetDimensions(PetscDS, PetscInt *[]);
PETSC_EXTERN PetscErrorCode PetscDSGetComponents(PetscDS, PetscInt *[]);
PETSC_EXTERN PetscErrorCode PetscDSGetComponentOffset(PetscDS, PetscInt, PetscInt *);
PETSC_EXTERN PetscErrorCode PetscDSGetComponentOffsets(PetscDS, PetscInt *[]);
PETSC_EXTERN PetscErrorCode PetscDSGetComponentDerivativeOffsets(PetscDS, PetscInt *[]);

PETSC_EXTERN PetscErrorCode PetscDSGetDiscretization(PetscDS, PetscInt, PetscObject *);
PETSC_EXTERN PetscErrorCode PetscDSSetDiscretization(PetscDS, PetscInt, PetscObject);
PETSC_EXTERN PetscErrorCode PetscDSAddDiscretization(PetscDS, PetscObject);
PETSC_EXTERN PetscErrorCode PetscDSGetImplicit(PetscDS, PetscInt, PetscBool*);
PETSC_EXTERN PetscErrorCode PetscDSSetImplicit(PetscDS, PetscInt, PetscBool);
PETSC_EXTERN PetscErrorCode PetscDSGetAdjacency(PetscDS, PetscInt, PetscBool*, PetscBool*);
PETSC_EXTERN PetscErrorCode PetscDSSetAdjacency(PetscDS, PetscInt, PetscBool,  PetscBool);
PETSC_EXTERN PetscErrorCode PetscDSGetConstants(PetscDS, PetscInt *, const PetscScalar *[]);
PETSC_EXTERN PetscErrorCode PetscDSSetConstants(PetscDS, PetscInt, PetscScalar[]);
PETSC_EXTERN PetscErrorCode PetscDSGetObjective(PetscDS, PetscInt,
                                                void (**)(PetscInt, PetscInt, PetscInt,
                                                          const PetscInt[], const PetscInt[], const PetscScalar[], const PetscScalar[], const PetscScalar[],
                                                          const PetscInt[], const PetscInt[], const PetscScalar[], const PetscScalar[], const PetscScalar[],
                                                          PetscReal, const PetscReal[], PetscInt, const PetscScalar[], PetscScalar[]));
PETSC_EXTERN PetscErrorCode PetscDSSetObjective(PetscDS, PetscInt,
                                                void (*)(PetscInt, PetscInt, PetscInt,
                                                         const PetscInt[], const PetscInt[], const PetscScalar[], const PetscScalar[], const PetscScalar[],
                                                         const PetscInt[], const PetscInt[], const PetscScalar[], const PetscScalar[], const PetscScalar[],
                                                         PetscReal, const PetscReal[], PetscInt, const PetscScalar[], PetscScalar[]));
PETSC_EXTERN PetscErrorCode PetscDSGetResidual(PetscDS, PetscInt,
                                               void (**)(PetscInt, PetscInt, PetscInt,
                                                         const PetscInt[], const PetscInt[], const PetscScalar[], const PetscScalar[], const PetscScalar[],
                                                         const PetscInt[], const PetscInt[], const PetscScalar[], const PetscScalar[], const PetscScalar[],
                                                         PetscReal, const PetscReal[], PetscInt, const PetscScalar[], PetscScalar[]),
                                               void (**)(PetscInt, PetscInt, PetscInt,
                                                         const PetscInt[], const PetscInt[], const PetscScalar[], const PetscScalar[], const PetscScalar[],
                                                         const PetscInt[], const PetscInt[], const PetscScalar[], const PetscScalar[], const PetscScalar[],
                                                         PetscReal, const PetscReal[], PetscInt, const PetscScalar[], PetscScalar[]));
PETSC_EXTERN PetscErrorCode PetscDSSetResidual(PetscDS, PetscInt,
                                               void (*)(PetscInt, PetscInt, PetscInt,
                                                        const PetscInt[], const PetscInt[], const PetscScalar[], const PetscScalar[], const PetscScalar[],
                                                        const PetscInt[], const PetscInt[], const PetscScalar[], const PetscScalar[], const PetscScalar[],
                                                        PetscReal, const PetscReal[], PetscInt, const PetscScalar[], PetscScalar[]),
                                               void (*)(PetscInt, PetscInt, PetscInt,
                                                        const PetscInt[], const PetscInt[], const PetscScalar[], const PetscScalar[], const PetscScalar[],
                                                        const PetscInt[], const PetscInt[], const PetscScalar[], const PetscScalar[], const PetscScalar[],
                                                        PetscReal, const PetscReal[], PetscInt, const PetscScalar[], PetscScalar[]));
PETSC_EXTERN PetscErrorCode PetscDSHasJacobian(PetscDS, PetscBool *);
PETSC_EXTERN PetscErrorCode PetscDSGetJacobian(PetscDS, PetscInt, PetscInt,
                                               void (**)(PetscInt, PetscInt, PetscInt,
                                                         const PetscInt[], const PetscInt[], const PetscScalar[], const PetscScalar[], const PetscScalar[],
                                                         const PetscInt[], const PetscInt[], const PetscScalar[], const PetscScalar[], const PetscScalar[],
                                                         PetscReal, PetscReal, const PetscReal[], PetscInt, const PetscScalar[], PetscScalar[]),
                                               void (**)(PetscInt, PetscInt, PetscInt,
                                                         const PetscInt[], const PetscInt[], const PetscScalar[], const PetscScalar[], const PetscScalar[],
                                                         const PetscInt[], const PetscInt[], const PetscScalar[], const PetscScalar[], const PetscScalar[],
                                                         PetscReal, PetscReal, const PetscReal[], PetscInt, const PetscScalar[], PetscScalar[]),
                                               void (**)(PetscInt, PetscInt, PetscInt,
                                                         const PetscInt[], const PetscInt[], const PetscScalar[], const PetscScalar[], const PetscScalar[],
                                                         const PetscInt[], const PetscInt[], const PetscScalar[], const PetscScalar[], const PetscScalar[],
                                                         PetscReal, PetscReal, const PetscReal[], PetscInt, const PetscScalar[], PetscScalar[]),
                                               void (**)(PetscInt, PetscInt, PetscInt,
                                                         const PetscInt[], const PetscInt[], const PetscScalar[], const PetscScalar[], const PetscScalar[],
                                                         const PetscInt[], const PetscInt[], const PetscScalar[], const PetscScalar[], const PetscScalar[],
                                                         PetscReal, PetscReal, const PetscReal[], PetscInt, const PetscScalar[], PetscScalar[]));
PETSC_EXTERN PetscErrorCode PetscDSSetJacobian(PetscDS, PetscInt, PetscInt,
                                               void (*)(PetscInt, PetscInt, PetscInt,
                                                        const PetscInt[], const PetscInt[], const PetscScalar[], const PetscScalar[], const PetscScalar[],
                                                        const PetscInt[], const PetscInt[], const PetscScalar[], const PetscScalar[], const PetscScalar[],
                                                        PetscReal, PetscReal, const PetscReal[], PetscInt, const PetscScalar[], PetscScalar[]),
                                               void (*)(PetscInt, PetscInt, PetscInt,
                                                        const PetscInt[], const PetscInt[], const PetscScalar[], const PetscScalar[], const PetscScalar[],
                                                        const PetscInt[], const PetscInt[], const PetscScalar[], const PetscScalar[], const PetscScalar[],
                                                        PetscReal, PetscReal, const PetscReal[], PetscInt, const PetscScalar[], PetscScalar[]),
                                               void (*)(PetscInt, PetscInt, PetscInt,
                                                        const PetscInt[], const PetscInt[], const PetscScalar[], const PetscScalar[], const PetscScalar[],
                                                        const PetscInt[], const PetscInt[], const PetscScalar[], const PetscScalar[], const PetscScalar[],
                                                        PetscReal, PetscReal, const PetscReal[], PetscInt, const PetscScalar[], PetscScalar[]),
                                               void (*)(PetscInt, PetscInt, PetscInt,
                                                        const PetscInt[], const PetscInt[], const PetscScalar[], const PetscScalar[], const PetscScalar[],
                                                        const PetscInt[], const PetscInt[], const PetscScalar[], const PetscScalar[], const PetscScalar[],
                                                        PetscReal, PetscReal, const PetscReal[], PetscInt, const PetscScalar[], PetscScalar[]));
PETSC_EXTERN PetscErrorCode PetscDSUseJacobianPreconditioner(PetscDS, PetscBool);
PETSC_EXTERN PetscErrorCode PetscDSHasJacobianPreconditioner(PetscDS, PetscBool *);
PETSC_EXTERN PetscErrorCode PetscDSGetJacobianPreconditioner(PetscDS, PetscInt, PetscInt,
                                               void (**)(PetscInt, PetscInt, PetscInt,
                                                         const PetscInt[], const PetscInt[], const PetscScalar[], const PetscScalar[], const PetscScalar[],
                                                         const PetscInt[], const PetscInt[], const PetscScalar[], const PetscScalar[], const PetscScalar[],
                                                         PetscReal, PetscReal, const PetscReal[], PetscInt, const PetscScalar[], PetscScalar[]),
                                               void (**)(PetscInt, PetscInt, PetscInt,
                                                         const PetscInt[], const PetscInt[], const PetscScalar[], const PetscScalar[], const PetscScalar[],
                                                         const PetscInt[], const PetscInt[], const PetscScalar[], const PetscScalar[], const PetscScalar[],
                                                         PetscReal, PetscReal, const PetscReal[], PetscInt, const PetscScalar[], PetscScalar[]),
                                               void (**)(PetscInt, PetscInt, PetscInt,
                                                         const PetscInt[], const PetscInt[], const PetscScalar[], const PetscScalar[], const PetscScalar[],
                                                         const PetscInt[], const PetscInt[], const PetscScalar[], const PetscScalar[], const PetscScalar[],
                                                         PetscReal, PetscReal, const PetscReal[], PetscInt, const PetscScalar[], PetscScalar[]),
                                               void (**)(PetscInt, PetscInt, PetscInt,
                                                         const PetscInt[], const PetscInt[], const PetscScalar[], const PetscScalar[], const PetscScalar[],
                                                         const PetscInt[], const PetscInt[], const PetscScalar[], const PetscScalar[], const PetscScalar[],
                                                         PetscReal, PetscReal, const PetscReal[], PetscInt, const PetscScalar[], PetscScalar[]));
PETSC_EXTERN PetscErrorCode PetscDSSetJacobianPreconditioner(PetscDS, PetscInt, PetscInt,
                                               void (*)(PetscInt, PetscInt, PetscInt,
                                                        const PetscInt[], const PetscInt[], const PetscScalar[], const PetscScalar[], const PetscScalar[],
                                                        const PetscInt[], const PetscInt[], const PetscScalar[], const PetscScalar[], const PetscScalar[],
                                                        PetscReal, PetscReal, const PetscReal[], PetscInt, const PetscScalar[], PetscScalar[]),
                                               void (*)(PetscInt, PetscInt, PetscInt,
                                                        const PetscInt[], const PetscInt[], const PetscScalar[], const PetscScalar[], const PetscScalar[],
                                                        const PetscInt[], const PetscInt[], const PetscScalar[], const PetscScalar[], const PetscScalar[],
                                                        PetscReal, PetscReal, const PetscReal[], PetscInt, const PetscScalar[], PetscScalar[]),
                                               void (*)(PetscInt, PetscInt, PetscInt,
                                                        const PetscInt[], const PetscInt[], const PetscScalar[], const PetscScalar[], const PetscScalar[],
                                                        const PetscInt[], const PetscInt[], const PetscScalar[], const PetscScalar[], const PetscScalar[],
                                                        PetscReal, PetscReal, const PetscReal[], PetscInt, const PetscScalar[], PetscScalar[]),
                                               void (*)(PetscInt, PetscInt, PetscInt,
                                                        const PetscInt[], const PetscInt[], const PetscScalar[], const PetscScalar[], const PetscScalar[],
                                                        const PetscInt[], const PetscInt[], const PetscScalar[], const PetscScalar[], const PetscScalar[],
                                                        PetscReal, PetscReal, const PetscReal[], PetscInt, const PetscScalar[], PetscScalar[]));
PETSC_EXTERN PetscErrorCode PetscDSHasDynamicJacobian(PetscDS, PetscBool *);
PETSC_EXTERN PetscErrorCode PetscDSGetDynamicJacobian(PetscDS, PetscInt, PetscInt,
                                                      void (**)(PetscInt, PetscInt, PetscInt,
                                                                const PetscInt[], const PetscInt[], const PetscScalar[], const PetscScalar[], const PetscScalar[],
                                                                const PetscInt[], const PetscInt[], const PetscScalar[], const PetscScalar[], const PetscScalar[],
                                                                PetscReal, PetscReal, const PetscReal[], PetscInt, const PetscScalar[], PetscScalar[]),
                                                      void (**)(PetscInt, PetscInt, PetscInt,
                                                                const PetscInt[], const PetscInt[], const PetscScalar[], const PetscScalar[], const PetscScalar[],
                                                                const PetscInt[], const PetscInt[], const PetscScalar[], const PetscScalar[], const PetscScalar[],
                                                                PetscReal, PetscReal, const PetscReal[], PetscInt, const PetscScalar[], PetscScalar[]),
                                                      void (**)(PetscInt, PetscInt, PetscInt,
                                                                const PetscInt[], const PetscInt[], const PetscScalar[], const PetscScalar[], const PetscScalar[],
                                                                const PetscInt[], const PetscInt[], const PetscScalar[], const PetscScalar[], const PetscScalar[],
                                                                PetscReal, PetscReal, const PetscReal[], PetscInt, const PetscScalar[], PetscScalar[]),
                                                      void (**)(PetscInt, PetscInt, PetscInt,
                                                                const PetscInt[], const PetscInt[], const PetscScalar[], const PetscScalar[], const PetscScalar[],
                                                                const PetscInt[], const PetscInt[], const PetscScalar[], const PetscScalar[], const PetscScalar[],
                                                                PetscReal, PetscReal, const PetscReal[], PetscInt, const PetscScalar[], PetscScalar[]));
PETSC_EXTERN PetscErrorCode PetscDSSetDynamicJacobian(PetscDS, PetscInt, PetscInt,
                                                      void (*)(PetscInt, PetscInt, PetscInt,
                                                               const PetscInt[], const PetscInt[], const PetscScalar[], const PetscScalar[], const PetscScalar[],
                                                               const PetscInt[], const PetscInt[], const PetscScalar[], const PetscScalar[], const PetscScalar[],
                                                               PetscReal, PetscReal, const PetscReal[], PetscInt, const PetscScalar[], PetscScalar[]),
                                                      void (*)(PetscInt, PetscInt, PetscInt,
                                                               const PetscInt[], const PetscInt[], const PetscScalar[], const PetscScalar[], const PetscScalar[],
                                                               const PetscInt[], const PetscInt[], const PetscScalar[], const PetscScalar[], const PetscScalar[],
                                                               PetscReal, PetscReal, const PetscReal[], PetscInt, const PetscScalar[], PetscScalar[]),
                                                      void (*)(PetscInt, PetscInt, PetscInt,
                                                               const PetscInt[], const PetscInt[], const PetscScalar[], const PetscScalar[], const PetscScalar[],
                                                               const PetscInt[], const PetscInt[], const PetscScalar[], const PetscScalar[], const PetscScalar[],
                                                               PetscReal, PetscReal, const PetscReal[], PetscInt, const PetscScalar[], PetscScalar[]),
                                                      void (*)(PetscInt, PetscInt, PetscInt,
                                                               const PetscInt[], const PetscInt[], const PetscScalar[], const PetscScalar[], const PetscScalar[],
                                                               const PetscInt[], const PetscInt[], const PetscScalar[], const PetscScalar[], const PetscScalar[],
                                                               PetscReal, PetscReal, const PetscReal[], PetscInt, const PetscScalar[], PetscScalar[]));
PETSC_EXTERN PetscErrorCode PetscDSGetRiemannSolver(PetscDS, PetscInt,
                                                    void (**)(PetscInt, PetscInt, const PetscReal[], const PetscReal[], const PetscScalar[], const PetscScalar[], PetscInt, const PetscScalar[], PetscScalar[], void *));
PETSC_EXTERN PetscErrorCode PetscDSSetRiemannSolver(PetscDS, PetscInt,
                                                    void (*)(PetscInt, PetscInt, const PetscReal[], const PetscReal[], const PetscScalar[], const PetscScalar[], PetscInt, const PetscScalar[], PetscScalar[], void *));
PETSC_EXTERN PetscErrorCode PetscDSGetUpdate(PetscDS, PetscInt,
                                             void (**)(PetscInt, PetscInt, PetscInt,
                                                       const PetscInt[], const PetscInt[], const PetscScalar[], const PetscScalar[], const PetscScalar[],
                                                       const PetscInt[], const PetscInt[], const PetscScalar[], const PetscScalar[], const PetscScalar[],
                                                       PetscReal, const PetscReal[], PetscInt, const PetscScalar[], PetscScalar[]));
PETSC_EXTERN PetscErrorCode PetscDSSetUpdate(PetscDS, PetscInt,
                                             void (*)(PetscInt, PetscInt, PetscInt,
                                                      const PetscInt[], const PetscInt[], const PetscScalar[], const PetscScalar[], const PetscScalar[],
                                                      const PetscInt[], const PetscInt[], const PetscScalar[], const PetscScalar[], const PetscScalar[],
                                                      PetscReal, const PetscReal[], PetscInt, const PetscScalar[], PetscScalar[]));
PETSC_EXTERN PetscErrorCode PetscDSGetContext(PetscDS, PetscInt, void **);
PETSC_EXTERN PetscErrorCode PetscDSSetContext(PetscDS, PetscInt, void *);
PETSC_EXTERN PetscErrorCode PetscDSGetBdResidual(PetscDS, PetscInt,
                                                 void (**)(PetscInt, PetscInt, PetscInt,
                                                           const PetscInt[], const PetscInt[], const PetscScalar[], const PetscScalar[], const PetscScalar[],
                                                           const PetscInt[], const PetscInt[], const PetscScalar[], const PetscScalar[], const PetscScalar[],
                                                           PetscReal, const PetscReal[], const PetscReal[], PetscInt, const PetscScalar[], PetscScalar[]),
                                                 void (**)(PetscInt, PetscInt, PetscInt,
                                                           const PetscInt[], const PetscInt[], const PetscScalar[], const PetscScalar[], const PetscScalar[],
                                                           const PetscInt[], const PetscInt[], const PetscScalar[], const PetscScalar[], const PetscScalar[],
                                                           PetscReal, const PetscReal[], const PetscReal[], PetscInt, const PetscScalar[], PetscScalar[]));
PETSC_EXTERN PetscErrorCode PetscDSSetBdResidual(PetscDS, PetscInt,
                                                 void (*)(PetscInt, PetscInt, PetscInt,
                                                          const PetscInt[], const PetscInt[], const PetscScalar[], const PetscScalar[], const PetscScalar[],
                                                          const PetscInt[], const PetscInt[], const PetscScalar[], const PetscScalar[], const PetscScalar[],
                                                          PetscReal, const PetscReal[], const PetscReal[], PetscInt, const PetscScalar[], PetscScalar[]),
                                                 void (*)(PetscInt, PetscInt, PetscInt,
                                                          const PetscInt[], const PetscInt[], const PetscScalar[], const PetscScalar[], const PetscScalar[],
                                                          const PetscInt[], const PetscInt[], const PetscScalar[], const PetscScalar[], const PetscScalar[],
                                                          PetscReal, const PetscReal[], const PetscReal[], PetscInt, const PetscScalar[], PetscScalar[]));
PETSC_EXTERN PetscErrorCode PetscDSGetBdJacobian(PetscDS, PetscInt, PetscInt,
                                                 void (**)(PetscInt, PetscInt, PetscInt,
                                                           const PetscInt[], const PetscInt[], const PetscScalar[], const PetscScalar[], const PetscScalar[],
                                                           const PetscInt[], const PetscInt[], const PetscScalar[], const PetscScalar[], const PetscScalar[],
                                                           PetscReal, PetscReal, const PetscReal[], const PetscReal[], PetscInt, const PetscScalar[], PetscScalar[]),
                                                 void (**)(PetscInt, PetscInt, PetscInt,
                                                           const PetscInt[], const PetscInt[], const PetscScalar[], const PetscScalar[], const PetscScalar[],
                                                           const PetscInt[], const PetscInt[], const PetscScalar[], const PetscScalar[], const PetscScalar[],
                                                           PetscReal, PetscReal, const PetscReal[], const PetscReal[], PetscInt, const PetscScalar[], PetscScalar[]),
                                                 void (**)(PetscInt, PetscInt, PetscInt,
                                                           const PetscInt[], const PetscInt[], const PetscScalar[], const PetscScalar[], const PetscScalar[],
                                                           const PetscInt[], const PetscInt[], const PetscScalar[], const PetscScalar[], const PetscScalar[],
                                                           PetscReal, PetscReal, const PetscReal[], const PetscReal[], PetscInt, const PetscScalar[], PetscScalar[]),
                                                 void (**)(PetscInt, PetscInt, PetscInt,
                                                           const PetscInt[], const PetscInt[], const PetscScalar[], const PetscScalar[], const PetscScalar[],
                                                           const PetscInt[], const PetscInt[], const PetscScalar[], const PetscScalar[], const PetscScalar[],
                                                           PetscReal, PetscReal, const PetscReal[], const PetscReal[], PetscInt, const PetscScalar[], PetscScalar[]));
PETSC_EXTERN PetscErrorCode PetscDSSetBdJacobian(PetscDS, PetscInt, PetscInt,
                                                 void (*)(PetscInt, PetscInt, PetscInt,
                                                          const PetscInt[], const PetscInt[], const PetscScalar[], const PetscScalar[], const PetscScalar[],
                                                          const PetscInt[], const PetscInt[], const PetscScalar[], const PetscScalar[], const PetscScalar[],
                                                          PetscReal, PetscReal, const PetscReal[], const PetscReal[], PetscInt, const PetscScalar[], PetscScalar[]),
                                                 void (*)(PetscInt, PetscInt, PetscInt,
                                                          const PetscInt[], const PetscInt[], const PetscScalar[], const PetscScalar[], const PetscScalar[],
                                                          const PetscInt[], const PetscInt[], const PetscScalar[], const PetscScalar[], const PetscScalar[],
                                                          PetscReal, PetscReal, const PetscReal[], const PetscReal[], PetscInt, const PetscScalar[], PetscScalar[]),
                                                 void (*)(PetscInt, PetscInt, PetscInt,
                                                          const PetscInt[], const PetscInt[], const PetscScalar[], const PetscScalar[], const PetscScalar[],
                                                          const PetscInt[], const PetscInt[], const PetscScalar[], const PetscScalar[], const PetscScalar[],
                                                          PetscReal, PetscReal, const PetscReal[], const PetscReal[], PetscInt, const PetscScalar[], PetscScalar[]),
                                                 void (*)(PetscInt, PetscInt, PetscInt,
                                                          const PetscInt[], const PetscInt[], const PetscScalar[], const PetscScalar[], const PetscScalar[],
                                                          const PetscInt[], const PetscInt[], const PetscScalar[], const PetscScalar[], const PetscScalar[],
                                                          PetscReal, PetscReal, const PetscReal[], const PetscReal[], PetscInt, const PetscScalar[], PetscScalar[]));
PETSC_EXTERN PetscErrorCode PetscDSGetExactSolution(PetscDS, PetscInt, PetscErrorCode (**)(PetscInt, PetscReal, const PetscReal[], PetscInt, PetscScalar[], void *));
PETSC_EXTERN PetscErrorCode PetscDSSetExactSolution(PetscDS, PetscInt, PetscErrorCode (*)(PetscInt, PetscReal, const PetscReal[], PetscInt, PetscScalar[], void *));
PETSC_EXTERN PetscErrorCode PetscDSGetTabulation(PetscDS, PetscReal ***, PetscReal ***);
PETSC_EXTERN PetscErrorCode PetscDSGetFaceTabulation(PetscDS, PetscReal ***, PetscReal ***);
PETSC_EXTERN PetscErrorCode PetscDSGetEvaluationArrays(PetscDS, PetscScalar **, PetscScalar **, PetscScalar **);
PETSC_EXTERN PetscErrorCode PetscDSGetWeakFormArrays(PetscDS, PetscScalar **, PetscScalar **, PetscScalar **, PetscScalar **, PetscScalar **, PetscScalar **);
PETSC_EXTERN PetscErrorCode PetscDSGetRefCoordArrays(PetscDS, PetscReal **, PetscScalar **);
PETSC_EXTERN PetscErrorCode PetscDSCopyConstants(PetscDS, PetscDS);
PETSC_EXTERN PetscErrorCode PetscDSCopyEquations(PetscDS, PetscDS);
PETSC_EXTERN PetscErrorCode PetscDSSelectEquations(PetscDS, PetscInt, const PetscInt[], PetscDS);
PETSC_EXTERN PetscErrorCode PetscDSAddBoundary(PetscDS, DMBoundaryConditionType, const char[], const char[], PetscInt, PetscInt, const PetscInt *, void (*)(void), PetscInt, const PetscInt *, void *);
PETSC_EXTERN PetscErrorCode PetscDSGetNumBoundary(PetscDS, PetscInt *);
PETSC_EXTERN PetscErrorCode PetscDSGetBoundary(PetscDS, PetscInt, DMBoundaryConditionType *, const char **, const char **, PetscInt *, PetscInt *, const PetscInt **, void (**)(void), PetscInt *, const PetscInt **, void **);
PETSC_EXTERN PetscErrorCode PetscDSCopyBoundary(PetscDS, PetscDS);

#endif
