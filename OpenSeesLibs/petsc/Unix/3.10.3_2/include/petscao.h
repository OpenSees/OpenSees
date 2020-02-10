/*
  An application ordering is mapping between an application-centric
  ordering (the ordering that is "natural" for the application) and
  the parallel ordering that PETSc uses.
*/
#if !defined(__PETSCAO_H)
#define __PETSCAO_H
#include <petscis.h>

/*S
     AO - Abstract PETSc object that manages mapping between different global numbering

   Level: intermediate

  Concepts: global numbering

.seealso:  AOCreateBasic(), AOCreateBasicIS(), AOPetscToApplication(), AOView(), AOApplicationToPetsc()
S*/
typedef struct _p_AO* AO;

/*J
    AOType - String with the name of a PETSc application ordering or the creation function
       with an optional dynamic library name.

   Level: beginner

.seealso: AOSetType(), AO
J*/
typedef const char* AOType;
#define AOBASIC               "basic"
#define AOADVANCED            "advanced"
#define AOMAPPING             "mapping"
#define AOMEMORYSCALABLE      "memoryscalable"

/* Logging support */
PETSC_EXTERN PetscClassId AO_CLASSID;

PETSC_EXTERN PetscErrorCode AOInitializePackage(void);

PETSC_EXTERN PetscErrorCode AOCreate(MPI_Comm,AO*);
PETSC_EXTERN PetscErrorCode AOSetIS(AO,IS,IS);
PETSC_EXTERN PetscErrorCode AOSetFromOptions(AO);

PETSC_EXTERN PetscErrorCode AOCreateBasic(MPI_Comm,PetscInt,const PetscInt[],const PetscInt[],AO*);
PETSC_EXTERN PetscErrorCode AOCreateBasicIS(IS,IS,AO*);
PETSC_EXTERN PetscErrorCode AOCreateMemoryScalable(MPI_Comm,PetscInt,const PetscInt[],const PetscInt[],AO*);
PETSC_EXTERN PetscErrorCode AOCreateMemoryScalableIS(IS,IS,AO*);
PETSC_EXTERN PetscErrorCode AOCreateMapping(MPI_Comm,PetscInt,const PetscInt[],const PetscInt[],AO*);
PETSC_EXTERN PetscErrorCode AOCreateMappingIS(IS,IS,AO*);

PETSC_EXTERN PetscErrorCode AOView(AO,PetscViewer);
PETSC_STATIC_INLINE PetscErrorCode AOViewFromOptions(AO A,PetscObject obj,const char name[]) {return PetscObjectViewFromOptions((PetscObject)A,obj,name);}
PETSC_EXTERN PetscErrorCode AODestroy(AO*);

/* Dynamic creation and loading functions */
PETSC_EXTERN PetscFunctionList AOList;
PETSC_EXTERN PetscErrorCode AOSetType(AO, AOType);
PETSC_EXTERN PetscErrorCode AOGetType(AO, AOType *);

PETSC_EXTERN PetscErrorCode AORegister(const char [], PetscErrorCode (*)(AO));

PETSC_EXTERN PetscErrorCode AOPetscToApplication(AO,PetscInt,PetscInt[]);
PETSC_EXTERN PetscErrorCode AOApplicationToPetsc(AO,PetscInt,PetscInt[]);
PETSC_EXTERN PetscErrorCode AOPetscToApplicationIS(AO,IS);
PETSC_EXTERN PetscErrorCode AOApplicationToPetscIS(AO,IS);

PETSC_EXTERN PetscErrorCode AOPetscToApplicationPermuteInt(AO, PetscInt, PetscInt[]);
PETSC_EXTERN PetscErrorCode AOApplicationToPetscPermuteInt(AO, PetscInt, PetscInt[]);
PETSC_EXTERN PetscErrorCode AOPetscToApplicationPermuteReal(AO, PetscInt, PetscReal[]);
PETSC_EXTERN PetscErrorCode AOApplicationToPetscPermuteReal(AO, PetscInt, PetscReal[]);

PETSC_EXTERN PetscErrorCode AOMappingHasApplicationIndex(AO, PetscInt, PetscBool  *);
PETSC_EXTERN PetscErrorCode AOMappingHasPetscIndex(AO, PetscInt, PetscBool  *);

/* ----------------------------------------------------*/
#endif
