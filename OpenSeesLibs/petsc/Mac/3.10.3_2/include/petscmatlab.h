/*
    Defines an interface to the MATLAB Engine from PETSc
*/

#if !defined(__PETSCMATLAB_H)
#define __PETSCMATLAB_H

PETSC_EXTERN PetscClassId MATLABENGINE_CLASSID;

/*S
     PetscMatlabEngine - Object used to communicate with MATLAB

   Level: intermediate

.seealso:  PetscMatlabEngineCreate(), PetscMatlabEngineDestroy(), PetscMatlabEngineEvaluate(),
           PetscMatlabEngineGetOutput(), PetscMatlabEnginePut(), PetscMatlabEngineGet(),
           PetscMatlabEnginePrintOutput(), PetscMatlabEnginePutArray(), PetscMatlabEngineGetArray(),
           PETSC_MATLAB_ENGINE_(), PETSC_MATLAB_ENGINE_SELF, PETSC_MATLAB_ENGINE_WORLD
S*/
typedef struct _p_PetscMatlabEngine* PetscMatlabEngine;

PETSC_EXTERN PetscErrorCode PetscMatlabEngineCreate(MPI_Comm,const char[],PetscMatlabEngine*);
PETSC_EXTERN PetscErrorCode PetscMatlabEngineDestroy(PetscMatlabEngine*);
PETSC_EXTERN PetscErrorCode PetscMatlabEngineEvaluate(PetscMatlabEngine,const char[],...);
PETSC_EXTERN PetscErrorCode PetscMatlabEngineGetOutput(PetscMatlabEngine,char **);
PETSC_EXTERN PetscErrorCode PetscMatlabEnginePrintOutput(PetscMatlabEngine,FILE*);
PETSC_EXTERN PetscErrorCode PetscMatlabEnginePut(PetscMatlabEngine,PetscObject);
PETSC_EXTERN PetscErrorCode PetscMatlabEngineGet(PetscMatlabEngine,PetscObject);
PETSC_EXTERN PetscErrorCode PetscMatlabEnginePutArray(PetscMatlabEngine,int,int,const PetscScalar*,const char[]);
PETSC_EXTERN PetscErrorCode PetscMatlabEngineGetArray(PetscMatlabEngine,int,int,PetscScalar*,const char[]);

PETSC_EXTERN PetscMatlabEngine  PETSC_MATLAB_ENGINE_(MPI_Comm);

/*MC
  PETSC_MATLAB_ENGINE_WORLD - same as PETSC_MATLAB_ENGINE_(PETSC_COMM_WORLD)

  Level: developer
M*/
#define PETSC_MATLAB_ENGINE_WORLD PETSC_MATLAB_ENGINE_(PETSC_COMM_WORLD)

/*MC
  PETSC_MATLAB_ENGINE_SELF - same as PETSC_MATLAB_ENGINE_(PETSC_COMM_SELF)

  Level: developer
M*/
#define PETSC_MATLAB_ENGINE_SELF  PETSC_MATLAB_ENGINE_(PETSC_COMM_SELF)

#endif
