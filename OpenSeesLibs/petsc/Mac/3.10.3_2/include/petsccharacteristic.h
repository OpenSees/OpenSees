/*
   Defines the interface functions for the method of characteristics solvers
*/
#ifndef __PETSCCHARACTERISTICS_H
#define __PETSCCHARACTERISTICS_H

#include <petscvec.h>
#include <petscdmdatypes.h>

PETSC_EXTERN PetscErrorCode CharacteristicInitializePackage(void);

/*S
     Characteristic - Abstract PETSc object that manages method of characteristics solves

   Level: beginner

  Concepts: Method of characteristics

.seealso:  CharacteristicCreate(), CharacteristicSetType(), CharacteristicType, SNES, TS, PC, KSP
S*/
typedef struct _p_Characteristic *Characteristic;

/*J
    CharacteristicType - String with the name of a characteristics method.

   Level: beginner

.seealso: CharacteristicSetType(), Characteristic
J*/
#define CHARACTERISTICDA "da"
typedef const char* CharacteristicType;

PETSC_EXTERN PetscErrorCode CharacteristicCreate(MPI_Comm, Characteristic *);
PETSC_EXTERN PetscErrorCode CharacteristicSetType(Characteristic, CharacteristicType);
PETSC_EXTERN PetscErrorCode CharacteristicSetUp(Characteristic);
PETSC_EXTERN PetscErrorCode CharacteristicSetVelocityInterpolation(Characteristic, DM, Vec, Vec, PetscInt, PetscInt[], PetscErrorCode (*)(Vec, PetscReal[], PetscInt, PetscInt[], PetscScalar[], void *), void *);
PETSC_EXTERN PetscErrorCode CharacteristicSetVelocityInterpolationLocal(Characteristic, DM, Vec, Vec, PetscInt, PetscInt[], PetscErrorCode (*)(void *, PetscReal[], PetscInt, PetscInt[], PetscScalar[], void *), void *);
PETSC_EXTERN PetscErrorCode CharacteristicSetFieldInterpolation(Characteristic, DM, Vec, PetscInt, PetscInt[], PetscErrorCode (*)(Vec, PetscReal[], PetscInt, PetscInt[], PetscScalar[], void *), void *);
PETSC_EXTERN PetscErrorCode CharacteristicSetFieldInterpolationLocal(Characteristic, DM, Vec, PetscInt, PetscInt[], PetscErrorCode (*)(void *, PetscReal[], PetscInt, PetscInt[], PetscScalar[], void *), void *);
PETSC_EXTERN PetscErrorCode CharacteristicSolve(Characteristic, PetscReal, Vec);
PETSC_EXTERN PetscErrorCode CharacteristicDestroy(Characteristic*);

PETSC_EXTERN PetscFunctionList CharacteristicList;

PETSC_EXTERN PetscErrorCode CharacteristicRegister(const char[],PetscErrorCode (*)(Characteristic));

#endif /*__PETSCCHARACTERISTICS_H*/
