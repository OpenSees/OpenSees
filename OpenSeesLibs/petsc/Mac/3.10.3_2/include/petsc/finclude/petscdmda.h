
!
!  Include file for Fortran use of the DMDA (distributed array) package in PETSc
!
#if !defined (__PETSCDMDADEF_H)
#define __PETSCDMDADEF_H

#include "petsc/finclude/petscis.h"
#include "petsc/finclude/petscvec.h"
#include "petsc/finclude/petscmat.h"
#include "petsc/finclude/petscdm.h"

#define DMDAStencilType PetscEnum
#define DMDADirection PetscEnum

#define DMDALocalInfo PetscInt
!
!   DMDA_LOCAL_INFO_SIZE is one large than the size incase the DA is larger than an integer (on 64 bit systems).
!   non-int fields are not accessiable from fortran.
!
#define DMDA_LOCAL_INFO_SIZE 25
#define DMDA_LOCAL_INFO_DIM 1
#define DMDA_LOCAL_INFO_DOF 2
#define DMDA_LOCAL_INFO_MX 4
#define DMDA_LOCAL_INFO_MY 5
#define DMDA_LOCAL_INFO_MZ 6
#define DMDA_LOCAL_INFO_XS 7
#define DMDA_LOCAL_INFO_YS 8
#define DMDA_LOCAL_INFO_ZS 9
#define DMDA_LOCAL_INFO_XM 10
#define DMDA_LOCAL_INFO_YM 11
#define DMDA_LOCAL_INFO_ZM 12
#define DMDA_LOCAL_INFO_GXS 13
#define DMDA_LOCAL_INFO_GYS 14
#define DMDA_LOCAL_INFO_GZS 15
#define DMDA_LOCAL_INFO_GXM 16
#define DMDA_LOCAL_INFO_GYM 17
#define DMDA_LOCAL_INFO_GZM 18

#define XG_RANGE in(DMDA_LOCAL_INFO_GXS)+1:in(DMDA_LOCAL_INFO_GXS)+in(DMDA_LOCAL_INFO_GXM)
#define YG_RANGE in(DMDA_LOCAL_INFO_GYS)+1:in(DMDA_LOCAL_INFO_GYS)+in(DMDA_LOCAL_INFO_GYM)
#define ZG_RANGE in(DMDA_LOCAL_INFO_GZS)+1:in(DMDA_LOCAL_INFO_GZS)+in(DMDA_LOCAL_INFO_GZM)
#define X_RANGE in(DMDA_LOCAL_INFO_XS)+1:in(DMDA_LOCAL_INFO_XS)+in(DMDA_LOCAL_INFO_XM)
#define Y_RANGE in(DMDA_LOCAL_INFO_YS)+1:in(DMDA_LOCAL_INFO_YS)+in(DMDA_LOCAL_INFO_YM)
#define Z_RANGE in(DMDA_LOCAL_INFO_ZS)+1:in(DMDA_LOCAL_INFO_ZS)+in(DMDA_LOCAL_INFO_ZM)

#define DMDAInterpolationType PetscEnum
#define DMDAElementType PetscEnum

#define PetscGLL   PetscFortranAddr

#endif
