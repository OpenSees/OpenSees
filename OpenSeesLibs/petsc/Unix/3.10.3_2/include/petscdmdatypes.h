#if !defined(_PETSCDMDATYPES_H)
#define _PETSCDMDATYPES_H

#include <petscdmtypes.h>

/*E
    DMDAStencilType - Determines if the stencil extends only along the coordinate directions, or also
      to the northeast, northwest etc

   Level: beginner

.seealso: DMDACreate1d(), DMDACreate2d(), DMDACreate3d(), DMDACreate(), DMDASetStencilType()
E*/
typedef enum { DMDA_STENCIL_STAR,DMDA_STENCIL_BOX } DMDAStencilType;

/*E
    DMDAInterpolationType - Defines the type of interpolation that will be returned by
       DMCreateInterpolation.

   Level: beginner

.seealso: DMDACreate1d(), DMDACreate2d(), DMDACreate3d(), DMCreateInterpolation(), DMDASetInterpolationType(), DMDACreate()
E*/
typedef enum { DMDA_Q0, DMDA_Q1 } DMDAInterpolationType;

/*E
    DMDAElementType - Defines the type of elements that will be returned by
       DMDAGetElements()

   Level: beginner

.seealso: DMDACreate1d(), DMDACreate2d(), DMDACreate3d(), DMCreateInterpolation(), DMDASetInterpolationType(),
          DMDASetElementType(), DMDAGetElements(), DMDARestoreElements(), DMDACreate()
E*/
typedef enum { DMDA_ELEMENT_P1, DMDA_ELEMENT_Q1 } DMDAElementType;

/*S
     DMDALocalInfo - C struct that contains information about a structured grid and a processors logical
              location in it.

   Level: beginner

  Concepts: distributed array

  Fortran Notes - This should be declared as
$    DMDALocalInfo :: info(DMDA_LOCAL_INFO_SIZE)
     and the entries accessed via
$    info(DMDA_LOCAL_INFO_DIM)
$    info(DMDA_LOCAL_INFO_DOF) etc.
   The entries bx,by,bz, st, and da are not accessible from Fortran.

.seealso:  DMDACreate1d(), DMDACreate2d(), DMDACreate3d(), DMDestroy(), DM, DMDAGetLocalInfo(), DMDAGetInfo()
S*/
typedef struct {
  PetscInt         dim,dof,sw;
  PetscInt         mx,my,mz;    /* global number of grid points in each direction */
  PetscInt         xs,ys,zs;    /* starting point of this processor, excluding ghosts */
  PetscInt         xm,ym,zm;    /* number of grid points on this processor, excluding ghosts */
  PetscInt         gxs,gys,gzs;    /* starting point of this processor including ghosts */
  PetscInt         gxm,gym,gzm;    /* number of grid points on this processor including ghosts */
  DMBoundaryType   bx,by,bz; /* type of ghost nodes at boundary */
  DMDAStencilType  st;
  DM               da;
} DMDALocalInfo;

#endif
