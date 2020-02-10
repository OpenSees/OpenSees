/*
   Distributed arrays - communication tools for parallel, rectangular grids.
*/

#if !defined(_DAIMPL_H)
#define _DAIMPL_H

#include <petscdmda.h>
#include <petsc/private/dmimpl.h>

typedef struct {
  PetscInt              M,N,P;                 /* array dimensions */
  PetscInt              m,n,p;                 /* processor layout */
  PetscInt              w;                     /* degrees of freedom per node */
  PetscInt              s;                     /* stencil width */
  PetscInt              xs,xe,ys,ye,zs,ze;     /* range of local values */
  PetscInt              Xs,Xe,Ys,Ye,Zs,Ze;     /* range including ghost values
                                                   values above already scaled by w */
  PetscInt              base;                  /* global number of 1st local node, includes the * w term */
  DMBoundaryType        bx,by,bz;              /* indicates type of ghost nodes at boundary */
  VecScatter            gtol,ltol;        /* scatters, see below for details */
  DMDAStencilType       stencil_type;          /* stencil, either box or star */
  DMDAInterpolationType interptype;

  PetscInt              nlocal,Nlocal;         /* local size of local vector and global vector, includes the * w term */

  PetscInt              xol,yol,zol;           /* overlap of local subdomains */
  PetscInt              xo,yo,zo;              /* offsets for the indices in x y and z */
  PetscInt              Mo,No,Po;              /* the size of the problem the offset is in to */
  PetscInt              Nsub;                  /* number of local subdomains to decompose into */
  PetscInt              nonxs,nonys,nonzs;     /* the nonoverlapping starts in the case of a subdomain da */
  PetscInt              nonxm,nonym,nonzm;     /* the nonoverlapping sizes in the case of a subdomain da */

  AO                    ao;                    /* application ordering context */
  AOType                aotype;                /* type of application ordering */

  char                  **fieldname;           /* names of individual components in vectors */
  char                  **coordinatename;      /* names of coordinate directions, for example, x, y, z */

  PetscInt              *lx,*ly,*lz;        /* number of nodes in each partition block along 3 axis */
  Vec                   natural;            /* global vector for storing items in natural order */
  VecScatter            gton;               /* vector scatter from global to natural */
  PetscMPIInt           *neighbors;         /* ranks of all neighbors and self */

  ISColoring            localcoloring;       /* set by DMCreateColoring() */
  ISColoring            ghostedcoloring;

  DMDAElementType       elementtype;
  PetscInt              ne;                  /* number of elements */
  PetscInt              nen;                 /* number of nodes per element */
  PetscInt              *e;                  /* the elements */
  IS                    ecorners;            /* corners of the subdomain */

  PetscInt              refine_x,refine_y,refine_z;    /* ratio used in refining */
  PetscInt              coarsen_x,coarsen_y,coarsen_z; /* ratio used for coarsening */
                        /* if the refinement is done differently on different levels */
  PetscInt              refine_x_hier_n,*refine_x_hier,refine_y_hier_n,*refine_y_hier,refine_z_hier_n,*refine_z_hier;

#define DMDA_MAX_WORK_ARRAYS 2 /* work arrays for holding work via DMDAGetArray() */
  void                  *arrayin[DMDA_MAX_WORK_ARRAYS],*arrayout[DMDA_MAX_WORK_ARRAYS];
  void                  *arrayghostedin[DMDA_MAX_WORK_ARRAYS],*arrayghostedout[DMDA_MAX_WORK_ARRAYS];
  void                  *startin[DMDA_MAX_WORK_ARRAYS],*startout[DMDA_MAX_WORK_ARRAYS];
  void                  *startghostedin[DMDA_MAX_WORK_ARRAYS],*startghostedout[DMDA_MAX_WORK_ARRAYS];

  PetscErrorCode (*lf)(DM, Vec, Vec, void *);
  PetscErrorCode (*lj)(DM, Vec, Vec, void *);

  /* used by DMDASetBlockFills() */
  PetscInt              *ofill,*dfill;
  PetscInt              *ofillcols;

  /* used by DMDASetMatPreallocateOnly() */
  PetscBool             prealloc_only;
  PetscInt              preallocCenterDim; /* Dimension of the points which connect adjacent points for preallocation */
} DM_DA;

/*
  Vectors:
     Global has on each processor the interior degrees of freedom and
         no ghost points. This vector is what the solvers usually see.
     Local has on each processor the ghost points as well. This is
          what code to calculate Jacobians, etc. usually sees.
  Vector scatters:
     gtol - Global representation to local
     ltol - Local representation to local representation, updates the
            ghostpoint values in the second vector from (correct) interior
            values in the first vector.  This is good for explicit
            nearest neighbor timestepping.
*/

PETSC_INTERN PetscErrorCode VecView_MPI_DA(Vec,PetscViewer);
PETSC_INTERN PetscErrorCode VecLoad_Default_DA(Vec, PetscViewer);
PETSC_INTERN PetscErrorCode DMView_DA_Matlab(DM,PetscViewer);
PETSC_INTERN PetscErrorCode DMView_DA_Binary(DM,PetscViewer);
PETSC_INTERN PetscErrorCode DMView_DA_VTK(DM,PetscViewer);
PETSC_INTERN PetscErrorCode DMView_DA_GLVis(DM,PetscViewer);
PETSC_EXTERN PetscErrorCode DMDAVTKWriteAll(PetscObject,PetscViewer);
PETSC_EXTERN PetscErrorCode DMDASelectFields(DM,PetscInt*,PetscInt**);

#endif
