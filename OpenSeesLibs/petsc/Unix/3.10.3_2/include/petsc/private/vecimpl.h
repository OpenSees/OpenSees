
/*
   This private file should not be included in users' code.
   Defines the fields shared by all vector implementations.

*/

#ifndef __VECIMPL_H
#define __VECIMPL_H

#include <petscvec.h>
#include <petsc/private/petscimpl.h>
#include <petscviewer.h>

PETSC_EXTERN PetscBool VecRegisterAllCalled;
PETSC_EXTERN PetscErrorCode VecRegisterAll(void);

/* ----------------------------------------------------------------------------*/

typedef struct _VecOps *VecOps;
struct _VecOps {
  PetscErrorCode (*duplicate)(Vec,Vec*);         /* get single vector */
  PetscErrorCode (*duplicatevecs)(Vec,PetscInt,Vec**);     /* get array of vectors */
  PetscErrorCode (*destroyvecs)(PetscInt,Vec[]);           /* free array of vectors */
  PetscErrorCode (*dot)(Vec,Vec,PetscScalar*);             /* z = x^H * y */
  PetscErrorCode (*mdot)(Vec,PetscInt,const Vec[],PetscScalar*); /* z[j] = x dot y[j] */
  PetscErrorCode (*norm)(Vec,NormType,PetscReal*);        /* z = sqrt(x^H * x) */
  PetscErrorCode (*tdot)(Vec,Vec,PetscScalar*);             /* x'*y */
  PetscErrorCode (*mtdot)(Vec,PetscInt,const Vec[],PetscScalar*);/* z[j] = x dot y[j] */
  PetscErrorCode (*scale)(Vec,PetscScalar);                 /* x = alpha * x   */
  PetscErrorCode (*copy)(Vec,Vec);                     /* y = x */
  PetscErrorCode (*set)(Vec,PetscScalar);                        /* y = alpha  */
  PetscErrorCode (*swap)(Vec,Vec);                               /* exchange x and y */
  PetscErrorCode (*axpy)(Vec,PetscScalar,Vec);                   /* y = y + alpha * x */
  PetscErrorCode (*axpby)(Vec,PetscScalar,PetscScalar,Vec);      /* y = alpha * x + beta * y*/
  PetscErrorCode (*maxpy)(Vec,PetscInt,const PetscScalar*,Vec*); /* y = y + alpha[j] x[j] */
  PetscErrorCode (*aypx)(Vec,PetscScalar,Vec);                   /* y = x + alpha * y */
  PetscErrorCode (*waxpy)(Vec,PetscScalar,Vec,Vec);         /* w = y + alpha * x */
  PetscErrorCode (*axpbypcz)(Vec,PetscScalar,PetscScalar,PetscScalar,Vec,Vec);   /* z = alpha * x + beta *y + gamma *z*/
  PetscErrorCode (*pointwisemult)(Vec,Vec,Vec);        /* w = x .* y */
  PetscErrorCode (*pointwisedivide)(Vec,Vec,Vec);      /* w = x ./ y */
  PetscErrorCode (*setvalues)(Vec,PetscInt,const PetscInt[],const PetscScalar[],InsertMode);
  PetscErrorCode (*assemblybegin)(Vec);                /* start global assembly */
  PetscErrorCode (*assemblyend)(Vec);                  /* end global assembly */
  PetscErrorCode (*getarray)(Vec,PetscScalar**);            /* get data array */
  PetscErrorCode (*getsize)(Vec,PetscInt*);
  PetscErrorCode (*getlocalsize)(Vec,PetscInt*);
  PetscErrorCode (*restorearray)(Vec,PetscScalar**);        /* restore data array */
  PetscErrorCode (*max)(Vec,PetscInt*,PetscReal*);      /* z = max(x); idx=index of max(x) */
  PetscErrorCode (*min)(Vec,PetscInt*,PetscReal*);      /* z = min(x); idx=index of min(x) */
  PetscErrorCode (*setrandom)(Vec,PetscRandom);         /* set y[j] = random numbers */
  PetscErrorCode (*setoption)(Vec,VecOption,PetscBool );
  PetscErrorCode (*setvaluesblocked)(Vec,PetscInt,const PetscInt[],const PetscScalar[],InsertMode);
  PetscErrorCode (*destroy)(Vec);
  PetscErrorCode (*view)(Vec,PetscViewer);
  PetscErrorCode (*placearray)(Vec,const PetscScalar*);     /* place data array */
  PetscErrorCode (*replacearray)(Vec,const PetscScalar*);     /* replace data array */
  PetscErrorCode (*dot_local)(Vec,Vec,PetscScalar*);
  PetscErrorCode (*tdot_local)(Vec,Vec,PetscScalar*);
  PetscErrorCode (*norm_local)(Vec,NormType,PetscReal*);
  PetscErrorCode (*mdot_local)(Vec,PetscInt,const Vec[],PetscScalar*);
  PetscErrorCode (*mtdot_local)(Vec,PetscInt,const Vec[],PetscScalar*);
  PetscErrorCode (*load)(Vec,PetscViewer);
  PetscErrorCode (*reciprocal)(Vec);
  PetscErrorCode (*conjugate)(Vec);
  PetscErrorCode (*setlocaltoglobalmapping)(Vec,ISLocalToGlobalMapping);
  PetscErrorCode (*setvalueslocal)(Vec,PetscInt,const PetscInt *,const PetscScalar *,InsertMode);
  PetscErrorCode (*resetarray)(Vec);      /* vector points to its original array, i.e. undoes any VecPlaceArray() */
  PetscErrorCode (*setfromoptions)(PetscOptionItems*,Vec);
  PetscErrorCode (*maxpointwisedivide)(Vec,Vec,PetscReal*);      /* m = max abs(x ./ y) */
  PetscErrorCode (*pointwisemax)(Vec,Vec,Vec);
  PetscErrorCode (*pointwisemaxabs)(Vec,Vec,Vec);
  PetscErrorCode (*pointwisemin)(Vec,Vec,Vec);
  PetscErrorCode (*getvalues)(Vec,PetscInt,const PetscInt[],PetscScalar[]);
  PetscErrorCode (*sqrt)(Vec);
  PetscErrorCode (*abs)(Vec);
  PetscErrorCode (*exp)(Vec);
  PetscErrorCode (*log)(Vec);
  PetscErrorCode (*shift)(Vec,PetscScalar);
  PetscErrorCode (*create)(Vec);
  PetscErrorCode (*stridegather)(Vec,PetscInt,Vec,InsertMode);
  PetscErrorCode (*stridescatter)(Vec,PetscInt,Vec,InsertMode);
  PetscErrorCode (*dotnorm2)(Vec,Vec,PetscScalar*,PetscScalar*);
  PetscErrorCode (*getsubvector)(Vec,IS,Vec*);
  PetscErrorCode (*restoresubvector)(Vec,IS,Vec*);
  PetscErrorCode (*getarrayread)(Vec,const PetscScalar**);
  PetscErrorCode (*restorearrayread)(Vec,const PetscScalar**);
  PetscErrorCode (*stridesubsetgather)(Vec,PetscInt,const PetscInt[],const PetscInt[],Vec,InsertMode);
  PetscErrorCode (*stridesubsetscatter)(Vec,PetscInt,const PetscInt[],const PetscInt[],Vec,InsertMode);
  PetscErrorCode (*viewnative)(Vec,PetscViewer);
  PetscErrorCode (*loadnative)(Vec,PetscViewer);
  PetscErrorCode (*getlocalvector)(Vec,Vec);
  PetscErrorCode (*restorelocalvector)(Vec,Vec);
  PetscErrorCode (*getlocalvectorread)(Vec,Vec);
  PetscErrorCode (*restorelocalvectorread)(Vec,Vec);
};

/*
    The stash is used to temporarily store inserted vec values that
  belong to another processor. During the assembly phase the stashed
  values are moved to the correct processor and
*/

typedef struct {
  PetscInt      nmax;                   /* maximum stash size */
  PetscInt      umax;                   /* max stash size user wants */
  PetscInt      oldnmax;                /* the nmax value used previously */
  PetscInt      n;                      /* stash size */
  PetscInt      bs;                     /* block size of the stash */
  PetscInt      reallocs;               /* preserve the no of mallocs invoked */
  PetscInt      *idx;                   /* global row numbers in stash */
  PetscScalar   *array;                 /* array to hold stashed values */
  /* The following variables are used for communication */
  MPI_Comm      comm;
  PetscMPIInt   size,rank;
  PetscMPIInt   tag1,tag2;
  MPI_Request   *send_waits;            /* array of send requests */
  MPI_Request   *recv_waits;            /* array of receive requests */
  MPI_Status    *send_status;           /* array of send status */
  PetscInt      nsends,nrecvs;          /* numbers of sends and receives */
  PetscScalar   *svalues,*rvalues;      /* sending and receiving data */
  PetscInt      *sindices,*rindices;
  PetscInt      rmax;                   /* maximum message length */
  PetscInt      *nprocs;                /* tmp data used both during scatterbegin and end */
  PetscInt      nprocessed;             /* number of messages already processed */
  PetscBool     donotstash;
  PetscBool     ignorenegidx;           /* ignore negative indices passed into VecSetValues/VetGetValues */
  InsertMode    insertmode;
  PetscInt      *bowners;
} VecStash;

struct _p_Vec {
  PETSCHEADER(struct _VecOps);
  PetscLayout            map;
  void                   *data;     /* implementation-specific data */
  PetscBool              array_gotten;
  VecStash               stash,bstash; /* used for storing off-proc values during assembly */
  PetscBool              petscnative;  /* means the ->data starts with VECHEADER and can use VecGetArrayFast()*/
  PetscInt               lock;   /* vector is locked to read only */
#if defined(PETSC_HAVE_VIENNACL) || defined(PETSC_HAVE_VECCUDA)
  PetscOffloadFlag       valid_GPU_array;    /* indicates where the most recently modified vector data is (GPU or CPU) */
  void                   *spptr; /* this is the special pointer to the array on the GPU */
#endif
};

PETSC_EXTERN PetscLogEvent VEC_SetRandom;
PETSC_EXTERN PetscLogEvent VEC_View;
PETSC_EXTERN PetscLogEvent VEC_Max;
PETSC_EXTERN PetscLogEvent VEC_Min;
PETSC_EXTERN PetscLogEvent VEC_Dot;
PETSC_EXTERN PetscLogEvent VEC_MDot;
PETSC_EXTERN PetscLogEvent VEC_TDot;
PETSC_EXTERN PetscLogEvent VEC_MTDot;
PETSC_EXTERN PetscLogEvent VEC_Norm;
PETSC_EXTERN PetscLogEvent VEC_Normalize;
PETSC_EXTERN PetscLogEvent VEC_Scale;
PETSC_EXTERN PetscLogEvent VEC_Copy;
PETSC_EXTERN PetscLogEvent VEC_Set;
PETSC_EXTERN PetscLogEvent VEC_AXPY;
PETSC_EXTERN PetscLogEvent VEC_AYPX;
PETSC_EXTERN PetscLogEvent VEC_WAXPY;
PETSC_EXTERN PetscLogEvent VEC_MAXPY;
PETSC_EXTERN PetscLogEvent VEC_AssemblyEnd;
PETSC_EXTERN PetscLogEvent VEC_PointwiseMult;
PETSC_EXTERN PetscLogEvent VEC_SetValues;
PETSC_EXTERN PetscLogEvent VEC_Load;
PETSC_EXTERN PetscLogEvent VEC_ScatterBegin;
PETSC_EXTERN PetscLogEvent VEC_ScatterEnd;
PETSC_EXTERN PetscLogEvent VEC_ReduceArithmetic;
PETSC_EXTERN PetscLogEvent VEC_ReduceCommunication;
PETSC_EXTERN PetscLogEvent VEC_ReduceBegin;
PETSC_EXTERN PetscLogEvent VEC_ReduceEnd;
PETSC_EXTERN PetscLogEvent VEC_Swap;
PETSC_EXTERN PetscLogEvent VEC_AssemblyBegin;
PETSC_EXTERN PetscLogEvent VEC_DotNorm2;
PETSC_EXTERN PetscLogEvent VEC_AXPBYPCZ;
PETSC_EXTERN PetscLogEvent VEC_Ops;
PETSC_EXTERN PetscLogEvent VEC_ViennaCLCopyToGPU;
PETSC_EXTERN PetscLogEvent VEC_ViennaCLCopyFromGPU;
PETSC_EXTERN PetscLogEvent VEC_CUDACopyToGPU;
PETSC_EXTERN PetscLogEvent VEC_CUDACopyFromGPU;
PETSC_EXTERN PetscLogEvent VEC_CUDACopyToGPUSome;
PETSC_EXTERN PetscLogEvent VEC_CUDACopyFromGPUSome;

PETSC_EXTERN PetscErrorCode VecView_Seq(Vec,PetscViewer);
#if defined(PETSC_HAVE_VIENNACL)
PETSC_EXTERN PetscErrorCode VecViennaCLAllocateCheckHost(Vec v);
PETSC_EXTERN PetscErrorCode VecViennaCLCopyFromGPU(Vec v);
#endif
#if defined(PETSC_HAVE_VECCUDA)
PETSC_EXTERN PetscErrorCode VecCUDAAllocateCheckHost(Vec v);
PETSC_EXTERN PetscErrorCode VecCUDACopyFromGPU(Vec v);
#endif


/*
     Common header shared by array based vectors,
   currently Vec_Seq and Vec_MPI
*/
#define VECHEADER                          \
  PetscScalar *array;                      \
  PetscScalar *array_allocated;                        /* if the array was allocated by PETSc this is its pointer */  \
  PetscScalar *unplacedarray;                           /* if one called VecPlaceArray(), this is where it stashed the original */

/* Default obtain and release vectors; can be used by any implementation */
PETSC_EXTERN PetscErrorCode VecDuplicateVecs_Default(Vec,PetscInt,Vec *[]);
PETSC_EXTERN PetscErrorCode VecDestroyVecs_Default(PetscInt,Vec []);
PETSC_INTERN PetscErrorCode VecLoad_Binary(Vec, PetscViewer);
PETSC_EXTERN PetscErrorCode VecLoad_Default(Vec, PetscViewer);

PETSC_EXTERN PetscInt  NormIds[7];  /* map from NormType to IDs used to cache/retreive values of norms */

PETSC_INTERN PetscErrorCode VecStashCreate_Private(MPI_Comm,PetscInt,VecStash*);
PETSC_INTERN PetscErrorCode VecStashDestroy_Private(VecStash*);
PETSC_EXTERN PetscErrorCode VecStashExpand_Private(VecStash*,PetscInt);
PETSC_INTERN PetscErrorCode VecStashScatterEnd_Private(VecStash*);
PETSC_INTERN PetscErrorCode VecStashSetInitialSize_Private(VecStash*,PetscInt);
PETSC_INTERN PetscErrorCode VecStashGetInfo_Private(VecStash*,PetscInt*,PetscInt*);
PETSC_INTERN PetscErrorCode VecStashScatterBegin_Private(VecStash*,PetscInt*);
PETSC_INTERN PetscErrorCode VecStashScatterGetMesg_Private(VecStash*,PetscMPIInt*,PetscInt**,PetscScalar**,PetscInt*);
PETSC_INTERN PetscErrorCode VecStashSortCompress_Private(VecStash*);
PETSC_INTERN PetscErrorCode VecStashGetOwnerList_Private(VecStash*,PetscLayout,PetscMPIInt*,PetscMPIInt**);

/*
  VecStashValue_Private - inserts a single value into the stash.

  Input Parameters:
  stash  - the stash
  idx    - the global of the inserted value
  values - the value inserted
*/
PETSC_STATIC_INLINE PetscErrorCode VecStashValue_Private(VecStash *stash,PetscInt row,PetscScalar value)
{
  PetscErrorCode ierr;
  /* Check and see if we have sufficient memory */
  if (((stash)->n + 1) > (stash)->nmax) {
    ierr = VecStashExpand_Private(stash,1);CHKERRQ(ierr);
  }
  (stash)->idx[(stash)->n]   = row;
  (stash)->array[(stash)->n] = value;
  (stash)->n++;
  return 0;
}

/*
  VecStashValuesBlocked_Private - inserts 1 block of values into the stash.

  Input Parameters:
  stash  - the stash
  idx    - the global block index
  values - the values inserted
*/
PETSC_STATIC_INLINE PetscErrorCode VecStashValuesBlocked_Private(VecStash *stash,PetscInt row,PetscScalar *values)
{
  PetscInt       jj,stash_bs=(stash)->bs;
  PetscScalar    *array;
  PetscErrorCode ierr;
  if (((stash)->n+1) > (stash)->nmax) {
    ierr = VecStashExpand_Private(stash,1);CHKERRQ(ierr);
  }
  array = (stash)->array + stash_bs*(stash)->n;
  (stash)->idx[(stash)->n]   = row;
  for (jj=0; jj<stash_bs; jj++) { array[jj] = values[jj];}
  (stash)->n++;
  return 0;
}

PETSC_INTERN PetscErrorCode VecStrideGather_Default(Vec,PetscInt,Vec,InsertMode);
PETSC_INTERN PetscErrorCode VecStrideScatter_Default(Vec,PetscInt,Vec,InsertMode);
PETSC_INTERN PetscErrorCode VecReciprocal_Default(Vec);
PETSC_INTERN PetscErrorCode VecStrideSubSetGather_Default(Vec,PetscInt,const PetscInt[],const PetscInt[],Vec,InsertMode);
PETSC_INTERN PetscErrorCode VecStrideSubSetScatter_Default(Vec,PetscInt,const PetscInt[],const PetscInt[],Vec,InsertMode);

#if defined(PETSC_HAVE_MATLAB_ENGINE)
PETSC_EXTERN PetscErrorCode VecMatlabEnginePut_Default(PetscObject,void*);
PETSC_EXTERN PetscErrorCode VecMatlabEngineGet_Default(PetscObject,void*);
#endif

PETSC_EXTERN PetscErrorCode PetscSectionGetField_Internal(PetscSection, PetscSection, Vec, PetscInt, PetscInt, PetscInt, IS *, Vec *);
PETSC_EXTERN PetscErrorCode PetscSectionRestoreField_Internal(PetscSection, PetscSection, Vec, PetscInt, PetscInt, PetscInt, IS *, Vec *);

#define VecCheckSameLocalSize(x,ar1,y,ar2)                                 \
  if ((x)->map->n != (y)->map->n) SETERRQ4(PETSC_COMM_SELF,PETSC_ERR_ARG_INCOMP,"Incompatible vector local lengths parameter # %d local size %D != parameter # %d local size %D", ar1,(x)->map->n, ar2,(y)->map->n);

#define VecCheckSameSize(x,ar1,y,ar2)                                      \
  if ((x)->map->N != (y)->map->N) SETERRQ4(PetscObjectComm((PetscObject)x),PETSC_ERR_ARG_INCOMP,"Incompatible vector global lengths parameter # %d global size %D != paramter # %d global size %D", ar1,(x)->map->N, ar2,(y)->map->N);\
  VecCheckSameLocalSize(x,ar1,y,ar2);
  
#define VecCheckLocalSize(x,ar1,n)                                         \
  if ((x)->map->n != n) SETERRQ3(PETSC_COMM_SELF,PETSC_ERR_ARG_INCOMP,"Incorrect vector local size: parameter # %d local size %D != %D",ar1,(x)->map->n,n);
  
#define VecCheckSize(x,ar1,n,N)                                            \
  if ((x)->map->N != N) SETERRQ(PetscObjectComm((PetscObject)x),PETSC_ERR_ARG_INCOMP,"Incorrect vector global size: parameter # %d global size %D != %D",ar1,(x)->map->N, N);\
  VecCheckLocalSize(x,ar1,n);

typedef struct _VecTaggerOps *VecTaggerOps;
struct _VecTaggerOps {
  PetscErrorCode (*create) (VecTagger);
  PetscErrorCode (*destroy) (VecTagger);
  PetscErrorCode (*setfromoptions) (PetscOptionItems*,VecTagger);
  PetscErrorCode (*setup) (VecTagger);
  PetscErrorCode (*view) (VecTagger,PetscViewer);
  PetscErrorCode (*computeboxes) (VecTagger,Vec,PetscInt *,VecTaggerBox **);
  PetscErrorCode (*computeis) (VecTagger,Vec,IS *);
};
struct _p_VecTagger {
  PETSCHEADER(struct _VecTaggerOps);
  void      *data;
  PetscInt  blocksize;
  PetscBool invert;
  PetscBool setupcalled;
};

PETSC_EXTERN PetscBool      VecTaggerRegisterAllCalled;
PETSC_EXTERN PetscErrorCode VecTaggerRegisterAll(void);
PETSC_EXTERN PetscErrorCode VecTaggerComputeIS_FromBoxes(VecTagger,Vec,IS*);
PETSC_EXTERN PetscMPIInt Petsc_Reduction_keyval;

#endif
