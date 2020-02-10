
/*
   This private file should not be included in users' code.
   Defines the fields shared by all vector scatter implementations.

*/

#ifndef __VECSCATTERIMPL_H
#define __VECSCATTERIMPL_H

#include <petscvec.h>
#include <petsc/private/petscimpl.h>
#include <petsc/private/vecimpl.h>
#include <petscviewer.h>

PETSC_EXTERN PetscBool VecScatterRegisterAllCalled;
PETSC_EXTERN PetscErrorCode VecScatterRegisterAll(void);

typedef enum { VEC_SCATTER_SEQ_GENERAL,VEC_SCATTER_SEQ_STRIDE,
               VEC_SCATTER_MPI_GENERAL,VEC_SCATTER_MPI_TOALL,
               VEC_SCATTER_MPI_TOONE} VecScatterFormat;

#define VECSCATTER_IMPL_HEADER \
      VecScatterFormat format;

typedef struct {
  VECSCATTER_IMPL_HEADER
} VecScatter_Common;

/* A plan to optimize individual memory copies (e.g., pack/unpack to/from send/recv buffers, or local scatters)
   in VecScatter. Currently, a scatter to a neighbor processor may be transformed into 1) multiple (including one)
   contiguous memory copies, e.g., memcpy; OR 2) one strided memory copies.

    For brevity, we call them memory copies. In reality, the optimization applies not only to INSERT_VALUES, but also to ADD_VALUES, etc.
 */
typedef struct {
  PetscInt  n;                /* number of processors */
  PetscBool *optimized;       /* [n] is the scatter to procs[i] optimized? */
  PetscInt  *copy_offsets;    /* [n+1] we number all copies. Scatter to procs[i] is optimized into copies in [copy_offsets[i],copy_offsets[i+1]) */
  PetscInt  *copy_starts;     /* [*] j-th copy starts at index copy_starts[j] of the vector */
  PetscInt  *copy_lengths;    /* [*] with length copy_lengths[j] in bytes */
  PetscInt  *stride_first;    /* [n] if optimized[i] is TRUE but copy_offsets[i] = copy_offsets[i+1], then scatter to procs[i] is strided. The first */
  PetscInt  *stride_step;     /* [n]   index is stride_first[i], step is stride_step[i], */
  PetscInt  *stride_n;        /* [n]   and total stride_n[i] steps */
  PetscBool same_copy_starts; /* used only by VecScatterMemcpyPlanCreate_SGToSG(). If true, to's copy_starts[] values
                                 are as same as from's. Used to quickly test if we are doing a self-copy */
} VecScatterMemcpyPlan;

/*
   These scatters are for the purely local case.
*/
typedef struct {
  VECSCATTER_IMPL_HEADER
  PetscInt       n;                    /* number of components to scatter */
  PetscInt       *vslots;              /* locations of components */
  /*
       The next three fields are used in parallel scatters, they contain
       optimization in the special case that the "to" vector and the "from"
       vector are the same, so one only needs copy components that truly
       copies instead of just y[idx[i]] = y[jdx[i]] where idx[i] == jdx[i].
  */
  PetscBool      nonmatching_computed;
  PetscInt       n_nonmatching;        /* number of "from"s  != "to"s */
  PetscInt       *slots_nonmatching;   /* locations of "from"s  != "to"s */
  VecScatterMemcpyPlan memcpy_plan;    /* a plan to optimize pack/unpack with memcpy */
} VecScatter_Seq_General;

typedef struct {
  VECSCATTER_IMPL_HEADER
  PetscInt       n;
  PetscInt       first;
  PetscInt       step;
} VecScatter_Seq_Stride;

/*
   This scatter is for a global vector copied (completely) to each processor (or all to one)
*/
typedef struct {
  VECSCATTER_IMPL_HEADER
  PetscMPIInt    *count;        /* elements of vector on each processor */
  PetscMPIInt    *displx;
  PetscScalar    *work1;
  PetscScalar    *work2;
} VecScatter_MPI_ToAll;

/*
   This is the general parallel scatter
*/
typedef struct {
  VECSCATTER_IMPL_HEADER
  PetscInt               n;        /* number of processors to send/receive */
  PetscInt               *starts;  /* starting point in indices and values for each proc*/
  PetscInt               *indices; /* list of all components sent or received */
  PetscMPIInt            *procs;   /* processors we are communicating with in scatter */
  VecScatterMemcpyPlan   memcpy_plan; /* a plan to optimize pack/unpack/scatter */
  MPI_Request            *requests,*rev_requests;
  PetscScalar            *values;  /* buffer for all sends or receives */
  VecScatter_Seq_General local;    /* any part that happens to be local */
  MPI_Status             *sstatus,*rstatus;
  PetscInt               bs;
  PetscBool              contiq;
#if defined(PETSC_HAVE_MPI_WIN_CREATE_FEATURE)      /* these uses windows for communication only within each node */
  PetscMPIInt            msize,sharedcnt;           /* total to entries that are going to processes with the same shared memory space */
  PetscScalar            *sharedspace;              /* space each process puts data to be read from other processes; allocated by MPI */
  PetscScalar            **sharedspaces;            /* [msize] space other processes put data to be read from this processes. */
  PetscInt               *sharedspacesoffset;       /* [msize] offset into sharedspaces, that I will read from */
  PetscInt               *sharedspacestarts;        /* [msize+1] for each shared memory partner this maps to the part of sharedspaceindices of that partner */
  PetscInt               *sharedspaceindices;       /* [] for each shared memory partner contains indices where values are to be copied to */
  MPI_Win                sharedwin;                 /* Window that owns sharedspace */
  PetscInt               notdone;                   /* used by VecScatterEndMPI3Node() */
#endif
} VecScatter_MPI_General;

/* Routines to create, copy, destroy or execute a memcpy plan */

/* Create a memcpy plan based on a list of indices */
PETSC_INTERN PetscErrorCode VecScatterMemcpyPlanCreate_Index(PetscInt,const PetscInt*,const PetscInt*,PetscInt,VecScatterMemcpyPlan*);
/* Create a memcpy plan for a SG (sequential general vector) to SG scatter */
PETSC_INTERN PetscErrorCode VecScatterMemcpyPlanCreate_SGToSG(PetscInt,VecScatter_Seq_General*,VecScatter_Seq_General*);
PETSC_INTERN PetscErrorCode VecScatterMemcpyPlanCopy(const VecScatterMemcpyPlan*,VecScatterMemcpyPlan*);
PETSC_INTERN PetscErrorCode VecScatterMemcpyPlanDestroy(VecScatterMemcpyPlan*);
/* Create/copy/destroy a memcpy plan for a P (parallel vector) to P scatter */
PETSC_INTERN PetscErrorCode VecScatterMemcpyPlanCreate_PtoP(VecScatter_MPI_General*,VecScatter_MPI_General*);
PETSC_INTERN PetscErrorCode VecScatterMemcpyPlanCopy_PtoP(const VecScatter_MPI_General*,const VecScatter_MPI_General*,VecScatter_MPI_General*,VecScatter_MPI_General*);
PETSC_INTERN PetscErrorCode VecScatterMemcpyPlanDestroy_PtoP(VecScatter_MPI_General*,VecScatter_MPI_General*);

/* Pack data from x to y according to the i-th memcpy plan in xplan */
PETSC_STATIC_INLINE PetscErrorCode VecScatterMemcpyPlanExecute_Pack(PetscInt i,const PetscScalar *PETSC_RESTRICT x,const VecScatterMemcpyPlan *xplan,PetscScalar *PETSC_RESTRICT y,InsertMode addv,PetscInt bs)
{
  PetscErrorCode    ierr;
  PetscInt          j,k,len,step,n;
  const PetscScalar *xv;
  PetscBool         strided;

  PetscFunctionBegin;
  strided = (xplan->copy_offsets[i] == xplan->copy_offsets[i+1]) ? PETSC_TRUE : PETSC_FALSE;
  if (strided) {
    xv   = x+xplan->stride_first[i];
    step = xplan->stride_step[i];
    n    = xplan->stride_n[i];
  }

  if (addv == INSERT_VALUES) {
    if (strided) {
      for (j=0; j<n; j++)
        for (k=0; k<bs; k++) y[j*bs+k] = xv[j*step+k];
    } else {
      for (j=xplan->copy_offsets[i]; j<xplan->copy_offsets[i+1]; j++) {
        len  = xplan->copy_lengths[j];
        ierr = PetscMemcpy(y,x+xplan->copy_starts[j],len);CHKERRQ(ierr);
        y    = (PetscScalar*)((PetscChar*)y + len);
      }
    }
  } else if (addv == ADD_VALUES) {
    if (strided) {
      for (j=0; j<n; j++)
        for (k=0; k<bs; k++) y[j*bs+k] += xv[j*step+k];
    } else {
      for (j=xplan->copy_offsets[i]; j<xplan->copy_offsets[i+1]; j++) {
        len  = xplan->copy_lengths[j]/sizeof(PetscScalar);
        xv   = x+xplan->copy_starts[i];
        for (k=0; k<len; k++) y[k] += xv[k];
        y   += len;
      }
    }
  }
#if !defined(PETSC_USE_COMPLEX)
  else if (addv == MAX_VALUES) {
    if (strided) {
      for (j=0; j<n; j++)
        for (k=0; k<bs; k++) y[j*bs+k] = PetscMax(y[j*bs+k],xv[j*step+k]);
    } else {
      for (j=xplan->copy_offsets[i]; j<xplan->copy_offsets[i+1]; j++) {
        len  = xplan->copy_lengths[j]/sizeof(PetscScalar);
        xv   = x+xplan->copy_starts[i];
        for (k=0; k<len; k++) y[k] = PetscMax(y[k],xv[k]);
        y   += len;
      }
    }
  }
#endif
  else {
    SETERRQ1(PETSC_COMM_SELF, PETSC_ERR_ARG_WRONG, "Cannot handle insert mode %d in packing",addv);
  }
  PetscFunctionReturn(0);
}

/* Unpack data from x to y according to the i-th memcpy plan in yplan */
PETSC_STATIC_INLINE PetscErrorCode VecScatterMemcpyPlanExecute_Unpack(PetscInt i,const PetscScalar *PETSC_RESTRICT x,PetscScalar *PETSC_RESTRICT y,const VecScatterMemcpyPlan *yplan,InsertMode addv,PetscInt bs)
{
  PetscErrorCode ierr;
  PetscInt       j,k,len,step,n;
  PetscScalar    *yv;
  PetscBool      strided;

  PetscFunctionBegin;
  strided = (yplan->copy_offsets[i] == yplan->copy_offsets[i+1]) ? PETSC_TRUE : PETSC_FALSE;
  if (strided) {
    yv   = y+yplan->stride_first[i];
    step = yplan->stride_step[i];
    n    = yplan->stride_n[i];
  }

  if (addv == INSERT_VALUES) {
    if (strided) {
      for (j=0; j<n; j++)
        for (k=0; k<bs; k++) yv[j*step+k] = x[j*bs+k];
    } else {
      for (j=yplan->copy_offsets[i]; j<yplan->copy_offsets[i+1]; j++) {
        len  = yplan->copy_lengths[j];
        ierr = PetscMemcpy(y+yplan->copy_starts[j],x,len);CHKERRQ(ierr);
        x    = (PetscScalar*)((PetscChar*)x + len);
      }
    }
  } else if (addv == ADD_VALUES) {
    if (strided) {
      for (j=0; j<n; j++)
        for (k=0; k<bs; k++) yv[j*step+k] += x[j*bs+k];
    } else {
      for (j=yplan->copy_offsets[i]; j<yplan->copy_offsets[i+1]; j++) {
        len  = yplan->copy_lengths[j]/sizeof(PetscScalar);
        yv   = y+yplan->copy_starts[j];
        for (k=0; k<len; k++) yv[k] += x[k];
        x   += len;
      }
    }
  }
#if !defined(PETSC_USE_COMPLEX)
  else if (addv == MAX_VALUES) {
    if (strided) {
      for (j=0; j<n; j++)
        for (k=0; k<bs; k++) yv[j*step+k] = PetscMax(yv[j*step+k],x[j*bs+k]);
    } else {
      for (j=yplan->copy_offsets[i]; j<yplan->copy_offsets[i+1]; j++) {
        len  = yplan->copy_lengths[j]/sizeof(PetscScalar);
        yv   = y+yplan->copy_starts[j];
        for (k=0; k<len; k++) yv[k] = PetscMax(yv[k],x[k]);
        x   += len;
      }
    }
  }
#endif
  else {
    SETERRQ1(PETSC_COMM_SELF, PETSC_ERR_ARG_WRONG, "Cannot handle insert mode %d in unpacking",addv);
  }
  PetscFunctionReturn(0);
}

/* Scatter data from piece-wise contiguous x to (conforming) piece-wise contiguous y according to the i-th memcpy plan in xplan and yplan respectively */
PETSC_STATIC_INLINE PetscErrorCode VecScatterMemcpyPlanExecute_Scatter(PetscInt i,const PetscScalar *PETSC_RESTRICT x,const VecScatterMemcpyPlan *xplan,PetscScalar *PETSC_RESTRICT y,const VecScatterMemcpyPlan *yplan,InsertMode addv)
{
  PetscErrorCode    ierr;
  PetscInt          j,k,len;
  const PetscScalar *xv;
  PetscScalar       *yv;

  PetscFunctionBegin;
  if (addv == INSERT_VALUES) {
    for (j=xplan->copy_offsets[i]; j<xplan->copy_offsets[i+1]; j++) {
      ierr = PetscMemcpy(y+yplan->copy_starts[j],x+xplan->copy_starts[j],xplan->copy_lengths[j]);CHKERRQ(ierr);
    }
  } else if (addv == ADD_VALUES) {
    for (j=xplan->copy_offsets[i]; j<xplan->copy_offsets[i+1]; j++) {
      len = xplan->copy_lengths[j]/sizeof(PetscScalar);
      xv  = x+xplan->copy_starts[j];
      yv  = y+yplan->copy_starts[j];
      for (k=0; k<len; k++) yv[k] += xv[k];
    }
  }
#if !defined(PETSC_USE_COMPLEX)
  else if (addv == MAX_VALUES) {
    for (j=xplan->copy_offsets[i]; j<xplan->copy_offsets[i+1]; j++) {
      len = xplan->copy_lengths[j]/sizeof(PetscScalar);
      xv  = x+xplan->copy_starts[j];
      yv  = y+yplan->copy_starts[j];
      for (k=0; k<len; k++) yv[k] = PetscMax(yv[k],xv[k]);
    }
  }
#endif
  else {
    SETERRQ1(PETSC_COMM_SELF, PETSC_ERR_ARG_WRONG, "Cannot handle insert mode %d in scattering",addv);
  }
  PetscFunctionReturn(0);
}

PETSC_INTERN PetscErrorCode VecScatterGetTypes_Private(VecScatter,VecScatterFormat*,VecScatterFormat*);
PETSC_INTERN PetscErrorCode VecScatterIsSequential_Private(VecScatter_Common*,PetscBool*);

typedef struct _VecScatterOps *VecScatterOps;
struct _VecScatterOps {
  PetscErrorCode (*begin)(VecScatter,Vec,Vec,InsertMode,ScatterMode);
  PetscErrorCode (*end)(VecScatter,Vec,Vec,InsertMode,ScatterMode);
  PetscErrorCode (*copy)(VecScatter,VecScatter);
  PetscErrorCode (*destroy)(VecScatter);
  PetscErrorCode (*view)(VecScatter,PetscViewer);
  PetscErrorCode (*viewfromoptions)(VecScatter,const char prefix[],const char name[]);
  PetscErrorCode (*remap)(VecScatter,PetscInt *,PetscInt*);
  PetscErrorCode (*getmerged)(VecScatter,PetscBool *);
};

struct _p_VecScatter {
  PETSCHEADER(struct _VecScatterOps);
  PetscInt       to_n,from_n;
  PetscBool      inuse;                /* prevents corruption from mixing two scatters */
  PetscBool      beginandendtogether;  /* indicates that the scatter begin and end  function are called together, VecScatterEnd() is then treated as a nop */
  void           *fromdata,*todata;
  void           *spptr;
  PetscBool      is_duplicate;         /* IS has duplicate indices, would cause writing error in the case StoP of VecScatterEndMPI3Node */
  Vec            to_v,from_v;          /* used in VecScatterCreate() */
  IS             to_is,from_is;        /* used in VecScatterCreate() */
};

PETSC_INTERN PetscErrorCode VecScatterCreate_Seq(VecScatter);
PETSC_INTERN PetscErrorCode VecScatterCreate_MPI1(VecScatter);
PETSC_INTERN PetscErrorCode VecScatterCreate_MPI3(VecScatter);
PETSC_INTERN PetscErrorCode VecScatterCreate_MPI3Node(VecScatter);

PETSC_EXTERN PetscMPIInt Petsc_Reduction_keyval;

#endif
