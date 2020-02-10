#if !defined(_PETSCSFIMPL_H)
#define _PETSCSFIMPL_H

#include <petscsf.h>
#include <petsc/private/petscimpl.h>
#include <petscviewer.h>

PETSC_EXTERN PetscLogEvent PETSCSF_SetGraph;
PETSC_EXTERN PetscLogEvent PETSCSF_SetUp;
PETSC_EXTERN PetscLogEvent PETSCSF_BcastBegin;
PETSC_EXTERN PetscLogEvent PETSCSF_BcastEnd;
PETSC_EXTERN PetscLogEvent PETSCSF_ReduceBegin;
PETSC_EXTERN PetscLogEvent PETSCSF_ReduceEnd;
PETSC_EXTERN PetscLogEvent PETSCSF_FetchAndOpBegin;
PETSC_EXTERN PetscLogEvent PETSCSF_FetchAndOpEnd;

struct _PetscSFOps {
  PetscErrorCode (*Reset)(PetscSF);
  PetscErrorCode (*Destroy)(PetscSF);
  PetscErrorCode (*SetUp)(PetscSF);
  PetscErrorCode (*SetFromOptions)(PetscOptionItems*,PetscSF);
  PetscErrorCode (*View)(PetscSF,PetscViewer);
  PetscErrorCode (*Duplicate)(PetscSF,PetscSFDuplicateOption,PetscSF);
  PetscErrorCode (*BcastBegin)(PetscSF,MPI_Datatype,const void*,void*);
  PetscErrorCode (*BcastEnd)(PetscSF,MPI_Datatype,const void*,void*);
  PetscErrorCode (*ReduceBegin)(PetscSF,MPI_Datatype,const void*,void*,MPI_Op);
  PetscErrorCode (*ReduceEnd)(PetscSF,MPI_Datatype,const void*,void*,MPI_Op);
  PetscErrorCode (*FetchAndOpBegin)(PetscSF,MPI_Datatype,void*,const void*,void*,MPI_Op);
  PetscErrorCode (*FetchAndOpEnd)(PetscSF,MPI_Datatype,void*,const void *,void *,MPI_Op);
};

struct _p_PetscSF {
  PETSCHEADER(struct _PetscSFOps);
  PetscInt        nroots;       /* Number of root vertices on current process (candidates for incoming edges) */
  PetscInt        nleaves;      /* Number of leaf vertices on current process (this process specifies a root for each leaf) */
  PetscInt        *mine;        /* Location of leaves in leafdata arrays provided to the communication routines */
  PetscInt        *mine_alloc;
  PetscInt        minleaf,maxleaf;
  PetscSFNode     *remote;      /* Remote references to roots for each local leaf */
  PetscSFNode     *remote_alloc;
  PetscInt        nranks;       /* Number of ranks owning roots connected to my leaves */
  PetscInt        ndranks;      /* Number of ranks in distinguished group holding roots connected to my leaves */
  PetscMPIInt     *ranks;       /* List of ranks referenced by "remote" */
  PetscInt        *roffset;     /* Array of length nranks+1, offset in rmine/rremote for each rank */
  PetscInt        *rmine;       /* Concatenated array holding local indices referencing each remote rank */
  PetscInt        *rremote;     /* Concatenated array holding remote indices referenced for each remote rank */
  PetscBool       degreeknown;  /* The degree is currently known, do not have to recompute */
  PetscInt        *degree;      /* Degree of each of my root vertices */
  PetscInt        *degreetmp;   /* Temporary local array for computing degree */
  PetscBool       rankorder;    /* Sort ranks for gather and scatter operations */
  MPI_Group       ingroup;      /* Group of processes connected to my roots */
  MPI_Group       outgroup;     /* Group of processes connected to my leaves */
  PetscSF         multi;        /* Internal graph used to implement gather and scatter operations */
  PetscBool       graphset;     /* Flag indicating that the graph has been set, required before calling communication routines */
  PetscBool       setupcalled;  /* Type and communication structures have been set up */

  void *data;                   /* Pointer to implementation */
};

PETSC_EXTERN PetscBool PetscSFRegisterAllCalled;
PETSC_EXTERN PetscErrorCode PetscSFRegisterAll(void);

PETSC_EXTERN PetscErrorCode MPIPetsc_Type_unwrap(MPI_Datatype,MPI_Datatype*,PetscBool*);
PETSC_EXTERN PetscErrorCode MPIPetsc_Type_compare(MPI_Datatype,MPI_Datatype,PetscBool*);
PETSC_EXTERN PetscErrorCode MPIPetsc_Type_compare_contig(MPI_Datatype,MPI_Datatype,PetscInt*);

#endif
