/*
   A star forest (SF) describes a communication pattern
*/
#if !defined(__PETSCSF_H)
#define __PETSCSF_H
#include <petscsys.h>
#include <petscsftypes.h>

PETSC_EXTERN PetscClassId PETSCSF_CLASSID;

/*J
    PetscSFType - String with the name of a PetscSF method or the creation function
       with an optional dynamic library name, for example
       http://www.mcs.anl.gov/petsc/lib.so:mysfcreate()

   Level: beginner

   Notes:
    The two approaches provided are
$     PETSCSFBASIC which uses MPI 1 message passing to perform the communication and
$     PETSCSFWINDOW which uses MPI 2 one-sided operations to perform the communication, this may be more efficient,
$                   but may not be available for all MPI distributions. In particular OpenMPI has bugs in its one-sided
$                   operations that prevent its use.

.seealso: PetscSFSetType(), PetscSF
J*/
typedef const char *PetscSFType;
#define PETSCSFBASIC  "basic"
#define PETSCSFWINDOW "window"

/*E
    PetscSFWindowSyncType - Type of synchronization for PETSCSFWINDOW

$  PETSCSF_WINDOW_SYNC_FENCE - simplest model, synchronizing across communicator
$  PETSCSF_WINDOW_SYNC_LOCK - passive model, less synchronous, requires less setup than PETSCSF_WINDOW_SYNC_ACTIVE, but may require more handshakes
$  PETSCSF_WINDOW_SYNC_ACTIVE - active model, provides most information to MPI implementation, needs to construct 2-way process groups (more setup than PETSCSF_WINDOW_SYNC_LOCK)

   Level: advanced

.seealso: PetscSFWindowSetSyncType(), PetscSFWindowGetSyncType()
E*/
typedef enum {PETSCSF_WINDOW_SYNC_FENCE,PETSCSF_WINDOW_SYNC_LOCK,PETSCSF_WINDOW_SYNC_ACTIVE} PetscSFWindowSyncType;
PETSC_EXTERN const char *const PetscSFWindowSyncTypes[];

/*E
    PetscSFDuplicateOption - Aspects to preserve when duplicating a PetscSF

$  PETSCSF_DUPLICATE_CONFONLY - configuration only, user must call PetscSFSetGraph()
$  PETSCSF_DUPLICATE_RANKS - communication ranks preserved, but different graph (allows simpler setup after calling PetscSFSetGraph())
$  PETSCSF_DUPLICATE_GRAPH - entire graph duplicated

   Level: beginner

.seealso: PetscSFDuplicate()
E*/
typedef enum {PETSCSF_DUPLICATE_CONFONLY,PETSCSF_DUPLICATE_RANKS,PETSCSF_DUPLICATE_GRAPH} PetscSFDuplicateOption;
PETSC_EXTERN const char *const PetscSFDuplicateOptions[];

PETSC_EXTERN PetscFunctionList PetscSFList;
PETSC_EXTERN PetscErrorCode PetscSFRegister(const char[],PetscErrorCode (*)(PetscSF));

PETSC_EXTERN PetscErrorCode PetscSFInitializePackage(void);
PETSC_EXTERN PetscErrorCode PetscSFFinalizePackage(void);
PETSC_EXTERN PetscErrorCode PetscSFCreate(MPI_Comm,PetscSF*);
PETSC_EXTERN PetscErrorCode PetscSFDestroy(PetscSF*);
PETSC_EXTERN PetscErrorCode PetscSFSetType(PetscSF,PetscSFType);
PETSC_EXTERN PetscErrorCode PetscSFGetType(PetscSF,PetscSFType*);
PETSC_EXTERN PetscErrorCode PetscSFView(PetscSF,PetscViewer);
PETSC_STATIC_INLINE PetscErrorCode PetscSFViewFromOptions(PetscSF A,PetscObject obj,const char name[]) {return PetscObjectViewFromOptions((PetscObject)A,obj,name);}
PETSC_EXTERN PetscErrorCode PetscSFSetUp(PetscSF);
PETSC_EXTERN PetscErrorCode PetscSFSetFromOptions(PetscSF);
PETSC_EXTERN PetscErrorCode PetscSFDuplicate(PetscSF,PetscSFDuplicateOption,PetscSF*);
PETSC_EXTERN PetscErrorCode PetscSFWindowSetSyncType(PetscSF,PetscSFWindowSyncType);
PETSC_EXTERN PetscErrorCode PetscSFWindowGetSyncType(PetscSF,PetscSFWindowSyncType*);
PETSC_EXTERN PetscErrorCode PetscSFSetRankOrder(PetscSF,PetscBool);
PETSC_EXTERN PetscErrorCode PetscSFSetGraph(PetscSF,PetscInt,PetscInt,const PetscInt*,PetscCopyMode,const PetscSFNode*,PetscCopyMode);
PETSC_EXTERN PetscErrorCode PetscSFGetGraph(PetscSF,PetscInt*,PetscInt*,const PetscInt**,const PetscSFNode**);
PETSC_EXTERN PetscErrorCode PetscSFGetLeafRange(PetscSF,PetscInt*,PetscInt*);
PETSC_EXTERN PetscErrorCode PetscSFCreateEmbeddedSF(PetscSF,PetscInt,const PetscInt*,PetscSF*);
PETSC_EXTERN PetscErrorCode PetscSFCreateEmbeddedLeafSF(PetscSF,PetscInt,const PetscInt *, PetscSF *);
PETSC_EXTERN PetscErrorCode PetscSFReset(PetscSF);
PETSC_EXTERN PetscErrorCode PetscSFSetUpRanks(PetscSF,MPI_Group);
PETSC_EXTERN PetscErrorCode PetscSFGetRanks(PetscSF,PetscInt*,const PetscMPIInt**,const PetscInt**,const PetscInt**,const PetscInt**);
PETSC_EXTERN PetscErrorCode PetscSFGetGroups(PetscSF,MPI_Group*,MPI_Group*);
PETSC_EXTERN PetscErrorCode PetscSFGetMultiSF(PetscSF,PetscSF*);
PETSC_EXTERN PetscErrorCode PetscSFCreateInverseSF(PetscSF,PetscSF*);

/* broadcasts rootdata to leafdata */
PETSC_EXTERN PetscErrorCode PetscSFBcastBegin(PetscSF,MPI_Datatype,const void*,void*)
  PetscAttrMPIPointerWithType(3,2) PetscAttrMPIPointerWithType(4,2);
PETSC_EXTERN PetscErrorCode PetscSFBcastEnd(PetscSF,MPI_Datatype,const void*,void*)
  PetscAttrMPIPointerWithType(3,2) PetscAttrMPIPointerWithType(4,2);
/* Reduce leafdata into rootdata using provided operation */
PETSC_EXTERN PetscErrorCode PetscSFReduceBegin(PetscSF,MPI_Datatype,const void*,void *,MPI_Op)
  PetscAttrMPIPointerWithType(3,2) PetscAttrMPIPointerWithType(4,2);
PETSC_EXTERN PetscErrorCode PetscSFReduceEnd(PetscSF,MPI_Datatype,const void*,void*,MPI_Op)
  PetscAttrMPIPointerWithType(3,2) PetscAttrMPIPointerWithType(4,2);
/* Atomically modifies (using provided operation) rootdata using leafdata from each leaf, value at root at time of modification is returned in leafupdate. */
PETSC_EXTERN PetscErrorCode PetscSFFetchAndOpBegin(PetscSF,MPI_Datatype,void*,const void*,void*,MPI_Op)
  PetscAttrMPIPointerWithType(3,2) PetscAttrMPIPointerWithType(4,2) PetscAttrMPIPointerWithType(5,2);
PETSC_EXTERN PetscErrorCode PetscSFFetchAndOpEnd(PetscSF,MPI_Datatype,void*,const void*,void*,MPI_Op)
  PetscAttrMPIPointerWithType(3,2) PetscAttrMPIPointerWithType(4,2) PetscAttrMPIPointerWithType(5,2);
/* Compute the degree of every root vertex (number of leaves in its star) */
PETSC_EXTERN PetscErrorCode PetscSFComputeDegreeBegin(PetscSF,const PetscInt**);
PETSC_EXTERN PetscErrorCode PetscSFComputeDegreeEnd(PetscSF,const PetscInt**);
/* Concatenate data from all leaves into roots */
PETSC_EXTERN PetscErrorCode PetscSFGatherBegin(PetscSF,MPI_Datatype,const void*,void*)
  PetscAttrMPIPointerWithType(3,2) PetscAttrMPIPointerWithType(4,2);
PETSC_EXTERN PetscErrorCode PetscSFGatherEnd(PetscSF,MPI_Datatype,const void*,void*)
  PetscAttrMPIPointerWithType(3,2) PetscAttrMPIPointerWithType(4,2);
/* Distribute distinct values to each leaf from roots */
PETSC_EXTERN PetscErrorCode PetscSFScatterBegin(PetscSF,MPI_Datatype,const void*,void*)
  PetscAttrMPIPointerWithType(3,2) PetscAttrMPIPointerWithType(4,2);
PETSC_EXTERN PetscErrorCode PetscSFScatterEnd(PetscSF,MPI_Datatype,const void*,void*)
  PetscAttrMPIPointerWithType(3,2) PetscAttrMPIPointerWithType(4,2);

PETSC_EXTERN PetscErrorCode PetscSFCompose(PetscSF,PetscSF,PetscSF*);

#if defined(MPI_REPLACE)
#  define MPIU_REPLACE MPI_REPLACE
#else
/* When using an old MPI such that MPI_REPLACE is not defined, we do not pass MPI_REPLACE to MPI at all.  Instead, we
 * use it as a flag for our own reducer in the PETSCSFBASIC implementation.  This could be any unique value unlikely to
 * collide with another MPI_Op so we'll just use the value that has been used by every version of MPICH since
 * MPICH2-1.0.6. */
#  define MPIU_REPLACE (MPI_Op)(0x5800000d)
#endif

#endif
