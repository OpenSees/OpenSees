/*
    Defines profile/logging in PETSc.
*/

#if !defined(__PetscLog_H)
#define __PetscLog_H
#include <petscsys.h>

/* General logging of information; different from event logging */
PETSC_EXTERN PetscErrorCode PetscInfo_Private(const char[],void*,const char[],...);
#if defined(PETSC_USE_INFO)
#define PetscInfo(A,S)                       PetscInfo_Private(PETSC_FUNCTION_NAME,A,S)
#define PetscInfo1(A,S,a1)                   PetscInfo_Private(PETSC_FUNCTION_NAME,A,S,a1)
#define PetscInfo2(A,S,a1,a2)                PetscInfo_Private(PETSC_FUNCTION_NAME,A,S,a1,a2)
#define PetscInfo3(A,S,a1,a2,a3)             PetscInfo_Private(PETSC_FUNCTION_NAME,A,S,a1,a2,a3)
#define PetscInfo4(A,S,a1,a2,a3,a4)          PetscInfo_Private(PETSC_FUNCTION_NAME,A,S,a1,a2,a3,a4)
#define PetscInfo5(A,S,a1,a2,a3,a4,a5)       PetscInfo_Private(PETSC_FUNCTION_NAME,A,S,a1,a2,a3,a4,a5)
#define PetscInfo6(A,S,a1,a2,a3,a4,a5,a6)    PetscInfo_Private(PETSC_FUNCTION_NAME,A,S,a1,a2,a3,a4,a5,a6)
#define PetscInfo7(A,S,a1,a2,a3,a4,a5,a6,a7) PetscInfo_Private(PETSC_FUNCTION_NAME,A,S,a1,a2,a3,a4,a5,a6,a7)
#else
#define PetscInfo(A,S)                       0
#define PetscInfo1(A,S,a1)                   0
#define PetscInfo2(A,S,a1,a2)                0
#define PetscInfo3(A,S,a1,a2,a3)             0
#define PetscInfo4(A,S,a1,a2,a3,a4)          0
#define PetscInfo5(A,S,a1,a2,a3,a4,a5)       0
#define PetscInfo6(A,S,a1,a2,a3,a4,a5,a6)    0
#define PetscInfo7(A,S,a1,a2,a3,a4,a5,a6,a7) 0
#endif
PETSC_EXTERN PetscErrorCode PetscInfoDeactivateClass(PetscClassId);
PETSC_EXTERN PetscErrorCode PetscInfoActivateClass(PetscClassId);
PETSC_EXTERN PetscBool PetscLogPrintInfo;  /* if true, indicates PetscInfo() is turned on */

/*MC
    PetscLogEvent - id used to identify PETSc or user events which timed portions (blocks of executable)
     code.

    Level: intermediate

.seealso: PetscLogEventRegister(), PetscLogEventBegin(), PetscLogEventEnd(), PetscLogStage
M*/
typedef int PetscLogEvent;

/*MC
    PetscLogStage - id used to identify user stages (phases, sections) of runs - for logging

    Level: intermediate

.seealso: PetscLogStageRegister(), PetscLogStagePush(), PetscLogStagePop(), PetscLogEvent
M*/
typedef int PetscLogStage;

#define PETSC_EVENT  1311311
PETSC_EXTERN PetscLogEvent PETSC_LARGEST_EVENT;

/* Global flop counter */
PETSC_EXTERN PetscLogDouble petsc_TotalFlops;
PETSC_EXTERN PetscLogDouble petsc_tmp_flops;

/* We must make the following structures available to access the event
     activation flags in the PetscLogEventBegin/End() macros. These are not part of the PETSc public
     API and are not intended to be used by other parts of PETSc or by users.

     The code that manipulates these structures is in src/sys/logging/utils.
*/
typedef struct _n_PetscIntStack *PetscIntStack;

/* -----------------------------------------------------------------------------------------------------*/
/*
    PetscClassRegInfo, PetscClassPerfInfo - Each class has two data structures associated with it. The first has
       static information about it, the second collects statistics on how many objects of the class are created,
       how much memory they use, etc.

    PetscClassRegLog, PetscClassPerfLog - arrays of the PetscClassRegInfo and PetscClassPerfInfo for all classes.
*/
typedef struct  {
  char           *name;   /* The class name */
  PetscClassId   classid; /* The integer identifying this class */
} PetscClassRegInfo;

typedef struct {
  PetscClassId   id;           /* The integer identifying this class */
  int            creations;    /* The number of objects of this class created */
  int            destructions; /* The number of objects of this class destroyed */
  PetscLogDouble mem;          /* The total memory allocated by objects of this class */
  PetscLogDouble descMem;      /* The total memory allocated by descendents of these objects */
} PetscClassPerfInfo;

typedef struct _n_PetscClassRegLog *PetscClassRegLog;
struct _n_PetscClassRegLog {
  int               numClasses; /* The number of classes registered */
  int               maxClasses; /* The maximum number of classes */
  PetscClassRegInfo *classInfo; /* The structure for class information (classids are monotonicly increasing) */
};

typedef struct _n_PetscClassPerfLog *PetscClassPerfLog;
struct _n_PetscClassPerfLog {
  int                numClasses; /* The number of logging classes */
  int                maxClasses; /* The maximum number of classes */
  PetscClassPerfInfo *classInfo; /* The structure for class information (classids are monotonicly increasing) */
};
/* -----------------------------------------------------------------------------------------------------*/
/*
    PetscEventRegInfo, PetscEventPerfInfo - Each event has two data structures associated with it. The first has
       static information about it, the second collects statistics on how many times the event is used, how
       much time it takes, etc.

    PetscEventRegLog, PetscEventPerfLog - an array of all PetscEventRegInfo and PetscEventPerfInfo for all events. There is one
      of these for each stage.

*/
typedef struct {
  char         *name;         /* The name of this event */
  PetscClassId classid;       /* The class the event is associated with */
  PetscBool    collective;    /* Flag this event as collective */
#if defined (PETSC_HAVE_MPE)
  int          mpe_id_begin;  /* MPE IDs that define the event */
  int          mpe_id_end;
#endif
} PetscEventRegInfo;

typedef struct {
  int            id;            /* The integer identifying this event */
  PetscBool      active;        /* The flag to activate logging */
  PetscBool      visible;       /* The flag to print info in summary */
  int            depth;         /* The nesting depth of the event call */
  int            count;         /* The number of times this event was executed */
  PetscLogDouble flops, flops2, flopsTmp; /* The flops and flops^2 used in this event */
  PetscLogDouble time, time2, timeTmp;    /* The time and time^2 taken for this event */
  PetscLogDouble syncTime;                /* The synchronization barrier time */
  PetscLogDouble numMessages;   /* The number of messages in this event */
  PetscLogDouble messageLength; /* The total message lengths in this event */
  PetscLogDouble numReductions; /* The number of reductions in this event */
} PetscEventPerfInfo;

typedef struct _n_PetscEventRegLog *PetscEventRegLog;
struct _n_PetscEventRegLog {
  int               numEvents;  /* The number of registered events */
  int               maxEvents;  /* The maximum number of events */
  PetscEventRegInfo *eventInfo; /* The registration information for each event */
};

typedef struct _n_PetscEventPerfLog *PetscEventPerfLog;
struct _n_PetscEventPerfLog {
  int                numEvents;  /* The number of logging events */
  int                maxEvents;  /* The maximum number of events */
  PetscEventPerfInfo *eventInfo; /* The performance information for each event */
};
/* ------------------------------------------------------------------------------------------------------------*/
/*
   PetscStageInfo - Contains all the information about a particular stage.

   PetscStageLog - An array of PetscStageInfo for each registered stage. There is a single one of these in the code.
*/
typedef struct _PetscStageInfo {
  char               *name;     /* The stage name */
  PetscBool          used;      /* The stage was pushed on this processor */
  PetscEventPerfInfo perfInfo;  /* The stage performance information */
  PetscEventPerfLog  eventLog;  /* The event information for this stage */
  PetscClassPerfLog  classLog;  /* The class information for this stage */
} PetscStageInfo;

typedef struct _n_PetscStageLog *PetscStageLog;
struct _n_PetscStageLog {
  int              numStages;   /* The number of registered stages */
  int              maxStages;   /* The maximum number of stages */
  PetscIntStack    stack;       /* The stack for active stages */
  int              curStage;    /* The current stage (only used in macros so we don't call PetscIntStackTop) */
  PetscStageInfo   *stageInfo;  /* The information for each stage */
  PetscEventRegLog eventLog;    /* The registered events */
  PetscClassRegLog classLog;    /* The registered classes */
};
/* -----------------------------------------------------------------------------------------------------*/

PETSC_EXTERN PetscErrorCode PetscLogObjectParent(PetscObject,PetscObject);
PETSC_EXTERN PetscErrorCode PetscLogObjectMemory(PetscObject,PetscLogDouble);


#if defined(PETSC_USE_LOG)  /* --- Logging is turned on --------------------------------*/
PETSC_EXTERN PetscStageLog petsc_stageLog;
PETSC_EXTERN PetscErrorCode PetscLogGetStageLog(PetscStageLog*);
PETSC_EXTERN PetscErrorCode PetscStageLogGetCurrent(PetscStageLog,int*);
PETSC_EXTERN PetscErrorCode PetscStageLogGetEventPerfLog(PetscStageLog,int,PetscEventPerfLog*);

/*
   Flop counting:  We count each arithmetic operation (e.g., addition, multiplication) separately.

   For the complex numbers version, note that
       1 complex addition = 2 flops
       1 complex multiplication = 6 flops,
   where we define 1 flop as that for a double precision scalar.  We roughly approximate
   flop counting for complex numbers by multiplying the total flops by 4; this corresponds
   to the assumption that we're counting mostly additions and multiplications -- and
   roughly the same number of each.  More accurate counting could be done by distinguishing
   among the various arithmetic operations.
 */

#if defined(PETSC_USE_COMPLEX)
#define PETSC_FLOPS_PER_OP 4.0
#else
#define PETSC_FLOPS_PER_OP 1.0
#endif

PETSC_STATIC_INLINE PetscErrorCode PetscLogFlops(PetscLogDouble n)
{
  PetscFunctionBegin;
#if defined(PETSC_USE_DEBUG)
  if (n < 0) SETERRQ(PETSC_COMM_SELF,PETSC_ERR_ARG_OUTOFRANGE,"Cannot log negative flops");
#endif
  petsc_TotalFlops += PETSC_FLOPS_PER_OP*n;
  PetscFunctionReturn(0);
}

PETSC_EXTERN PetscErrorCode PetscGetFlops(PetscLogDouble *);

#if defined (PETSC_HAVE_MPE)
PETSC_EXTERN PetscErrorCode PetscLogMPEBegin(void);
PETSC_EXTERN PetscErrorCode PetscLogMPEDump(const char[]);
#endif

PETSC_EXTERN PetscErrorCode (*PetscLogPLB)(PetscLogEvent,int,PetscObject,PetscObject,PetscObject,PetscObject);
PETSC_EXTERN PetscErrorCode (*PetscLogPLE)(PetscLogEvent,int,PetscObject,PetscObject,PetscObject,PetscObject);
PETSC_EXTERN PetscErrorCode (*PetscLogPHC)(PetscObject);
PETSC_EXTERN PetscErrorCode (*PetscLogPHD)(PetscObject);

#define PetscLogObjectParents(p,n,d)  0;{int _i; for (_i=0; _i<(n); _i++) {ierr = PetscLogObjectParent((PetscObject)(p),(PetscObject)(d)[_i]);CHKERRQ(ierr);}}
#define PetscLogObjectCreate(h)      ((PetscLogPHC) ? (*PetscLogPHC)((PetscObject)(h)) : 0)
#define PetscLogObjectDestroy(h)     ((PetscLogPHD) ? (*PetscLogPHD)((PetscObject)(h)) : 0)
PETSC_EXTERN PetscErrorCode PetscLogObjectState(PetscObject, const char[], ...);

/* Initialization functions */
PETSC_EXTERN PetscErrorCode PetscLogDefaultBegin(void);
PETSC_EXTERN PetscErrorCode PetscLogAllBegin(void);
PETSC_EXTERN PetscErrorCode PetscLogNestedBegin(void);
PETSC_EXTERN PetscErrorCode PetscLogTraceBegin(FILE *);
PETSC_EXTERN PetscErrorCode PetscLogActions(PetscBool);
PETSC_EXTERN PetscErrorCode PetscLogObjects(PetscBool);
PETSC_EXTERN PetscErrorCode PetscLogSetThreshold(PetscLogDouble,PetscLogDouble*);
PETSC_EXTERN PetscErrorCode PetscLogSet(PetscErrorCode (*)(int, int, PetscObject, PetscObject, PetscObject, PetscObject),
                                        PetscErrorCode (*)(int, int, PetscObject, PetscObject, PetscObject, PetscObject));

/* Output functions */
PETSC_EXTERN PetscErrorCode PetscLogView(PetscViewer);
PETSC_EXTERN PetscErrorCode PetscLogViewFromOptions(void);
PETSC_EXTERN PetscErrorCode PetscLogDump(const char[]);

/* Stage functions */
PETSC_EXTERN PetscErrorCode PetscLogStageRegister(const char[],PetscLogStage*);
PETSC_EXTERN PetscErrorCode PetscLogStagePush(PetscLogStage);
PETSC_EXTERN PetscErrorCode PetscLogStagePop(void);
PETSC_EXTERN PetscErrorCode PetscLogStageSetActive(PetscLogStage,PetscBool);
PETSC_EXTERN PetscErrorCode PetscLogStageGetActive(PetscLogStage,PetscBool*);
PETSC_EXTERN PetscErrorCode PetscLogStageSetVisible(PetscLogStage,PetscBool);
PETSC_EXTERN PetscErrorCode PetscLogStageGetVisible(PetscLogStage,PetscBool*);
PETSC_EXTERN PetscErrorCode PetscLogStageGetId(const char[],PetscLogStage*);

/* Event functions */
PETSC_EXTERN PetscErrorCode PetscLogEventRegister(const char[],PetscClassId,PetscLogEvent*);
PETSC_EXTERN PetscErrorCode PetscLogEventSetCollective(PetscLogEvent,PetscBool);
PETSC_EXTERN PetscErrorCode PetscLogEventIncludeClass(PetscClassId);
PETSC_EXTERN PetscErrorCode PetscLogEventExcludeClass(PetscClassId);
PETSC_EXTERN PetscErrorCode PetscLogEventActivate(PetscLogEvent);
PETSC_EXTERN PetscErrorCode PetscLogEventDeactivate(PetscLogEvent);
PETSC_EXTERN PetscErrorCode PetscLogEventSetActiveAll(PetscLogEvent,PetscBool);
PETSC_EXTERN PetscErrorCode PetscLogEventActivateClass(PetscClassId);
PETSC_EXTERN PetscErrorCode PetscLogEventDeactivateClass(PetscClassId);
PETSC_EXTERN PetscErrorCode PetscLogEventGetId(const char[],PetscLogEvent*);
PETSC_EXTERN PetscErrorCode PetscLogEventGetPerfInfo(int,PetscLogEvent,PetscEventPerfInfo*);

/* Global counters */
PETSC_EXTERN PetscLogDouble petsc_irecv_ct;
PETSC_EXTERN PetscLogDouble petsc_isend_ct;
PETSC_EXTERN PetscLogDouble petsc_recv_ct;
PETSC_EXTERN PetscLogDouble petsc_send_ct;
PETSC_EXTERN PetscLogDouble petsc_irecv_len;
PETSC_EXTERN PetscLogDouble petsc_isend_len;
PETSC_EXTERN PetscLogDouble petsc_recv_len;
PETSC_EXTERN PetscLogDouble petsc_send_len;
PETSC_EXTERN PetscLogDouble petsc_allreduce_ct;
PETSC_EXTERN PetscLogDouble petsc_gather_ct;
PETSC_EXTERN PetscLogDouble petsc_scatter_ct;
PETSC_EXTERN PetscLogDouble petsc_wait_ct;
PETSC_EXTERN PetscLogDouble petsc_wait_any_ct;
PETSC_EXTERN PetscLogDouble petsc_wait_all_ct;
PETSC_EXTERN PetscLogDouble petsc_sum_of_waits_ct;

PETSC_EXTERN PetscBool PetscLogSyncOn;  /* true if logging synchronization is enabled */

#define PetscLogEventBegin(e,o1,o2,o3,o4) \
  (((PetscLogPLB && petsc_stageLog->stageInfo[petsc_stageLog->curStage].perfInfo.active && petsc_stageLog->stageInfo[petsc_stageLog->curStage].eventLog->eventInfo[e].active) ? \
    (*PetscLogPLB)((e),0,(PetscObject)(o1),(PetscObject)(o2),(PetscObject)(o3),(PetscObject)(o4)) : 0 ))

#define PetscLogEventEnd(e,o1,o2,o3,o4) \
  (((PetscLogPLE && petsc_stageLog->stageInfo[petsc_stageLog->curStage].perfInfo.active && petsc_stageLog->stageInfo[petsc_stageLog->curStage].eventLog->eventInfo[e].active) ? \
    (*PetscLogPLE)((e),0,(PetscObject)(o1),(PetscObject)(o2),(PetscObject)(o3),(PetscObject)(o4)) : 0 ))

PETSC_EXTERN PetscErrorCode PetscLogEventGetFlops(PetscLogEvent,PetscLogDouble*);
PETSC_EXTERN PetscErrorCode PetscLogEventZeroFlops(PetscLogEvent);

/*
     These are used internally in the PETSc routines to keep a count of MPI messages and
   their sizes.

     This does not work for MPI-Uni because our include/petsc/mpiuni/mpi.h file
   uses macros to defined the MPI operations.

     It does not work correctly from HP-UX because it processes the
   macros in a way that sometimes it double counts, hence
   PETSC_HAVE_BROKEN_RECURSIVE_MACRO

     It does not work with Windows because winmpich lacks MPI_Type_size()
*/
#if !defined(__MPIUNI_H) && !defined(PETSC_HAVE_BROKEN_RECURSIVE_MACRO) && !defined (PETSC_HAVE_MPI_MISSING_TYPESIZE)
/*
   Logging of MPI activities
*/
PETSC_STATIC_INLINE PetscErrorCode PetscMPITypeSize(PetscLogDouble *buff,PetscMPIInt count,MPI_Datatype type)
{
  PetscMPIInt mysize;
  PetscErrorCode _myierr;
  if (type == MPI_DATATYPE_NULL) return 0;
  _myierr = MPI_Type_size(type,&mysize);CHKERRQ(_myierr);
  *buff += (PetscLogDouble) (count*mysize);
  return 0;
}

PETSC_STATIC_INLINE PetscErrorCode PetscMPITypeSizeComm(MPI_Comm comm, PetscLogDouble *buff,PetscMPIInt *counts,MPI_Datatype type)
{
  PetscMPIInt mysize, commsize, p;
  PetscErrorCode _myierr;

  if (type == MPI_DATATYPE_NULL) return 0;
  _myierr = MPI_Comm_size(comm,&commsize);CHKERRQ(_myierr);
  _myierr = MPI_Type_size(type,&mysize);CHKERRQ(_myierr);
  for (p = 0; p < commsize; ++p) {
    *buff += (PetscLogDouble) (counts[p]*mysize);
  }
  return 0;
}

/*
    Returns 1 if the communicator is parallel else zero
*/
PETSC_STATIC_INLINE int PetscMPIParallelComm(MPI_Comm comm)
{
  PetscMPIInt size; MPI_Comm_size(comm,&size); return size > 1;
}

#define MPI_Irecv(buf,count,datatype,source,tag,comm,request) \
  ((petsc_irecv_ct++,0) || PetscMPITypeSize(&(petsc_irecv_len),(count),(datatype)) || MPI_Irecv((buf),(count),(datatype),(source),(tag),(comm),(request)))

#define MPI_Isend(buf,count,datatype,dest,tag,comm,request) \
  ((petsc_isend_ct++,0) || PetscMPITypeSize(&(petsc_isend_len),(count),(datatype)) || MPI_Isend((buf),(count),(datatype),(dest),(tag),(comm),(request)))

#define MPI_Startall_irecv(count,number,requests) \
  ((petsc_irecv_ct += (PetscLogDouble)(number),0) || PetscMPITypeSize(&(petsc_irecv_len),(count),(MPIU_SCALAR)) || MPI_Startall((number),(requests)))

#define MPI_Startall_isend(count,number,requests) \
  ((petsc_isend_ct += (PetscLogDouble)(number),0) || PetscMPITypeSize(&(petsc_isend_len),(count),(MPIU_SCALAR)) || MPI_Startall((number),(requests)))

#define MPI_Start_isend(count,requests) \
  ((petsc_isend_ct++,0) || PetscMPITypeSize((&petsc_isend_len),(count),(MPIU_SCALAR)) || MPI_Start((requests)))

#define MPI_Recv(buf,count,datatype,source,tag,comm,status) \
  ((petsc_recv_ct++,0) || PetscMPITypeSize((&petsc_recv_len),(count),(datatype)) || MPI_Recv((buf),(count),(datatype),(source),(tag),(comm),(status)))

#define MPI_Send(buf,count,datatype,dest,tag,comm) \
  ((petsc_send_ct++,0) || PetscMPITypeSize((&petsc_send_len),(count),(datatype)) || MPI_Send((buf),(count),(datatype),(dest),(tag),(comm)))

#define MPI_Wait(request,status) \
  ((petsc_wait_ct++,petsc_sum_of_waits_ct++,0) || MPI_Wait((request),(status)))

#define MPI_Waitany(a,b,c,d) \
  ((petsc_wait_any_ct++,petsc_sum_of_waits_ct++,0) || MPI_Waitany((a),(b),(c),(d)))

#define MPI_Waitall(count,array_of_requests,array_of_statuses) \
  ((petsc_wait_all_ct++,petsc_sum_of_waits_ct += (PetscLogDouble) (count),0) || MPI_Waitall((count),(array_of_requests),(array_of_statuses)))

#define MPI_Allreduce(sendbuf,recvbuf,count,datatype,op,comm) \
  ((petsc_allreduce_ct += PetscMPIParallelComm((comm)),0) || MPI_Allreduce((sendbuf),(recvbuf),(count),(datatype),(op),(comm)))

#define MPI_Bcast(buffer,count,datatype,root,comm) \
  ((petsc_allreduce_ct += PetscMPIParallelComm((comm)),0) || MPI_Bcast((buffer),(count),(datatype),(root),(comm)))

#define MPI_Reduce_scatter_block(sendbuf,recvbuf,recvcount,datatype,op,comm) \
  ((petsc_allreduce_ct += PetscMPIParallelComm((comm)),0) || MPI_Reduce_scatter_block((sendbuf),(recvbuf),(recvcount),(datatype),(op),(comm)))

#define MPI_Alltoall(sendbuf,sendcount,sendtype,recvbuf,recvcount,recvtype,comm) \
  ((petsc_allreduce_ct += PetscMPIParallelComm((comm)),0) || PetscMPITypeSize((&petsc_send_len),(sendcount),(sendtype)) || MPI_Alltoall((sendbuf),(sendcount),(sendtype),(recvbuf),(recvcount),(recvtype),(comm)))

#define MPI_Alltoallv(sendbuf,sendcnts,sdispls,sendtype,recvbuf,recvcnts,rdispls,recvtype,comm) \
  ((petsc_allreduce_ct += PetscMPIParallelComm((comm)),0) || PetscMPITypeSizeComm((comm),(&petsc_send_len),(sendcnts),(sendtype)) || MPI_Alltoallv((sendbuf),(sendcnts),(sdispls),(sendtype),(recvbuf),(recvcnts),(rdispls),(recvtype),(comm)))

#define MPI_Allgather(sendbuf,sendcount,sendtype,recvbuf,recvcount,recvtype,comm) \
  ((petsc_gather_ct += PetscMPIParallelComm((comm)),0) || MPI_Allgather((sendbuf),(sendcount),(sendtype),(recvbuf),(recvcount),(recvtype),(comm)))

#define MPI_Allgatherv(sendbuf,sendcount,sendtype,recvbuf,recvcount,displs,recvtype,comm) \
  ((petsc_gather_ct += PetscMPIParallelComm((comm)),0) || MPI_Allgatherv((sendbuf),(sendcount),(sendtype),(recvbuf),(recvcount),(displs),(recvtype),(comm)))

#define MPI_Gather(sendbuf,sendcount,sendtype,recvbuf,recvcount,recvtype,root,comm) \
  ((petsc_gather_ct++,0) || PetscMPITypeSize((&petsc_send_len),(sendcount),(sendtype)) || MPI_Gather((sendbuf),(sendcount),(sendtype),(recvbuf),(recvcount),(recvtype),(root),(comm)))

#define MPI_Gatherv(sendbuf,sendcount,sendtype,recvbuf,recvcount,displs,recvtype,root,comm) \
  ((petsc_gather_ct++,0) || PetscMPITypeSize((&petsc_send_len),(sendcount),(sendtype)) || MPI_Gatherv((sendbuf),(sendcount),(sendtype),(recvbuf),(recvcount),(displs),(recvtype),(root),(comm)))

#define MPI_Scatter(sendbuf,sendcount,sendtype,recvbuf,recvcount,recvtype,root,comm) \
  ((petsc_scatter_ct++,0) || PetscMPITypeSize((&petsc_recv_len),(recvcount),(recvtype)) || MPI_Scatter((sendbuf),(sendcount),(sendtype),(recvbuf),(recvcount),(recvtype),(root),(comm)))

#define MPI_Scatterv(sendbuf,sendcount,displs,sendtype,recvbuf,recvcount,recvtype,root,comm) \
  ((petsc_scatter_ct++,0) || PetscMPITypeSize((&petsc_recv_len),(recvcount),(recvtype)) || MPI_Scatterv((sendbuf),(sendcount),(displs),(sendtype),(recvbuf),(recvcount),(recvtype),(root),(comm)))

#else

#define MPI_Startall_irecv(count,number,requests) \
  (MPI_Startall((number),(requests)))

#define MPI_Startall_isend(count,number,requests) \
  (MPI_Startall((number),(requests)))

#define MPI_Start_isend(count,requests) \
  (MPI_Start((requests)))

#endif /* !__MPIUNI_H && ! PETSC_HAVE_BROKEN_RECURSIVE_MACRO */

#else  /* ---Logging is turned off --------------------------------------------*/

#define PetscLogFlops(n)                   0
#define PetscGetFlops(a)                   (*(a) = 0.0,0)

#define PetscLogStageRegister(a,b)         0
#define PetscLogStagePush(a)               0
#define PetscLogStagePop()                 0
#define PetscLogStageSetActive(a,b)        0
#define PetscLogStageGetActive(a,b)        0
#define PetscLogStageGetVisible(a,b)       0
#define PetscLogStageSetVisible(a,b)       0
#define PetscLogStageGetId(a,b)            (*(b)=0,0)

#define PetscLogEventRegister(a,b,c)       0
#define PetscLogEventSetCollective(a,b)    0
#define PetscLogEventIncludeClass(a)       0
#define PetscLogEventExcludeClass(a)       0
#define PetscLogEventActivate(a)           0
#define PetscLogEventDeactivate(a)         0
#define PetscLogEventActivateClass(a)      0
#define PetscLogEventDeactivateClass(a)    0
#define PetscLogEventSetActiveAll(a,b)     0
#define PetscLogEventGetId(a,b)            (*(b)=0,0)
#define PetscLogEventGetPerfInfo(a,b,c)    0

#define PetscLogPLB                        0
#define PetscLogPLE                        0
#define PetscLogPHC                        0
#define PetscLogPHD                        0

#define PetscLogObjectParents(p,n,c)       0
#define PetscLogObjectCreate(h)            0
#define PetscLogObjectDestroy(h)           0
PETSC_EXTERN PetscErrorCode PetscLogObjectState(PetscObject,const char[],...);

#define PetscLogDefaultBegin()             0
#define PetscLogAllBegin()                 0
#define PetscLogNestedBegin()              0
#define PetscLogTraceBegin(file)           0
#define PetscLogActions(a)                 0
#define PetscLogObjects(a)                 0
#define PetscLogSetThreshold(a,b)          0
#define PetscLogSet(lb,le)                 0

#define PetscLogView(viewer)               0
#define PetscLogViewFromOptions()          0
#define PetscLogDump(c)                    0

#define PetscLogEventBegin(e,o1,o2,o3,o4)  0
#define PetscLogEventEnd(e,o1,o2,o3,o4)    0

/* If PETSC_USE_LOG is NOT defined, these still need to be! */
#define MPI_Startall_irecv(count,number,requests) MPI_Startall(number,requests)
#define MPI_Startall_isend(count,number,requests) MPI_Startall(number,requests)
#define MPI_Start_isend(count,requests)           MPI_Start(requests)

#endif   /* PETSC_USE_LOG */

#define PetscPreLoadBegin(flag,name) \
do {\
  PetscBool      PetscPreLoading = flag;\
  int            PetscPreLoadMax,PetscPreLoadIt;\
  PetscLogStage  _stageNum;\
  PetscErrorCode _3_ierr; \
  _3_ierr = PetscOptionsGetBool(NULL,NULL,"-preload",&PetscPreLoading,NULL);CHKERRQ(_3_ierr); \
  PetscPreLoadMax = (int)(PetscPreLoading);\
  PetscPreLoadingUsed = PetscPreLoading ? PETSC_TRUE : PetscPreLoadingUsed;\
  for (PetscPreLoadIt=0; PetscPreLoadIt<=PetscPreLoadMax; PetscPreLoadIt++) {\
    PetscPreLoadingOn = PetscPreLoading;\
    _3_ierr = PetscBarrier(NULL);CHKERRQ(_3_ierr);\
    if (PetscPreLoadIt>0) {\
      _3_ierr = PetscLogStageGetId(name,&_stageNum);CHKERRQ(_3_ierr);\
    } else {\
      _3_ierr = PetscLogStageRegister(name,&_stageNum);CHKERRQ(_3_ierr); \
    }\
    _3_ierr = PetscLogStageSetActive(_stageNum,(PetscBool)(!PetscPreLoadMax || PetscPreLoadIt));\
    _3_ierr = PetscLogStagePush(_stageNum);CHKERRQ(_3_ierr);

#define PetscPreLoadEnd() \
    _3_ierr = PetscLogStagePop();CHKERRQ(_3_ierr);\
    PetscPreLoading = PETSC_FALSE;\
  }\
} while (0)

#define PetscPreLoadStage(name) do {                                         \
    _3_ierr = PetscLogStagePop();CHKERRQ(_3_ierr);                      \
    if (PetscPreLoadIt>0) {                                                  \
      _3_ierr = PetscLogStageGetId(name,&_stageNum);CHKERRQ(_3_ierr);   \
    } else {                                                            \
      _3_ierr = PetscLogStageRegister(name,&_stageNum);CHKERRQ(_3_ierr); \
    }                                                                   \
    _3_ierr = PetscLogStageSetActive(_stageNum,(PetscBool)(!PetscPreLoadMax || PetscPreLoadIt)); \
    _3_ierr = PetscLogStagePush(_stageNum);CHKERRQ(_3_ierr);            \
  } while (0)

/* some vars for logging */
PETSC_EXTERN PetscBool PetscPreLoadingUsed;       /* true if we are or have done preloading */
PETSC_EXTERN PetscBool PetscPreLoadingOn;         /* true if we are currently in a preloading calculation */

#endif
