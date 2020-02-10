#include <petsc/private/petscimpl.h>
#include <petsctime.h>

/* A simple stack */
struct _n_PetscIntStack {
  int  top;   /* The top of the stack */
  int  max;   /* The maximum stack size */
  int *stack; /* The storage */
};

/* The structure for action logging */
#define CREATE      0
#define DESTROY     1
#define ACTIONBEGIN 2
#define ACTIONEND   3
typedef struct _Action {
  int            action;        /* The type of execution */
  PetscLogEvent  event;         /* The event number */
  PetscClassId   classid;        /* The event class id */
  PetscLogDouble time;          /* The time of occurence */
  PetscLogDouble flops;         /* The cumlative flops */
  PetscLogDouble mem;           /* The current memory usage */
  PetscLogDouble maxmem;        /* The maximum memory usage */
  int            id1, id2, id3; /* The ids of associated objects */
} Action;

/* The structure for object logging */
typedef struct _Object {
  PetscObject    obj;      /* The associated PetscObject */
  int            parent;   /* The parent id */
  PetscLogDouble mem;      /* The memory associated with the object */
  char           name[64]; /* The object name */
  char           info[64]; /* The information string */
} Object;

/* Action and object logging variables */
PETSC_EXTERN Action    *petsc_actions;
PETSC_EXTERN Object    *petsc_objects;
PETSC_EXTERN PetscBool  petsc_logActions;
PETSC_EXTERN PetscBool  petsc_logObjects;
PETSC_EXTERN int        petsc_numActions;
PETSC_EXTERN int        petsc_maxActions;
PETSC_EXTERN int        petsc_numObjects;
PETSC_EXTERN int        petsc_maxObjects;
PETSC_EXTERN int        petsc_numObjectsDestroyed;

PETSC_EXTERN FILE          *petsc_tracefile;
PETSC_EXTERN int            petsc_tracelevel;
PETSC_EXTERN const char    *petsc_traceblanks;
PETSC_EXTERN char           petsc_tracespace[128];
PETSC_EXTERN PetscLogDouble petsc_tracetime;

#ifdef PETSC_USE_LOG

PETSC_EXTERN PetscErrorCode PetscIntStackCreate(PetscIntStack *);
PETSC_EXTERN PetscErrorCode PetscIntStackDestroy(PetscIntStack);
PETSC_EXTERN PetscErrorCode PetscIntStackPush(PetscIntStack, int);
PETSC_EXTERN PetscErrorCode PetscIntStackPop(PetscIntStack, int *);
PETSC_EXTERN PetscErrorCode PetscIntStackTop(PetscIntStack, int *);
PETSC_EXTERN PetscErrorCode PetscIntStackEmpty(PetscIntStack, PetscBool  *);

/* Creation and destruction functions */
PETSC_EXTERN PetscErrorCode PetscEventRegLogCreate(PetscEventRegLog *);
PETSC_EXTERN PetscErrorCode PetscEventRegLogDestroy(PetscEventRegLog);
PETSC_EXTERN PetscErrorCode PetscEventPerfLogCreate(PetscEventPerfLog *);
PETSC_EXTERN PetscErrorCode PetscEventPerfLogDestroy(PetscEventPerfLog);
/* General functions */
PETSC_EXTERN PetscErrorCode PetscEventPerfLogEnsureSize(PetscEventPerfLog, int);
PETSC_EXTERN PetscErrorCode PetscEventPerfInfoClear(PetscEventPerfInfo *);
PETSC_EXTERN PetscErrorCode PetscEventPerfInfoCopy(PetscEventPerfInfo *, PetscEventPerfInfo *);
/* Registration functions */
PETSC_EXTERN PetscErrorCode PetscEventRegLogRegister(PetscEventRegLog, const char [], PetscClassId, PetscLogEvent *);
/* Query functions */
PETSC_EXTERN PetscErrorCode PetscEventPerfLogSetVisible(PetscEventPerfLog, PetscLogEvent, PetscBool );
PETSC_EXTERN PetscErrorCode PetscEventPerfLogGetVisible(PetscEventPerfLog, PetscLogEvent, PetscBool  *);
/* Activaton functions */
PETSC_EXTERN PetscErrorCode PetscEventPerfLogActivate(PetscEventPerfLog, PetscLogEvent);
PETSC_EXTERN PetscErrorCode PetscEventPerfLogDeactivate(PetscEventPerfLog, PetscLogEvent);
PETSC_EXTERN PetscErrorCode PetscEventPerfLogActivateClass(PetscEventPerfLog, PetscEventRegLog, PetscClassId);
PETSC_EXTERN PetscErrorCode PetscEventPerfLogDeactivateClass(PetscEventPerfLog, PetscEventRegLog, PetscClassId);

/* Logging functions */
PETSC_EXTERN PetscErrorCode PetscLogEventSynchronize(PetscLogEvent, MPI_Comm);
PETSC_EXTERN PetscErrorCode PetscLogEventBeginDefault(PetscLogEvent, int, PetscObject, PetscObject, PetscObject, PetscObject);
PETSC_EXTERN PetscErrorCode PetscLogEventEndDefault(PetscLogEvent, int, PetscObject, PetscObject, PetscObject, PetscObject);
PETSC_EXTERN PetscErrorCode PetscLogEventBeginComplete(PetscLogEvent, int, PetscObject, PetscObject, PetscObject, PetscObject);
PETSC_EXTERN PetscErrorCode PetscLogEventEndComplete(PetscLogEvent, int, PetscObject, PetscObject, PetscObject, PetscObject);
PETSC_EXTERN PetscErrorCode PetscLogEventBeginTrace(PetscLogEvent, int, PetscObject, PetscObject, PetscObject, PetscObject);
PETSC_EXTERN PetscErrorCode PetscLogEventEndTrace(PetscLogEvent, int, PetscObject, PetscObject, PetscObject, PetscObject);

/* Creation and destruction functions */
PETSC_EXTERN PetscErrorCode PetscClassRegLogCreate(PetscClassRegLog *);
PETSC_EXTERN PetscErrorCode PetscClassRegLogDestroy(PetscClassRegLog);
PETSC_EXTERN PetscErrorCode PetscClassPerfLogCreate(PetscClassPerfLog *);
PETSC_EXTERN PetscErrorCode PetscClassPerfLogDestroy(PetscClassPerfLog);
PETSC_EXTERN PetscErrorCode PetscClassRegInfoDestroy(PetscClassRegInfo *);
/* General functions */
PETSC_EXTERN PetscErrorCode PetscClassPerfLogEnsureSize(PetscClassPerfLog, int);
PETSC_EXTERN PetscErrorCode PetscClassPerfInfoClear(PetscClassPerfInfo *);
/* Registration functions */
PETSC_EXTERN PetscErrorCode PetscClassRegLogRegister(PetscClassRegLog, const char [], PetscClassId);
/* Query functions */
PETSC_EXTERN PetscErrorCode PetscClassRegLogGetClass(PetscClassRegLog, PetscClassId, int *);
/* Logging functions */
PETSC_EXTERN PetscErrorCode PetscLogObjCreateDefault(PetscObject);
PETSC_EXTERN PetscErrorCode PetscLogObjDestroyDefault(PetscObject);

/* Creation and destruction functions */
PETSC_EXTERN PetscErrorCode PetscStageLogCreate(PetscStageLog *);
PETSC_EXTERN PetscErrorCode PetscStageLogDestroy(PetscStageLog);
/* Registration functions */
PETSC_EXTERN PetscErrorCode PetscStageLogRegister(PetscStageLog, const char [], int *);
/* Runtime functions */
PETSC_EXTERN PetscErrorCode PetscStageLogPush(PetscStageLog, int);
PETSC_EXTERN PetscErrorCode PetscStageLogPop(PetscStageLog);
PETSC_EXTERN PetscErrorCode PetscStageLogSetActive(PetscStageLog, int, PetscBool );
PETSC_EXTERN PetscErrorCode PetscStageLogGetActive(PetscStageLog, int, PetscBool  *);
PETSC_EXTERN PetscErrorCode PetscStageLogSetVisible(PetscStageLog, int, PetscBool );
PETSC_EXTERN PetscErrorCode PetscStageLogGetVisible(PetscStageLog, int, PetscBool  *);
PETSC_EXTERN PetscErrorCode PetscStageLogGetStage(PetscStageLog, const char [], PetscLogStage *);
PETSC_EXTERN PetscErrorCode PetscStageLogGetClassRegLog(PetscStageLog, PetscClassRegLog *);
PETSC_EXTERN PetscErrorCode PetscStageLogGetEventRegLog(PetscStageLog, PetscEventRegLog *);
PETSC_EXTERN PetscErrorCode PetscStageLogGetClassPerfLog(PetscStageLog, int, PetscClassPerfLog *);

PETSC_EXTERN PetscErrorCode PetscEventRegLogGetEvent(PetscEventRegLog, const char [], PetscLogEvent *);

PETSC_INTERN PetscErrorCode PetscLogView_Nested(PetscViewer);
PETSC_INTERN PetscErrorCode PetscLogNestedEnd(void);

#endif /* PETSC_USE_LOG */
