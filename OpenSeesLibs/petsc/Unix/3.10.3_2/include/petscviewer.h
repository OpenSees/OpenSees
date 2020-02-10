/*
     PetscViewers are objects where other objects can be looked at or stored.
*/

#if !defined(__PETSCVIEWER_H)
#define __PETSCVIEWER_H

#include <petscsys.h>
#include <petscviewertypes.h>
#include <petscdrawtypes.h>

PETSC_EXTERN PetscClassId PETSC_VIEWER_CLASSID;

/*J
    PetscViewerType - String with the name of a PETSc PETScViewer

   Level: beginner

.seealso: PetscViewerSetType(), PetscViewer, PetscViewerRegister(), PetscViewerCreate()
J*/
typedef const char* PetscViewerType;
#define PETSCVIEWERSOCKET       "socket"
#define PETSCVIEWERASCII        "ascii"
#define PETSCVIEWERBINARY       "binary"
#define PETSCVIEWERSTRING       "string"
#define PETSCVIEWERDRAW         "draw"
#define PETSCVIEWERVU           "vu"
#define PETSCVIEWERMATHEMATICA  "mathematica"
#define PETSCVIEWERHDF5         "hdf5"
#define PETSCVIEWERVTK          "vtk"
#define PETSCVIEWERMATLAB       "matlab"
#define PETSCVIEWERSAWS         "saws"
#define PETSCVIEWERGLVIS        "glvis"
#define PETSCVIEWERADIOS        "adios"
#define PETSCVIEWERADIOS2       "adios2"

PETSC_EXTERN PetscFunctionList PetscViewerList;
PETSC_EXTERN PetscErrorCode PetscViewerInitializePackage(void);

PETSC_EXTERN PetscErrorCode PetscViewerRegister(const char[],PetscErrorCode (*)(PetscViewer));

PETSC_EXTERN PetscErrorCode PetscViewerCreate(MPI_Comm,PetscViewer*);
PETSC_EXTERN PetscErrorCode PetscViewerSetFromOptions(PetscViewer);
PETSC_EXTERN PetscErrorCode PetscViewerASCIIOpenWithFILE(MPI_Comm,FILE*,PetscViewer*);

PETSC_EXTERN PetscErrorCode PetscViewerASCIIOpen(MPI_Comm,const char[],PetscViewer*);
PETSC_EXTERN PetscErrorCode PetscViewerASCIISetFILE(PetscViewer,FILE*);
PETSC_EXTERN PetscErrorCode PetscViewerBinaryOpen(MPI_Comm,const char[],PetscFileMode,PetscViewer*);
PETSC_EXTERN PetscErrorCode PetscViewerADIOSOpen(MPI_Comm,const char[],PetscFileMode,PetscViewer*);
PETSC_EXTERN PetscErrorCode PetscViewerADIOS2Open(MPI_Comm,const char[],PetscFileMode,PetscViewer*);
PETSC_EXTERN PetscErrorCode PetscViewerBinaryGetFlowControl(PetscViewer,PetscInt*);
PETSC_EXTERN PetscErrorCode PetscViewerBinarySetFlowControl(PetscViewer,PetscInt);
PETSC_EXTERN PetscErrorCode PetscViewerBinarySetUseMPIIO(PetscViewer,PetscBool);
PETSC_EXTERN PetscErrorCode PetscViewerBinaryGetUseMPIIO(PetscViewer,PetscBool *);
#if defined(PETSC_HAVE_MPIIO)
PETSC_EXTERN PetscErrorCode PetscViewerBinaryGetMPIIODescriptor(PetscViewer,MPI_File*);
PETSC_EXTERN PetscErrorCode PetscViewerBinaryGetMPIIOOffset(PetscViewer,MPI_Offset*);
PETSC_EXTERN PetscErrorCode PetscViewerBinaryAddMPIIOOffset(PetscViewer,MPI_Offset);
#endif

PETSC_EXTERN PetscErrorCode PetscViewerSocketOpen(MPI_Comm,const char[],int,PetscViewer*);
PETSC_EXTERN PetscErrorCode PetscViewerStringOpen(MPI_Comm,char[],size_t,PetscViewer*);
PETSC_EXTERN PetscErrorCode PetscViewerDrawOpen(MPI_Comm,const char[],const char[],int,int,int,int,PetscViewer*);
PETSC_EXTERN PetscErrorCode PetscViewerDrawSetDrawType(PetscViewer,PetscDrawType);
PETSC_EXTERN PetscErrorCode PetscViewerDrawGetDrawType(PetscViewer,PetscDrawType*);
PETSC_EXTERN PetscErrorCode PetscViewerDrawSetTitle(PetscViewer,const char[]);
PETSC_EXTERN PetscErrorCode PetscViewerDrawGetTitle(PetscViewer,const char*[]);
PETSC_EXTERN PetscErrorCode PetscViewerDrawGetDraw(PetscViewer,PetscInt,PetscDraw*);
PETSC_EXTERN PetscErrorCode PetscViewerDrawBaseAdd(PetscViewer,PetscInt);
PETSC_EXTERN PetscErrorCode PetscViewerDrawBaseSet(PetscViewer,PetscInt);
PETSC_EXTERN PetscErrorCode PetscViewerDrawGetDrawLG(PetscViewer,PetscInt,PetscDrawLG*);
PETSC_EXTERN PetscErrorCode PetscViewerDrawGetDrawAxis(PetscViewer,PetscInt,PetscDrawAxis*);

PETSC_EXTERN PetscErrorCode PetscViewerMathematicaOpen(MPI_Comm, int, const char[], const char[], PetscViewer *);
PETSC_EXTERN PetscErrorCode PetscViewerSiloOpen(MPI_Comm, const char[], PetscViewer *);
PETSC_EXTERN PetscErrorCode PetscViewerMatlabOpen(MPI_Comm,const char[],PetscFileMode,PetscViewer*);

/*E
    PetscViewerGLVisType - indicates what type of GLVis viewer to use

    Level: beginner

.seealso: PetscViewerGLVisOpen()
E*/
typedef enum {PETSC_VIEWER_GLVIS_DUMP, PETSC_VIEWER_GLVIS_SOCKET} PetscViewerGLVisType;
PETSC_EXTERN PetscErrorCode PetscViewerGLVisOpen(MPI_Comm,PetscViewerGLVisType,const char*,PetscInt,PetscViewer*);
PETSC_EXTERN PetscErrorCode PetscViewerGLVisSetPrecision(PetscViewer,PetscInt);
PETSC_EXTERN PetscErrorCode PetscViewerGLVisSetSnapId(PetscViewer,PetscInt);
PETSC_EXTERN PetscErrorCode PetscViewerGLVisSetFields(PetscViewer,PetscInt,const char*[],PetscInt[],PetscErrorCode(*)(PetscObject,PetscInt,PetscObject[],void*),PetscObject[],void*,PetscErrorCode(*)(void*));

PETSC_EXTERN PetscErrorCode PetscViewerGetType(PetscViewer,PetscViewerType*);
PETSC_EXTERN PetscErrorCode PetscViewerSetType(PetscViewer,PetscViewerType);
PETSC_EXTERN PetscErrorCode PetscViewerDestroy(PetscViewer*);
PETSC_EXTERN PetscErrorCode PetscViewerGetSubViewer(PetscViewer,MPI_Comm,PetscViewer*);
PETSC_EXTERN PetscErrorCode PetscViewerRestoreSubViewer(PetscViewer,MPI_Comm,PetscViewer*);

PETSC_EXTERN PetscErrorCode PetscViewerSetUp(PetscViewer);
PETSC_EXTERN PetscErrorCode PetscViewerView(PetscViewer,PetscViewer);
PETSC_STATIC_INLINE PetscErrorCode PetscViewerViewFromOptions(PetscViewer A,PetscObject obj,const char name[]) {return PetscObjectViewFromOptions((PetscObject)A,obj,name);}

PETSC_EXTERN PetscErrorCode PetscViewerSetOptionsPrefix(PetscViewer,const char[]);
PETSC_EXTERN PetscErrorCode PetscViewerAppendOptionsPrefix(PetscViewer,const char[]);
PETSC_EXTERN PetscErrorCode PetscViewerGetOptionsPrefix(PetscViewer,const char*[]);

/*E
    PetscViewerFormat - Way a viewer presents the object

   Level: beginner

.seealso: PetscViewer, PetscViewerType, PetscViewerPushFormat(), PetscViewerPopFormat()
E*/
typedef enum {
  PETSC_VIEWER_DEFAULT,
  PETSC_VIEWER_ASCII_MATLAB,
  PETSC_VIEWER_ASCII_MATHEMATICA,
  PETSC_VIEWER_ASCII_IMPL,
  PETSC_VIEWER_ASCII_INFO,
  PETSC_VIEWER_ASCII_INFO_DETAIL,
  PETSC_VIEWER_ASCII_COMMON,
  PETSC_VIEWER_ASCII_SYMMODU,
  PETSC_VIEWER_ASCII_INDEX,
  PETSC_VIEWER_ASCII_DENSE,
  PETSC_VIEWER_ASCII_MATRIXMARKET,
  PETSC_VIEWER_ASCII_VTK,
  PETSC_VIEWER_ASCII_VTK_CELL,
  PETSC_VIEWER_ASCII_VTK_COORDS,
  PETSC_VIEWER_ASCII_PCICE,
  PETSC_VIEWER_ASCII_PYTHON,
  PETSC_VIEWER_ASCII_FACTOR_INFO,
  PETSC_VIEWER_ASCII_LATEX,
  PETSC_VIEWER_ASCII_XML,
  PETSC_VIEWER_ASCII_GLVIS,
  PETSC_VIEWER_DRAW_BASIC,
  PETSC_VIEWER_DRAW_LG,
  PETSC_VIEWER_DRAW_CONTOUR,
  PETSC_VIEWER_DRAW_PORTS,
  PETSC_VIEWER_VTK_VTS,
  PETSC_VIEWER_VTK_VTR,
  PETSC_VIEWER_VTK_VTU,
  PETSC_VIEWER_BINARY_MATLAB,
  PETSC_VIEWER_NATIVE,
  PETSC_VIEWER_HDF5_PETSC,
  PETSC_VIEWER_HDF5_VIZ,
  PETSC_VIEWER_HDF5_XDMF,
  PETSC_VIEWER_NOFORMAT,
  PETSC_VIEWER_LOAD_BALANCE
  } PetscViewerFormat;
PETSC_EXTERN const char *const PetscViewerFormats[];

PETSC_EXTERN PETSC_DEPRECATED("Use PetscViewerPushFormat()/PetscViewerPopFormat()") PetscErrorCode PetscViewerSetFormat(PetscViewer,PetscViewerFormat);
PETSC_EXTERN PetscErrorCode PetscViewerPushFormat(PetscViewer,PetscViewerFormat);
PETSC_EXTERN PetscErrorCode PetscViewerPopFormat(PetscViewer);
PETSC_EXTERN PetscErrorCode PetscViewerGetFormat(PetscViewer,PetscViewerFormat*);
PETSC_EXTERN PetscErrorCode PetscViewerFlush(PetscViewer);

PETSC_EXTERN PetscErrorCode PetscOptionsPushGetViewerOff(PetscBool);
PETSC_EXTERN PetscErrorCode PetscOptionsPopGetViewerOff(void);
PETSC_EXTERN PetscErrorCode PetscOptionsGetViewerOff(PetscBool*);
PETSC_EXTERN PetscErrorCode PetscOptionsGetViewer(MPI_Comm,const char[],const char[],PetscViewer*,PetscViewerFormat*,PetscBool*);
#define PetscOptionsViewer(a,b,c,d,e,f) PetscOptionsViewer_Private(PetscOptionsObject,a,b,c,d,e,f);
PETSC_EXTERN PetscErrorCode PetscOptionsViewer_Private(PetscOptionItems*,const char[],const char[],const char[],PetscViewer*,PetscViewerFormat *,PetscBool *);

typedef struct {PetscViewer viewer;PetscViewerFormat format;} PetscViewerAndFormat;
PETSC_EXTERN PetscErrorCode  PetscViewerAndFormatCreate(PetscViewer,PetscViewerFormat,PetscViewerAndFormat**);
PETSC_EXTERN PetscErrorCode  PetscViewerAndFormatDestroy(PetscViewerAndFormat**);

/*
   Operations explicit to a particular class of viewers
*/

PETSC_EXTERN PetscErrorCode PetscViewerASCIIGetPointer(PetscViewer,FILE**);
PETSC_EXTERN PetscErrorCode PetscViewerFileGetMode(PetscViewer,PetscFileMode*);
PETSC_EXTERN PetscErrorCode PetscViewerFileSetMode(PetscViewer,PetscFileMode);
PETSC_EXTERN PetscErrorCode PetscViewerRead(PetscViewer,void*,PetscInt,PetscInt*,PetscDataType);
PETSC_EXTERN PetscErrorCode PetscViewerASCIIPrintf(PetscViewer,const char[],...);
PETSC_EXTERN PetscErrorCode PetscViewerASCIISynchronizedPrintf(PetscViewer,const char[],...);
PETSC_EXTERN PetscErrorCode PetscViewerASCIIPushSynchronized(PetscViewer);
PETSC_EXTERN PetscErrorCode PetscViewerASCIIPopSynchronized(PetscViewer);
PETSC_EXTERN PetscErrorCode PetscViewerASCIIPushTab(PetscViewer);
PETSC_EXTERN PetscErrorCode PetscViewerASCIIPopTab(PetscViewer);
PETSC_EXTERN PetscErrorCode PetscViewerASCIIUseTabs(PetscViewer,PetscBool );
PETSC_EXTERN PetscErrorCode PetscViewerASCIISetTab(PetscViewer,PetscInt);
PETSC_EXTERN PetscErrorCode PetscViewerASCIIGetTab(PetscViewer,PetscInt*);
PETSC_EXTERN PetscErrorCode PetscViewerASCIIAddTab(PetscViewer,PetscInt);
PETSC_EXTERN PetscErrorCode PetscViewerASCIISubtractTab(PetscViewer,PetscInt);
PETSC_EXTERN PetscErrorCode PetscViewerASCIIRead(PetscViewer,void *,PetscInt,PetscInt*,PetscDataType);
PETSC_EXTERN PetscErrorCode PetscViewerBinaryGetDescriptor(PetscViewer,int*);
PETSC_EXTERN PetscErrorCode PetscViewerBinaryGetInfoPointer(PetscViewer,FILE **);
PETSC_EXTERN PetscErrorCode PetscViewerBinaryRead(PetscViewer,void*,PetscInt,PetscInt*,PetscDataType);
PETSC_EXTERN PetscErrorCode PetscViewerBinaryWrite(PetscViewer,void*,PetscInt,PetscDataType,PetscBool );
PETSC_EXTERN PetscErrorCode PetscViewerStringSPrintf(PetscViewer,const char[],...);
PETSC_EXTERN PetscErrorCode PetscViewerStringSetString(PetscViewer,char[],PetscInt);
PETSC_EXTERN PetscErrorCode PetscViewerDrawClear(PetscViewer);
PETSC_EXTERN PetscErrorCode PetscViewerDrawSetHold(PetscViewer,PetscBool);
PETSC_EXTERN PetscErrorCode PetscViewerDrawGetHold(PetscViewer,PetscBool*);
PETSC_EXTERN PetscErrorCode PetscViewerDrawSetPause(PetscViewer,PetscReal);
PETSC_EXTERN PetscErrorCode PetscViewerDrawGetPause(PetscViewer,PetscReal*);
PETSC_EXTERN PetscErrorCode PetscViewerDrawSetInfo(PetscViewer,const char[],const char[],int,int,int,int);
PETSC_EXTERN PetscErrorCode PetscViewerDrawResize(PetscViewer,int,int);
PETSC_EXTERN PetscErrorCode PetscViewerDrawSetBounds(PetscViewer,PetscInt,const PetscReal*);
PETSC_EXTERN PetscErrorCode PetscViewerDrawGetBounds(PetscViewer,PetscInt*,const PetscReal**);
PETSC_EXTERN PetscErrorCode PetscViewerSocketSetConnection(PetscViewer,const char[],int);
PETSC_EXTERN PetscErrorCode PetscViewerBinarySkipInfo(PetscViewer);
PETSC_EXTERN PetscErrorCode PetscViewerBinarySetSkipInfo(PetscViewer,PetscBool);
PETSC_EXTERN PetscErrorCode PetscViewerBinaryGetSkipInfo(PetscViewer,PetscBool*);
PETSC_EXTERN PetscErrorCode PetscViewerBinarySetSkipOptions(PetscViewer,PetscBool );
PETSC_EXTERN PetscErrorCode PetscViewerBinaryGetSkipOptions(PetscViewer,PetscBool *);
PETSC_EXTERN PetscErrorCode PetscViewerBinarySetSkipHeader(PetscViewer,PetscBool);
PETSC_EXTERN PetscErrorCode PetscViewerBinaryGetSkipHeader(PetscViewer,PetscBool*);
PETSC_EXTERN PetscErrorCode PetscViewerBinaryReadStringArray(PetscViewer,char***);
PETSC_EXTERN PetscErrorCode PetscViewerBinaryWriteStringArray(PetscViewer,const char *const*);

PETSC_EXTERN PetscErrorCode PetscViewerFileSetName(PetscViewer,const char[]);
PETSC_EXTERN PetscErrorCode PetscViewerFileGetName(PetscViewer,const char**);

PETSC_EXTERN PetscErrorCode PetscViewerVUGetPointer(PetscViewer, FILE**);
PETSC_EXTERN PetscErrorCode PetscViewerVUSetVecSeen(PetscViewer, PetscBool );
PETSC_EXTERN PetscErrorCode PetscViewerVUGetVecSeen(PetscViewer, PetscBool  *);
PETSC_EXTERN PetscErrorCode PetscViewerVUPrintDeferred(PetscViewer, const char [], ...);
PETSC_EXTERN PetscErrorCode PetscViewerVUFlushDeferred(PetscViewer);

PETSC_EXTERN PetscErrorCode PetscViewerMathematicaInitializePackage(void);
PETSC_EXTERN PetscErrorCode PetscViewerMathematicaFinalizePackage(void);
PETSC_EXTERN PetscErrorCode PetscViewerMathematicaGetName(PetscViewer, const char **);
PETSC_EXTERN PetscErrorCode PetscViewerMathematicaSetName(PetscViewer, const char []);
PETSC_EXTERN PetscErrorCode PetscViewerMathematicaClearName(PetscViewer);
PETSC_EXTERN PetscErrorCode PetscViewerMathematicaSkipPackets(PetscViewer, int);

PETSC_EXTERN PetscErrorCode PetscViewerSiloGetName(PetscViewer, char **);
PETSC_EXTERN PetscErrorCode PetscViewerSiloSetName(PetscViewer, const char []);
PETSC_EXTERN PetscErrorCode PetscViewerSiloClearName(PetscViewer);
PETSC_EXTERN PetscErrorCode PetscViewerSiloGetMeshName(PetscViewer, char **);
PETSC_EXTERN PetscErrorCode PetscViewerSiloSetMeshName(PetscViewer, const char []);
PETSC_EXTERN PetscErrorCode PetscViewerSiloClearMeshName(PetscViewer);

typedef enum {PETSC_VTK_POINT_FIELD, PETSC_VTK_POINT_VECTOR_FIELD, PETSC_VTK_CELL_FIELD, PETSC_VTK_CELL_VECTOR_FIELD} PetscViewerVTKFieldType;
PETSC_EXTERN PetscErrorCode PetscViewerVTKAddField(PetscViewer,PetscObject,PetscErrorCode (*PetscViewerVTKWriteFunction)(PetscObject,PetscViewer),PetscViewerVTKFieldType,PetscBool,PetscObject);
PETSC_EXTERN PetscErrorCode PetscViewerVTKGetDM(PetscViewer,PetscObject*);
PETSC_EXTERN PetscErrorCode PetscViewerVTKOpen(MPI_Comm,const char[],PetscFileMode,PetscViewer*);

/*
     These are all the default viewers that do not have to be explicitly opened
*/
PETSC_EXTERN PetscViewer    PETSC_VIEWER_STDOUT_(MPI_Comm);
PETSC_EXTERN PetscErrorCode PetscViewerASCIIGetStdout(MPI_Comm,PetscViewer*);
PETSC_EXTERN PetscViewer    PETSC_VIEWER_STDERR_(MPI_Comm);
PETSC_EXTERN PetscErrorCode PetscViewerASCIIGetStderr(MPI_Comm,PetscViewer*);
PETSC_EXTERN PetscViewer    PETSC_VIEWER_DRAW_(MPI_Comm);
PETSC_EXTERN PetscViewer    PETSC_VIEWER_SOCKET_(MPI_Comm);
PETSC_EXTERN PetscViewer    PETSC_VIEWER_BINARY_(MPI_Comm);
PETSC_EXTERN PetscViewer    PETSC_VIEWER_MATLAB_(MPI_Comm);
PETSC_EXTERN PetscViewer    PETSC_VIEWER_HDF5_(MPI_Comm);
PETSC_EXTERN PetscViewer    PETSC_VIEWER_GLVIS_(MPI_Comm);
PETSC_EXTERN PetscViewer    PETSC_VIEWER_MATHEMATICA_WORLD_PRIVATE;

#define PETSC_VIEWER_STDERR_SELF  PETSC_VIEWER_STDERR_(PETSC_COMM_SELF)
#define PETSC_VIEWER_STDERR_WORLD PETSC_VIEWER_STDERR_(PETSC_COMM_WORLD)

/*MC
  PETSC_VIEWER_STDOUT_WORLD  - same as PETSC_VIEWER_STDOUT_(PETSC_COMM_WORLD)

  Level: beginner
M*/
#define PETSC_VIEWER_STDOUT_WORLD PETSC_VIEWER_STDOUT_(PETSC_COMM_WORLD)

/*MC
  PETSC_VIEWER_STDOUT_SELF  - same as PETSC_VIEWER_STDOUT_(PETSC_COMM_SELF)

  Level: beginner
M*/
#define PETSC_VIEWER_STDOUT_SELF  PETSC_VIEWER_STDOUT_(PETSC_COMM_SELF)

/*MC
  PETSC_VIEWER_DRAW_WORLD  - same as PETSC_VIEWER_DRAW_(PETSC_COMM_WORLD)

  Level: intermediate
M*/
#define PETSC_VIEWER_DRAW_WORLD   PETSC_VIEWER_DRAW_(PETSC_COMM_WORLD)

/*MC
  PETSC_VIEWER_DRAW_SELF  - same as PETSC_VIEWER_DRAW_(PETSC_COMM_SELF)

  Level: intermediate
M*/
#define PETSC_VIEWER_DRAW_SELF    PETSC_VIEWER_DRAW_(PETSC_COMM_SELF)

/*MC
  PETSC_VIEWER_SOCKET_WORLD  - same as PETSC_VIEWER_SOCKET_(PETSC_COMM_WORLD)

  Level: intermediate
M*/
#define PETSC_VIEWER_SOCKET_WORLD PETSC_VIEWER_SOCKET_(PETSC_COMM_WORLD)

/*MC
  PETSC_VIEWER_SOCKET_SELF  - same as PETSC_VIEWER_SOCKET_(PETSC_COMM_SELF)

  Level: intermediate
M*/
#define PETSC_VIEWER_SOCKET_SELF  PETSC_VIEWER_SOCKET_(PETSC_COMM_SELF)

/*MC
  PETSC_VIEWER_BINARY_WORLD  - same as PETSC_VIEWER_BINARY_(PETSC_COMM_WORLD)

  Level: intermediate
M*/
#define PETSC_VIEWER_BINARY_WORLD PETSC_VIEWER_BINARY_(PETSC_COMM_WORLD)

/*MC
  PETSC_VIEWER_BINARY_SELF  - same as PETSC_VIEWER_BINARY_(PETSC_COMM_SELF)

  Level: intermediate
M*/
#define PETSC_VIEWER_BINARY_SELF  PETSC_VIEWER_BINARY_(PETSC_COMM_SELF)

/*MC
  PETSC_VIEWER_MATLAB_WORLD  - same as PETSC_VIEWER_MATLAB_(PETSC_COMM_WORLD)

  Level: intermediate
M*/
#define PETSC_VIEWER_MATLAB_WORLD PETSC_VIEWER_MATLAB_(PETSC_COMM_WORLD)

/*MC
  PETSC_VIEWER_MATLAB_SELF  - same as PETSC_VIEWER_MATLAB_(PETSC_COMM_SELF)

  Level: intermediate
M*/
#define PETSC_VIEWER_MATLAB_SELF  PETSC_VIEWER_MATLAB_(PETSC_COMM_SELF)

#define PETSC_VIEWER_MATHEMATICA_WORLD (PetscViewerInitializeMathematicaWorld_Private(),PETSC_VIEWER_MATHEMATICA_WORLD_PRIVATE)

PETSC_EXTERN PetscErrorCode PetscViewerFlowControlStart(PetscViewer,PetscInt*,PetscInt*);
PETSC_EXTERN PetscErrorCode PetscViewerFlowControlStepMaster(PetscViewer,PetscInt,PetscInt*,PetscInt);
PETSC_EXTERN PetscErrorCode PetscViewerFlowControlEndMaster(PetscViewer,PetscInt*);
PETSC_EXTERN PetscErrorCode PetscViewerFlowControlStepWorker(PetscViewer,PetscMPIInt,PetscInt*);
PETSC_EXTERN PetscErrorCode PetscViewerFlowControlEndWorker(PetscViewer,PetscInt*);

/*
   PetscViewer writes to MATLAB .mat file
*/
PETSC_EXTERN PetscErrorCode PetscViewerMatlabPutArray(PetscViewer,int,int,const PetscScalar*,const char*);
PETSC_EXTERN PetscErrorCode PetscViewerMatlabGetArray(PetscViewer,int,int,PetscScalar*,const char*);
PETSC_EXTERN PetscErrorCode PetscViewerMatlabPutVariable(PetscViewer,const char*,void*);

#if defined(PETSC_HAVE_SAWS)
PETSC_EXTERN PetscErrorCode PetscObjectViewSAWs(PetscObject,PetscViewer);
#endif

/*S
     PetscViewers - Abstract collection of PetscViewers. It is just an expandable array of viewers.

   Level: intermediate

  Concepts: viewing

.seealso:  PetscViewerCreate(), PetscViewerSetType(), PetscViewerType, PetscViewer, PetscViewersCreate(),
           PetscViewersGetViewer()
S*/
typedef struct _n_PetscViewers* PetscViewers;
PETSC_EXTERN PetscErrorCode PetscViewersCreate(MPI_Comm,PetscViewers*);
PETSC_EXTERN PetscErrorCode PetscViewersDestroy(PetscViewers*);
PETSC_EXTERN PetscErrorCode PetscViewersGetViewer(PetscViewers,PetscInt,PetscViewer*);

#endif
