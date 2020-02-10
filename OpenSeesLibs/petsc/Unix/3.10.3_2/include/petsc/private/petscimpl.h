
/*
    Defines the basic header of all PETSc objects.
*/

#if !defined(_PETSCHEAD_H)
#define _PETSCHEAD_H
#include <petscsys.h>

/* These are used internally by PETSc ASCII IO routines*/
#include <stdarg.h>
PETSC_EXTERN PetscErrorCode PetscVFPrintfDefault(FILE*,const char[],va_list);

#if defined(PETSC_HAVE_CLOSURES)
PETSC_EXTERN PetscErrorCode PetscVFPrintfSetClosure(int (^)(const char*));
#endif


#if defined(PETSC_HAVE_CUDA)
#include <cuda.h>
#include <cublas_v2.h>
#endif

/*
   All major PETSc data structures have a common core; this is defined
   below by PETSCHEADER.

   PetscHeaderCreate() should be used whenever creating a PETSc structure.
*/

/*
   PetscOps: structure of core operations that all PETSc objects support.

      getcomm()         - Gets the object's communicator.
      view()            - Is the routine for viewing the entire PETSc object; for
                          example, MatView() is the general matrix viewing routine.
                          This is used by PetscObjectView((PetscObject)obj) to allow
                          viewing any PETSc object.
      destroy()         - Is the routine for destroying the entire PETSc object;
                          for example,MatDestroy() is the general matrix
                          destruction routine.
                          This is used by PetscObjectDestroy((PetscObject*)&obj) to allow
                          destroying any PETSc object.
      compose()         - Associates a PETSc object with another PETSc object with a name
      query()           - Returns a different PETSc object that has been associated
                          with the first object using a name.
      composefunction() - Attaches an a function to a PETSc object with a name.
      queryfunction()   - Requests a registered function that has been attached to a PETSc object.
*/

typedef struct {
   PetscErrorCode (*getcomm)(PetscObject,MPI_Comm *);
   PetscErrorCode (*view)(PetscObject,PetscViewer);
   PetscErrorCode (*destroy)(PetscObject*);
   PetscErrorCode (*compose)(PetscObject,const char[],PetscObject);
   PetscErrorCode (*query)(PetscObject,const char[],PetscObject *);
   PetscErrorCode (*composefunction)(PetscObject,const char[],void (*)(void));
   PetscErrorCode (*queryfunction)(PetscObject,const char[],void (**)(void));
} PetscOps;

typedef enum {PETSC_FORTRAN_CALLBACK_CLASS,PETSC_FORTRAN_CALLBACK_SUBTYPE,PETSC_FORTRAN_CALLBACK_MAXTYPE} PetscFortranCallbackType;
typedef int PetscFortranCallbackId;
#define PETSC_SMALLEST_FORTRAN_CALLBACK ((PetscFortranCallbackId)1000)
PETSC_EXTERN PetscErrorCode PetscFortranCallbackRegister(PetscClassId,const char*,PetscFortranCallbackId*);
PETSC_EXTERN PetscErrorCode PetscFortranCallbackGetSizes(PetscClassId,PetscInt*,PetscInt*);

typedef struct {
  void (*func)(void);
  void *ctx;
} PetscFortranCallback;

/*
   All PETSc objects begin with the fields defined in PETSCHEADER.
   The PetscObject is a way of examining these fields regardless of
   the specific object. In C++ this could be a base abstract class
   from which all objects are derived.
*/
#define PETSC_MAX_OPTIONS_HANDLER 5
typedef struct _p_PetscObject {
  PetscClassId         classid;
  PetscOps             bops[1];
  MPI_Comm             comm;
  PetscInt             type;
  PetscLogDouble       flops,time,mem,memchildren;
  PetscObjectId        id;
  PetscInt             refct;
  PetscMPIInt          tag;
  PetscFunctionList    qlist;
  PetscObjectList      olist;
  char                 *class_name;    /*  for example, "Vec" */
  char                 *description;
  char                 *mansec;
  char                 *type_name;     /*  this is the subclass, for example VECSEQ which equals "seq" */
  PetscObject          parent;
  PetscObjectId        parentid;
  char*                name;
  char                 *prefix;
  PetscInt             tablevel;
  void                 *cpp;
  PetscObjectState     state;
  PetscInt             int_idmax,        intstar_idmax;
  PetscObjectState     *intcomposedstate,*intstarcomposedstate;
  PetscInt             *intcomposeddata, **intstarcomposeddata;
  PetscInt             real_idmax,        realstar_idmax;
  PetscObjectState     *realcomposedstate,*realstarcomposedstate;
  PetscReal            *realcomposeddata, **realstarcomposeddata;
  PetscInt             scalar_idmax,        scalarstar_idmax;
  PetscObjectState     *scalarcomposedstate,*scalarstarcomposedstate;
  PetscScalar          *scalarcomposeddata, **scalarstarcomposeddata;
  void                 (**fortran_func_pointers)(void);                  /* used by Fortran interface functions to stash user provided Fortran functions */
  PetscInt             num_fortran_func_pointers;                        /* number of Fortran function pointers allocated */
  PetscFortranCallback *fortrancallback[PETSC_FORTRAN_CALLBACK_MAXTYPE];
  PetscInt             num_fortrancallback[PETSC_FORTRAN_CALLBACK_MAXTYPE];
  void                 *python_context;
  PetscErrorCode       (*python_destroy)(void*);

  PetscInt             noptionhandler;
  PetscErrorCode       (*optionhandler[PETSC_MAX_OPTIONS_HANDLER])(PetscOptionItems*,PetscObject,void*);
  PetscErrorCode       (*optiondestroy[PETSC_MAX_OPTIONS_HANDLER])(PetscObject,void*);
  void                 *optionctx[PETSC_MAX_OPTIONS_HANDLER];
  PetscBool            optionsprinted;
#if defined(PETSC_HAVE_SAWS)
  PetscBool            amsmem;          /* if PETSC_TRUE then this object is registered with SAWs and visible to clients */
  PetscBool            amspublishblock; /* if PETSC_TRUE and publishing objects then will block at PetscObjectSAWsBlock() */
#endif
  PetscOptions         options;         /* options database used, NULL means default */
  PetscBool            donotPetscObjectPrintClassNamePrefixType;
} _p_PetscObject;

#define PETSCHEADER(ObjectOps) \
  _p_PetscObject hdr;          \
  ObjectOps      ops[1]

#define  PETSCFREEDHEADER -1

PETSC_EXTERN_TYPEDEF typedef PetscErrorCode (*PetscObjectDestroyFunction)(PetscObject*); /* force cast in next macro to NEVER use extern "C" style */
PETSC_EXTERN_TYPEDEF typedef PetscErrorCode (*PetscObjectViewFunction)(PetscObject,PetscViewer);

/*@C
    PetscHeaderCreate - Creates a PETSc object of a particular class

    Input Parameters:
+   classid - the classid associated with this object (for example VEC_CLASSID)
.   class_name - string name of class; should be static (for example "Vec")
.   descr - string containing short description; should be static (for example "Vector")
.   mansec - string indicating section in manual pages; should be static (for example "Vec")
.   comm - the MPI Communicator
.   destroy - the destroy routine for this object (for example VecDestroy())
-   view - the view routine for this object (for example VecView())

    Output Parameter:
.   h - the newly created object

    Level: developer

.seealso: PetscHeaderDestroy(), PetscClassIdRegister()

@*/
#define PetscHeaderCreate(h,classid,class_name,descr,mansec,comm,destroy,view) \
  (PetscNew(&(h)) || \
   PetscHeaderCreate_Private((PetscObject)h,classid,class_name,descr,mansec,comm,(PetscObjectDestroyFunction)destroy,(PetscObjectViewFunction)view) || \
   PetscLogObjectCreate(h) || \
   PetscLogObjectMemory((PetscObject)h,sizeof(*(h))))

PETSC_EXTERN PetscErrorCode PetscComposedQuantitiesDestroy(PetscObject obj);
PETSC_EXTERN PetscErrorCode PetscHeaderCreate_Private(PetscObject,PetscClassId,const char[],const char[],const char[],MPI_Comm,PetscObjectDestroyFunction,PetscObjectViewFunction);

/*@C
    PetscHeaderDestroy - Final step in destroying a PetscObject

    Input Parameters:
.   h - the header created with PetscHeaderCreate()

    Level: developer

.seealso: PetscHeaderCreate()
@*/
#define PetscHeaderDestroy(h) (PetscHeaderDestroy_Private((PetscObject)(*h)) || PetscFree(*h))

PETSC_EXTERN PetscErrorCode PetscHeaderDestroy_Private(PetscObject);
PETSC_EXTERN PetscErrorCode PetscObjectCopyFortranFunctionPointers(PetscObject,PetscObject);
PETSC_EXTERN PetscErrorCode PetscObjectSetFortranCallback(PetscObject,PetscFortranCallbackType,PetscFortranCallbackId*,void(*)(void),void *ctx);
PETSC_EXTERN PetscErrorCode PetscObjectGetFortranCallback(PetscObject,PetscFortranCallbackType,PetscFortranCallbackId,void(**)(void),void **ctx);

PETSC_INTERN PetscErrorCode PetscCitationsInitialize(void);
PETSC_INTERN PetscErrorCode PetscFreeMPIResources(void);



PETSC_EXTERN PetscBool PetscCheckPointer(const void*,PetscDataType);

/*
    Macros to test if a PETSc object is valid and if pointers are valid
*/
#if !defined(PETSC_USE_DEBUG)

#define PetscValidHeaderSpecific(h,ck,arg) do {} while (0)
#define PetscValidHeaderSpecificType(h,ck,arg,t) do {} while (0)
#define PetscValidHeader(h,arg) do {} while (0)
#define PetscValidPointer(h,arg) do {} while (0)
#define PetscValidCharPointer(h,arg) do {} while (0)
#define PetscValidIntPointer(h,arg) do {} while (0)
#define PetscValidScalarPointer(h,arg) do {} while (0)
#define PetscValidRealPointer(h,arg) do {} while (0)
#define PetscValidFunction(h,arg) do {} while (0)

#else

/*  This check is for subtype methods such as DMDAGetCorners() that do not use the PetscTryMethod() or PetscUseMethod() paradigm */
#define PetscValidHeaderSpecificType(h,ck,arg,t) \
  do {   \
    PetscErrorCode __ierr; \
    PetscBool      same; \
    PetscValidHeaderSpecific(h,ck,arg); \
    __ierr = PetscObjectTypeCompare((PetscObject)h,t,&same);CHKERRQ(__ierr);      \
    if (!same) SETERRQ3(PETSC_COMM_SELF,PETSC_ERR_ARG_WRONG,"Wrong subtype object:Parameter # %d must have implementation %s it is %s",arg,t,((PetscObject)h)->type_name); \
  } while (0)

#define PetscValidHeaderSpecific(h,ck,arg)                              \
  do {                                                                  \
    if (!h) SETERRQ1(PETSC_COMM_SELF,PETSC_ERR_ARG_NULL,"Null Object: Parameter # %d",arg); \
    if (!PetscCheckPointer(h,PETSC_OBJECT)) SETERRQ1(PETSC_COMM_SELF,PETSC_ERR_ARG_CORRUPT,"Invalid Pointer to Object: Parameter # %d",arg); \
    if (((PetscObject)(h))->classid != ck) {                            \
      if (((PetscObject)(h))->classid == PETSCFREEDHEADER) SETERRQ1(PETSC_COMM_SELF,PETSC_ERR_ARG_CORRUPT,"Object already free: Parameter # %d",arg); \
      else SETERRQ1(PETSC_COMM_SELF,PETSC_ERR_ARG_WRONG,"Wrong type of object: Parameter # %d",arg); \
    }                                                                   \
  } while (0)

#define PetscValidHeader(h,arg)                                         \
  do {                                                                  \
    if (!h) SETERRQ1(PETSC_COMM_SELF,PETSC_ERR_ARG_NULL,"Null Object: Parameter # %d",arg); \
    if (!PetscCheckPointer(h,PETSC_OBJECT)) SETERRQ1(PETSC_COMM_SELF,PETSC_ERR_ARG_CORRUPT,"Invalid Pointer to Object: Parameter # %d",arg); \
    if (((PetscObject)(h))->classid == PETSCFREEDHEADER) SETERRQ1(PETSC_COMM_SELF,PETSC_ERR_ARG_CORRUPT,"Object already free: Parameter # %d",arg); \
    else if (((PetscObject)(h))->classid < PETSC_SMALLEST_CLASSID || ((PetscObject)(h))->classid > PETSC_LARGEST_CLASSID) SETERRQ1(PETSC_COMM_SELF,PETSC_ERR_ARG_CORRUPT,"Invalid type of object: Parameter # %d",arg); \
  } while (0)

#define PetscValidPointer(h,arg)                                        \
  do {                                                                  \
    if (!h) SETERRQ1(PETSC_COMM_SELF,PETSC_ERR_ARG_NULL,"Null Pointer: Parameter # %d",arg); \
    if (!PetscCheckPointer(h,PETSC_CHAR)) SETERRQ1(PETSC_COMM_SELF,PETSC_ERR_ARG_BADPTR,"Invalid Pointer: Parameter # %d",arg); \
  } while (0)

#define PetscValidCharPointer(h,arg)                                    \
  do {                                                                  \
    if (!h) SETERRQ1(PETSC_COMM_SELF,PETSC_ERR_ARG_NULL,"Null Pointer: Parameter # %d",arg);\
    if (!PetscCheckPointer(h,PETSC_CHAR)) SETERRQ1(PETSC_COMM_SELF,PETSC_ERR_ARG_BADPTR,"Invalid Pointer to char: Parameter # %d",arg); \
  } while (0)

#define PetscValidIntPointer(h,arg)                                     \
  do {                                                                  \
    if (!h) SETERRQ1(PETSC_COMM_SELF,PETSC_ERR_ARG_BADPTR,"Null Pointer: Parameter # %d",arg); \
    if (!PetscCheckPointer(h,PETSC_INT)) SETERRQ1(PETSC_COMM_SELF,PETSC_ERR_ARG_BADPTR,"Invalid Pointer to Int: Parameter # %d",arg); \
  } while (0)

#define PetscValidScalarPointer(h,arg)                                  \
  do {                                                                  \
    if (!h) SETERRQ1(PETSC_COMM_SELF,PETSC_ERR_ARG_NULL,"Null Pointer: Parameter # %d",arg); \
    if (!PetscCheckPointer(h,PETSC_SCALAR)) SETERRQ1(PETSC_COMM_SELF,PETSC_ERR_ARG_BADPTR,"Invalid Pointer to PetscScalar: Parameter # %d",arg); \
  } while (0)

#define PetscValidRealPointer(h,arg)                                  \
  do {                                                                  \
    if (!h) SETERRQ1(PETSC_COMM_SELF,PETSC_ERR_ARG_NULL,"Null Pointer: Parameter # %d",arg); \
    if (!PetscCheckPointer(h,PETSC_REAL)) SETERRQ1(PETSC_COMM_SELF,PETSC_ERR_ARG_BADPTR,"Invalid Pointer to PetscReal: Parameter # %d",arg); \
  } while (0)

#define PetscValidFunction(f,arg)                                       \
  do {                                                                  \
    if (!f) SETERRQ1(PETSC_COMM_SELF,PETSC_ERR_ARG_NULL,"Null Function Pointer: Parameter # %d",arg); \
  } while (0)

#endif

#if !defined(PETSC_USE_DEBUG)

#define PetscCheckSameType(a,arga,b,argb) do {} while (0)
#define PetscCheckTypeName(a,type) do {} while (0)
#define PetscCheckTypeNames(a,type1,type2) do {} while (0)
#define PetscValidType(a,arg) do {} while (0)
#define PetscCheckSameComm(a,arga,b,argb) do {} while (0)
#define PetscCheckSameTypeAndComm(a,arga,b,argb) do {} while (0)
#define PetscValidLogicalCollectiveScalar(a,b,c) do {} while (0)
#define PetscValidLogicalCollectiveReal(a,b,c) do {} while (0)
#define PetscValidLogicalCollectiveInt(a,b,c) do {} while (0)
#define PetscValidLogicalCollectiveMPIInt(a,b,c) do {} while (0)
#define PetscValidLogicalCollectiveBool(a,b,c) do {} while (0)
#define PetscValidLogicalCollectiveEnum(a,b,c) do {} while (0)

#else

/*
    For example, in the dot product between two vectors,
  both vectors must be either Seq or MPI, not one of each
*/
#define PetscCheckSameType(a,arga,b,argb) \
  if (((PetscObject)a)->type != ((PetscObject)b)->type) SETERRQ2(PETSC_COMM_SELF,PETSC_ERR_ARG_NOTSAMETYPE,"Objects not of same type: Argument # %d and %d",arga,argb);
/*
    Check type_name
*/
#define PetscCheckTypeName(a,type)                                      \
  do {                                                                  \
    PetscBool      __match;                                             \
    PetscErrorCode _7_ierr;                                             \
    _7_ierr = PetscObjectTypeCompare(((PetscObject)a),(type),&__match);CHKERRQ(_7_ierr); \
    if (!__match) SETERRQ2(PETSC_COMM_SELF,PETSC_ERR_ARG_WRONG,"Object (%s) is not %s",(char*)(((PetscObject)a)->type_name),type); \
  } while (0)

#define PetscCheckTypeNames(a,type1,type2)                              \
  do {                                                                  \
    PetscBool       __match;                                            \
    PetscErrorCode _7_ierr;                                             \
    _7_ierr = PetscObjectTypeCompareAny(((PetscObject)a),&__match,(type1),(type2),"");CHKERRQ(_7_ierr); \
    if (!__match) SETERRQ3(PETSC_COMM_SELF,PETSC_ERR_ARG_WRONG,"Object (%s) is not %s or %s",(char*)(((PetscObject)a)->type_name),type1,type2); \
  } while (0)
/*
   Use this macro to check if the type is set
*/
#define PetscValidType(a,arg) \
  if (!((PetscObject)a)->type_name) SETERRQ2(PETSC_COMM_SELF,PETSC_ERR_ARG_WRONGSTATE,"%s object's type is not set: Argument # %d",((PetscObject)a)->class_name,arg);
/*
   Sometimes object must live on same communicator to inter-operate
*/
#define PetscCheckSameComm(a,arga,b,argb)                               \
  do {                                                                  \
    PetscErrorCode _6_ierr,__flag;                                      \
    _6_ierr = MPI_Comm_compare(PetscObjectComm((PetscObject)a),PetscObjectComm((PetscObject)b),&__flag);CHKERRQ(_6_ierr);                                                   \
    if (__flag != MPI_CONGRUENT && __flag != MPI_IDENT) SETERRQ3(PETSC_COMM_SELF,PETSC_ERR_ARG_NOTSAMECOMM,"Different communicators in the two objects: Argument # %d and %d flag %d",arga,argb,__flag); \
  } while (0)

#define PetscCheckSameTypeAndComm(a,arga,b,argb)        \
  do {                                                  \
    PetscCheckSameType(a,arga,b,argb);                  \
    PetscCheckSameComm(a,arga,b,argb);                  \
  } while (0)

#define PetscValidLogicalCollectiveScalar(a,b,c)                        \
  do {                                                                  \
    PetscErrorCode _7_ierr;                                             \
    PetscReal b1[5],b2[5];                                              \
    if (PetscIsNanScalar(b)) {b1[4] = 1;} else {b1[4] = 0;};            \
    b1[0] = -PetscRealPart(b); b1[1] = PetscRealPart(b);b1[2] = -PetscImaginaryPart(b); b1[3] = PetscImaginaryPart(b);         \
    _7_ierr = MPI_Allreduce(b1,b2,5,MPIU_REAL,MPIU_MAX,PetscObjectComm((PetscObject)a));CHKERRQ(_7_ierr); \
    if (!(b2[4] > 0) && !(PetscEqualReal(-b2[0],b2[1]) && PetscEqualReal(-b2[2],b2[3]))) SETERRQ1(PetscObjectComm((PetscObject)a),PETSC_ERR_ARG_WRONG,"Scalar value must be same on all processes, argument # %d",c); \
  } while (0)

#define PetscValidLogicalCollectiveReal(a,b,c)                          \
  do {                                                                  \
    PetscErrorCode _7_ierr;                                             \
    PetscReal b1[3],b2[3];                                              \
    if (PetscIsNanReal(b)) {b1[2] = 1;} else {b1[2] = 0;};              \
    b1[0] = -b; b1[1] = b;                                              \
    _7_ierr = MPI_Allreduce(b1,b2,3,MPIU_REAL,MPIU_MAX,PetscObjectComm((PetscObject)a));CHKERRQ(_7_ierr); \
    if (!(b2[2] > 0) && !PetscEqualReal(-b2[0],b2[1])) SETERRQ1(PetscObjectComm((PetscObject)a),PETSC_ERR_ARG_WRONG,"Real value must be same on all processes, argument # %d",c); \
  } while (0)

#define PetscValidLogicalCollectiveInt(a,b,c)                           \
  do {                                                                  \
    PetscErrorCode _7_ierr;                                             \
    PetscInt b1[2],b2[2];                                               \
    b1[0] = -b; b1[1] = b;                                              \
    _7_ierr = MPIU_Allreduce(b1,b2,2,MPIU_INT,MPI_MAX,PetscObjectComm((PetscObject)a));CHKERRQ(_7_ierr); \
    if (-b2[0] != b2[1]) SETERRQ1(PetscObjectComm((PetscObject)a),PETSC_ERR_ARG_WRONG,"Int value must be same on all processes, argument # %d",c); \
  } while (0)

#define PetscValidLogicalCollectiveMPIInt(a,b,c) do {} while (0)

#define PetscValidLogicalCollectiveBool(a,b,c)                          \
  do {                                                                  \
    PetscErrorCode _7_ierr;                                             \
    PetscMPIInt b1[2],b2[2];                                            \
    b1[0] = -(PetscMPIInt)b; b1[1] = (PetscMPIInt)b;                    \
    _7_ierr = MPIU_Allreduce(b1,b2,2,MPI_INT,MPI_MAX,PetscObjectComm((PetscObject)a));CHKERRQ(_7_ierr); \
    if (-b2[0] != b2[1]) SETERRQ1(PetscObjectComm((PetscObject)a),PETSC_ERR_ARG_WRONG,"Bool value must be same on all processes, argument # %d",c); \
  } while (0)

#define PetscValidLogicalCollectiveEnum(a,b,c)                          \
  do {                                                                  \
    PetscErrorCode _7_ierr;                                             \
    PetscMPIInt b1[2],b2[2];                                            \
    b1[0] = -(PetscMPIInt)b; b1[1] = (PetscMPIInt)b;                    \
    _7_ierr = MPIU_Allreduce(b1,b2,2,MPI_INT,MPI_MAX,PetscObjectComm((PetscObject)a));CHKERRQ(_7_ierr); \
    if (-b2[0] != b2[1]) SETERRQ1(PetscObjectComm((PetscObject)a),PETSC_ERR_ARG_WRONG,"Enum value must be same on all processes, argument # %d",c); \
  } while (0)

#endif

/*
   PetscTryMethod - Queries an object for a method, if it exists then calls it.
              These are intended to be used only inside PETSc functions.

   Level: developer

.seealso: PetscUseMethod()
*/
#define  PetscTryMethod(obj,A,B,C) \
  0;{ PetscErrorCode (*f)B, __ierr; \
    __ierr = PetscObjectQueryFunction((PetscObject)obj,A,&f);CHKERRQ(__ierr); \
    if (f) {__ierr = (*f)C;CHKERRQ(__ierr);}\
  }

/*
   PetscUseMethod - Queries an object for a method, if it exists then calls it, otherwise generates an error.
              These are intended to be used only inside PETSc functions.

   Level: developer

.seealso: PetscTryMethod()
*/
#define  PetscUseMethod(obj,A,B,C) \
  0;{ PetscErrorCode (*f)B, __ierr; \
    __ierr = PetscObjectQueryFunction((PetscObject)obj,A,&f);CHKERRQ(__ierr); \
    if (f) {__ierr = (*f)C;CHKERRQ(__ierr);}\
    else SETERRQ1(PetscObjectComm((PetscObject)obj),PETSC_ERR_SUP,"Cannot locate function %s in object",A); \
  }

/*MC
   PetscObjectStateIncrease - Increases the state of any PetscObject

   Synopsis:
   #include "petsc/private/petscimpl.h"
   PetscErrorCode PetscObjectStateIncrease(PetscObject obj)

   Logically Collective

   Input Parameter:
.  obj - any PETSc object, for example a Vec, Mat or KSP. This must be
         cast with a (PetscObject), for example,
         PetscObjectStateIncrease((PetscObject)mat);

   Notes:
    object state is an integer which gets increased every time
   the object is changed internally. By saving and later querying the object state
   one can determine whether information about the object is still current.
   Currently, state is maintained for Vec and Mat objects.

   This routine is mostly for internal use by PETSc; a developer need only
   call it after explicit access to an object's internals. Routines such
   as VecSet() or MatScale() already call this routine. It is also called, as a
   precaution, in VecRestoreArray(), MatRestoreRow(), MatDenseRestoreArray().

   This routine is logically collective because state equality comparison needs to be possible without communication.

   Level: developer

   seealso: PetscObjectStateGet()

   Concepts: state

M*/
#define PetscObjectStateIncrease(obj) ((obj)->state++,0)

PETSC_EXTERN PetscErrorCode PetscObjectStateGet(PetscObject,PetscObjectState*);
PETSC_EXTERN PetscErrorCode PetscObjectStateSet(PetscObject,PetscObjectState);
PETSC_EXTERN PetscErrorCode PetscObjectComposedDataRegister(PetscInt*);
PETSC_EXTERN PetscErrorCode PetscObjectComposedDataIncreaseInt(PetscObject);
PETSC_EXTERN PetscErrorCode PetscObjectComposedDataIncreaseIntstar(PetscObject);
PETSC_EXTERN PetscErrorCode PetscObjectComposedDataIncreaseReal(PetscObject);
PETSC_EXTERN PetscErrorCode PetscObjectComposedDataIncreaseRealstar(PetscObject);
PETSC_EXTERN PetscErrorCode PetscObjectComposedDataIncreaseScalar(PetscObject);
PETSC_EXTERN PetscErrorCode PetscObjectComposedDataIncreaseScalarstar(PetscObject);
PETSC_EXTERN PetscInt         PetscObjectComposedDataMax;
/*MC
   PetscObjectComposedDataSetInt - attach integer data to a PetscObject

   Synopsis:
   #include "petsc/private/petscimpl.h"
   PetscErrorCode PetscObjectComposedDataSetInt(PetscObject obj,int id,int data)

   Not collective

   Input parameters:
+  obj - the object to which data is to be attached
.  id - the identifier for the data
-  data - the data to  be attached

   Notes
   The data identifier can best be created through a call to  PetscObjectComposedDataRegister()

   Level: developer
M*/
#define PetscObjectComposedDataSetInt(obj,id,data)                                      \
  ((((obj)->int_idmax < PetscObjectComposedDataMax) && PetscObjectComposedDataIncreaseInt(obj)) ||  \
   ((obj)->intcomposeddata[id] = data,(obj)->intcomposedstate[id] = (obj)->state, 0))

/*MC
   PetscObjectComposedDataGetInt - retrieve integer data attached to an object

   Synopsis:
   #include "petsc/private/petscimpl.h"
   PetscErrorCode PetscObjectComposedDataGetInt(PetscObject obj,int id,int data,PetscBool  flag)

   Not collective

   Input parameters:
+  obj - the object from which data is to be retrieved
-  id - the identifier for the data

   Output parameters
+  data - the data to be retrieved
-  flag - PETSC_TRUE if the data item exists and is valid, PETSC_FALSE otherwise

   The 'data' and 'flag' variables are inlined, so they are not pointers.

   Level: developer
M*/
#define PetscObjectComposedDataGetInt(obj,id,data,flag)                            \
  ((((obj)->intcomposedstate && ((obj)->intcomposedstate[id] == (obj)->state)) ?   \
   (data = (obj)->intcomposeddata[id],flag = PETSC_TRUE) : (flag = PETSC_FALSE)),0)

/*MC
   PetscObjectComposedDataSetIntstar - attach integer array data to a PetscObject

   Synopsis:
   #include "petsc/private/petscimpl.h"
   PetscErrorCode PetscObjectComposedDataSetIntstar(PetscObject obj,int id,int *data)

   Not collective

   Input parameters:
+  obj - the object to which data is to be attached
.  id - the identifier for the data
-  data - the data to  be attached

   Notes
   The data identifier can best be determined through a call to
   PetscObjectComposedDataRegister()

   Level: developer
M*/
#define PetscObjectComposedDataSetIntstar(obj,id,data)                                          \
  ((((obj)->intstar_idmax < PetscObjectComposedDataMax) && PetscObjectComposedDataIncreaseIntstar(obj)) ||  \
   ((obj)->intstarcomposeddata[id] = data,(obj)->intstarcomposedstate[id] = (obj)->state, 0))

/*MC
   PetscObjectComposedDataGetIntstar - retrieve integer array data
   attached to an object

   Synopsis:
   #include "petsc/private/petscimpl.h"
   PetscErrorCode PetscObjectComposedDataGetIntstar(PetscObject obj,int id,int *data,PetscBool  flag)

   Not collective

   Input parameters:
+  obj - the object from which data is to be retrieved
-  id - the identifier for the data

   Output parameters
+  data - the data to be retrieved
-  flag - PETSC_TRUE if the data item exists and is valid, PETSC_FALSE otherwise

   The 'data' and 'flag' variables are inlined, so they are not pointers.

   Level: developer
M*/
#define PetscObjectComposedDataGetIntstar(obj,id,data,flag)                               \
  ((((obj)->intstarcomposedstate && ((obj)->intstarcomposedstate[id] == (obj)->state)) ?  \
   (data = (obj)->intstarcomposeddata[id],flag = PETSC_TRUE) : (flag = PETSC_FALSE)),0)

/*MC
   PetscObjectComposedDataSetReal - attach real data to a PetscObject

   Synopsis:
   #include "petsc/private/petscimpl.h"
   PetscErrorCode PetscObjectComposedDataSetReal(PetscObject obj,int id,PetscReal data)

   Not collective

   Input parameters:
+  obj - the object to which data is to be attached
.  id - the identifier for the data
-  data - the data to  be attached

   Notes
   The data identifier can best be determined through a call to
   PetscObjectComposedDataRegister()

   Level: developer
M*/
#define PetscObjectComposedDataSetReal(obj,id,data)                                       \
  ((((obj)->real_idmax < PetscObjectComposedDataMax) && PetscObjectComposedDataIncreaseReal(obj)) ||  \
   ((obj)->realcomposeddata[id] = data,(obj)->realcomposedstate[id] = (obj)->state, 0))

/*MC
   PetscObjectComposedDataGetReal - retrieve real data attached to an object

   Synopsis:
   #include "petsc/private/petscimpl.h"
   PetscErrorCode PetscObjectComposedDataGetReal(PetscObject obj,int id,PetscReal data,PetscBool  flag)

   Not collective

   Input parameters:
+  obj - the object from which data is to be retrieved
-  id - the identifier for the data

   Output parameters
+  data - the data to be retrieved
-  flag - PETSC_TRUE if the data item exists and is valid, PETSC_FALSE otherwise

   The 'data' and 'flag' variables are inlined, so they are not pointers.

   Level: developer
M*/
#define PetscObjectComposedDataGetReal(obj,id,data,flag)                            \
  ((((obj)->realcomposedstate && ((obj)->realcomposedstate[id] == (obj)->state)) ?  \
   (data = (obj)->realcomposeddata[id],flag = PETSC_TRUE) : (flag = PETSC_FALSE)),0)

/*MC
   PetscObjectComposedDataSetRealstar - attach real array data to a PetscObject

   Synopsis:
   #include "petsc/private/petscimpl.h"
   PetscErrorCode PetscObjectComposedDataSetRealstar(PetscObject obj,int id,PetscReal *data)

   Not collective

   Input parameters:
+  obj - the object to which data is to be attached
.  id - the identifier for the data
-  data - the data to  be attached

   Notes
   The data identifier can best be determined through a call to
   PetscObjectComposedDataRegister()

   Level: developer
M*/
#define PetscObjectComposedDataSetRealstar(obj,id,data)                                           \
  ((((obj)->realstar_idmax < PetscObjectComposedDataMax) && PetscObjectComposedDataIncreaseRealstar(obj)) ||  \
   ((obj)->realstarcomposeddata[id] = data, (obj)->realstarcomposedstate[id] = (obj)->state, 0))

/*MC
   PetscObjectComposedDataGetRealstar - retrieve real array data
   attached to an object

   Synopsis:
   #include "petsc/private/petscimpl.h"
   PetscErrorCode PetscObjectComposedDataGetRealstar(PetscObject obj,int id,PetscReal *data,PetscBool  flag)

   Not collective

   Input parameters:
+  obj - the object from which data is to be retrieved
-  id - the identifier for the data

   Output parameters
+  data - the data to be retrieved
-  flag - PETSC_TRUE if the data item exists and is valid, PETSC_FALSE otherwise

   The 'data' and 'flag' variables are inlined, so they are not pointers.

   Level: developer
M*/
#define PetscObjectComposedDataGetRealstar(obj,id,data,flag)                                \
  ((((obj)->realstarcomposedstate && ((obj)->realstarcomposedstate[id] == (obj)->state)) ?  \
   (data = (obj)->realstarcomposeddata[id],flag = PETSC_TRUE) : (flag = PETSC_FALSE)),0)

/*MC
   PetscObjectComposedDataSetScalar - attach scalar data to a PetscObject

   Synopsis:
   #include "petsc/private/petscimpl.h"
   PetscErrorCode PetscObjectComposedDataSetScalar(PetscObject obj,int id,PetscScalar data)

   Not collective

   Input parameters:
+  obj - the object to which data is to be attached
.  id - the identifier for the data
-  data - the data to  be attached

   Notes
   The data identifier can best be determined through a call to
   PetscObjectComposedDataRegister()

   Level: developer
M*/
#if defined(PETSC_USE_COMPLEX)
#define PetscObjectComposedDataSetScalar(obj,id,data)                                        \
  ((((obj)->scalar_idmax < PetscObjectComposedDataMax) && PetscObjectComposedDataIncreaseScalar(obj)) || \
   ((obj)->scalarcomposeddata[id] = data,(obj)->scalarcomposedstate[id] = (obj)->state, 0))
#else
#define PetscObjectComposedDataSetScalar(obj,id,data) \
        PetscObjectComposedDataSetReal(obj,id,data)
#endif
/*MC
   PetscObjectComposedDataGetScalar - retrieve scalar data attached to an object

   Synopsis:
   #include "petsc/private/petscimpl.h"
   PetscErrorCode PetscObjectComposedDataGetScalar(PetscObject obj,int id,PetscScalar data,PetscBool  flag)

   Not collective

   Input parameters:
+  obj - the object from which data is to be retrieved
-  id - the identifier for the data

   Output parameters
+  data - the data to be retrieved
-  flag - PETSC_TRUE if the data item exists and is valid, PETSC_FALSE otherwise

   The 'data' and 'flag' variables are inlined, so they are not pointers.

   Level: developer
M*/
#if defined(PETSC_USE_COMPLEX)
#define PetscObjectComposedDataGetScalar(obj,id,data,flag)                              \
  ((((obj)->scalarcomposedstate && ((obj)->scalarcomposedstate[id] == (obj)->state) ) ? \
   (data = (obj)->scalarcomposeddata[id],flag = PETSC_TRUE) : (flag = PETSC_FALSE)),0)
#else
#define PetscObjectComposedDataGetScalar(obj,id,data,flag)                             \
        PetscObjectComposedDataGetReal(obj,id,data,flag)
#endif

/*MC
   PetscObjectComposedDataSetScalarstar - attach scalar array data to a PetscObject

   Synopsis:
   #include "petsc/private/petscimpl.h"
   PetscErrorCode PetscObjectComposedDataSetScalarstar(PetscObject obj,int id,PetscScalar *data)

   Not collective

   Input parameters:
+  obj - the object to which data is to be attached
.  id - the identifier for the data
-  data - the data to  be attached

   Notes
   The data identifier can best be determined through a call to
   PetscObjectComposedDataRegister()

   Level: developer
M*/
#if defined(PETSC_USE_COMPLEX)
#define PetscObjectComposedDataSetScalarstar(obj,id,data)                                             \
  ((((obj)->scalarstar_idmax < PetscObjectComposedDataMax) && PetscObjectComposedDataIncreaseScalarstar(obj)) ||  \
   ((obj)->scalarstarcomposeddata[id] = data,(obj)->scalarstarcomposedstate[id] = (obj)->state, 0))
#else
#define PetscObjectComposedDataSetScalarstar(obj,id,data) \
        PetscObjectComposedDataSetRealstar(obj,id,data)
#endif
/*MC
   PetscObjectComposedDataGetScalarstar - retrieve scalar array data
   attached to an object

   Synopsis:
   #include "petsc/private/petscimpl.h"
   PetscErrorCode PetscObjectComposedDataGetScalarstar(PetscObject obj,int id,PetscScalar *data,PetscBool  flag)

   Not collective

   Input parameters:
+  obj - the object from which data is to be retrieved
-  id - the identifier for the data

   Output parameters
+  data - the data to be retrieved
-  flag - PETSC_TRUE if the data item exists and is valid, PETSC_FALSE otherwise

   The 'data' and 'flag' variables are inlined, so they are not pointers.

   Level: developer
M*/
#if defined(PETSC_USE_COMPLEX)
#define PetscObjectComposedDataGetScalarstar(obj,id,data,flag)                                 \
  ((((obj)->scalarstarcomposedstate && ((obj)->scalarstarcomposedstate[id] == (obj)->state)) ? \
       (data = (obj)->scalarstarcomposeddata[id],flag = PETSC_TRUE) : (flag = PETSC_FALSE)),0)
#else
#define PetscObjectComposedDataGetScalarstar(obj,id,data,flag)         \
        PetscObjectComposedDataGetRealstar(obj,id,data,flag)
#endif

PETSC_EXTERN PetscErrorCode PetscObjectGetId(PetscObject,PetscObjectId*);

PETSC_EXTERN PetscErrorCode PetscMonitorCompare(PetscErrorCode (*)(void),void *,PetscErrorCode (*)(void**),PetscErrorCode (*)(void),void *,PetscErrorCode (*)(void**),PetscBool *);

PETSC_EXTERN PetscMPIInt Petsc_Counter_keyval;
PETSC_EXTERN PetscMPIInt Petsc_InnerComm_keyval;
PETSC_EXTERN PetscMPIInt Petsc_OuterComm_keyval;
PETSC_EXTERN PetscMPIInt Petsc_Seq_keyval;
PETSC_EXTERN PetscMPIInt Petsc_ShmComm_keyval;

/*
  PETSc communicators have this attribute, see
  PetscCommDuplicate(), PetscCommDestroy(), PetscCommGetNewTag(), PetscObjectGetName()
*/
typedef struct {
  PetscMPIInt tag;              /* next free tag value */
  PetscInt    refcount;         /* number of references, communicator can be freed when this reaches 0 */
  PetscInt    namecount;        /* used to generate the next name, as in Vec_0, Mat_1, ... */
} PetscCommCounter;

/*E
    PetscOffloadFlag - indicates which memory (CPU, GPU, or none contains valid vector

   PETSC_OFFLOAD_UNALLOCATED  - no memory contains valid matrix entries; NEVER used for vectors
   PETSC_OFFLOAD_GPU - GPU has valid vector/matrix entries
   PETSC_OFFLOAD_CPU - CPU has valid vector/matrix entries
   PETSC_OFFLOAD_BOTH - Both GPU and CPU have valid vector/matrix entries and they match

   Level: developer
E*/
typedef enum {PETSC_OFFLOAD_UNALLOCATED,PETSC_OFFLOAD_GPU,PETSC_OFFLOAD_CPU,PETSC_OFFLOAD_BOTH} PetscOffloadFlag;

typedef enum {STATE_BEGIN, STATE_PENDING, STATE_END} SRState;

typedef enum {PETSC_SR_REDUCE_SUM=0,PETSC_SR_REDUCE_MAX=1,PETSC_SR_REDUCE_MIN=2} PetscSRReductionType;

typedef struct {
  MPI_Comm    comm;
  MPI_Request request;
  PetscBool   async;
  PetscScalar *lvalues;     /* this are the reduced values before call to MPI_Allreduce() */
  PetscScalar *gvalues;     /* values after call to MPI_Allreduce() */
  void        **invecs;     /* for debugging only, vector/memory used with each op */
  PetscInt    *reducetype;  /* is particular value to be summed or maxed? */
  SRState     state;        /* are we calling xxxBegin() or xxxEnd()? */
  PetscInt    maxops;       /* total amount of space we have for requests */
  PetscInt    numopsbegin;  /* number of requests that have been queued in */
  PetscInt    numopsend;    /* number of requests that have been gotten by user */
} PetscSplitReduction;

PETSC_EXTERN PetscErrorCode PetscSplitReductionGet(MPI_Comm,PetscSplitReduction**);
PETSC_EXTERN PetscErrorCode PetscSplitReductionEnd(PetscSplitReduction*);
PETSC_EXTERN PetscErrorCode PetscSplitReductionExtend(PetscSplitReduction*);

#if !defined(PETSC_SKIP_SPINLOCK)
#if defined(PETSC_HAVE_THREADSAFETY)
#  if defined(PETSC_HAVE_CONCURRENCYKIT)
#if defined(__cplusplus)
/*  CK does not have extern "C" protection in their include files */
extern "C" {
#endif
#include <ck_spinlock.h>
#if defined(__cplusplus)
}
#endif
typedef ck_spinlock_t PetscSpinlock;
PETSC_STATIC_INLINE PetscErrorCode PetscSpinlockCreate(PetscSpinlock *ck_spinlock)
{
  ck_spinlock_init(ck_spinlock);
  return 0;
}
PETSC_STATIC_INLINE PetscErrorCode PetscSpinlockLock(PetscSpinlock *ck_spinlock)
{
  ck_spinlock_lock(ck_spinlock);
  return 0;
}
PETSC_STATIC_INLINE PetscErrorCode PetscSpinlockUnlock(PetscSpinlock *ck_spinlock)
{
  ck_spinlock_unlock(ck_spinlock);
  return 0;
}
PETSC_STATIC_INLINE PetscErrorCode PetscSpinlockDestroy(PetscSpinlock *ck_spinlock)
{
  return 0;
}
#  elif defined(PETSC_HAVE_OPENMP)

#include <omp.h>
typedef omp_lock_t PetscSpinlock;
PETSC_STATIC_INLINE PetscErrorCode PetscSpinlockCreate(PetscSpinlock *omp_lock)
{
  omp_init_lock(omp_lock);
  return 0;
}
PETSC_STATIC_INLINE PetscErrorCode PetscSpinlockLock(PetscSpinlock *omp_lock)
{
  omp_set_lock(omp_lock);
  return 0;
}
PETSC_STATIC_INLINE PetscErrorCode PetscSpinlockUnlock(PetscSpinlock *omp_lock)
{
  omp_unset_lock(omp_lock);
  return 0;
}
PETSC_STATIC_INLINE PetscErrorCode PetscSpinlockDestroy(PetscSpinlock *omp_lock)
{
  omp_destroy_lock(omp_lock);
  return 0;
}
#else
Thread safety requires either --with-openmp or --download-concurrencykit
#endif

#else
typedef int PetscSpinlock;
#define PetscSpinlockCreate(a)  0
#define PetscSpinlockLock(a)    0
#define PetscSpinlockUnlock(a)  0
#define PetscSpinlockDestroy(a) 0
#endif

#if defined(PETSC_HAVE_THREADSAFETY)
PETSC_INTERN PetscSpinlock PetscViewerASCIISpinLockOpen;
PETSC_INTERN PetscSpinlock PetscViewerASCIISpinLockStdout;
PETSC_INTERN PetscSpinlock PetscViewerASCIISpinLockStderr;
PETSC_INTERN PetscSpinlock PetscCommSpinLock;
#endif
#endif

PETSC_EXTERN PetscLogEvent PETSC_Barrier;
PETSC_EXTERN PetscLogEvent PETSC_BuildTwoSided;
PETSC_EXTERN PetscLogEvent PETSC_BuildTwoSidedF;

#if defined(PETSC_HAVE_ADIOS)
PETSC_EXTERN int64_t Petsc_adios_group;
#endif

#endif /* _PETSCHEAD_H */
