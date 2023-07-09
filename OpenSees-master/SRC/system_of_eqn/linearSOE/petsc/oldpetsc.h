/* ****************************************************************** **
**    OpenSees - Open System for Earthquake Engineering Simulation    **
**          Pacific Earthquake Engineering Research Center            **
**                                                                    **
**                                                                    **
** (C) Copyright 1999, The Regents of the University of California    **
** All Rights Reserved.                                               **
**                                                                    **
** Commercial use of this program without express permission of the   **
** University of California, Berkeley, is strictly prohibited.  See   **
** file 'COPYRIGHT'  in main directory for information on usage and   **
** redistribution,  and for a DISCLAIMER OF ALL WARRANTIES.           **
**                                                                    **
** Developed by:                                                      **
**   Frank McKenna (fmckenna@ce.berkeley.edu)                         **
**   Gregory L. Fenves (fenves@ce.berkeley.edu)                       **
**   Filip C. Filippou (filippou@ce.berkeley.edu)                     **
**                                                                    **
** ****************************************************************** */
                                                                        
// $Revision: 1.2 $
// $Date: 2003-02-14 23:02:02 $
// $Source: /usr/local/cvs/OpenSees/SRC/system_of_eqn/linearSOE/petsc/oldpetsc.h,v $
                                                                        
                                                                        
/* $Id: oldpetsc.h,v 1.2 2003-02-14 23:02:02 fmk Exp $ */
/*
   This is the main PETSc include file (for C and C++).  It is included by all
   other PETSc include files, so it almost never has to be specifically included.
*/
#if !defined(__PETSC_PACKAGE)
#define __PETSC_PACKAGE

/* 
   Current PETSc Version 
*/
#define PETSC_VERSION_NUMBER "PETSc Version 2.0.22, Released April 28, 1998."

#define PETSC_VERSION_MAJOR    2
#define PETSC_VERSION_MINOR    0
#define PETSC_VERSION_SUBMINOR 22
#define PETSC_VERSION_DATE     "April 29, 1998"
#define PETSC_AUTHOR_INFO      "The PETSc Team:\
 Satish Balay, Bill Gropp, Lois Curfman McInnes, Barry Smith\n\
 Bug reports, questions: petsc-maint@mcs.anl.gov\n\
 Web page: http://www.mcs.anl.gov/petsc/\n"

/* ========================================================================== */
/* Before anything else, include the PETSc configuration file.  This 
   contains various definitions that handle portability issues and the 
   presence of important features. 

   petscconf.h is contained in bmake/${PETSC_ARCH}/petscconf.h  
*/
#include <petscconf.h>


/* ========================================================================== */

#include <stdio.h>
/*
    Defines the interface to MPI allowing the use of all MPI functions.
*/
#include "mpi.h"

/*
    Defines some elementary mathematics functions and constants.
*/
#include "petscmath.h"
/*
    This shouuld be in petscmath.h?
*/
#if defined(USE_POINTER_CONVERSION)
#define PetscFortranAddr   int
#else
#define PetscFortranAddr   long
#endif

extern MPI_Comm PETSC_COMM_WORLD;
extern MPI_Comm PETSC_COMM_SELF;
extern int      PetscInitializedCalled;
extern int      PetscSetCommWorld(MPI_Comm);

/*
    Defines the malloc employed by PETSc. Users may use these routines as well. 
*/
#define PetscMalloc(a)       (*PetscTrMalloc)(a,__LINE__,__FUNC__,__FILE__,__SDIR__)
#define PetscNew(A)          (A*) PetscMalloc(sizeof(A))
#define PetscFree(a)         (*PetscTrFree)(a,__LINE__,__FUNC__,__FILE__,__SDIR__)
extern void *(*PetscTrMalloc)(unsigned int,int,char*,char*,char*);
extern int  (*PetscTrFree)(void *,int,char*,char*,char*);
extern int  PetscSetMalloc(void *(*)(unsigned int,int,char*,char*,char*),
                           int (*)(void *,int,char*,char*,char*));
extern int  PetscClearMalloc(void);

extern int   PetscTrDump(FILE *);
extern int   PetscTrSpace( PLogDouble *, PLogDouble *,PLogDouble *);
extern int   PetscTrValid(int,char *,char *,char *);
extern int   PetscTrDebugLevel(int);
extern int   PetscTrLog(void);
extern int   PetscTrLogDump(FILE *);
extern int   PetscGetResidentSetSize(PLogDouble *);

#include <src/inline/bitarray.h>

typedef enum {PETSC_INT = 0, PETSC_DOUBLE = 1, PETSC_SHORT = 2, PETSC_FLOAT = 3,
              PETSC_COMPLEX = 4, PETSC_CHAR = 5, PETSC_LOGICAL = 6} PetscDataType;
#if defined(USE_PETSC_COMPLEX)
#define PETSC_SCALAR PETSC_COMPLEX
#else
#define PETSC_SCALAR PETSC_DOUBLE
#endif

typedef enum {PETSC_INT_SIZE = sizeof(int), PETSC_DOUBLE_SIZE = sizeof(double),
              PETSC_SCALAR_SIZE = sizeof(Scalar), PETSC_COMPLEX_SIZE = sizeof(double),
              PETSC_CHAR_SIZE = sizeof(char), PETSC_LOGICAL_SIZE = 1} PetscDataTypeSize;
extern int PetscDataTypeToMPIDataType(PetscDataType,MPI_Datatype*);
extern int PetscDataTypeGetSize(PetscDataType,int*);
extern int PetscDataTypeGetName(PetscDataType,char**);

/*
    Basic memory and string operations
*/
extern int   PetscMemcpy(void *,void *,int);
extern int   PetscBitMemcpy(void*,int,void*,int,int,PetscDataType);
extern int   PetscMemmove(void *,void *,int);
extern int   PetscMemzero(void *,int);
extern int   PetscMemcmp(void*, void*, int);
extern int   PetscStrlen(char *);
extern int   PetscStrcmp(char *,char *);
extern int   PetscStrcasecmp(char *,char *);
extern int   PetscStrncmp(char *,char *,int );
extern int   PetscStrcpy(char *,char *);
extern int   PetscStrcat(char *,char *);
extern int   PetscStrncat(char *,char *,int);
extern int   PetscStrncpy(char *,char *,int);
extern char* PetscStrchr(char *,char);
extern char* PetscStrrchr(char *,char);
extern char* PetscStrstr(char*,char*);
extern char* PetscStrtok(char*,char*);
extern char* PetscStrrtok(char*,char*);

/*
       Basic PETSc constants
*/
typedef enum { PETSC_FALSE, PETSC_TRUE } PetscTruth;
#define PETSC_NULL            0
#define PETSC_DECIDE         -1
#define PETSC_DETERMINE      PETSC_DECIDE
#define PETSC_DEFAULT        -2

/*
    Each PETSc object class has it's own cookie (internal integer in the 
  data structure used for error checking). These are all defined by an offset 
  from the lowest one, PETSC_COOKIE. If you increase these you must 
  increase the field sizes in petsc/src/plog/src/plog.c
*/
#define PETSC_COOKIE                    1211211
#define LARGEST_PETSC_COOKIE_PREDEFINED PETSC_COOKIE + 30
#define LARGEST_PETSC_COOKIE_ALLOWED    PETSC_COOKIE + 50
extern int LARGEST_PETSC_COOKIE;

#include "viewer.h"
#include "options.h"

/*
    Defines basic graphics available from PETSc.
*/
#include "draw.h"

extern int PetscGetTime(PLogDouble*);
extern int PetscGetCPUTime(PLogDouble*);
extern int PetscSleep(int);

extern int  PetscInitialize(int*,char***,char*,char*);
extern int  PetscInitializeNoArguments(void);
extern int  PetscFinalize(void);
extern void PetscInitializeFortran(void);

/*
    Functions that can act on any PETSc object.
*/
typedef struct _p_PetscObject* PetscObject;
extern int PetscObjectDestroy(PetscObject);
extern int PetscObjectExists(PetscObject,int*);
extern int PetscObjectGetComm(PetscObject,MPI_Comm *comm);
extern int PetscObjectGetCookie(PetscObject,int *cookie);
extern int PetscObjectGetType(PetscObject,int *type);
extern int PetscObjectSetName(PetscObject,char*);
extern int PetscObjectGetName(PetscObject,char**);
extern int PetscObjectReference(PetscObject);
extern int PetscObjectGetReference(PetscObject,int*);
extern int PetscObjectDereference(PetscObject);
extern int PetscObjectGetNewTag(PetscObject,int *);
extern int PetscObjectRestoreNewTag(PetscObject,int *);
extern int PetscObjectView(PetscObject,Viewer);

extern int PetscObjectCompose(PetscObject,char *,PetscObject);
extern int PetscObjectQuery(PetscObject,char *,PetscObject *);
#if defined(USE_DYNAMIC_LIBRARIES)
#define PetscObjectComposeFunction(a,b,c,d) PetscObjectComposeFunction_Private(a,b,c,0)
#else
#define PetscObjectComposeFunction(a,b,c,d) PetscObjectComposeFunction_Private(a,b,c,d)
#endif
extern int PetscObjectComposeFunction_Private(PetscObject,char *,char *,void *);
extern int PetscObjectQueryFunction(PetscObject,char *,void **);


/*
    Defines PETSc error handling.
*/
#include "petsopserror.h"

/*
    Mechanism for managing lists of objects attached (composed) with 
   a PETSc object.
*/
typedef struct _OList *OList;
extern int OListDestroy(OList *);
extern int OListFind(OList,char *,PetscObject*);
extern int OListAdd(OList *,char *,PetscObject);
extern int OListRemove(OList *,char *);
extern int OListDuplicate(OList,OList *);

/*
    Dynamic library lists. Lists of names of routines in dynamic 
  link libraries that will be loaded as needed.
*/
typedef struct _DLList *DLList;
extern int    DLRegister_Private(DLList*,char*,char*,int (*)(void *));
extern int    DLRegisterCreate(DLList *);
extern int    DLRegisterDestroy(DLList);
extern int    DLRegisterFind(MPI_Comm,DLList,char*,int (**)(void*));
extern int    DLRegisterPrintTypes(MPI_Comm,FILE*,char*,char *,DLList);
#if defined(USE_DYNAMIC_LIBRARIES)
#define       DLRegister(a,b,p,c) DLRegister_Private(a,b,p,0)
#else
#define       DLRegister(a,b,p,c) DLRegister_Private(a,b,p,(int (*)(void *))c)
#endif

typedef struct _DLLibraryList *DLLibraryList;
extern DLLibraryList DLLibrariesLoaded;
extern int DLLibraryOpen(MPI_Comm,char *,void **);
extern int DLLibrarySym(MPI_Comm,DLLibraryList *,char*,char *, void **);
extern int DLLibraryAppend(MPI_Comm,DLLibraryList *,char *);
extern int DLLibraryPrepend(MPI_Comm,DLLibraryList *,char *);
extern int DLLibraryClose(DLLibraryList);


#include "petschead.h"

/*
     Defines PETSc profiling.
*/
#include "petsclog.h"

extern int  PetscSequentialPhaseBegin(MPI_Comm,int);
extern int  PetscSequentialPhaseEnd(MPI_Comm,int);

/*M 
    PetscBarrier - Blocks until this routine is executed by all
                   processors owning the object A.

   Input Parameters:
.  A - PETSc object  ( Mat, Vec, IS, SNES etc...)

   Synopsis:
   void PetscBarrier(PetscObject obj)

  Notes: 
  This routine calls MPI_Barrier with the communicator
  of the PETSc Object "A". 

.keywords: barrier, petscobject
M*/

#define PetscBarrier(A) \
  { \
    PetscValidHeader(A); \
    PLogEventBegin(Petsc_Barrier,A,0,0,0); \
    MPI_Barrier(((PetscObject)A)->comm); \
    PLogEventEnd(Petsc_Barrier,A,0,0,0); \
  }

extern int PetscMPIDump(FILE *);

/*
      This code allows one to pass a PETSc object in C
  to a Fortran routine, where (like all PETSc objects in 
  Fortran) it is treated as an integer.
*/
extern int  PetscCObjectToFortranObject(void *,PetscFortranAddr *);
extern int  PetscFortranObjectToCObject(PetscFortranAddr,void *);
extern int  MPICCommToFortranComm(MPI_Comm,int *);
extern int  MPIFortranCommToCComm(int,MPI_Comm*);

/*
      Simple PETSc parallel IO for ASCII printing
*/
extern int  PetscFixFilename(char*);
extern FILE *PetscFOpen(MPI_Comm,char *,char *);
extern int  PetscFClose(MPI_Comm,FILE*);
extern int  PetscFPrintf(MPI_Comm,FILE*,char *,...);
extern int  PetscPrintf(MPI_Comm,char *,...);
extern int  (*PetscErrorPrintf)(char *,...);
extern int  (*PetscHelpPrintf)(MPI_Comm,char *,...);

extern int  PetscSynchronizedPrintf(MPI_Comm,char *,...);
extern int  PetscSynchronizedFPrintf(MPI_Comm,FILE*,char *,...);
extern int  PetscSynchronizedFlush(MPI_Comm);


typedef struct _p_PetscObjectContainer*  PetscObjectContainer;
extern int PetscObjectContainerGetPointer(PetscObjectContainer,void **);
extern int PetscObjectContainerSetPointer(PetscObjectContainer,void *);
extern int PetscObjectContainerDestroy(PetscObjectContainer);
extern int PetscObjectContainerCreate(MPI_Comm comm,PetscObjectContainer *);

/*
    C code optimization is often enhanced by telling the compiler 
  that certain pointer arguments to functions are not aliased to 
  to other arguments. This is not yet ANSI C standard so we define 
  the macro "restrict" to indicate that the variable is not aliased 
  to any other argument.
*/
#if defined(HAVE_RESTRICT) && !defined(__cplusplus)
#define restrict _Restrict
#else
#define restrict
#endif

/*
   For incremental debugging
*/
extern int PetscCompare;
extern int PetscCompareDouble(double);
extern int PetscCompareScalar(Scalar);
extern int PetscCompareInt(int);

/*
   For use in debuggers 
*/
extern int PetscGlobalRank,PetscGlobalSize;
extern int PetscIntView(int,int*,Viewer);
extern int PetscDoubleView(int,double *,Viewer);

#include <mat.h>

#endif
