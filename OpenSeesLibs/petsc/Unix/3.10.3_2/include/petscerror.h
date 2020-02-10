/*
    Contains all error handling interfaces for PETSc.
*/
#if !defined(__PETSCERROR_H)
#define __PETSCERROR_H

/*
     These are the generic error codes. These error codes are used
     many different places in the PETSc source code. The string versions are
     at src/sys/error/err.c any changes here must also be made there
     These are also define in include/petsc/finclude/petscerror.h any CHANGES here
     must be also made there.

*/
#define PETSC_ERR_MIN_VALUE        54   /* should always be one less then the smallest value */

#define PETSC_ERR_MEM              55   /* unable to allocate requested memory */
#define PETSC_ERR_SUP              56   /* no support for requested operation */
#define PETSC_ERR_SUP_SYS          57   /* no support for requested operation on this computer system */
#define PETSC_ERR_ORDER            58   /* operation done in wrong order */
#define PETSC_ERR_SIG              59   /* signal received */
#define PETSC_ERR_FP               72   /* floating point exception */
#define PETSC_ERR_COR              74   /* corrupted PETSc object */
#define PETSC_ERR_LIB              76   /* error in library called by PETSc */
#define PETSC_ERR_PLIB             77   /* PETSc library generated inconsistent data */
#define PETSC_ERR_MEMC             78   /* memory corruption */
#define PETSC_ERR_CONV_FAILED      82   /* iterative method (KSP or SNES) failed */
#define PETSC_ERR_USER             83   /* user has not provided needed function */
#define PETSC_ERR_SYS              88   /* error in system call */
#define PETSC_ERR_POINTER          70   /* pointer does not point to valid address */
#define PETSC_ERR_MPI_LIB_INCOMP   87   /* MPI library at runtime is not compatible with MPI user compiled with */

#define PETSC_ERR_ARG_SIZ          60   /* nonconforming object sizes used in operation */
#define PETSC_ERR_ARG_IDN          61   /* two arguments not allowed to be the same */
#define PETSC_ERR_ARG_WRONG        62   /* wrong argument (but object probably ok) */
#define PETSC_ERR_ARG_CORRUPT      64   /* null or corrupted PETSc object as argument */
#define PETSC_ERR_ARG_OUTOFRANGE   63   /* input argument, out of range */
#define PETSC_ERR_ARG_BADPTR       68   /* invalid pointer argument */
#define PETSC_ERR_ARG_NOTSAMETYPE  69   /* two args must be same object type */
#define PETSC_ERR_ARG_NOTSAMECOMM  80   /* two args must be same communicators */
#define PETSC_ERR_ARG_WRONGSTATE   73   /* object in argument is in wrong state, e.g. unassembled mat */
#define PETSC_ERR_ARG_TYPENOTSET   89   /* the type of the object has not yet been set */
#define PETSC_ERR_ARG_INCOMP       75   /* two arguments are incompatible */
#define PETSC_ERR_ARG_NULL         85   /* argument is null that should not be */
#define PETSC_ERR_ARG_UNKNOWN_TYPE 86   /* type name doesn't match any registered type */

#define PETSC_ERR_FILE_OPEN        65   /* unable to open file */
#define PETSC_ERR_FILE_READ        66   /* unable to read from file */
#define PETSC_ERR_FILE_WRITE       67   /* unable to write to file */
#define PETSC_ERR_FILE_UNEXPECTED  79   /* unexpected data in file */

#define PETSC_ERR_MAT_LU_ZRPVT     71   /* detected a zero pivot during LU factorization */
#define PETSC_ERR_MAT_CH_ZRPVT     81   /* detected a zero pivot during Cholesky factorization */

#define PETSC_ERR_INT_OVERFLOW     84

#define PETSC_ERR_FLOP_COUNT       90
#define PETSC_ERR_NOT_CONVERGED    91  /* solver did not converge */
#define PETSC_ERR_MISSING_FACTOR   92  /* MatGetFactor() failed */
#define PETSC_ERR_OPT_OVERWRITE    93  /* attempted to over wrote options which should not be changed */

#define PETSC_ERR_MAX_VALUE        94  /* this is always the one more than the largest error code */

#define PetscStringizeArg(a) #a
#define PetscStringize(a) PetscStringizeArg(a)


/*MC
   SETERRQ - Macro to be called when an error has been detected,

   Synopsis:
   #include <petscsys.h>
   PetscErrorCode SETERRQ(MPI_Comm comm,PetscErrorCode ierr,char *message)

   Collective on MPI_Comm

   Input Parameters:
+  comm - A communicator, use PETSC_COMM_SELF unless you know all ranks of another communicator will detect the error
.  ierr - nonzero error code, see the list of standard error codes in include/petscerror.h
-  message - error message

  Level: beginner

   Notes:
    Once the error handler is called the calling function is then returned from with the given error code.

    See SETERRQ1(), SETERRQ2(), SETERRQ3() for versions that take arguments

    Experienced users can set the error handler with PetscPushErrorHandler().

   Concepts: error^setting condition

.seealso: PetscTraceBackErrorHandler(), PetscPushErrorHandler(), PetscError(), CHKERRQ(), CHKMEMQ, SETERRQ1(), SETERRQ2(), SETERRQ3()
M*/
#define SETERRQ(comm,ierr,s) return PetscError(comm,__LINE__,PETSC_FUNCTION_NAME,__FILE__,ierr,PETSC_ERROR_INITIAL,s)

/*MC
   SETERRMPI - Macro to be called when an error has been detected within an MPI callback function

   Synopsis:
   #include <petscsys.h>
   PetscErrorCode SETERRMPI(MPI_Comm comm,PetscErrorCode ierr,char *message)

   Collective on MPI_Comm

   Input Parameters:
+  comm - A communicator, use PETSC_COMM_SELF unless you know all ranks of another communicator will detect the error
.  ierr - nonzero error code, see the list of standard error codes in include/petscerror.h
-  message - error message

  Level: developer

   Notes:
    This macro is FOR USE IN MPI CALLBACK FUNCTIONS ONLY, such as those passed to MPI_Comm_create_keyval(). It always returns the error code PETSC_MPI_ERROR_CODE
    which is registered with MPI_Add_error_code() when PETSc is initialized.

   Concepts: error^setting condition

.seealso: SETERRQ(), CHKERRQ(), CHKERRMPI(), PetscTraceBackErrorHandler(), PetscPushErrorHandler(), PetscError(), CHKERRQ(), CHKMEMQ, SETERRQ1(), SETERRQ2(), SETERRQ3()
M*/
#define SETERRMPI(comm,ierr,s) return (PetscError(comm,__LINE__,PETSC_FUNCTION_NAME,__FILE__,ierr,PETSC_ERROR_INITIAL,s),PETSC_MPI_ERROR_CODE)

/*MC
   SETERRQ1 - Macro that is called when an error has been detected,

   Synopsis:
   #include <petscsys.h>
   PetscErrorCode SETERRQ1(MPI_Comm comm,PetscErrorCode ierr,char *formatmessage,arg)

   Collective on MPI_Comm

   Input Parameters:
+  comm - A communicator, so that the error can be collective
.  ierr - nonzero error code, see the list of standard error codes in include/petscerror.h
.  message - error message in the printf format
-  arg - argument (for example an integer, string or double)

  Level: beginner

   Notes:
    Once the error handler is called the calling function is then returned from with the given error code.

   Experienced users can set the error handler with PetscPushErrorHandler().

   Concepts: error^setting condition

.seealso: PetscTraceBackErrorHandler(), PetscPushErrorHandler(), PetscError(), CHKERRQ(), CHKMEMQ, SETERRQ(), SETERRQ2(), SETERRQ3()
M*/
#define SETERRQ1(comm,ierr,s,a1) return PetscError(comm,__LINE__,PETSC_FUNCTION_NAME,__FILE__,ierr,PETSC_ERROR_INITIAL,s,a1)

/*MC
   SETERRQ2 - Macro that is called when an error has been detected,

   Synopsis:
   #include <petscsys.h>
   PetscErrorCode SETERRQ2(MPI_Comm comm,PetscErrorCode ierr,char *formatmessage,arg1,arg2)

   Collective on MPI_Comm

   Input Parameters:
+  comm - A communicator, so that the error can be collective
.  ierr - nonzero error code, see the list of standard error codes in include/petscerror.h
.  message - error message in the printf format
.  arg1 - argument (for example an integer, string or double)
-  arg2 - argument (for example an integer, string or double)

  Level: beginner

   Notes:
    Once the error handler is called the calling function is then returned from with the given error code.

   Experienced users can set the error handler with PetscPushErrorHandler().

   Concepts: error^setting condition

.seealso: PetscTraceBackErrorHandler(), PetscPushErrorHandler(), PetscError(), CHKERRQ(), CHKMEMQ, SETERRQ1(), SETERRQ3()
M*/
#define SETERRQ2(comm,ierr,s,a1,a2) return PetscError(comm,__LINE__,PETSC_FUNCTION_NAME,__FILE__,ierr,PETSC_ERROR_INITIAL,s,a1,a2)

/*MC
   SETERRQ3 - Macro that is called when an error has been detected,

   Synopsis:
   #include <petscsys.h>
   PetscErrorCode SETERRQ3(MPI_Comm comm,PetscErrorCode ierr,char *formatmessage,arg1,arg2,arg3)

   Collective on MPI_Comm

   Input Parameters:
+  comm - A communicator, so that the error can be collective
.  ierr - nonzero error code, see the list of standard error codes in include/petscerror.h
.  message - error message in the printf format
.  arg1 - argument (for example an integer, string or double)
.  arg2 - argument (for example an integer, string or double)
-  arg3 - argument (for example an integer, string or double)

  Level: beginner

   Notes:
    Once the error handler is called the calling function is then returned from with the given error code.

    There are also versions for 4, 5, 6 and 7 arguments.

   Experienced users can set the error handler with PetscPushErrorHandler().

   Concepts: error^setting condition

.seealso: PetscTraceBackErrorHandler(), PetscPushErrorHandler(), PetscError(), CHKERRQ(), CHKMEMQ, SETERRQ1(), SETERRQ2()
M*/
#define SETERRQ3(comm,ierr,s,a1,a2,a3) return PetscError(comm,__LINE__,PETSC_FUNCTION_NAME,__FILE__,ierr,PETSC_ERROR_INITIAL,s,a1,a2,a3)

/*MC
   SETERRQ4 - Macro that is called when an error has been detected,

   Synopsis:
   #include <petscsys.h>
   PetscErrorCode SETERRQ4(MPI_Comm comm,PetscErrorCode ierr,char *formatmessage,arg1,arg2,arg3)

   Collective on MPI_Comm

   Input Parameters:
+  comm - A communicator, so that the error can be collective
.  ierr - nonzero error code, see the list of standard error codes in include/petscerror.h
.  message - error message in the printf format
.  arg1 - argument (for example an integer, string or double)
.  arg2 - argument (for example an integer, string or double)
.  arg3 - argument (for example an integer, string or double)
-  arg4 - argument (for example an integer, string or double)

  Level: beginner

   Notes:
    Once the error handler is called the calling function is then returned from with the given error code.

    There are also versions for 4, 5, 6 and 7 arguments.

   Experienced users can set the error handler with PetscPushErrorHandler().

   Concepts: error^setting condition

.seealso: PetscTraceBackErrorHandler(), PetscPushErrorHandler(), PetscError(), CHKERRQ(), CHKMEMQ, SETERRQ1(), SETERRQ2()
M*/
#define SETERRQ4(comm,ierr,s,a1,a2,a3,a4) return PetscError(comm,__LINE__,PETSC_FUNCTION_NAME,__FILE__,ierr,PETSC_ERROR_INITIAL,s,a1,a2,a3,a4)

/*MC
   SETERRQ5 - Macro that is called when an error has been detected,

   Synopsis:
   #include <petscsys.h>
   PetscErrorCode SETERRQ5(MPI_Comm comm,PetscErrorCode ierr,char *formatmessage,arg1,arg2,arg3)

   Collective on MPI_COmm

   Input Parameters:
+  comm - A communicator, so that the error can be collective
.  ierr - nonzero error code, see the list of standard error codes in include/petscerror.h
.  message - error message in the printf format
.  arg1 - argument (for example an integer, string or double)
.  arg2 - argument (for example an integer, string or double)
.  arg3 - argument (for example an integer, string or double)
.  arg4 - argument (for example an integer, string or double)
-  arg5 - argument (for example an integer, string or double)

  Level: beginner

   Notes:
    Once the error handler is called the calling function is then returned from with the given error code.

    There are also versions for 4, 5, 6 and 7 arguments.

   Experienced users can set the error handler with PetscPushErrorHandler().

   Concepts: error^setting condition

.seealso: PetscTraceBackErrorHandler(), PetscPushErrorHandler(), PetscError(), CHKERRQ(), CHKMEMQ, SETERRQ1(), SETERRQ2()
M*/
#define SETERRQ5(comm,ierr,s,a1,a2,a3,a4,a5) return PetscError(comm,__LINE__,PETSC_FUNCTION_NAME,__FILE__,ierr,PETSC_ERROR_INITIAL,s,a1,a2,a3,a4,a5)

/*MC
   SETERRQ6 - Macro that is called when an error has been detected,

   Synopsis:
   #include <petscsys.h>
   PetscErrorCode SETERRQ6(MPI_Comm comm,PetscErrorCode ierr,char *formatmessage,arg1,arg2,arg3)

   Collective on MPI_Comm

   Input Parameters:
+  comm - A communicator, so that the error can be collective
.  ierr - nonzero error code, see the list of standard error codes in include/petscerror.h
.  message - error message in the printf format
.  arg1 - argument (for example an integer, string or double)
.  arg2 - argument (for example an integer, string or double)
.  arg3 - argument (for example an integer, string or double)
.  arg4 - argument (for example an integer, string or double)
.  arg5 - argument (for example an integer, string or double)
-  arg6 - argument (for example an integer, string or double)

  Level: beginner

   Notes:
    Once the error handler is called the calling function is then returned from with the given error code.

    There are also versions for 4, 5, 6 and 7 arguments.

   Experienced users can set the error handler with PetscPushErrorHandler().

   Concepts: error^setting condition

.seealso: PetscTraceBackErrorHandler(), PetscPushErrorHandler(), PetscError(), CHKERRQ(), CHKMEMQ, SETERRQ1(), SETERRQ2()
M*/
#define SETERRQ6(comm,ierr,s,a1,a2,a3,a4,a5,a6) return PetscError(comm,__LINE__,PETSC_FUNCTION_NAME,__FILE__,ierr,PETSC_ERROR_INITIAL,s,a1,a2,a3,a4,a5,a6)

/*MC
   SETERRQ7 - Macro that is called when an error has been detected,

   Synopsis:
   #include <petscsys.h>
   PetscErrorCode SETERRQ7(MPI_Comm comm,PetscErrorCode ierr,char *formatmessage,arg1,arg2,arg3)

   Collective on MPI_Comm

   Input Parameters:
+  comm - A communicator, so that the error can be collective
.  ierr - nonzero error code, see the list of standard error codes in include/petscerror.h
.  message - error message in the printf format
.  arg1 - argument (for example an integer, string or double)
.  arg2 - argument (for example an integer, string or double)
.  arg3 - argument (for example an integer, string or double)
.  arg4 - argument (for example an integer, string or double)
.  arg5 - argument (for example an integer, string or double)
.  arg6 - argument (for example an integer, string or double)
-  arg7 - argument (for example an integer, string or double)

  Level: beginner

   Notes:
    Once the error handler is called the calling function is then returned from with the given error code.

    There are also versions for 4, 5, 6 and 7 arguments.

   Experienced users can set the error handler with PetscPushErrorHandler().

   Concepts: error^setting condition

.seealso: PetscTraceBackErrorHandler(), PetscPushErrorHandler(), PetscError(), CHKERRQ(), CHKMEMQ, SETERRQ1(), SETERRQ2()
M*/
#define SETERRQ7(comm,ierr,s,a1,a2,a3,a4,a5,a6,a7) return PetscError(comm,__LINE__,PETSC_FUNCTION_NAME,__FILE__,ierr,PETSC_ERROR_INITIAL,s,a1,a2,a3,a4,a5,a6,a7)

/*MC
   SETERRQ8 - Macro that is called when an error has been detected,

   Synopsis:
   #include <petscsys.h>
   PetscErrorCode SETERRQ8(MPI_Comm comm,PetscErrorCode ierr,char *formatmessage,arg1,arg2,arg3)

   Collective on MPI_Comm

   Input Parameters:
+  comm - A communicator, so that the error can be collective
.  ierr - nonzero error code, see the list of standard error codes in include/petscerror.h
.  message - error message in the printf format
.  arg1 - argument (for example an integer, string or double)
.  arg2 - argument (for example an integer, string or double)
.  arg3 - argument (for example an integer, string or double)
.  arg4 - argument (for example an integer, string or double)
.  arg5 - argument (for example an integer, string or double)
.  arg6 - argument (for example an integer, string or double)
.  arg7 - argument (for example an integer, string or double)
-  arg8 - argument (for example an integer, string or double)

  Level: beginner

   Notes:
    Once the error handler is called the calling function is then returned from with the given error code.

    There are also versions for 4, 5, 6 and 7 arguments.

   Experienced users can set the error handler with PetscPushErrorHandler().

   Concepts: error^setting condition

.seealso: PetscTraceBackErrorHandler(), PetscPushErrorHandler(), PetscError(), CHKERRQ(), CHKMEMQ, SETERRQ1(), SETERRQ2()
M*/
#define SETERRQ8(comm,ierr,s,a1,a2,a3,a4,a5,a6,a7,a8) return PetscError(comm,__LINE__,PETSC_FUNCTION_NAME,__FILE__,ierr,PETSC_ERROR_INITIAL,s,a1,a2,a3,a4,a5,a6,a7,a8)

/*MC
   SETERRABORT - Macro that can be called when an error has been detected,

   Synopsis:
   #include <petscsys.h>
   PetscErrorCode SETERRABORT(MPI_Comm comm,PetscErrorCode ierr,char *message)

   Collective on MPI_Comm

   Input Parameters:
+  comm - A communicator, so that the error can be collective
.  ierr - nonzero error code, see the list of standard error codes in include/petscerror.h
-  message - error message in the printf format

  Level: beginner

   Notes:
    This function just calls MPI_Abort().

   Concepts: error^setting condition

.seealso: PetscTraceBackErrorHandler(), PetscPushErrorHandler(), PetscError(), CHKERRQ(), CHKMEMQ, SETERRQ1(), SETERRQ2()
M*/
#define SETERRABORT(comm,ierr,s) do {PetscError(comm,__LINE__,PETSC_FUNCTION_NAME,__FILE__,ierr,PETSC_ERROR_INITIAL,s);MPI_Abort(comm,ierr);} while (0)

/*MC
   CHKERRQ - Checks error code, if non-zero it calls the error handler and then returns

   Synopsis:
   #include <petscsys.h>
   PetscErrorCode CHKERRQ(PetscErrorCode ierr)

   Not Collective

   Input Parameters:
.  ierr - nonzero error code, see the list of standard error codes in include/petscerror.h

  Level: beginner

   Notes:
    Once the error handler is called the calling function is then returned from with the given error code.

    Experienced users can set the error handler with PetscPushErrorHandler().

    CHKERRQ(ierr) is fundamentally a macro replacement for
         if (ierr) return(PetscError(...,ierr,...));

    Although typical usage resembles "void CHKERRQ(PetscErrorCode)" as described above, for certain uses it is
    highly inappropriate to use it in this manner as it invokes return(PetscErrorCode). In particular,
    it cannot be used in functions which return(void) or any other datatype.  In these types of functions,
    you can use CHKERRV() which returns without an error code (bad idea since the error is ignored or
         if (ierr) {PetscError(....); return(YourReturnType);}
    where you may pass back a NULL to indicate an error. You can also call CHKERRABORT(comm,n) to have
    MPI_Abort() returned immediately.

   Concepts: error^setting condition

.seealso: PetscTraceBackErrorHandler(), PetscPushErrorHandler(), PetscError(), SETERRQ(), CHKMEMQ, SETERRQ1(), SETERRQ2(), SETERRQ2()
M*/
#define CHKERRQ(ierr)          do {if (PetscUnlikely(ierr)) return PetscError(PETSC_COMM_SELF,__LINE__,PETSC_FUNCTION_NAME,__FILE__,ierr,PETSC_ERROR_REPEAT," ");} while (0)
#define CHKERRV(ierr)          do {if (PetscUnlikely(ierr)) {ierr = PetscError(PETSC_COMM_SELF,__LINE__,PETSC_FUNCTION_NAME,__FILE__,ierr,PETSC_ERROR_REPEAT," ");return;}} while(0)
#define CHKERRABORT(comm,ierr) do {if (PetscUnlikely(ierr)) {PetscError(PETSC_COMM_SELF,__LINE__,PETSC_FUNCTION_NAME,__FILE__,ierr,PETSC_ERROR_REPEAT," ");MPI_Abort(comm,ierr);}} while (0)
#define CHKERRCONTINUE(ierr)   do {if (PetscUnlikely(ierr)) {PetscError(PETSC_COMM_SELF,__LINE__,PETSC_FUNCTION_NAME,__FILE__,ierr,PETSC_ERROR_REPEAT," ");}} while (0)


/*MC
   CHKERRMPI - Checks error code, if non-zero it calls the error handler and then returns

   Synopsis:
   #include <petscsys.h>
   PetscErrorCode CHKERRMPI(PetscErrorCode ierr)

   Not Collective

   Input Parameters:
.  ierr - nonzero error code, see the list of standard error codes in include/petscerror.h

  Level: developer

   Notes:
    This macro is FOR USE IN MPI CALLBACK FUNCTIONS ONLY, such as those passed to MPI_Comm_create_keyval(). It always returns the error code PETSC_MPI_ERROR_CODE
    which is registered with MPI_Add_error_code() when PETSc is initialized.

   Concepts: error^setting condition

.seealso: CHKERRQ(), PetscTraceBackErrorHandler(), PetscPushErrorHandler(), PetscError(), SETERRQ(), CHKMEMQ, SETERRQ1(), SETERRQ2(), SETERRQ2()
M*/
#define CHKERRMPI(ierr)        do {if (PetscUnlikely(ierr)) return (PetscError(PETSC_COMM_SELF,__LINE__,PETSC_FUNCTION_NAME,__FILE__,ierr,PETSC_ERROR_REPEAT," "),PETSC_MPI_ERROR_CODE);} while (0)

#ifdef PETSC_CLANGUAGE_CXX

/*MC
   CHKERRXX - Checks error code, if non-zero it calls the C++ error handler which throws an exception

   Synopsis:
   #include <petscsys.h>
   void CHKERRXX(PetscErrorCode ierr)

   Not Collective

   Input Parameters:
.  ierr - nonzero error code, see the list of standard error codes in include/petscerror.h

  Level: beginner

   Notes:
    Once the error handler throws a ??? exception.

    You can use CHKERRV() which returns without an error code (bad idea since the error is ignored)
    or CHKERRABORT(comm,n) to have MPI_Abort() returned immediately.

   Concepts: error^setting condition

.seealso: PetscTraceBackErrorHandler(), PetscPushErrorHandler(), PetscError(), SETERRQ(), CHKERRQ(), CHKMEMQ
M*/
#define CHKERRXX(ierr)  do {if (PetscUnlikely(ierr)) {PetscError(PETSC_COMM_SELF,__LINE__,PETSC_FUNCTION_NAME,__FILE__,ierr,PETSC_ERROR_IN_CXX,0);}} while(0)

#endif

#define CHKERRCUDA(err)   do {if (PetscUnlikely(err)) SETERRQ1(PETSC_COMM_SELF,PETSC_ERR_LIB,"CUDA error %d",err);} while(0)
#define CHKERRCUBLAS(err) do {if (PetscUnlikely(err)) SETERRQ1(PETSC_COMM_SELF,PETSC_ERR_LIB,"CUBLAS error %d",err);} while(0)

/*MC
   CHKMEMQ - Checks the memory for corruption, calls error handler if any is detected

   Synopsis:
   #include <petscsys.h>
   CHKMEMQ;

   Not Collective

  Level: beginner

   Notes:
    We highly recommend using valgrind http://www.mcs.anl.gov/petsc/documentation/faq.html#valgrind for finding memory problems. This is useful
    on systems that do not have valgrind, but much much less useful.

    Must run with the option -malloc_debug to enable this option

    Once the error handler is called the calling function is then returned from with the given error code.

    By defaults prints location where memory that is corrupted was allocated.

    Use CHKMEMA for functions that return void

   Concepts: memory corruption

.seealso: PetscTraceBackErrorHandler(), PetscPushErrorHandler(), PetscError(), SETERRQ(), CHKMEMQ, SETERRQ1(), SETERRQ2(), SETERRQ3(),
          PetscMallocValidate()
M*/
#define CHKMEMQ do {PetscErrorCode _7_ierr = PetscMallocValidate(__LINE__,PETSC_FUNCTION_NAME,__FILE__);CHKERRQ(_7_ierr);} while(0)

#define CHKMEMA PetscMallocValidate(__LINE__,PETSC_FUNCTION_NAME,__FILE__)

/*E
  PetscErrorType - passed to the PETSc error handling routines indicating if this is the first or a later call to the error handlers

  Level: advanced

  PETSC_ERROR_IN_CXX indicates the error was detected in C++ and an exception should be generated

  Developer Notes:
    This is currently used to decide when to print the detailed information about the run in PetscTraceBackErrorHandler()

.seealso: PetscError(), SETERRXX()
E*/
typedef enum {PETSC_ERROR_INITIAL=0,PETSC_ERROR_REPEAT=1,PETSC_ERROR_IN_CXX = 2} PetscErrorType;

#if defined(__clang_analyzer__)
__attribute__((analyzer_noreturn))
#endif
PETSC_EXTERN PetscErrorCode PetscError(MPI_Comm,int,const char*,const char*,PetscErrorCode,PetscErrorType,const char*,...);

PETSC_EXTERN PetscErrorCode PetscErrorPrintfInitialize(void);
PETSC_EXTERN PetscErrorCode PetscErrorMessage(int,const char*[],char **);
PETSC_EXTERN PetscErrorCode PetscTraceBackErrorHandler(MPI_Comm,int,const char*,const char*,PetscErrorCode,PetscErrorType,const char*,void*);
PETSC_EXTERN PetscErrorCode PetscIgnoreErrorHandler(MPI_Comm,int,const char*,const char*,PetscErrorCode,PetscErrorType,const char*,void*);
PETSC_EXTERN PetscErrorCode PetscEmacsClientErrorHandler(MPI_Comm,int,const char*,const char*,PetscErrorCode,PetscErrorType,const char*,void*);
PETSC_EXTERN PetscErrorCode PetscMPIAbortErrorHandler(MPI_Comm,int,const char*,const char*,PetscErrorCode,PetscErrorType,const char*,void*);
PETSC_EXTERN PetscErrorCode PetscAbortErrorHandler(MPI_Comm,int,const char*,const char*,PetscErrorCode,PetscErrorType,const char*,void*);
PETSC_EXTERN PetscErrorCode PetscAttachDebuggerErrorHandler(MPI_Comm,int,const char*,const char*,PetscErrorCode,PetscErrorType,const char*,void*);
PETSC_EXTERN PetscErrorCode PetscReturnErrorHandler(MPI_Comm,int,const char*,const char*,PetscErrorCode,PetscErrorType,const char*,void*);
PETSC_EXTERN PetscErrorCode PetscPushErrorHandler(PetscErrorCode (*handler)(MPI_Comm,int,const char*,const char*,PetscErrorCode,PetscErrorType,const char*,void*),void*);
PETSC_EXTERN PetscErrorCode PetscPopErrorHandler(void);
PETSC_EXTERN PetscErrorCode PetscSignalHandlerDefault(int,void*);
PETSC_EXTERN PetscErrorCode PetscPushSignalHandler(PetscErrorCode (*)(int,void *),void*);
PETSC_EXTERN PetscErrorCode PetscPopSignalHandler(void);
PETSC_EXTERN PetscErrorCode PetscCheckPointerSetIntensity(PetscInt);

/*MC
    PetscErrorPrintf - Prints error messages.

   Synopsis:
    #include <petscsys.h>
     PetscErrorCode (*PetscErrorPrintf)(const char format[],...);

    Not Collective

    Input Parameters:
.   format - the usual printf() format string

   Options Database Keys:
+    -error_output_stdout - cause error messages to be printed to stdout instead of the  (default) stderr
-    -error_output_none - to turn off all printing of error messages (does not change the way the error is handled.)

   Notes:
    Use
$     PetscErrorPrintf = PetscErrorPrintfNone; to turn off all printing of error messages (does not change the way the
$                        error is handled.) and
$     PetscErrorPrintf = PetscErrorPrintfDefault; to turn it back on or you can use your own function

          Use
     PETSC_STDERR = FILE* obtained from a file open etc. to have stderr printed to the file.
     PETSC_STDOUT = FILE* obtained from a file open etc. to have stdout printed to the file.

          Use
      PetscPushErrorHandler() to provide your own error handler that determines what kind of messages to print

   Level: developer

    Fortran Note:
    This routine is not supported in Fortran.

    Concepts: error messages^printing
    Concepts: printing^error messages

.seealso: PetscFPrintf(), PetscSynchronizedPrintf(), PetscHelpPrintf(), PetscPrintf(), PetscPushErrorHandler(), PetscVFPrintf(), PetscHelpPrintf()
M*/
PETSC_EXTERN PetscErrorCode (*PetscErrorPrintf)(const char[],...);

typedef enum {PETSC_FP_TRAP_OFF=0,PETSC_FP_TRAP_ON=1} PetscFPTrap;
PETSC_EXTERN PetscErrorCode PetscSetFPTrap(PetscFPTrap);
PETSC_EXTERN PetscErrorCode PetscFPTrapPush(PetscFPTrap);
PETSC_EXTERN PetscErrorCode PetscFPTrapPop(void);

/*
      Allows the code to build a stack frame as it runs
*/

#define PETSCSTACKSIZE 64

typedef struct  {
  const char      *function[PETSCSTACKSIZE];
  const char      *file[PETSCSTACKSIZE];
        int       line[PETSCSTACKSIZE];
        PetscBool petscroutine[PETSCSTACKSIZE];
        int       currentsize;
        int       hotdepth;
} PetscStack;

PETSC_EXTERN PetscStack *petscstack;

PetscErrorCode  PetscStackCopy(PetscStack*,PetscStack*);
PetscErrorCode  PetscStackPrint(PetscStack *,FILE*);
#if defined(PETSC_SERIALIZE_FUNCTIONS)
#include <petsc/private/petscfptimpl.h>
/*
   Registers the current function into the global function pointer to function name table

   Have to fix this to handle errors but cannot return error since used in PETSC_VIEWER_DRAW_() etc
*/
#define PetscRegister__FUNCT__() do { \
  static PetscBool __chked = PETSC_FALSE; \
  if (!__chked) {\
  void *ptr; PetscDLSym(NULL,PETSC_FUNCTION_NAME,&ptr);\
  __chked = PETSC_TRUE;\
  }} while (0)
#else
#define PetscRegister__FUNCT__()
#endif

#if defined(PETSC_USE_DEBUG)
PETSC_STATIC_INLINE PetscBool PetscStackActive(void)
{
  return(petscstack ? PETSC_TRUE : PETSC_FALSE);
}

/* Stack handling is based on the following two "NoCheck" macros.  These should only be called directly by other error
 * handling macros.  We record the line of the call, which may or may not be the location of the definition.  But is at
 * least more useful than "unknown" because it can distinguish multiple calls from the same function.
 */

#define PetscStackPushNoCheck(funct,petsc_routine,hot)                        \
  do {                                                                        \
    PetscStackSAWsTakeAccess();                                                \
    if (petscstack && (petscstack->currentsize < PETSCSTACKSIZE)) {         \
      petscstack->function[petscstack->currentsize]  = funct;               \
      petscstack->file[petscstack->currentsize]      = __FILE__;            \
      petscstack->line[petscstack->currentsize]      = __LINE__;            \
      petscstack->petscroutine[petscstack->currentsize] = petsc_routine;    \
      petscstack->currentsize++;                                             \
    }                                                                         \
    if (petscstack) {                                                        \
      petscstack->hotdepth += (hot || petscstack->hotdepth);                \
    }                                                                         \
    PetscStackSAWsGrantAccess();                                               \
  } while (0)

#define PetscStackPopNoCheck                                            \
  do {                                                                  \
    PetscStackSAWsTakeAccess();                                          \
    if (petscstack && petscstack->currentsize > 0) {                  \
      petscstack->currentsize--;                                       \
      petscstack->function[petscstack->currentsize]  = 0;             \
      petscstack->file[petscstack->currentsize]      = 0;             \
      petscstack->line[petscstack->currentsize]      = 0;             \
      petscstack->petscroutine[petscstack->currentsize] = PETSC_FALSE;\
    }                                                                   \
    if (petscstack) {                                                  \
      petscstack->hotdepth = PetscMax(petscstack->hotdepth-1,0);      \
    }                                                                   \
    PetscStackSAWsGrantAccess();                                         \
  } while (0)

/*MC
   PetscFunctionBegin - First executable line of each PETSc function,  used for error handling. Final
      line of PETSc functions should be PetscFunctionReturn(0);

   Synopsis:
   #include <petscsys.h>
   void PetscFunctionBegin;

   Not Collective

   Usage:
.vb
     int something;

     PetscFunctionBegin;
.ve

   Notes:
     Use PetscFunctionBeginUser for application codes.

     Not available in Fortran

   Level: developer

.seealso: PetscFunctionReturn(), PetscFunctionBeginHot(), PetscFunctionBeginUser()

.keywords: traceback, error handling
M*/
#define PetscFunctionBegin do {                                        \
    PetscStackPushNoCheck(PETSC_FUNCTION_NAME,PETSC_TRUE,PETSC_FALSE); \
    PetscRegister__FUNCT__();                                          \
  } while (0)

/*MC
   PetscFunctionBeginHot - Substitute for PetscFunctionBegin to be used in functions that are called in
   performance-critical circumstances.  Use of this function allows for lighter profiling by default.

   Synopsis:
   #include <petscsys.h>
   void PetscFunctionBeginHot;

   Not Collective

   Usage:
.vb
     int something;

     PetscFunctionBeginHot;
.ve

   Notes:
     Not available in Fortran

   Level: developer

.seealso: PetscFunctionBegin, PetscFunctionReturn()

.keywords: traceback, error handling
M*/
#define PetscFunctionBeginHot do {                                     \
    PetscStackPushNoCheck(PETSC_FUNCTION_NAME,PETSC_TRUE,PETSC_TRUE);  \
    PetscRegister__FUNCT__();                                          \
  } while (0)

/*MC
   PetscFunctionBeginUser - First executable line of user provided PETSc routine

   Synopsis:
   #include <petscsys.h>
   void PetscFunctionBeginUser;

   Not Collective

   Usage:
.vb
     int something;

     PetscFunctionBeginUser;
.ve

   Notes:
      Final line of PETSc functions should be PetscFunctionReturn(0) except for main().

      Not available in Fortran

      This is identical to PetscFunctionBegin except it labels the routine as a user
      routine instead of as a PETSc library routine.

   Level: intermediate

.seealso: PetscFunctionReturn(), PetscFunctionBegin, PetscFunctionBeginHot

.keywords: traceback, error handling
M*/
#define PetscFunctionBeginUser                                          \
  do {                                                                  \
    PetscStackPushNoCheck(PETSC_FUNCTION_NAME,PETSC_FALSE,PETSC_FALSE); \
    PetscRegister__FUNCT__();                                           \
  } while (0)


#define PetscStackPush(n) \
  do {                                                                  \
    PetscStackPushNoCheck(n,PETSC_FALSE,PETSC_FALSE);                   \
    CHKMEMQ;                                                            \
  } while (0)

#define PetscStackPop                           \
    do {                                        \
      CHKMEMQ;                                  \
      PetscStackPopNoCheck;                     \
    } while (0)

/*MC
   PetscFunctionReturn - Last executable line of each PETSc function
        used for error handling. Replaces return()

   Synopsis:
   #include <petscsys.h>
   void PetscFunctionReturn(0);

   Not Collective

   Usage:
.vb
    ....
     PetscFunctionReturn(0);
   }
.ve

   Notes:
     Not available in Fortran

   Level: developer

.seealso: PetscFunctionBegin()

.keywords: traceback, error handling
M*/
#define PetscFunctionReturn(a) \
  do {                                                                \
    PetscStackPopNoCheck;                                             \
    return(a);} while (0)

#define PetscFunctionReturnVoid() \
  do {                                                                \
    PetscStackPopNoCheck;                                             \
    return;} while (0)

#else

PETSC_STATIC_INLINE PetscBool PetscStackActive(void) {return PETSC_FALSE;}
#define PetscStackPushNoCheck(funct,petsc_routine,hot) do {} while (0)
#define PetscStackPopNoCheck                           do {} while (0)
#define PetscFunctionBegin
#define PetscFunctionBeginUser
#define PetscFunctionBeginHot
#define PetscFunctionReturn(a)    return(a)
#define PetscFunctionReturnVoid() return
#define PetscStackPop             CHKMEMQ
#define PetscStackPush(f)         CHKMEMQ

#endif

/*
    PetscStackCall - Calls an external library routine or user function after pushing the name of the routine on the stack.

   Input Parameters:
+   name - string that gives the name of the function being called
-   routine - actual call to the routine, including ierr = and CHKERRQ(ierr);

   Note: Often one should use PetscStackCallStandard() instead. This routine is intended for external library routines that DO NOT return error codes

   Developer Note: this is so that when a user or external library routine results in a crash or corrupts memory, they get blamed instead of PETSc.



*/
#define PetscStackCall(name,routine) do { PetscStackPush(name);routine;PetscStackPop; } while(0)

/*
    PetscStackCallStandard - Calls an external library routine after pushing the name of the routine on the stack.

   Input Parameters:
+   func-  name of the routine
-   args - arguments to the routine surrounded by ()

   Notes:
    This is intended for external package routines that return error codes. Use PetscStackCall() for those that do not.

   Developer Note: this is so that when an external packge routine results in a crash or corrupts memory, they get blamed instead of PETSc.

*/
#define PetscStackCallStandard(func,args) do {                        \
    PetscStackPush(#func);ierr = func args;PetscStackPop; if (ierr) SETERRQ2(PETSC_COMM_SELF,PETSC_ERR_LIB,"Error in %s(): error code %d",#func,(int)ierr); \
  } while (0)

PETSC_EXTERN PetscErrorCode PetscStackCreate(void);
PETSC_EXTERN PetscErrorCode PetscStackView(FILE*);
PETSC_EXTERN PetscErrorCode PetscStackDestroy(void);

#endif
