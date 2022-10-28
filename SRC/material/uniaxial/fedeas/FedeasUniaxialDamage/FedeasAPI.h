#pragma once
// #include "ArgumentFFI.h"
#include <tcl.h>

#ifdef __cplusplus
  extern "C" {
#endif

#ifdef ISW_COMMIT
# undef ISW_COMMIT
#endif
#ifdef ISW_DELETE
# undef ISW_DELETE
#endif

enum ISW_Task {
  ISW_COMMIT = 1<<0,
  ISW_DELETE = 1<<1,
  ISW_UPDATE = 1<<2,
  ISW_MALLOC = 1<<3,
  ISW_CREATE = 1<<4,
  ISW_COPY   = 1<<5,

  ISW_RETURN_RESIDUAL = 1<<6,
  ISW_RETURN_TANGENT  = 1<<7,
  ISW_RETURN_OTHER    = 1<<8,
  
  ISW_UPDATE_InitialTangent = 1<<9,
  ISW_UPDATE_CurrentTangent = 1<<10
};

const int ISW_ACTION = ISW_COMMIT 
                     | ISW_UPDATE 
                     | ISW_MALLOC 
                     | ISW_DELETE 
                     | ISW_CREATE
                     | ISW_COPY;

const int ISW_MODIFY = ISW_UPDATE_InitialTangent
                     | ISW_UPDATE_CurrentTangent;

typedef int StateRoutine(
    struct StateOperator* routine,
    Tcl_Interp* interp,
    const int action,
    const int argc, const char** argv,
    /* Double input data */
    const int argi, const double input[],
    /* Double output data, optional int flags */
    const int argo, double outputs[], int intv[]
);

typedef int StateObjectRoutine(
    struct StateOperator* routine,
    Tcl_Interp* interp,
    const int action,
    const int argc, Tcl_Obj* const* argv,
    /* Double input data */
    const int argi, const double input[],
    /* Double output data, optional int flags */
    const int argo, double outputs[], int intv[]
);

struct StateOperator {
  StateRoutine* call;
  StateObjectRoutine* callobj;

  void   *runtime;

  void   *data,
         *instance;
  size_t data_size;
};

#ifdef __cplusplus
}
#endif
