//===----------------------------------------------------------------------===//
//
//                                   xara
//
//===----------------------------------------------------------------------===//
//                              https://xara.so
//===----------------------------------------------------------------------===//
// These commands expect a Domain* as their clientData.
// 
// Written: cmp
//
#include <assert.h>
#include <Domain.h>
#include <tcl.h>
#include <Parsing.h>
#include <Logging.h>

int
TclCommand_setLoadConst(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char ** const argv)
{
  assert(clientData != nullptr);
  Domain* domain = (Domain*)clientData;
  
  domain->setLoadConstant();
  if (argc == 3) {
    if (strcmp(argv[1], "-time") == 0) {
      double newTime;
      if (Tcl_GetDouble(interp, argv[2], &newTime) != TCL_OK) {
        opserr << "WARNING readingvalue - loadConst -time value \n";
        return TCL_ERROR;
      } else {
        domain->setCurrentTime(newTime);
        domain->setCommittedTime(newTime);
      }
    }
  }
  return TCL_OK;
}

int
TclCommand_setCreep(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char ** const argv)
{
  assert(clientData != nullptr);
  Domain* the_domain = (Domain*)clientData;

  if (argc < 2) {
    opserr << "WARNING illegal command - setCreep value? \n";
    return TCL_ERROR;
  }
  int newFlag;
  if (Tcl_GetInt(interp, argv[1], &newFlag) != TCL_OK) {
    opserr << "WARNING reading creep value - setCreep newFlag? \n";
    return TCL_ERROR;
  } else {
    the_domain->setCreep(newFlag);
  }
  return TCL_OK;
}

int
TclCommand_setTime(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char ** const argv)
{
  assert(clientData != nullptr);
  Domain* domain = (Domain*)clientData;

  if (argc < 2) {
    opserr << "WARNING illegal command - time pseudoTime? \n";
    return TCL_ERROR;
  }
  double newTime;
  if (Tcl_GetDouble(interp, argv[1], &newTime) != TCL_OK) {
    opserr << "WARNING reading time value - time pseudoTime? \n";
    return TCL_ERROR;
  } else {
    domain->setCurrentTime(newTime);
    domain->setCommittedTime(newTime);
  }
  return TCL_OK;
}

int
TclCommand_getTime(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char ** const argv)
{
  assert(clientData != nullptr);
  Domain* domain = static_cast<Domain*>(clientData);

  Tcl_SetObjResult(interp, Tcl_NewDoubleObj(domain->getCurrentTime()));
  return TCL_OK;
}

