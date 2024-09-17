//===----------------------------------------------------------------------===//
//
//        OpenSees - Open System for Earthquake Engineering Simulation
//
//===----------------------------------------------------------------------===//
//
//
// Description: This file contains the function invoked when the user invokes
// the MySQL command in the interpreter.
//
// Written: fmk
//
#include <OPS_Globals.h>
#include <stdlib.h>
#include <string.h>
#include <tcl.h>

#include <BerkeleyDbDatastore.h>

#ifdef _USRDLL
#define DllExport _declspec(dllexport)
#else
#define DllExport
#endif

extern "C" DllExport int
TclCommand_BerkeleyDB(ClientData clientData, Tcl_Interp *interp, int argc,
                      TCL_Char ** const argv, Domain *theDomain,
                      FEM_ObjectBroker *theBroker, FE_Datastore **theDatabase)
{

  // delete the old database
  if (*theDatabase != 0)
    delete (*theDatabase);

  (*theDatabase) = new BerkeleyDbDatastore(argv[2], *theDomain, *theBroker);

  if (*theDatabase == 0) {
    opserr << "WARNING database MySql dabaseName? - out of memory\n";
    return TCL_ERROR;
  }

  return TCL_OK;
}

