//===----------------------------------------------------------------------===//
//
//        OpenSees - Open System for Earthquake Engineering Simulation
//
//===----------------------------------------------------------------------===//
#include <tcl.h>
#include <OPS_Globals.h>
#include <Channel.h>
#include <MachineBroker.h>

static int getPID(ClientData,  Tcl_Interp *, int, TCL_Char ** const argv);
static int getNP( ClientData,  Tcl_Interp *, int, TCL_Char ** const argv);


void Init_MachineRuntime(Tcl_Interp* interp, MachineBroker* theMachineBroker)
{
  Tcl_CreateCommand(interp, "getNP",     &getNP,   (ClientData)theMachineBroker, (Tcl_CmdDeleteProc *)NULL);
  Tcl_CreateCommand(interp, "getPID",    &getPID,  (ClientData)theMachineBroker, (Tcl_CmdDeleteProc *)NULL);
}

static int
getPID(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char ** const argv)
{
  int pid = 0;

  MachineBroker* theMachineBroker = (MachineBroker*)clientData;

  if (theMachineBroker != nullptr)
    pid = theMachineBroker->getPID();

  // now we copy the value to the tcl string that is returned
  Tcl_SetObjResult(interp, Tcl_NewIntObj(pid));

  return TCL_OK;
}

static int
getNP(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char ** const argv)
{
  int np = 1;
  MachineBroker* theMachineBroker = (MachineBroker*)clientData;

  if (theMachineBroker != nullptr)
    np = theMachineBroker->getNP();

  // now we copy the value to the tcl string that is returned
  Tcl_SetObjResult(interp, Tcl_NewIntObj(np));

  return TCL_OK;
}


