//===----------------------------------------------------------------------===//
//
//                                   xara
//
//===----------------------------------------------------------------------===//
//                              https://xara.so
//===----------------------------------------------------------------------===//
// Description: This unit contains implementations of parallel commands
// for use in non-parallel interpreters
//
#include <tcl.h>
#define TCL_Char CONST84 char

Tcl_CmdProc getPIDSequential;
Tcl_CmdProc getNPSequential;
Tcl_CmdProc opsBarrierSequential;
Tcl_CmdProc opsSendSequential;
Tcl_CmdProc opsRecvSequential;
Tcl_CmdProc opsPartitionSequential;

void G3_InitTclSequentialAPI(Tcl_Interp* interp)
{
  Tcl_CreateCommand(interp, "getNP",     &getNPSequential, (ClientData)NULL, (Tcl_CmdDeleteProc *)NULL);
  Tcl_CreateCommand(interp, "getPID",    &getPIDSequential, (ClientData)NULL, (Tcl_CmdDeleteProc *)NULL);
  Tcl_CreateCommand(interp, "barrier",   &opsBarrierSequential, (ClientData)NULL, (Tcl_CmdDeleteProc *)NULL);
//Tcl_CreateCommand(interp, "send",      &opsSendSequential, (ClientData)NULL, (Tcl_CmdDeleteProc *)NULL);
  Tcl_CreateCommand(interp, "recv",      &opsRecvSequential, (ClientData)NULL, (Tcl_CmdDeleteProc *)NULL);
  Tcl_CreateCommand(interp, "partition", &opsPartitionSequential, (ClientData)NULL, (Tcl_CmdDeleteProc *)NULL);
}


int
getPIDSequential(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char ** const argv)
{
  int pid = 0;
  char buffer[30];
  sprintf(buffer, "%d", pid);
  Tcl_SetResult(interp, buffer, TCL_VOLATILE);

  return TCL_OK;
}

int
getNPSequential(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char ** const argv)
{
  int np = 1;
  char buffer[30];
  sprintf(buffer, "%d", np);
  Tcl_SetResult(interp, buffer, TCL_VOLATILE);
  return TCL_OK;
}

int
opsPartitionSequential(ClientData clientData, Tcl_Interp *interp, int argc,
             TCL_Char ** const argv)
{
  return TCL_OK;
}


int
opsBarrierSequential(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char ** const argv)
{
  return TCL_OK;
}

int
opsRecvSequential(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char ** const argv)
{

  return TCL_OK;
}
