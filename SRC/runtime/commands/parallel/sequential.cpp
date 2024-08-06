//===----------------------------------------------------------------------===//
//
//        OpenSees - Open System for Earthquake Engineering Simulation
//
//===----------------------------------------------------------------------===//
//
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


#if 0
int
opsSendSequential(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char ** const argv)
{
  if (argc < 2)
    return TCL_OK;

  int otherPID = -1;
  int myPID = theMachineBroker->getPID();
  int np = theMachineBroker->getNP();
  const char *dataToSend = argv[argc - 1];
  int msgLength = strlen(dataToSend) + 1;

  const char *gMsg = dataToSend;
  //  strcpy(gMsg, dataToSend);

  if (strcmp(argv[1], "-pid") == 0 && argc > 3) {

    if (Tcl_GetInt(interp, argv[2], &otherPID) != TCL_OK) {
      opserr << "send -pid pid? data? - pid: " << argv[2] << " invalid\n";
      return TCL_ERROR;
    }

    if (otherPID > -1 && otherPID != myPID && otherPID < np) {

      MPI_Send((void *)(&msgLength), 1, MPI_INT, otherPID, 0, MPI_COMM_WORLD);
      MPI_Send((void *)gMsg, msgLength, MPI_CHAR, otherPID, 1, MPI_COMM_WORLD);

    } else {
      opserr << "send -pid pid? data? - pid: " << otherPID << " invalid\n";
      return TCL_ERROR;
    }

  } else {
    if (myPID == 0) {
      MPI_Bcast((void *)(&msgLength), 1, MPI_INT, 0, MPI_COMM_WORLD);
      MPI_Bcast((void *)gMsg, msgLength, MPI_CHAR, 0, MPI_COMM_WORLD);
    } else {
      opserr << "send data - only process 0 can do a broadcast - you may need "
                "to kill the application";
      return TCL_ERROR;
    }
  }
  return TCL_OK;
}

#endif
