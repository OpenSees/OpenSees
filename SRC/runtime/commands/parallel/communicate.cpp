
#include <tcl.h>
#include <mpi.h>
#include <string.h>
#include <Logging.h>
#include <Parsing.h>
#include <MachineBroker.h>


static int opsBarrier(ClientData, Tcl_Interp *, int, TCL_Char ** const argv);
static int opsSend(ClientData, Tcl_Interp *, int, TCL_Char ** const argv);
static int opsRecv(ClientData, Tcl_Interp *, int,TCL_Char ** const argv);

void Init_Communication(Tcl_Interp* interp, MachineBroker* theMachineBroker)
{
  Tcl_CreateCommand(interp, "send",      &opsSend, (ClientData)theMachineBroker, (Tcl_CmdDeleteProc *)NULL);
  Tcl_CreateCommand(interp, "recv",      &opsRecv, (ClientData)theMachineBroker, (Tcl_CmdDeleteProc *)NULL);
  Tcl_CreateCommand(interp, "barrier",   &opsBarrier, (ClientData)NULL, (Tcl_CmdDeleteProc *)NULL);
}


static int
opsSend(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char ** const argv)
{
  if (argc < 2)
    return TCL_OK;

  int otherPID = -1;
  MachineBroker* theMachineBroker = (MachineBroker*)clientData;
  int myPID = theMachineBroker->getPID();
  int np    = theMachineBroker->getNP();
  const char *dataToSend = argv[argc - 1];
  int msgLength = strlen(dataToSend) + 1;

  const char *gMsg = dataToSend;

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

static int
opsRecv(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char ** const argv)
{
  MachineBroker* theMachineBroker = (MachineBroker*)clientData;
  if (argc < 2)
    return TCL_OK;

  int otherPID = 0;
  int myPID = theMachineBroker->getPID();
  int np = theMachineBroker->getNP();
  TCL_Char *varToSet = argv[argc - 1];

  int msgLength = 0;
  char *gMsg = 0;

  if (strcmp(argv[1], "-pid") == 0 && argc > 3) {

    bool fromAny = false;

    if ((strcmp(argv[2], "ANY") == 0) || (strcmp(argv[2], "ANY_SOURCE") == 0) ||
        (strcmp(argv[2], "MPI_ANY_SOURCE") == 0)) {
      fromAny = true;
    } else {
      if (Tcl_GetInt(interp, argv[2], &otherPID) != TCL_OK) {
        opserr << "recv -pid pid? data? - pid: " << argv[2] << " invalid\n";
        return TCL_ERROR;
      }
    }

    if (otherPID > -1 && otherPID < np) {
      MPI_Status status;

      if (fromAny == false)
        if (myPID != otherPID)
          MPI_Recv((void *)(&msgLength), 1, MPI_INT, otherPID, 0,
                   MPI_COMM_WORLD, &status);
        else {
          opserr << "recv -pid pid? data? - " << otherPID
                 << " cant receive from self!\n";
          return TCL_ERROR;
        }
      else {
        MPI_Recv((void *)(&msgLength), 1, MPI_INT, MPI_ANY_SOURCE, 0,
                 MPI_COMM_WORLD, &status);
        otherPID = status.MPI_SOURCE;
      }

      if (msgLength > 0) {
        gMsg = new char[msgLength];

        MPI_Recv((void *)gMsg, msgLength, MPI_CHAR, otherPID, 1, MPI_COMM_WORLD,
                 &status);

        Tcl_SetVar(interp, varToSet, gMsg, TCL_LEAVE_ERR_MSG);
      }

    } else {
      opserr << "recv -pid pid? data? - " << otherPID << " invalid\n";
      return TCL_ERROR;
    }
  } else {
    if (myPID != 0) {
      MPI_Bcast((void *)(&msgLength), 1, MPI_INT, 0, MPI_COMM_WORLD);

      if (msgLength > 0) {
        gMsg = new char[msgLength];

        MPI_Bcast((void *)gMsg, msgLength, MPI_CHAR, 0, MPI_COMM_WORLD);

        Tcl_SetVar(interp, varToSet, gMsg, TCL_LEAVE_ERR_MSG);
      }

    } else {
      opserr << "recv data - only process 0 can do a broadcast - you may need "
                "to kill the application";
      return TCL_ERROR;
    }
  }
  return TCL_OK;
}

static int
opsBarrier(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char ** const argv)
{
  return MPI_Barrier(MPI_COMM_WORLD);
}

