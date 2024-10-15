//===----------------------------------------------------------------------===//
//
//        OpenSees - Open System for Earthquake Engineering Simulation
//
//===----------------------------------------------------------------------===//
//
//                                   OpenSeesMP                               
//
//===----------------------------------------------------------------------===//
//
#ifndef OPENSEESRT_VERSION
#  define OPENSEESRT_VERSION "0.0.0"
#endif

extern "C" {
#include <tcl.h>
}


// #include "commands.h"
#include <ID.h>
#include <stdio.h>

#include "G3_Runtime.h"
#include <MPI_MachineBroker.h>
#include <TclPackageClassBroker.h>

#include <Channel.h>
#include <Message.h>

int Init_OpenSees(Tcl_Interp *interp);

void Init_MachineRuntime(Tcl_Interp* interp, MachineBroker* theMachineBroker);

void Init_Communication(Tcl_Interp* interp, MachineBroker* theMachineBroker);

static int 
doNothing(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char **argv)
{
  return TCL_OK;
}


extern "C" int 
#ifdef _WIN32
__declspec(dllexport)
#endif
Libopenseesmp_Init(Tcl_Interp* interp)
{
  if (Tcl_InitStubs(interp, TCL_VERSION, 0) == NULL)
    return TCL_ERROR;

  if (Tcl_PkgProvide(interp, "OpenSeesMP", OPENSEESRT_VERSION) == TCL_ERROR)
    return TCL_ERROR;

  // Create a runtime instance, and store it with the interpreter
  G3_Runtime *rt = new G3_Runtime{interp};
  Tcl_SetAssocData(interp, "G3_Runtime", NULL, (ClientData)rt);

  int argc = 0; 
  char **argv = nullptr;

  FEM_ObjectBroker* theBroker = new TclPackageClassBroker();
  MachineBroker* theMachineBroker = new MPI_MachineBroker(theBroker, argc, argv);

  
  // Initialize process runtime


  int pid = theMachineBroker->getPID();
  int np  = theMachineBroker->getNP();

  // TODO: These need to be stored so they can be passed
  // to some SOE constructors
  Channel **theChannels = nullptr;
  int numChannels;

  if (pid == 0) {
    theChannels = new Channel *[np-1];
    numChannels = np-1;
  } else {
    theChannels = new Channel *[1];
    numChannels = 1;
  }

  //
  // if rank 0 we send all args
  //

  int numArg = 0;
  int sizeArg = 0;
  char **args = 0;
  char *dataArgs = 0;

  if (pid == 0) {

    for (int i=0; i<argc; i++)
      if (argv[i] == NULL) {
        i = argc+1;
      } else {
        numArg++;
        sizeArg += strlen(argv[i])+1;
      }

    static ID data(2);
    data(0) = numArg;
    data(1) = sizeArg;

    dataArgs = new char[sizeArg];
    int loc = 0;
    for (int i=0; i<numArg; i++) {
      int lengthArg = strlen(argv[i]);
      strncpy(&dataArgs[loc], argv[i],lengthArg);
      loc += lengthArg;
      dataArgs[loc] = '\0';
      loc++; 
    }

    Message msgChar(dataArgs, sizeArg);

    for (int j=0; j<np-1; j++) {
      Channel *otherChannel = theMachineBroker->getRemoteProcess();
      theChannels[j] = otherChannel;
      otherChannel->sendID(0,0,data);
      otherChannel->sendMsg(0,0,msgChar);
    }
      
  } else {

    static ID data(2);    

    Channel *myChannel = theMachineBroker->getMyChannel();
    theChannels[0] = myChannel;

    myChannel->recvID(0,0,data);
    numArg = data(0);
    sizeArg = data(1);
    dataArgs = new char[sizeArg];
    Message msgChar(dataArgs, sizeArg);
    
    myChannel->recvMsg(0,0,msgChar);

  }

  args = new char *[numArg];
  args[0] = dataArgs;
  int argCount = 1;
  for (int j=1; j<sizeArg-1; j++)
    if (argCount < numArg && dataArgs[j] == '\0') {
      args[argCount] = &dataArgs[j+1];
      argCount++;
    }


  Init_OpenSees(interp);

  // Add machine commands (getPID, getNP, etc);
  Init_MachineRuntime(interp, theMachineBroker);
  Init_Communication(interp, theMachineBroker);

  Tcl_CreateCommand(interp, "partition", &doNothing, (ClientData)NULL, (Tcl_CmdDeleteProc *)NULL);


  return 0;
}

