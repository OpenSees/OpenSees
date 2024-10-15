//===----------------------------------------------------------------------===//
//
//        OpenSees - Open System for Earthquake Engineering Simulation
//
//===----------------------------------------------------------------------===//
//
//                                   OpenSeesSP                               
//
//===----------------------------------------------------------------------===//
//
#ifndef OPENSEESRT_VERSION
#  define OPENSEESRT_VERSION "0.0.0"
#endif

extern "C" {
#include <tcl.h>
}

#include <stdio.h>
#include <string.h>
#include "G3_Runtime.h"

#include <PartitionedDomain.h>
#include <MPI_MachineBroker.h>
#include <ShadowSubdomain.h>
#include <ActorSubdomain.h>
#include <TclPackageClassBroker.h>
#include <DomainPartitioner.h>

#include <mpi.h>

int Init_OpenSees(Tcl_Interp *interp);

void Init_MachineRuntime(Tcl_Interp* interp, MachineBroker* theMachineBroker);
void Init_PartitionRuntime(Tcl_Interp* interp, MachineBroker*, FEM_ObjectBroker*);

static int 
doNothing(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char **argv)
{
  return TCL_OK;
}


extern "C" int 
#ifdef _WIN32
__declspec(dllexport)
#endif
Libopenseessp_Init(Tcl_Interp* interp)
{
  if (Tcl_InitStubs(interp, TCL_VERSION, 0) == NULL)
    return TCL_ERROR;

  if (Tcl_PkgProvide(interp, "OpenSeesSP", OPENSEESRT_VERSION) == TCL_ERROR)
    return TCL_ERROR;

  // Create a runtime instance, and store it with the interpreter
  G3_Runtime *rt = new G3_Runtime{interp};
  Tcl_SetAssocData(interp, "G3_Runtime", NULL, (ClientData)rt);

  int argc = 0; 
  char **argv = nullptr;

  MachineBroker* theMachineBroker = new MPI_MachineBroker(0, argc, argv);
  FEM_ObjectBroker* theBroker = new TclPackageClassBroker();
  theMachineBroker->setObjectBroker(theBroker);

  int pid = theMachineBroker->getPID();
  int np = theMachineBroker->getNP();

  //
  // depending on rank we do something
  //
  if (pid != 0) {

    //
    // on secondary processes we spin waiting to create & run actors
    //
    fprintf(stderr, "Secondary Process Running %d\n", pid);
    theMachineBroker->runActors();

  } else {

    //
    // on process 0 we create some ShadowSubdomains & then start the OpenSees interpreter
    //
    fprintf(stderr, "Primary Process Running OpenSees Interpreter %d\n", pid);   


    // Add machine commands (getPID, getNP, etc);
    Init_MachineRuntime(interp, theMachineBroker);
    Init_OpenSees(interp);
    Init_PartitionRuntime(interp, theMachineBroker, theBroker);

    Tcl_CreateCommand(interp, "barrier",   &doNothing, (ClientData)NULL, (Tcl_CmdDeleteProc *)NULL);

    // TODO
//  // some clean up to shut the remotes down if still running
//  theDomain.clearAll();
//  
//  // shutdown the remote machines
//  theMachineBroker->shutdown();
  }
  return 0;
}


