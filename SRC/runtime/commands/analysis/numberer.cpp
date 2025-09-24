//===----------------------------------------------------------------------===//
//
//        OpenSees - Open System for Earthquake Engineering Simulation
//
//===----------------------------------------------------------------------===//
// 
// Description: This file implements the selection of a Numberer object,
// which is used to optimally number the degrees of freedom of a problem.
//
#include <tcl.h>
#include <assert.h>
#include <BasicAnalysisBuilder.h>
#include <PlainNumberer.h>
#include <DOF_Numberer.h>
#include <RCM.h>
#include <AMDNumberer.h>

#if defined(_PARALLEL_PROCESSING) || defined(_PARALLEL_INTERPRETERS)
#  include <ParallelNumberer.h>
#endif

//
// command that sets the Numberer object
//
int
TclCommand_setNumberer(ClientData clientData, Tcl_Interp* interp, int argc, TCL_Char ** const argv)
{
  BasicAnalysisBuilder *builder = (BasicAnalysisBuilder *)clientData;
  assert(builder != nullptr);


  DOF_Numberer *theNumberer = nullptr;

  // make sure at least one other argument to contain numberer
  if (argc < 2) {
    opserr << "WARNING need to specify a Numberer type \n";
    return TCL_ERROR;
  }

#if defined(_PARALLEL_PROCESSING)
  // check argv[1] for type of Numberer and create the object
  if (strcmp(argv[1], "Plain") == 0) {
    theNumberer = new ParallelNumberer();
  } else if (strcmp(argv[1], "RCM") == 0) {
    RCM *theRCM = new RCM(false);
    theNumberer = new ParallelNumberer(*theRCM);
  } else {
    opserr << "WARNING No Numberer type exists (Plain, RCM only) \n";
    return TCL_ERROR;
  }
#else

  // check argv[1] for type of Numberer and create the object
  if (strcmp(argv[1], "Plain") == 0) {
    theNumberer = new PlainNumberer();

  } else if (strcmp(argv[1], "RCM") == 0) {
    RCM *theRCM = new RCM(false);
    theNumberer = new DOF_Numberer(*theRCM);

  } else if (strcmp(argv[1], "AMD") == 0) {
    AMD *theAMD = new AMD();
    theNumberer = new DOF_Numberer(*theAMD);
  }

#  ifdef _PARALLEL_INTERPRETERS
  else if ((strcmp(argv[1], "ParallelPlain") == 0) ||
           (strcmp(argv[1], "Parallel") == 0)) {
    ParallelNumberer *theParallelNumberer = new ParallelNumberer;
    theNumberer = theParallelNumberer;
    theParallelNumberer->setProcessID(OPS_rank);
    theParallelNumberer->setChannels(numChannels, theChannels);

  } else if (strcmp(argv[1], "ParallelRCM") == 0) {
    RCM *theRCM = new RCM(false);
    ParallelNumberer *theParallelNumberer = new ParallelNumberer(*theRCM);
    theNumberer = theParallelNumberer;
    theParallelNumberer->setProcessID(OPS_rank);
    theParallelNumberer->setChannels(numChannels, theChannels);
  }
#  endif
#endif

  else {
    opserr << "WARNING No Numberer type exists (Plain, RCM only) \n";
    return TCL_ERROR;
  }

  if (theNumberer == nullptr)
    return TCL_ERROR;

  builder->set(theNumberer);
  return TCL_OK;
}

