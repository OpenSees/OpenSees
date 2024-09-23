//===----------------------------------------------------------------------===//
//
//        OpenSees - Open System for Earthquake Engineering Simulation
//
//===----------------------------------------------------------------------===//
//
#include <tcl.h>
#include <Integrator.h>
#include <BasicModelBuilder.h>


int
TclCommand_sensitivityAlgorithm(ClientData builder, Tcl_Interp* interp, int argc, TCL_Char**const argv)
{
    if (builder == nullptr)
      return 0;

    int analysisTypeTag = 1;

    Integrator* theIntegrator = nullptr;

    if (builder->getStaticIntegrator() != nullptr) {
    	theIntegrator = builder->getStaticIntegrator();

    } else if(builder->getTransientIntegrator() != nullptr) {
    	theIntegrator = builder->getTransientIntegrator();
    }


    // 1: compute at each step (default); 
    // 2: compute by command; 
    if (OPS_GetNumRemainingInputArgs() < 1) {
    	opserr << "ERROR: Wrong number of parameters to sensitivity algorithm." << "\n";
    	return -1;
    }
    if (theIntegrator == nullptr) {
    	opserr << "The integrator needs to be instantiated before " << "\n"
    	       << " setting  sensitivity algorithm." << "\n";
    	return -1;
    }

    const char* type = OPS_GetString();
    if (strcmp(type,"-computeAtEachStep") == 0)
    	analysisTypeTag = 1;
    else if (strcmp(type,"-computeByCommand") == 0)
    	analysisTypeTag = 2;
    else {
    	opserr << "Unknown sensitivity algorithm option: " << type << "\n";
    	return -1;
    }

    theIntegrator->setComputeType(analysisTypeTag);
    theIntegrator->activateSensitivityKey();
	
    return 0;
}

