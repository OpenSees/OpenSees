/* *****************************************************************************
Copyright (c) 2015-2017, The Regents of the University of California (Regents).
All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this
   list of conditions and the following disclaimer.
2. Redistributions in binary form must reproduce the above copyright notice,
   this list of conditions and the following disclaimer in the documentation
   and/or other materials provided with the distribution.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR
ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
(INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

The views and conclusions contained in the software and documentation are those
of the authors and should not be interpreted as representing official policies,
either expressed or implied, of the FreeBSD Project.

REGENTS SPECIFICALLY DISCLAIMS ANY WARRANTIES, INCLUDING, BUT NOT LIMITED TO,
THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
THE SOFTWARE AND ACCOMPANYING DOCUMENTATION, IF ANY, PROVIDED HEREUNDER IS
PROVIDED "AS IS". REGENTS HAS NO OBLIGATION TO PROVIDE MAINTENANCE, SUPPORT,
UPDATES, ENHANCEMENTS, OR MODIFICATIONS.

*************************************************************************** */


// Written: Minjie

// Description: all opensees APIs are defined or declared here
//

#include "OpenSeesCommands.h"
#include <OPS_Globals.h>
#include <elementAPI.h>
#include <UniaxialMaterial.h>
#include <NDMaterial.h>
#include <SectionForceDeformation.h>
#include <SectionRepres.h>
#include <TimeSeries.h>
#include <CrdTransf.h>
#include <BeamIntegration.h>
#include <NodalLoad.h>
#include <AnalysisModel.h>
#include <PlainHandler.h>
#include <RCM.h>
#include <AMDNumberer.h>
#include <LimitCurve.h>
#include <DamageModel.h>
#include <FrictionModel.h>
#include <HystereticBackbone.h>
#include <YieldSurface_BC.h>
#include <CyclicModel.h>
#include <FileStream.h>
#include <CTestNormUnbalance.h>
#include <NewtonRaphson.h>
#include <TransformationConstraintHandler.h>
#include <Newmark.h>
#include <ProfileSPDLinSolver.h>
#include <ProfileSPDLinDirectSolver.h>
#include <ProfileSPDLinSOE.h>
#include <SymBandEigenSolver.h>
#include <SymBandEigenSOE.h>
#include <FullGenEigenSolver.h>
#include <FullGenEigenSOE.h>
#include <ArpackSOE.h>
#include <LoadControl.h>
#include <CTestPFEM.h>
#include <PFEMIntegrator.h>
#include <TransientIntegrator.h>
#include <PFEMSolver.h>
#include <PFEMLinSOE.h>
#include <Accelerator.h>
#include <KrylovAccelerator.h>
#include <AcceleratedNewton.h>
#include <RaphsonAccelerator.h>
#include <SecantAccelerator2.h>
#include <PeriodicAccelerator.h>
#include <LineSearch.h>
#include <InitialInterpolatedLineSearch.h>
#include <BisectionLineSearch.h>
#include <SecantLineSearch.h>
#include <RegulaFalsiLineSearch.h>
#include <NewtonLineSearch.h>
#include <FileDatastore.h>
#include <Mesh.h>

#ifdef _PARALLEL_INTERPRETERS
#include <mpi.h>
#include <MPI_MachineBroker.h>
#endif


// active object
static OpenSeesCommands* cmds = 0;

OpenSeesCommands::OpenSeesCommands(DL_Interpreter* interp)
    :interpreter(interp), theDomain(0), ndf(0), ndm(0),
     theSOE(0), theEigenSOE(0), theNumberer(0), theHandler(0),
     theStaticIntegrator(0), theTransientIntegrator(0),
     theAlgorithm(0), theStaticAnalysis(0), theTransientAnalysis(0),
     thePFEMAnalysis(0),
     theAnalysisModel(0), theTest(0), numEigen(0), theDatabase(0),
     theBroker(), theTimer(), theSimulationInfo(), theMachineBroker(0),
     reliability(0)
{
#ifdef _PARALLEL_INTERPRETERS
    theMachineBroker = new MPI_MachineBroker(&theBroker, 0, 0);
#endif

    cmds = this;

    theDomain = new Domain;

    reliability = new OpenSeesReliabilityCommands(theDomain);
}

OpenSeesCommands::~OpenSeesCommands()
{
    if (reliability != 0) delete reliability;
    if (theDomain != 0) delete theDomain;
    if (theDatabase != 0) delete theDatabase;
    cmds = 0;

#ifdef _PARALLEL_INTERPRETERS
    if (theMachineBroker != 0) {
	opserr << "Process "<<theMachineBroker->getPID() << " Terminating\n";
	theMachineBroker->shutdown();
	delete theMachineBroker;
    }
#endif

}

DL_Interpreter*
OpenSeesCommands::getInterpreter()
{
    return interpreter;
}

Domain*
OpenSeesCommands::getDomain()
{
    return theDomain;
}

void
OpenSeesCommands::setSOE(LinearSOE* soe)
{
    // if not in analysis object, delete old one
    if (theStaticAnalysis==0 && theTransientAnalysis==0) {
	if (theSOE != 0) {
	    delete theSOE;
	    theSOE = 0;
	}
    }

    // set new one
    theSOE = soe;
    if (soe == 0) return;

    // set in analysis object
    if (theStaticAnalysis != 0) {
	theStaticAnalysis->setLinearSOE(*soe);
    }
    if (theTransientAnalysis != 0) {
	theTransientAnalysis->setLinearSOE(*soe);
    }
}

int
OpenSeesCommands::eigen(int typeSolver, double shift,
			bool generalizedAlgo, bool findSmallest)
{
    //
    // create a transient analysis if no analysis exists
    //
    bool newanalysis = false;
    if (theStaticAnalysis == 0 && theTransientAnalysis == 0) {
	if (theAnalysisModel == 0)
	    theAnalysisModel = new AnalysisModel();
	if (theTest == 0)
	    theTest = new CTestNormUnbalance(1.0e-6,25,0);
	if (theAlgorithm == 0) {
	    theAlgorithm = new NewtonRaphson(*theTest);
	}
	if (theHandler == 0) {
	    theHandler = new TransformationConstraintHandler();
	}
	if (theNumberer == 0) {
	    RCM *theRCM = new RCM(false);
	    theNumberer = new DOF_Numberer(*theRCM);
	}
	if (theTransientIntegrator == 0) {
	    theTransientIntegrator = new Newmark(0.5,0.25);
	}
	if (theSOE == 0) {
	    ProfileSPDLinSolver *theSolver;
	    theSolver = new ProfileSPDLinDirectSolver();
	    theSOE = new ProfileSPDLinSOE(*theSolver);
	}

	theTransientAnalysis = new DirectIntegrationAnalysis(*theDomain,
							     *theHandler,
							     *theNumberer,
							     *theAnalysisModel,
							     *theAlgorithm,
							     *theSOE,
							     *theTransientIntegrator,
							     theTest);
	newanalysis = true;
    }

    //
    // create a new eigen system and solver
    //
    if (theEigenSOE != 0) {
	if (theEigenSOE->getClassTag() != typeSolver) {
	    //	delete theEigenSOE;
	    theEigenSOE = 0;
	}
    }

    if (theEigenSOE == 0) {

	if (typeSolver == EigenSOE_TAGS_SymBandEigenSOE) {
	    SymBandEigenSolver *theEigenSolver = new SymBandEigenSolver();
	    theEigenSOE = new SymBandEigenSOE(*theEigenSolver, *theAnalysisModel);

	} else if (typeSolver == EigenSOE_TAGS_FullGenEigenSOE) {

	    FullGenEigenSolver *theEigenSolver = new FullGenEigenSolver();
	    theEigenSOE = new FullGenEigenSOE(*theEigenSolver, *theAnalysisModel);

	} else {

	    theEigenSOE = new ArpackSOE(shift);

	}

	//
	// set the eigen soe in the system
	//

	if (theStaticAnalysis != 0) {
	    theStaticAnalysis->setEigenSOE(*theEigenSOE);
	} else if (theTransientAnalysis != 0) {
	    theTransientAnalysis->setEigenSOE(*theEigenSOE);
	}

    } // theEigenSOE != 0


    // run analysis
    int result = 0;
    if (theStaticAnalysis != 0) {
	result = theStaticAnalysis->eigen(numEigen, generalizedAlgo, findSmallest);
    } else if (theTransientAnalysis != 0) {
	result = theTransientAnalysis->eigen(numEigen, generalizedAlgo, findSmallest);
    }
    if (newanalysis) {
	delete theTransientAnalysis;
	theTransientAnalysis = 0;
    }

    if (result == 0) {
	const Vector &eigenvalues = theDomain->getEigenvalues();
	double* data = new double[numEigen];
	for (int i=0; i<numEigen; i++) {
	    data[i] = eigenvalues(i);
	}
	OPS_SetDoubleOutput(&numEigen, data);
	delete [] data;
    }

    return result;
}

int* OPS_GetNumEigen()                                                          
{                                                                               
    static int numEigen = 0;                                                    
    if (cmds == 0) return 0;                                                    
    numEigen = cmds->getNumEigen();                                             
    int numdata = 1;                                                            
    if (OPS_SetIntOutput(&numdata, &numEigen) < 0) {                            
        opserr << "WARNING failed to set output\n";                             
        return 0;                                                               
    }                                                                           
                                                                                
    return &numEigen;                                                           
}

void
OpenSeesCommands::setNumberer(DOF_Numberer* numberer)
{
    // if not in analysis object, delete old one
    if (theStaticAnalysis==0 && theTransientAnalysis==0) {
	if (theNumberer != 0) {
	    delete theNumberer;
	    theNumberer = 0;
	}
    }

    // set new one
    theNumberer = numberer;
    if (numberer == 0) return;

    // set in analysis object
    if (theStaticAnalysis != 0) {
	theStaticAnalysis->setNumberer(*numberer);
    }
    if (theTransientAnalysis != 0) {
	theTransientAnalysis->setNumberer(*numberer);
    }
}

void
OpenSeesCommands::setHandler(ConstraintHandler* handler)
{
    // if not in analysis object, delete old one and set new one
    if (theStaticAnalysis==0 && theTransientAnalysis==0) {
	if (theHandler != 0) {
	    delete theHandler;
	    theHandler = 0;
	}
	theHandler = handler;
	return;
    }

    // if analysis object is created, not set new one
    if (handler != 0) {
	opserr<<"WARNING can't set handler after analysis is created\n";
	delete handler;
    }

}

void
OpenSeesCommands::setCTest(ConvergenceTest* test)
{
    // if not in analysis object, delete old one and set new one
    if (theStaticAnalysis==0 && theTransientAnalysis==0) {
	if (theTest != 0) {
	    delete theTest;
	    theTest = 0;
	}
	theTest = test;
	return;
    }

    // set new one
    theTest = test;
    if (test == 0) return;

    // set in analysis object
    if (theStaticAnalysis != 0) {
	theStaticAnalysis->setConvergenceTest(*test);
    }
    if (theTransientAnalysis != 0) {
	theTransientAnalysis->setConvergenceTest(*test);
    }
}

void
OpenSeesCommands::setStaticIntegrator(StaticIntegrator* integrator)
{
    // error in transient analysis
    if (theTransientAnalysis != 0) {
	opserr<<"WARNING can't set static integrator in transient analysis\n";
	if (integrator != 0) {
	    delete integrator;
	}
	return;
    }

    // if not in static analysis object, delete old one
    if (theStaticAnalysis==0) {
	if (theStaticIntegrator != 0) {
	    delete theStaticIntegrator;
	    theStaticIntegrator = 0;
	}
    }

    // set new one
    theStaticIntegrator = integrator;
    if (integrator == 0) return;

    // set in analysis object
    if (theStaticAnalysis != 0) {
	theStaticAnalysis->setIntegrator(*integrator);
    }
}

void
OpenSeesCommands::setTransientIntegrator(TransientIntegrator* integrator)
{
    // error in static analysis
    if (theStaticAnalysis != 0) {
	opserr<<"WARNING can't set transient integrator in static analysis\n";
	if (integrator != 0) {
	    delete integrator;
	}
	return;
    }

    // if not in transient analysis object, delete old one
    if (theTransientAnalysis==0) {
	if (theTransientIntegrator != 0) {
	    delete theTransientIntegrator;
	    theTransientIntegrator = 0;
	}
    }

    // set new one
    theTransientIntegrator = integrator;
    if (integrator == 0) return;

    // set in analysis object
    if (theTransientAnalysis != 0) {
	theTransientAnalysis->setIntegrator(*integrator);
    }
}

void
OpenSeesCommands::setAlgorithm(EquiSolnAlgo* algorithm)
{
    // if not in analysis object, delete old one
    if (theStaticAnalysis==0 && theTransientAnalysis==0) {
	if (theAlgorithm != 0) {
	    delete theAlgorithm;
	    theAlgorithm = 0;
	}
    }

    // set new one
    theAlgorithm = algorithm;
    if (algorithm == 0) return;

    // set in analysis object
    if (theStaticAnalysis != 0) {
	theStaticAnalysis->setAlgorithm(*algorithm);
	if (theTest != 0) {
	    algorithm->setConvergenceTest(theTest);
	}
    }
    if (theTransientAnalysis != 0) {
	theTransientAnalysis->setAlgorithm(*algorithm);
	if (theTest != 0) {
	    algorithm->setConvergenceTest(theTest);
	}
    }
}

void
OpenSeesCommands::setStaticAnalysis()
{
    // delete the old analysis
    if (theStaticAnalysis != 0) {
	delete theStaticAnalysis;
	theStaticAnalysis = 0;
    }
    if (theTransientAnalysis != 0) {
	delete theTransientAnalysis;
	theTransientAnalysis = 0;
    }

    // create static analysis
    if (theAnalysisModel == 0) {
	theAnalysisModel = new AnalysisModel();
    }
    if (theTest == 0) {
	theTest = new CTestNormUnbalance(1.0e-6,25,0);
    }
    if (theAlgorithm == 0) {
	opserr << "WARNING analysis Static - no Algorithm yet specified, \n";
	opserr << " NewtonRaphson default will be used\n";
	theAlgorithm = new NewtonRaphson(*theTest);
    }
    if (theHandler == 0) {
	opserr << "WARNING analysis Static - no ConstraintHandler yet specified, \n";
	opserr << " PlainHandler default will be used\n";
	theHandler = new PlainHandler();
    }
    if (theNumberer == 0) {
	opserr << "WARNING analysis Static - no Numberer specified, \n";
	opserr << " RCM default will be used\n";
	RCM* theRCM = new RCM(false);
	theNumberer = new DOF_Numberer(*theRCM);
    }
    if (theStaticIntegrator == 0) {
	opserr << "WARNING analysis Static - no Integrator specified, \n";
	opserr << " StaticIntegrator default will be used\n";
	theStaticIntegrator = new LoadControl(1, 1, 1, 1);
    }
    if (theSOE == 0) {
	opserr << "WARNING analysis Static - no LinearSOE specified, \n";
	opserr << " ProfileSPDLinSOE default will be used\n";
	ProfileSPDLinSolver *theSolver;
	theSolver = new ProfileSPDLinDirectSolver();
	theSOE = new ProfileSPDLinSOE(*theSolver);
    }

    theStaticAnalysis = new StaticAnalysis(*theDomain,
    					   *theHandler,
    					   *theNumberer,
    					   *theAnalysisModel,
    					   *theAlgorithm,
    					   *theSOE,
    					   *theStaticIntegrator,
    					   theTest);

    if (theEigenSOE != 0) {
	theStaticAnalysis->setEigenSOE(*theEigenSOE);
    }
}

int
OpenSeesCommands::setPFEMAnalysis()
{
    // delete the old analysis
    if (theStaticAnalysis != 0) {
	delete theStaticAnalysis;
	theStaticAnalysis = 0;
    }
    if (theTransientAnalysis != 0) {
	delete theTransientAnalysis;
	theTransientAnalysis = 0;
    }

    // create PFEM analysis
    if(OPS_GetNumRemainingInputArgs() < 3) {
	opserr<<"WARNING: wrong no of args -- analysis PFEM dtmax dtmin gravity <ratio>\n";
	return -1;
    }

    int numdata = 1;
    double dtmax, dtmin, gravity, ratio=0.5;
    if (OPS_GetDoubleInput(&numdata, &dtmax) < 0) {
	opserr<<"WARNING: invalid dtmax \n";
	return -1;
    }
    if (OPS_GetDoubleInput(&numdata, &dtmin) < 0) {
	opserr<<"WARNING: invalid dtmin \n";
	return -1;
    }
    if (OPS_GetDoubleInput(&numdata, &gravity) < 0) {
	opserr<<"WARNING: invalid gravity \n";
	return -1;
    }
    if(OPS_GetNumRemainingInputArgs() > 0) {
	if (OPS_GetDoubleInput(&numdata, &ratio) < 0) {
	    opserr<<"WARNING: invalid ratio \n";
	    return -1;
	}
    }

    if(theAnalysisModel == 0) {
	theAnalysisModel = new AnalysisModel();
    }
    if(theTest == 0) {
	//theTest = new CTestNormUnbalance(1e-2,10000,1,2,3);
	theTest = new CTestPFEM(1e-2,1e-2,1e-2,1e-2,1e-4,1e-3,10000,100,1,2);
    }
    if(theAlgorithm == 0) {
	theAlgorithm = new NewtonRaphson(*theTest);
    }
    if(theHandler == 0) {
	theHandler = new TransformationConstraintHandler();
    }
    if(theNumberer == 0) {
	RCM* theRCM = new RCM(false);
	theNumberer = new DOF_Numberer(*theRCM);
    }
    if(theTransientIntegrator == 0) {
	theTransientIntegrator = new PFEMIntegrator();
    }
    if(theSOE == 0) {
	PFEMSolver* theSolver = new PFEMSolver();
	theSOE = new PFEMLinSOE(*theSolver);
    }
    thePFEMAnalysis = new PFEMAnalysis(*theDomain,
				       *theHandler,
				       *theNumberer,
				       *theAnalysisModel,
				       *theAlgorithm,
				       *theSOE,
				       *theTransientIntegrator,
				       theTest,dtmax,dtmin,gravity,ratio);

    theTransientAnalysis = thePFEMAnalysis;

    if (theEigenSOE != 0) {
	theTransientAnalysis->setEigenSOE(*theEigenSOE);
    }

    return 0;
}

void
OpenSeesCommands::setVariableAnalysis()
{
    // delete the old analysis
    if (theStaticAnalysis != 0) {
	delete theStaticAnalysis;
	theStaticAnalysis = 0;
    }
    if (theTransientAnalysis != 0) {
	delete theTransientAnalysis;
	theTransientAnalysis = 0;
    }

    // make sure all the components have been built,
    // otherwise print a warning and use some defaults
    if (theAnalysisModel == 0) {
	theAnalysisModel = new AnalysisModel();
    }

    if (theTest == 0) {
	theTest = new CTestNormUnbalance(1.0e-6,25,0);
    }

    if (theAlgorithm == 0) {
	opserr << "WARNING analysis Transient - no Algorithm yet specified, \n";
	opserr << " NewtonRaphson default will be used\n";
	theAlgorithm = new NewtonRaphson(*theTest);
    }

    if (theHandler == 0) {
	opserr << "WARNING analysis Transient dt tFinal - no ConstraintHandler\n";
	opserr << " yet specified, PlainHandler default will be used\n";
	theHandler = new PlainHandler();
    }

    if (theNumberer == 0) {
	opserr << "WARNING analysis Transient dt tFinal - no Numberer specified, \n";
	opserr << " RCM default will be used\n";
	RCM *theRCM = new RCM(false);
	theNumberer = new DOF_Numberer(*theRCM);
    }

    if (theTransientIntegrator == 0) {
	opserr << "WARNING analysis Transient dt tFinal - no Integrator specified, \n";
	opserr << " Newmark(.5,.25) default will be used\n";
	theTransientIntegrator = new Newmark(0.5,0.25);
    }

    if (theSOE == 0) {
	opserr << "WARNING analysis Transient dt tFinal - no LinearSOE specified, \n";
	opserr << " ProfileSPDLinSOE default will be used\n";
	ProfileSPDLinSolver *theSolver;
	theSolver = new ProfileSPDLinDirectSolver();
	theSOE = new ProfileSPDLinSOE(*theSolver);
    }

    theVariableTimeStepTransientAnalysis = new VariableTimeStepDirectIntegrationAnalysis
	(*theDomain,
	 *theHandler,
	 *theNumberer,
	 *theAnalysisModel,
	 *theAlgorithm,
	 *theSOE,
	 *theTransientIntegrator,
	 theTest);

    // set the pointer for variabble time step analysis
    theTransientAnalysis = theVariableTimeStepTransientAnalysis;

    if (theEigenSOE != 0) {
	theTransientAnalysis->setEigenSOE(*theEigenSOE);
    }

}

void
OpenSeesCommands::setTransientAnalysis()
{
    // delete the old analysis
    if (theStaticAnalysis != 0) {
	delete theStaticAnalysis;
	theStaticAnalysis = 0;
    }
    if (theTransientAnalysis != 0) {
	delete theTransientAnalysis;
	theTransientAnalysis = 0;
    }

    // create transient analysis
    if (theAnalysisModel == 0) {
	theAnalysisModel = new AnalysisModel();
    }
    if (theTest == 0) {
	theTest = new CTestNormUnbalance(1.0e-6,25,0);
    }
    if (theAlgorithm == 0) {
	opserr << "WARNING analysis Transient - no Algorithm yet specified, \n";
	opserr << " NewtonRaphson default will be used\n";
	theAlgorithm = new NewtonRaphson(*theTest);
    }
    if (theHandler == 0) {
	opserr << "WARNING analysis Transient - no ConstraintHandler yet specified, \n";
	opserr << " PlainHandler default will be used\n";
	theHandler = new PlainHandler();
    }
    if (theNumberer == 0) {
	opserr << "WARNING analysis Transient - no Numberer specified, \n";
	opserr << " RCM default will be used\n";
	RCM* theRCM = new RCM(false);
	theNumberer = new DOF_Numberer(*theRCM);
    }
    if (theTransientIntegrator == 0) {
	opserr << "WARNING analysis Transient - no Integrator specified, \n";
	opserr << " TransientIntegrator default will be used\n";
	theTransientIntegrator = new Newmark(0.5,0.25);
    }
    if (theSOE == 0) {
	opserr << "WARNING analysis Transient - no LinearSOE specified, \n";
	opserr << " ProfileSPDLinSOE default will be used\n";
	ProfileSPDLinSolver *theSolver;
	theSolver = new ProfileSPDLinDirectSolver();
	theSOE = new ProfileSPDLinSOE(*theSolver);
    }

    theTransientAnalysis = new DirectIntegrationAnalysis(*theDomain,
							 *theHandler,
							 *theNumberer,
							 *theAnalysisModel,
							 *theAlgorithm,
							 *theSOE,
							 *theTransientIntegrator,
							 theTest);
    if (theEigenSOE != 0) {
	theTransientAnalysis->setEigenSOE(*theEigenSOE);
    }

}

void
OpenSeesCommands::wipeAnalysis()
{
    if (theStaticAnalysis==0 && theTransientAnalysis==0) {
	if (theSOE != 0) delete theSOE;
	if (theEigenSOE != 0) delete theEigenSOE;
	if (theNumberer != 0) delete theNumberer;
	if (theHandler != 0) delete theHandler;
	if (theStaticIntegrator != 0) delete theStaticIntegrator;
	if (theTransientIntegrator != 0) delete theTransientIntegrator;
	if (theAlgorithm != 0) delete theAlgorithm;
	if (theTest != 0) delete theTest;
    }

    if (theStaticAnalysis != 0) {
    	theStaticAnalysis->clearAll();
    	delete theStaticAnalysis;
    }
    if (theTransientAnalysis != 0) {
    	theTransientAnalysis->clearAll();
    	delete theTransientAnalysis;
    }

    theAlgorithm = 0;
    theHandler = 0;
    theNumberer = 0;
    theAnalysisModel = 0;
    theSOE = 0;
    theEigenSOE = 0;
    theStaticIntegrator = 0;
    theTransientIntegrator = 0;
    theStaticAnalysis = 0;
    theTransientAnalysis = 0;
    thePFEMAnalysis = 0;
    theTest = 0;

}

void
OpenSeesCommands::wipe()
{
    this->wipeAnalysis();

    // data base
    if (theDatabase != 0) {
	delete theDatabase;
	theDatabase = 0;
    }

    // wipe all meshes
    OPS_clearAllMesh();

    // wipe domain
    if (theDomain != 0) {
	theDomain->clearAll();
    }

    // time set to zero
    ops_Dt = 0.0;

    // wipe uniaxial material
    OPS_clearAllUniaxialMaterial();
    OPS_clearAllNDMaterial();

    // wipe sections
    OPS_clearAllSectionForceDeformation();
    OPS_clearAllSectionRepres();

    // wipe time series
    OPS_clearAllTimeSeries();

    // wipe GeomTransf
    OPS_clearAllCrdTransf();

    // wipe BeamIntegration
    OPS_clearAllBeamIntegrationRule();

    // wipe limit state curve
    OPS_clearAllLimitCurve();

    // wipe damages model
    OPS_clearAllDamageModel();

    // wipe friction model
    OPS_clearAllFrictionModel();

    // wipe HystereticBackbone
    OPS_clearAllHystereticBackbone();

    // wipe YieldSurface_BC
    OPS_clearAllYieldSurface_BC();

    // wipe CyclicModel
    OPS_clearAllCyclicModel();

}

void
OpenSeesCommands::setFileDatabase(const char* filename)
{
    if (theDatabase != 0) delete theDatabase;
    theDatabase = new FileDatastore(filename, *theDomain, theBroker);
    if (theDatabase == 0) {
	opserr << "WARNING ran out of memory - database File " << filename << endln;
    }
}

/////////////////////////////
//// OpenSees APIs  /// /////
/////////////////////////////
int OPS_GetNumRemainingInputArgs()
{
    if (cmds == 0) return 0;
    DL_Interpreter* interp = cmds->getInterpreter();
    return interp->getNumRemainingInputArgs();
}

int OPS_GetIntInput(int *numData, int*data)
{
    if (cmds == 0) return 0;
    DL_Interpreter* interp = cmds->getInterpreter();
    if (numData == 0 || data == 0) return -1;
    return interp->getInt(data, *numData);
}

int OPS_SetIntOutput(int *numData, int*data)
{
    if (cmds == 0) return 0;
    DL_Interpreter* interp = cmds->getInterpreter();
    if (numData == 0 || data == 0) return -1;
    return interp->setInt(data, *numData);
}

int OPS_GetDoubleInput(int *numData, double *data)
{
    if (cmds == 0) return 0;
    DL_Interpreter* interp = cmds->getInterpreter();
    if (numData == 0 || data == 0) return -1;
    return interp->getDouble(data, *numData);
}

int OPS_SetDoubleOutput(int *numData, double *data)
{
    if (cmds == 0) return 0;
    DL_Interpreter* interp = cmds->getInterpreter();
    if (numData == 0 || data == 0) return -1;
    return interp->setDouble(data, *numData);
}

const char * OPS_GetString(void)
{
    if (cmds == 0) return "Invalid String Input!";
    DL_Interpreter* interp = cmds->getInterpreter();
    const char* res = interp->getString();
    if (res == 0) {
	return "Invalid String Input!";
    }
    return res;
}

int OPS_SetString(const char* str)
{
    if (cmds == 0) return 0;
    DL_Interpreter* interp = cmds->getInterpreter();
    return interp->setString(str);
}

Domain* OPS_GetDomain(void)
{
    if (cmds == 0) return 0;
    return cmds->getDomain();
}

int OPS_GetNDF()
{
    if (cmds == 0) return 0;
    return cmds->getNDF();
}

int OPS_GetNDM()
{
    if (cmds == 0) return 0;
    return cmds->getNDM();
}

int OPS_ResetCurrentInputArg(int cArg)
{
    if (cArg == 0) {
	opserr << "WARNING can't reset to argv[0]\n";
	return -1;
    }
    if (cmds == 0) return 0;
    DL_Interpreter* interp = cmds->getInterpreter();
    interp->resetInput(cArg);
    return 0;
}

UniaxialMaterial *OPS_GetUniaxialMaterial(int matTag)
{
    return OPS_getUniaxialMaterial(matTag);
}

int OPS_wipe()
{
    // wipe
    if (cmds != 0) {
	cmds->wipe();
    }

    return 0;
}

int OPS_wipeAnalysis()
{
    // wipe analysis
    if (cmds != 0) {
	cmds->wipeAnalysis();
    }

    return 0;
}

int OPS_model()
{
    // num args
    if(OPS_GetNumRemainingInputArgs() < 3) {
	opserr<<"WARNING insufficient args: model -ndm ndm <-ndf ndf>\n";
	return -1;
    }

    // model type
    const char* modeltype = OPS_GetString();
    if (strcmp(modeltype,"basic")!=0 && strcmp(modeltype,"Basic")!=0 &&
	strcmp(modeltype,"BasicBuilder")!=0 && strcmp(modeltype,"basicBuilder")!=0) {
	opserr<<"WARNING only basic builder is available at this time\n";
	return -1;
    }

    // ndm
    const char* ndmopt = OPS_GetString();
    if (strcmp(ndmopt,"-ndm") != 0) {
	opserr<<"WARNING frist option must be -ndm\n";
	return -1;
    }
    int numdata = 1;
    int ndm = 0;
    if (OPS_GetIntInput(&numdata, &ndm) < 0) {
	opserr<<"WARNING failed to read ndm\n";
	return -1;
    }
    if (ndm!=1 && ndm!=2 && ndm!=3) {
	opserr<<"ERROR ndm msut be 1, 2 or 3\n";
	return -1;
    }

    // ndf
    int ndf = 0;
    if (OPS_GetNumRemainingInputArgs() > 1) {
	const char* ndfopt = OPS_GetString();
	if (strcmp(ndfopt,"-ndf") != 0) {
	    opserr<<"WARNING second option must be -ndf\n";
	    return -1;
	}
	if (OPS_GetIntInput(&numdata, &ndf) < 0) {
	    opserr<<"WARNING failed to read ndf\n";
	    return -1;
	}
    }

    // set ndf
    if(ndf <= 0) {
	if (ndm == 1)
	    ndf = 1;
	else if (ndm == 2)
	    ndf = 3;
	else if (ndm == 3)
	    ndf = 6;
    }

    // set ndm and ndf
    if (cmds != 0) {
	cmds->setNDF(ndf);
	cmds->setNDM(ndm);
    }

    return 0;
}

int OPS_System()
{
    if (OPS_GetNumRemainingInputArgs() < 1) {
    	opserr << "WARNING insufficient args: system type ...\n";
    	return -1;
    }

    const char* type = OPS_GetString();

    // create soe
    LinearSOE* theSOE = 0;

    if ((strcmp(type,"BandGeneral") == 0) || (strcmp(type,"BandGEN") == 0)
    	|| (strcmp(type,"BandGen") == 0)){
	// BAND GENERAL SOE & SOLVER
    	theSOE = (LinearSOE*)OPS_BandGenLinLapack();

    } else if (strcmp(type,"BandSPD") == 0) {
	// BAND SPD SOE & SOLVER
    	theSOE = (LinearSOE*)OPS_BandSPDLinLapack();

    } else if (strcmp(type,"Diagonal") == 0) {
	// Diagonal SOE & SOLVER
	theSOE = (LinearSOE*)OPS_DiagonalDirectSolver();


    } else if (strcmp(type,"MPIDiagonal") == 0) {
	// Diagonal SOE & SOLVER
	theSOE = (LinearSOE*)OPS_DiagonalDirectSolver();

    } else if (strcmp(type,"SProfileSPD") == 0) {
	// PROFILE SPD SOE * SOLVER
	// now must determine the type of solver to create from rest of args
	theSOE = (LinearSOE*)OPS_SProfileSPDLinSolver();

    } else if (strcmp(type, "ProfileSPD") == 0) {

	theSOE = (LinearSOE*)OPS_ProfileSPDLinDirectSolver();

    } else if (strcmp(type, "PFEM") == 0) {
	// PFEM SOE & SOLVER

	if(OPS_GetNumRemainingInputArgs() < 1) {
	    theSOE = (LinearSOE*)OPS_PFEMSolver_Umfpack();
	} else {

	    const char* type = OPS_GetString();

	    if(strcmp(type, "-compressible") == 0) {

		theSOE = (LinearSOE*)OPS_PFEMCompressibleSolver();

	    } else if (strcmp(type, "-laplace") == 0) {

	    	theSOE = (LinearSOE*)OPS_PFEMSolver_Laplace();
	    } else if (strcmp(type, "-lumpM") == 0) {
	    	theSOE = (LinearSOE*)OPS_PFEMSolver_LumpM();		

	    // } else if(strcmp(type, "-umfpack") == 0) {

	    // 	theSOE = (LinearSOE*)OPS_PFEMSolver_Umfpack();

	    // } else if(strcmp(type, "-diag") == 0) {

	    // 	theSOE = (LinearSOE*)OPS_PFEMDiaSolver();

	    // } else if(strcmp(type, "-egen") == 0) {

	    //  theSOE = (LinearSOE*)OPS_EgenSolver();

	    } else if (strcmp(type,"-mumps") ==0) {
// #ifdef _PARALLEL_INTERPRETERS
// 	    int relax = 20;
// 	    if (argc > 3) {
// 		if (Tcl_GetInt(interp, argv[3], &relax) != TCL_OK) {
// 		    opserr<<"WARNING: failed to read relax\n";
// 		    return TCL_ERROR;
// 		}
// 	    }
// 	    PFEMSolver_Mumps* theSolver = new PFEMSolver_Mumps(relax,0,0,0);
// 	    theSOE = new PFEMLinSOE(*theSolver);
// #endif
	    } else if (strcmp(type,"-quasi-mumps")==0) {
// #ifdef _PARALLEL_INTERPRETERS
// 	    int relax = 20;
// 	    if (argc > 3) {
// 		if (Tcl_GetInt(interp, argv[3], &relax) != TCL_OK) {
// 		    opserr<<"WARNING: failed to read relax\n";
// 		    return TCL_ERROR;
// 		}
// 	    }
// 	    PFEMCompressibleSolver_Mumps* theSolver = new PFEMCompressibleSolver_Mumps(relax,0,0);
// 	    theSOE = new PFEMCompressibleLinSOE(*theSolver);
// #endif

	    }
	}


    } else if ((strcmp(type,"SparseGeneral") == 0) ||
	       (strcmp(type,"SuperLU") == 0) ||
	       (strcmp(type,"SparseGEN") == 0)) {

	// SPARSE GENERAL SOE * SOLVER
	theSOE = (LinearSOE*)OPS_SuperLUSolver();


    } else if ((strcmp(type,"SparseSPD") == 0) || (strcmp(type,"SparseSYM") == 0)) {
	// now must determine the type of solver to create from rest of args
	theSOE = (LinearSOE*)OPS_SymSparseLinSolver();

    } else if (strcmp(type, "UmfPack") == 0 || strcmp(type, "Umfpack") == 0) {

	theSOE = (LinearSOE*)OPS_UmfpackGenLinSolver();

    } else if (strcmp(type,"FullGeneral") == 0) {
	// now must determine the type of solver to create from rest of args
	theSOE = (LinearSOE*)OPS_FullGenLinLapackSolver();

    } else if (strcmp(type,"Petsc") == 0) {

    } else if (strcmp(type,"Mumps") == 0) {


    } else {
    	opserr<<"WARNING unknown system type "<<type<<"\n";
    	return -1;
    }

    // set soe
    if (cmds != 0) {
	cmds->setSOE(theSOE);
    }

    return 0;
}

int OPS_Numberer()
{
    if (OPS_GetNumRemainingInputArgs() < 1) {
    	opserr << "WARNING insufficient args: numberer type ...\n";
    	return -1;
    }

    const char* type = OPS_GetString();

    // create numberer
    DOF_Numberer *theNumberer = 0;
    if (strcmp(type,"Plain") == 0) {

    	theNumberer = (DOF_Numberer*)OPS_PlainNumberer();

    } else if (strcmp(type,"RCM") == 0) {

    	RCM *theRCM = new RCM(false);
    	theNumberer = new DOF_Numberer(*theRCM);

    } else if (strcmp(type,"AMD") == 0) {

    	AMD *theAMD = new AMD();
    	theNumberer = new DOF_Numberer(*theAMD);

    } else {
    	opserr<<"WARNING unknown numberer type "<<type<<"\n";
    	return -1;
    }

    // set numberer
    if (cmds != 0) {
	cmds->setNumberer(theNumberer);
    }

    return 0;
}

int OPS_ConstraintHandler()
{
    if (OPS_GetNumRemainingInputArgs() < 1) {
    	opserr << "WARNING insufficient args: constraints type ...\n";
    	return -1;
    }

    const char* type = OPS_GetString();

    // create handler
    ConstraintHandler* theHandler = 0;
    if (strcmp(type,"Plain") == 0) {

    	theHandler = (ConstraintHandler*)OPS_PlainHandler();

    } else if (strcmp(type,"Penalty") == 0) {

    	theHandler = (ConstraintHandler*)OPS_PenaltyConstraintHandler();

    } else if (strcmp(type,"Lagrange") == 0) {

    	theHandler = (ConstraintHandler*)OPS_LagrangeConstraintHandler();

    } else if (strcmp(type,"Transformation") == 0) {
    	theHandler = (ConstraintHandler*)OPS_TransformationConstraintHandler();

    } else {
    	opserr<<"WARNING unknown ConstraintHandler type "<<type<<"\n";
    	return -1;
    }

    // set handler
    if (cmds != 0) {
	cmds->setHandler(theHandler);
    }

    return 0;
}

int OPS_CTest()
{
    if (OPS_GetNumRemainingInputArgs() < 1) {
    	opserr << "WARNING insufficient args: test type ...\n";
    	return -1;
    }

    const char* type = OPS_GetString();

    // create ctest
    ConvergenceTest* theTest = 0;
    if (strcmp(type,"NormDispAndUnbalance") == 0) {
	theTest = (ConvergenceTest*)OPS_NormDispAndUnbalance();

    } else if (strcmp(type,"NormDispOrUnbalance") == 0) {
	theTest = (ConvergenceTest*)OPS_NormDispOrUnbalance();

    } else if (strcmp(type,"PFEM") == 0) {
	theTest = (ConvergenceTest*)OPS_CTestPFEM();

    } else if (strcmp(type,"FixedNumIter") == 0) {
	theTest = (ConvergenceTest*)OPS_CTestFixedNumIter();

    } else if (strcmp(type,"NormUnbalance") == 0) {
	theTest = (ConvergenceTest*)OPS_CTestNormUnbalance();

    } else if (strcmp(type,"NormDispIncr") == 0) {
	theTest = (ConvergenceTest*)OPS_CTestNormDispIncr();

    } else if (strcmp(type,"EnergyIncr") == 0) {
	theTest = (ConvergenceTest*)OPS_CTestEnergyIncr();

    } else if (strcmp(type,"RelativeNormUnbalance") == 0) {
	theTest = (ConvergenceTest*)OPS_CTestRelativeNormUnbalance();


    } else if (strcmp(type,"RelativeNormDispIncr") == 0) {
	theTest = (ConvergenceTest*)OPS_CTestRelativeNormDispIncr();

    } else if (strcmp(type,"RelativeEnergyIncr") == 0) {
	theTest = (ConvergenceTest*)OPS_CTestRelativeEnergyIncr();

    } else if (strcmp(type,"RelativeTotalNormDispIncr") == 0) {
	theTest = (ConvergenceTest*)OPS_CTestRelativeTotalNormDispIncr();

    } else {

	opserr<<"WARNING unknown CTest type "<<type<<"\n";
    	return -1;
    }

    // set test
    if (cmds != 0) {
	cmds->setCTest(theTest);
    }

    return 0;
}

int OPS_Integrator()
{
    if (OPS_GetNumRemainingInputArgs() < 1) {
    	opserr << "WARNING insufficient args: integrator type ...\n";
    	return -1;
    }

    const char* type = OPS_GetString();

    // create integrator
    StaticIntegrator* si = 0;
    TransientIntegrator* ti = 0;
    if (strcmp(type,"LoadControl") == 0) {

	si = (StaticIntegrator*)OPS_LoadControlIntegrator();

    } else if (strcmp(type,"DisplacementControl") == 0) {

	si = (StaticIntegrator*)OPS_DisplacementControlIntegrator();

    } else if (strcmp(type,"ArcLength") == 0) {
	si = (StaticIntegrator*)OPS_ArcLength();

    } else if (strcmp(type,"ArcLength1") == 0) {
	si = (StaticIntegrator*)OPS_ArcLength1();

    } else if (strcmp(type,"HSConstraint") == 0) {
	si = (StaticIntegrator*)OPS_HSConstraint();

    } else if (strcmp(type,"MinUnbalDispNorm") == 0) {
	si = (StaticIntegrator*)OPS_MinUnbalDispNorm();

    } else if (strcmp(type,"Newmark") == 0) {
	ti = (TransientIntegrator*)OPS_Newmark();

    } else if (strcmp(type,"TRBDF2") == 0 || strcmp(type,"Bathe") == 0) {
	ti = (TransientIntegrator*)OPS_TRBDF2();

    } else if (strcmp(type,"TRBDF3") == 0 || strcmp(type,"Bathe3") == 0) {
	ti = (TransientIntegrator*)OPS_TRBDF3();

    } else if (strcmp(type,"Houbolt") == 0) {
	ti = (TransientIntegrator*)OPS_Houbolt();

    } else if (strcmp(type,"BackwardEuler") == 0) {
	ti = (TransientIntegrator*)OPS_BackwardEuler();

    } else if (strcmp(type,"PFEM") == 0) {
	ti = (TransientIntegrator*)OPS_PFEMIntegrator();

    } else if (strcmp(type,"NewmarkExplicit") == 0) {
	ti = (TransientIntegrator*)OPS_NewmarkExplicit();

    } else if (strcmp(type,"NewmarkHSIncrLimit") == 0) {
	ti = (TransientIntegrator*)OPS_NewmarkHSIncrLimit();

    } else if (strcmp(type,"NewmarkHSIncrReduct") == 0) {
	ti = (TransientIntegrator*)OPS_NewmarkHSIncrReduct();

    } else if (strcmp(type,"NewmarkHSFixedNumIter") == 0) {
	ti = (TransientIntegrator*)OPS_NewmarkHSFixedNumIter();

    } else if (strcmp(type,"HHT") == 0) {
	ti = (TransientIntegrator*)OPS_HHT();

    } else if (strcmp(type,"HHT_TP") == 0) {
	ti = (TransientIntegrator*)OPS_HHT_TP();

    } else if (strcmp(type,"HHTGeneralized") == 0) {
	ti = (TransientIntegrator*)OPS_HHTGeneralized();

    } else if (strcmp(type,"HHTGeneralized_TP") == 0) {
	ti = (TransientIntegrator*)OPS_HHTGeneralized_TP();

    } else if (strcmp(type,"HHTExplicit") == 0) {
	ti = (TransientIntegrator*)OPS_HHTExplicit();

    } else if (strcmp(type,"HHTExplicit_TP") == 0) {
	ti = (TransientIntegrator*)OPS_HHTExplicit_TP();

    } else if (strcmp(type,"HHTGeneralizedExplicit") == 0) {
	ti = (TransientIntegrator*)OPS_HHTGeneralizedExplicit();

    } else if (strcmp(type,"HHTGeneralizedExplicit_TP") == 0) {
	ti = (TransientIntegrator*)OPS_HHTGeneralizedExplicit_TP();

    } else if (strcmp(type,"HHTHSIncrLimit") == 0) {
	ti = (TransientIntegrator*)OPS_HHTHSIncrLimit();

    } else if (strcmp(type,"HHTHSIncrLimit_TP") == 0) {
	ti = (TransientIntegrator*)OPS_HHTHSIncrLimit_TP();

    } else if (strcmp(type,"HHTHSIncrReduct") == 0) {
	ti = (TransientIntegrator*)OPS_HHTHSIncrReduct();

    } else if (strcmp(type,"HHTHSIncrReduct_TP") == 0) {
	ti = (TransientIntegrator*)OPS_HHTHSIncrReduct_TP();

    } else if (strcmp(type,"HHTHSFixedNumIter") == 0) {
	ti = (TransientIntegrator*)OPS_HHTHSIncrReduct();

    } else if (strcmp(type,"HHTHSFixedNumIter_TP") == 0) {
	ti = (TransientIntegrator*)OPS_HHTHSIncrReduct_TP();

    } else if (strcmp(type,"GeneralizedAlpha") == 0) {
	ti = (TransientIntegrator*)OPS_GeneralizedAlpha();

    } else if (strcmp(type,"KRAlphaExplicit") == 0) {
	ti = (TransientIntegrator*)OPS_KRAlphaExplicit();

    } else if (strcmp(type,"KRAlphaExplicit_TP") == 0) {
	ti = (TransientIntegrator*)OPS_KRAlphaExplicit_TP();

    } else if (strcmp(type,"AlphaOS") == 0) {
	ti = (TransientIntegrator*)OPS_AlphaOS();

    } else if (strcmp(type,"AlphaOS_TP") == 0) {
	ti = (TransientIntegrator*)OPS_AlphaOS_TP();

    } else if (strcmp(type,"AlphaOSGeneralized") == 0) {
	ti = (TransientIntegrator*)OPS_AlphaOSGeneralized();

    } else if (strcmp(type,"AlphaOSGeneralized_TP") == 0) {
	ti = (TransientIntegrator*)OPS_AlphaOSGeneralized_TP();

    } else if (strcmp(type,"Collocation") == 0) {
	ti = (TransientIntegrator*)OPS_Collocation();

    } else if (strcmp(type,"CollocationHSIncrReduct") == 0) {
	ti = (TransientIntegrator*)OPS_CollocationHSIncrReduct();

    } else if (strcmp(type,"CollocationHSIncrLimit") == 0) {
	ti = (TransientIntegrator*)OPS_CollocationHSIncrLimit();

    } else if (strcmp(type,"CollocationHSFixedNumIter") == 0) {
	ti = (TransientIntegrator*)OPS_CollocationHSFixedNumIter();

    } else if (strcmp(type,"Newmark1") == 0) {
	ti = (TransientIntegrator*)OPS_Newmark1();

    } else if (strcmp(type,"WilsonTheta") == 0) {
	ti = (TransientIntegrator*)OPS_WilsonTheta();

    } else if (strcmp(type,"CentralDifference") == 0) {
	ti = (TransientIntegrator*)OPS_CentralDifference();

    } else if (strcmp(type,"CentralDifferenceAlternative") == 0) {
	ti = (TransientIntegrator*)OPS_CentralDifferenceAlternative();

    } else if (strcmp(type,"CentralDifferenceNoDamping") == 0) {
	ti = (TransientIntegrator*)OPS_CentralDifferenceNoDamping();

    } else {
	opserr<<"WARNING unknown integrator type "<<type<<"\n";
    }

    // set integrator
    if (si != 0) {
	if (cmds != 0) {
	    cmds->setStaticIntegrator(si);
	}
    } else if (ti != 0) {
	if (cmds != 0) {
	    cmds->setTransientIntegrator(ti);
	}
    }

    return 0;
}

int OPS_Algorithm()
{
    if (OPS_GetNumRemainingInputArgs() < 1) {
    	opserr << "WARNING insufficient args: algorithm type ...\n";
    	return -1;
    }

    const char* type = OPS_GetString();

    // create algorithm
    EquiSolnAlgo* theAlgo = 0;
    if (strcmp(type, "Linear") == 0) {
	theAlgo = (EquiSolnAlgo*) OPS_LinearAlgorithm();

    } else if (strcmp(type, "Newton") == 0) {
	theAlgo = (EquiSolnAlgo*) OPS_NewtonRaphsonAlgorithm();

    } else if (strcmp(type, "ModifiedNewton") == 0) {
	theAlgo = (EquiSolnAlgo*) OPS_ModifiedNewton();

    } else if (strcmp(type, "KrylovNewton") == 0) {
	theAlgo = (EquiSolnAlgo*) OPS_KrylovNewton();

    } else if (strcmp(type, "RaphsonNewton") == 0) {
	theAlgo = (EquiSolnAlgo*) OPS_RaphsonNewton();

    } else if (strcmp(type, "MillerNewton") == 0) {
	theAlgo = (EquiSolnAlgo*) OPS_MillerNewton();

    } else if (strcmp(type, "SecantNewton") == 0) {
	theAlgo = (EquiSolnAlgo*) OPS_SecantNewton();

    } else if (strcmp(type, "PeriodicNewton") == 0) {
	theAlgo = (EquiSolnAlgo*) OPS_PeriodicNewton();


    } else if (strcmp(type, "Broyden") == 0) {
	theAlgo = (EquiSolnAlgo*)OPS_Broyden();

    } else if (strcmp(type, "BFGS") == 0) {
	theAlgo = (EquiSolnAlgo*)OPS_BFGS();

    } else if (strcmp(type, "NewtonLineSearch") == 0) {
	theAlgo = (EquiSolnAlgo*)OPS_NewtonLineSearch();

    } else {
	opserr<<"WARNING unknown algorithm type "<<type<<"\n";
    }

    // set algorithm
    if (theAlgo != 0) {
	if (cmds != 0) {
	    cmds->setAlgorithm(theAlgo);
	}
    }

    return 0;
}

int OPS_Analysis()
{
    if (OPS_GetNumRemainingInputArgs() < 1) {
    	opserr << "WARNING insufficient args: analysis type ...\n";
    	return -1;
    }

    const char* type = OPS_GetString();

    // create analysis
    if (strcmp(type, "Static") == 0) {
	if (cmds != 0) {
	    cmds->setStaticAnalysis();
	}
    } else if (strcmp(type, "Transient") == 0) {
	if (cmds != 0) {
	    cmds->setTransientAnalysis();
	}
    } else if (strcmp(type, "PFEM") == 0) {
	if (cmds != 0) {
	    if (cmds->setPFEMAnalysis() < 0) {
		return -1;
	    }
	}
    } else if (strcmp(type, "VariableTimeStepTransient") == 0 ||
	       (strcmp(type,"TransientWithVariableTimeStep") == 0) ||
	       (strcmp(type,"VariableTransient") == 0)) {
	if (cmds != 0) {
	    cmds->setVariableAnalysis();
	}

    } else {
	opserr<<"WARNING unknown analysis type "<<type<<"\n";
    }

    return 0;
}

int OPS_analyze()
{
    if (cmds == 0) return 0;

    int result = 0;
    StaticAnalysis* theStaticAnalysis = cmds->getStaticAnalysis();
    TransientAnalysis* theTransientAnalysis = cmds->getTransientAnalysis();
    PFEMAnalysis* thePFEMAnalysis = cmds->getPFEMAnalysis();

    if (theStaticAnalysis != 0) {
	if (OPS_GetNumRemainingInputArgs() < 1) {
	    opserr << "WARNING insufficient args: analyze numIncr ...\n";
	    return -1;
	}
	int numIncr;
	int numdata = 1;
	if (OPS_GetIntInput(&numdata, &numIncr) < 0) return -1;
	result = theStaticAnalysis->analyze(numIncr);

    } else if (thePFEMAnalysis != 0) {

	result = thePFEMAnalysis->analyze();

    } else if (theTransientAnalysis != 0) {
	if (OPS_GetNumRemainingInputArgs() < 2) {
	    opserr << "WARNING insufficient args: analyze numIncr deltaT ...\n";
	    return -1;
	}
	int numIncr;
	int numdata = 1;
	if (OPS_GetIntInput(&numdata, &numIncr) < 0) return -1;

	double dt;
	if (OPS_GetDoubleInput(&numdata, &dt) < 0) return -1;
	ops_Dt = dt;

	result = theTransientAnalysis->analyze(numIncr, dt);
    } else {
	opserr << "WARNING No Analysis type has been specified \n";
	return -1;
    }

    if (result < 0) {
	opserr << "OpenSees > analyze failed, returned: " << result << " error flag\n";
    }

    int numdata = 1;
    if (OPS_SetIntOutput(&numdata, &result) < 0) {
	opserr<<"WARNING failed to set output\n";
	return -1;
    }

    return 0;
}

int OPS_eigenAnalysis()
{
    // make sure at least one other argument to contain type of system
    if (OPS_GetNumRemainingInputArgs() < 1) {
	opserr << "WARNING want - eigen <type> numModes?\n";
	return -1;
    }

    // 0 - frequency/generalized (default),1 - standard, 2 - buckling
    bool generalizedAlgo = true;


    int typeSolver = EigenSOE_TAGS_ArpackSOE;
    double shift = 0.0;
    bool findSmallest = true;

    // Check type of eigenvalue analysis
    while (OPS_GetNumRemainingInputArgs() > 1) {

	const char* type = OPS_GetString();

	if ((strcmp(type,"frequency") == 0) ||
	    (strcmp(type,"-frequency") == 0) ||
	    (strcmp(type,"generalized") == 0) ||
	    (strcmp(type,"-generalized") == 0))
	    generalizedAlgo = true;

	else if ((strcmp(type,"standard") == 0) ||
		 (strcmp(type,"-standard") == 0))
	    generalizedAlgo = false;

	else if ((strcmp(type,"-findLargest") == 0))
	    findSmallest = false;

	else if ((strcmp(type,"genBandArpack") == 0) ||
		 (strcmp(type,"-genBandArpack") == 0) ||
		 (strcmp(type,"genBandArpackEigen") == 0) ||
		 (strcmp(type,"-genBandArpackEigen") == 0))
	    typeSolver = EigenSOE_TAGS_ArpackSOE;

	else if ((strcmp(type,"symmBandLapack") == 0) ||
		 (strcmp(type,"-symmBandLapack") == 0) ||
		 (strcmp(type,"symmBandLapackEigen") == 0) ||
		 (strcmp(type,"-symmBandLapackEigen") == 0))
	    typeSolver = EigenSOE_TAGS_SymBandEigenSOE;

	else if ((strcmp(type,"fullGenLapack") == 0) ||
		 (strcmp(type,"-fullGenLapack") == 0) ||
		 (strcmp(type,"fullGenLapackEigen") == 0) ||
		 (strcmp(type,"-fullGenLapackEigen") == 0))
	    typeSolver = EigenSOE_TAGS_FullGenEigenSOE;

	else {
	    opserr << "eigen - unknown option specified " << type << endln;
	}

    }

    // check argv[loc] for number of modes
    int numEigen;
    int numdata = 1;
    if (OPS_GetIntInput(&numdata, &numEigen) < 0) {
	opserr << "WARNING eigen numModes?  - can't read numModes\n";
	return -1;
    }

    if (numEigen < 0) {
	opserr << "WARNING eigen numModes?  - illegal numModes\n";
	return -1;
    }
    cmds->setNumEigen(numEigen);

    // set eigen soe
    if (cmds->eigen(typeSolver,shift,generalizedAlgo,findSmallest) < 0) {
	opserr<<"WANRING failed to do eigen analysis\n";
	return -1;
    }


    return 0;
}

int OPS_resetModel()
{
    if (cmds == 0) return 0;
    Domain* theDomain = OPS_GetDomain();
    if (theDomain != 0) {
	theDomain->revertToStart();
    }
    TransientIntegrator* theTransientIntegrator = cmds->getTransientIntegrator();
    if (theTransientIntegrator != 0) {
	theTransientIntegrator->revertToStart();
    }
    return 0;
}

int OPS_initializeAnalysis()
{
    if (cmds == 0) return 0;
    DirectIntegrationAnalysis* theTransientAnalysis =
	cmds->getTransientAnalysis();

    StaticAnalysis* theStaticAnalysis =
	cmds->getStaticAnalysis();

    if (theTransientAnalysis != 0) {
	theTransientAnalysis->initialize();
    }else if (theStaticAnalysis != 0) {
	theStaticAnalysis->initialize();
    }

    Domain* theDomain = OPS_GetDomain();
    if (theDomain != 0) {
	theDomain->initialize();
    }

    return 0;
}

int OPS_printA()
{
    if (cmds == 0) return 0;
    FileStream outputFile;
    OPS_Stream *output = &opserr;

    bool ret = false;
    if (OPS_GetNumRemainingInputArgs() > 0) {
	const char* flag = OPS_GetString();

	if ((strcmp(flag,"file") == 0) || (strcmp(flag,"-file") == 0)) {

	    const char* filename = OPS_GetString();
	    if (outputFile.setFile(filename) != 0) {
		opserr << "print <filename> .. - failed to open file: " << filename << endln;
		return -1;
	    }
	    output = &outputFile;
	} else if((strcmp(flag,"ret") == 0) || (strcmp(flag,"-ret") == 0)) {
	    ret = true;
	}
    }

    LinearSOE* theSOE = cmds->getSOE();
    StaticIntegrator* theStaticIntegrator = cmds->getStaticIntegrator();
    TransientIntegrator* theTransientIntegrator = cmds->getTransientIntegrator();

    if (theSOE != 0) {
	if (theStaticIntegrator != 0) {
	    theStaticIntegrator->formTangent();
	} else if (theTransientIntegrator != 0) {
	    theTransientIntegrator->formTangent(0);
	}

	Matrix *A = const_cast<Matrix*>(theSOE->getA());
	if (A != 0) {
	    if (ret) {
		int size = A->noRows() * A->noCols();
		if (size >0) {
		    double& ptr = (*A)(0,0);
		    if (OPS_SetDoubleOutput(&size, &ptr) < 0) {
			opserr << "WARNING: printA - failed to set output\n";
			return -1;
		    }
		}
	    } else {
		*output << *A;
	    }
	}
    }

    // close the output file
    outputFile.close();

    return 0;
}

int OPS_printB()
{
    if (cmds == 0) return 0;
    FileStream outputFile;
    OPS_Stream *output = &opserr;

    LinearSOE* theSOE = cmds->getSOE();
    StaticIntegrator* theStaticIntegrator = cmds->getStaticIntegrator();
    TransientIntegrator* theTransientIntegrator = cmds->getTransientIntegrator();

    bool ret = false;
    if (OPS_GetNumRemainingInputArgs() > 0) {
	const char* flag = OPS_GetString();

	if ((strcmp(flag,"file") == 0) || (strcmp(flag,"-file") == 0)) {
	    const char* filename = OPS_GetString();

	    if (outputFile.setFile(filename) != 0) {
		opserr << "print <filename> .. - failed to open file: " << filename << endln;
		return -1;
	    }
	    output = &outputFile;
	} else if((strcmp(flag,"ret") == 0) || (strcmp(flag,"-ret") == 0)) {
	    ret = true;
	}
    }
    if (theSOE != 0) {
	if (theStaticIntegrator != 0) {
	    theStaticIntegrator->formTangent();
	} else if (theTransientIntegrator != 0) {
	    theTransientIntegrator->formTangent(0);
	}

	Vector &b = const_cast<Vector&>(theSOE->getB());
	if (ret) {
	    int size = b.Size();
	    if (size > 0) {
		double &ptr = b(0);
		if (OPS_SetDoubleOutput(&size, &ptr) < 0) {
		    opserr << "WARNING: printb - failed to set output\n";
		    return -1;
		}
	    }
	} else {
	    *output << b;
	}
    }

    // close the output file
    outputFile.close();

    return 0;
}

int printNode(OPS_Stream& output);
int printElement(OPS_Stream& output);
int printAlgorithm(OPS_Stream& output);
int printIntegrator(OPS_Stream& output);

int OPS_printModel()
{
    int res = 0;

    int flag = OPS_PRINT_CURRENTSTATE;

    FileStream outputFile;
    OPS_Stream *output = &opserr;
    bool done = false;

    Domain* theDomain = OPS_GetDomain();
    if (theDomain == 0) return -1;

    // if just 'print' then print out the entire domain
    if (OPS_GetNumRemainingInputArgs() < 1) {
        opserr << *theDomain;
        return 0;
    }

    while (done == false && OPS_GetNumRemainingInputArgs() > 0) {

        const char* arg = OPS_GetString();

        // if 'print ele i j k..' print out some elements
        if ((strcmp(arg, "-ele") == 0) || (strcmp(arg, "ele") == 0)) {
            res = printElement(*output);
            done = true;
        }
        // if 'print node i j k ..' print out some nodes
        else if ((strcmp(arg, "-node") == 0) || (strcmp(arg, "node") == 0)) {
            res = printNode(*output);
            done = true;
        }

        // if 'print integrator flag' print out the integrator
        else if ((strcmp(arg, "integrator") == 0) || (strcmp(arg, "-integrator") == 0)) {
            res = printIntegrator(*output);
            done = true;
        }

        // if 'print algorithm flag' print out the algorithm
        else if ((strcmp(arg, "algorithm") == 0) || (strcmp(arg, "-algorithm") == 0)) {
            res = printAlgorithm(*output);
            done = true;
        }

        // if 'print -JSON' print using JSON format
        else if ((strcmp(arg, "JSON") == 0) || (strcmp(arg, "-JSON") == 0)) {
            flag = OPS_PRINT_PRINTMODEL_JSON;
        }

        else {

            if ((strcmp(arg, "file") == 0) || (strcmp(arg, "-file") == 0)) {}

            if (OPS_GetNumRemainingInputArgs() < 1) break;
            const char* filename = OPS_GetString();

            openMode mode = APPEND;
            if (flag == OPS_PRINT_PRINTMODEL_JSON)
                mode = OVERWRITE;
            if (outputFile.setFile(filename, mode) != 0) {
                opserr << "print <filename> .. - failed to open file: " << filename << endln;
                return -1;
            }

            SimulationInformation* simulationInfo = cmds->getSimulationInformation();
            if (simulationInfo == 0) return -1;

            // if just 'print <filename>' then print out the entire domain to eof
            if (OPS_GetNumRemainingInputArgs() < 1) {
                if (flag == OPS_PRINT_PRINTMODEL_JSON)
                    simulationInfo->Print(outputFile, flag);
                theDomain->Print(outputFile, flag);
                return 0;
            }

            output = &outputFile;

        }
    }

    // close the output file
    outputFile.close();
    return res;
}

void* OPS_KrylovNewton()
{
    if (cmds == 0) return 0;
    int incrementTangent = CURRENT_TANGENT;
    int iterateTangent = CURRENT_TANGENT;
    int maxDim = 3;
    while (OPS_GetNumRemainingInputArgs() > 0) {
	const char* flag = OPS_GetString();

	if (strcmp(flag,"-iterate") == 0 && OPS_GetNumRemainingInputArgs()>0) {
	    const char* flag2 = OPS_GetString();

	    if (strcmp(flag2,"current") == 0) {
		iterateTangent = CURRENT_TANGENT;
	    }
	    if (strcmp(flag2,"initial") == 0) {
		iterateTangent = INITIAL_TANGENT;
	    }
	    if (strcmp(flag2,"noTangent") == 0) {
		iterateTangent = NO_TANGENT;
	    }
	} else if (strcmp(flag,"-increment") == 0 && OPS_GetNumRemainingInputArgs()>0) {
	    const char* flag2 = OPS_GetString();

	    if (strcmp(flag2,"current") == 0) {
		incrementTangent = CURRENT_TANGENT;
	    }
	    if (strcmp(flag2,"initial") == 0) {
		incrementTangent = INITIAL_TANGENT;
	    }
	    if (strcmp(flag2,"noTangent") == 0) {
		incrementTangent = NO_TANGENT;
	    }
	} else if (strcmp(flag,"-maxDim") == 0 && OPS_GetNumRemainingInputArgs()>0) {

	    maxDim = atoi(flag);
	    int numdata = 1;
	    if (OPS_GetIntInput(&numdata, &maxDim) < 0) {
		opserr<< "WARNING KrylovNewton failed to read maxDim\n";
		return 0;
	    }
	}
    }

    ConvergenceTest* theTest = cmds->getCTest();
    if (theTest == 0) {
      opserr << "ERROR: No ConvergenceTest yet specified\n";
      return 0;
    }

    Accelerator *theAccel;
    theAccel = new KrylovAccelerator(maxDim, iterateTangent);

    return new AcceleratedNewton(*theTest, theAccel, incrementTangent);
}

void* OPS_RaphsonNewton()
{
    if (cmds == 0) return 0;
    int incrementTangent = CURRENT_TANGENT;
    int iterateTangent = CURRENT_TANGENT;

    while (OPS_GetNumRemainingInputArgs() > 0) {
	const char* flag = OPS_GetString();

	if (strcmp(flag,"-iterate") == 0 && OPS_GetNumRemainingInputArgs()>0) {
	    const char* flag2 = OPS_GetString();

	    if (strcmp(flag2,"current") == 0) {
		iterateTangent = CURRENT_TANGENT;
	    }
	    if (strcmp(flag2,"initial") == 0) {
		iterateTangent = INITIAL_TANGENT;
	    }
	    if (strcmp(flag2,"noTangent") == 0) {
		iterateTangent = NO_TANGENT;
	    }
	} else if (strcmp(flag,"-increment") == 0 && OPS_GetNumRemainingInputArgs()>0) {
	    const char* flag2 = OPS_GetString();

	    if (strcmp(flag2,"current") == 0) {
		incrementTangent = CURRENT_TANGENT;
	    }
	    if (strcmp(flag2,"initial") == 0) {
		incrementTangent = INITIAL_TANGENT;
	    }
	    if (strcmp(flag2,"noTangent") == 0) {
		incrementTangent = NO_TANGENT;
	    }
	}
    }

    ConvergenceTest* theTest = cmds->getCTest();
    if (theTest == 0) {
      opserr << "ERROR: No ConvergenceTest yet specified\n";
      return 0;
    }

    Accelerator *theAccel;
    theAccel = new RaphsonAccelerator(iterateTangent);

    return new AcceleratedNewton(*theTest, theAccel, incrementTangent);
}

void* OPS_MillerNewton()
{
    if (cmds == 0) return 0;
    int incrementTangent = CURRENT_TANGENT;
    int iterateTangent = CURRENT_TANGENT;
    int maxDim = 3;
    while (OPS_GetNumRemainingInputArgs() > 0) {
	const char* flag = OPS_GetString();

	if (strcmp(flag,"-iterate") == 0 && OPS_GetNumRemainingInputArgs()>0) {
	    const char* flag2 = OPS_GetString();

	    if (strcmp(flag2,"current") == 0) {
		iterateTangent = CURRENT_TANGENT;
	    }
	    if (strcmp(flag2,"initial") == 0) {
		iterateTangent = INITIAL_TANGENT;
	    }
	    if (strcmp(flag2,"noTangent") == 0) {
		iterateTangent = NO_TANGENT;
	    }
	} else if (strcmp(flag,"-increment") == 0 && OPS_GetNumRemainingInputArgs()>0) {
	    const char* flag2 = OPS_GetString();

	    if (strcmp(flag2,"current") == 0) {
		incrementTangent = CURRENT_TANGENT;
	    }
	    if (strcmp(flag2,"initial") == 0) {
		incrementTangent = INITIAL_TANGENT;
	    }
	    if (strcmp(flag2,"noTangent") == 0) {
		incrementTangent = NO_TANGENT;
	    }
	} else if (strcmp(flag,"-maxDim") == 0 && OPS_GetNumRemainingInputArgs()>0) {

	    maxDim = atoi(flag);
	    int numdata = 1;
	    if (OPS_GetIntInput(&numdata, &maxDim) < 0) {
		opserr<< "WARNING KrylovNewton failed to read maxDim\n";
		return 0;
	    }
	}
    }

    ConvergenceTest* theTest = cmds->getCTest();
    if (theTest == 0) {
      opserr << "ERROR: No ConvergenceTest yet specified\n";
      return 0;
    }

    Accelerator *theAccel = 0;
    return new AcceleratedNewton(*theTest, theAccel, incrementTangent);
}

void* OPS_SecantNewton()
{
    if (cmds == 0) return 0;
    int incrementTangent = CURRENT_TANGENT;
    int iterateTangent = CURRENT_TANGENT;
    int maxDim = 3;
    while (OPS_GetNumRemainingInputArgs() > 0) {
	const char* flag = OPS_GetString();

	if (strcmp(flag,"-iterate") == 0 && OPS_GetNumRemainingInputArgs()>0) {
	    const char* flag2 = OPS_GetString();

	    if (strcmp(flag2,"current") == 0) {
		iterateTangent = CURRENT_TANGENT;
	    }
	    if (strcmp(flag2,"initial") == 0) {
		iterateTangent = INITIAL_TANGENT;
	    }
	    if (strcmp(flag2,"noTangent") == 0) {
		iterateTangent = NO_TANGENT;
	    }
	} else if (strcmp(flag,"-increment") == 0 && OPS_GetNumRemainingInputArgs()>0) {
	    const char* flag2 = OPS_GetString();

	    if (strcmp(flag2,"current") == 0) {
		incrementTangent = CURRENT_TANGENT;
	    }
	    if (strcmp(flag2,"initial") == 0) {
		incrementTangent = INITIAL_TANGENT;
	    }
	    if (strcmp(flag2,"noTangent") == 0) {
		incrementTangent = NO_TANGENT;
	    }
	} else if (strcmp(flag,"-maxDim") == 0 && OPS_GetNumRemainingInputArgs()>0) {

	    maxDim = atoi(flag);
	    int numdata = 1;
	    if (OPS_GetIntInput(&numdata, &maxDim) < 0) {
		opserr<< "WARNING KrylovNewton failed to read maxDim\n";
		return 0;
	    }
	}
    }

    ConvergenceTest* theTest = cmds->getCTest();
    if (theTest == 0) {
      opserr << "ERROR: No ConvergenceTest yet specified\n";
      return 0;
    }

    Accelerator *theAccel;
    theAccel = new SecantAccelerator2(maxDim, iterateTangent);

    return new AcceleratedNewton(*theTest, theAccel, incrementTangent);
}

void* OPS_PeriodicNewton()
{
    if (cmds == 0) return 0;
    int incrementTangent = CURRENT_TANGENT;
    int iterateTangent = CURRENT_TANGENT;
    int maxDim = 3;
    while (OPS_GetNumRemainingInputArgs() > 0) {
	const char* flag = OPS_GetString();

	if (strcmp(flag,"-iterate") == 0 && OPS_GetNumRemainingInputArgs()>0) {
	    const char* flag2 = OPS_GetString();

	    if (strcmp(flag2,"current") == 0) {
		iterateTangent = CURRENT_TANGENT;
	    }
	    if (strcmp(flag2,"initial") == 0) {
		iterateTangent = INITIAL_TANGENT;
	    }
	    if (strcmp(flag2,"noTangent") == 0) {
		iterateTangent = NO_TANGENT;
	    }
	} else if (strcmp(flag,"-increment") == 0 && OPS_GetNumRemainingInputArgs()>0) {
	    const char* flag2 = OPS_GetString();

	    if (strcmp(flag2,"current") == 0) {
		incrementTangent = CURRENT_TANGENT;
	    }
	    if (strcmp(flag2,"initial") == 0) {
		incrementTangent = INITIAL_TANGENT;
	    }
	    if (strcmp(flag2,"noTangent") == 0) {
		incrementTangent = NO_TANGENT;
	    }
	} else if (strcmp(flag,"-maxDim") == 0 && OPS_GetNumRemainingInputArgs()>0) {

	    maxDim = atoi(flag);
	    int numdata = 1;
	    if (OPS_GetIntInput(&numdata, &maxDim) < 0) {
		opserr<< "WARNING KrylovNewton failed to read maxDim\n";
		return 0;
	    }
	}
    }

    ConvergenceTest* theTest = cmds->getCTest();
    if (theTest == 0) {
      opserr << "ERROR: No ConvergenceTest yet specified\n";
      return 0;
    }

    Accelerator *theAccel;
    theAccel = new PeriodicAccelerator(maxDim, iterateTangent);

    return new AcceleratedNewton(*theTest, theAccel, incrementTangent);
}

void* OPS_NewtonLineSearch()
{
    if (cmds == 0) return 0;
    ConvergenceTest* theTest = cmds->getCTest();

    if (theTest == 0) {
	opserr << "ERROR: No ConvergenceTest yet specified\n";
	return 0;
    }

    // set some default variable
    double tol        = 0.8;
    int    maxIter    = 10;
    double maxEta     = 10.0;
    double minEta     = 0.1;
    int    pFlag      = 1;
    int    typeSearch = 0;

    int numdata = 1;

    while (OPS_GetNumRemainingInputArgs() > 0) {
	const char* flag = OPS_GetString();

	if (strcmp(flag, "-tol") == 0 && OPS_GetNumRemainingInputArgs()>0) {

	    if (OPS_GetDoubleInput(&numdata, &tol) < 0) {
		opserr << "WARNING NewtonLineSearch failed to read tol\n";
		return 0;
	    }

	} else if (strcmp(flag, "-maxIter") == 0 && OPS_GetNumRemainingInputArgs()>0) {

	    if (OPS_GetIntInput(&numdata, &maxIter) < 0) {
		opserr << "WARNING NewtonLineSearch failed to read maxIter\n";
		return 0;
	    }

	} else if (strcmp(flag, "-pFlag") == 0 && OPS_GetNumRemainingInputArgs()>0) {

	    if (OPS_GetIntInput(&numdata, &pFlag) < 0) {
		opserr << "WARNING NewtonLineSearch failed to read pFlag\n";
		return 0;
	    }

	} else if (strcmp(flag, "-minEta") == 0 && OPS_GetNumRemainingInputArgs()>0) {

	    if (OPS_GetDoubleInput(&numdata, &minEta) < 0) {
		opserr << "WARNING NewtonLineSearch failed to read minEta\n";
		return 0;
	    }

	} else if (strcmp(flag, "-maxEta") == 0 && OPS_GetNumRemainingInputArgs()>0) {

	    if (OPS_GetDoubleInput(&numdata, &maxEta) < 0) {
		opserr << "WARNING NewtonLineSearch failed to read maxEta\n";
		return 0;
	    }

	} else if (strcmp(flag, "-type") == 0 && OPS_GetNumRemainingInputArgs()>0) {
	    const char* flag2 = OPS_GetString();

	    if (strcmp(flag2, "Bisection") == 0)
		typeSearch = 1;
	    else if (strcmp(flag2, "Secant") == 0)
		typeSearch = 2;
	    else if (strcmp(flag2, "RegulaFalsi") == 0)
		typeSearch = 3;
	    else if (strcmp(flag2, "LinearInterpolated") == 0)
		typeSearch = 3;
	    else if (strcmp(flag2, "InitialInterpolated") == 0)
		typeSearch = 0;
	}
    }

    LineSearch *theLineSearch = 0;
    if (typeSearch == 0)
	theLineSearch = new InitialInterpolatedLineSearch(tol, maxIter, minEta, maxEta, pFlag);

    else if (typeSearch == 1)
	theLineSearch = new BisectionLineSearch(tol, maxIter, minEta, maxEta, pFlag);
    else if (typeSearch == 2)
	theLineSearch = new SecantLineSearch(tol, maxIter, minEta, maxEta, pFlag);
    else if (typeSearch == 3)
	theLineSearch = new RegulaFalsiLineSearch(tol, maxIter, minEta, maxEta, pFlag);

    return new NewtonLineSearch(*theTest, theLineSearch);
}

int OPS_getCTestNorms()
{
    if (cmds == 0) return 0;
    ConvergenceTest* theTest = cmds->getCTest();

    if (theTest != 0) {
	const Vector &norms = theTest->getNorms();
	int numdata = norms.Size();
	double* data = new double[numdata];

	for (int i=0; i<numdata; i++) {
	    data[i] = norms(i);
	}

	if (OPS_SetDoubleOutput(&numdata, data) < 0) {
	    opserr << "WARNING failed to set test norms\n";
	    delete [] data;
	    return -1;
	}
	delete [] data;
	return 0;
    }

    opserr << "ERROR testNorms - no convergence test!\n";
    return -1;
}

int OPS_getCTestIter()
{
    if (cmds == 0) return 0;
    ConvergenceTest* theTest = cmds->getCTest();

    if (theTest != 0) {
	int res = theTest->getNumTests();
	int numdata = 1;
	if (OPS_SetIntOutput(&numdata, &res) < 0) {
	    opserr << "WARNING failed to set test iter\n";
	    return -1;
	}

	return 0;
    }

    opserr << "ERROR testIter - no convergence test!\n";
    return -1;
}

int OPS_Database()
{
    if (cmds == 0) return 0;
    // make sure at least one other argument to contain integrator
    if (OPS_GetNumRemainingInputArgs() < 1) {
	opserr << "WARNING need to specify a Database type; valid type File, MySQL, BerkeleyDB \n";
	return -1;
    }

    //
    // check argv[1] for type of Database, parse in rest of arguments
    // needed for the type of Database, create the object and add to Domain
    //

    // a File Database
    const char* type = OPS_GetString();
    if (strcmp(type,"File") == 0) {
	if (OPS_GetNumRemainingInputArgs() < 1) {
	    opserr << "WARNING database File fileName? ";
	    return -1;
	}

	const char* filename = OPS_GetString();
	cmds->setFileDatabase(filename);

	return 0;
    }
    opserr << "WARNING No database type exists ";
    opserr << "for database of type:" << type << "valid database type File\n";

    return -1;
}

int OPS_save()
{
    if (cmds == 0) return 0;
    // make sure at least one other argument to contain type of system
    if (OPS_GetNumRemainingInputArgs() < 1) {
	opserr << "WARNING save no commit tag - want save commitTag?";
	return -1;
    }

    // check argv[1] for commitTag
    int commitTag;
    int numdata = 1;
    if (OPS_GetIntInput(&numdata, &commitTag) < 0) {
	opserr << "WARNING - save could not read commitTag " << endln;
	return -1;
    }

    FE_Datastore* theDatabase = cmds->getDatabase();
    if (theDatabase == 0) {
	opserr << "WARNING: save - no database has been constructed\n";
	return -1;
    }

    if (theDatabase->commitState(commitTag) < 0) {
	opserr << "WARNING - database failed to commitState \n";
	return -1;
    }

    return 0;
}

int OPS_restore()
{
    if (cmds == 0) return 0;
    // make sure at least one other argument to contain type of system
    if (OPS_GetNumRemainingInputArgs() < 1) {
	opserr << "WARNING restore no commit tag - want restore commitTag?";
	return -1;
    }

    // check argv[1] for commitTag
    int commitTag;
    int numdata = 1;
    if (OPS_GetIntInput(&numdata, &commitTag) < 0) {
	opserr << "WARNING - restore could not read commitTag " << endln;
	return -1;
    }

    FE_Datastore* theDatabase = cmds->getDatabase();
    if (theDatabase == 0) {
	opserr << "WARNING: save - no database has been constructed\n";
	return -1;
    }

    if (theDatabase->restoreState(commitTag) < 0) {
	opserr << "WARNING - database failed to restore state \n";
	return -1;
    }

    return 0;
}

int OPS_startTimer()
{
    if (cmds == 0) return 0;
    Timer* timer = cmds->getTimer();
    if (timer == 0) return -1;
    timer->start();
    return 0;
}

int OPS_stopTimer()
{
    if (cmds == 0) return 0;
    Timer* theTimer = cmds->getTimer();
    if (theTimer == 0) return -1;
    theTimer->pause();
    opserr << *theTimer;
    return 0;
}

int OPS_modalDamping()
{
    if (cmds == 0) return 0;
    if (OPS_GetNumRemainingInputArgs() < 1) {
	opserr << "WARNING modalDamping ?factor - not enough arguments to command\n";
	return -1;
    }

    int numEigen = cmds->getNumEigen();
    EigenSOE* theEigenSOE = cmds->getEigenSOE();

    if (numEigen == 0 || theEigenSOE == 0) {
	opserr << "WARINING - modalDmping - eigen command needs to be called first - NO MODAL DAMPING APPLIED\n ";
	return -1;
    }

    double factor;
    int numdata = 1;
    if (OPS_GetDoubleInput(&numdata, &factor) < 0) {
	opserr << "WARNING rayleigh alphaM? betaK? betaK0? betaKc? - could not read betaK? \n";
	return -1;
    }

    Vector modalDampingValues(numEigen);
    for (int i=0; i<numEigen; i++) {
	modalDampingValues(i) = factor;
    }

    Domain* theDomain = OPS_GetDomain();
    if (theDomain != 0) {
	theDomain->setModalDampingFactors(&modalDampingValues, true);
    }

    return 0;
}

int OPS_modalDampingQ()
{
    if (cmds == 0) return 0;
    if (OPS_GetNumRemainingInputArgs() < 1) {
	opserr << "WARNING modalDamping ?factor - not enough arguments to command\n";
	return -1;
    }

    int numEigen = cmds->getNumEigen();
    EigenSOE* theEigenSOE = cmds->getEigenSOE();

    if (numEigen == 0 || theEigenSOE == 0) {
	opserr << "WARINING - modalDmping - eigen command needs to be called first - NO MODAL DAMPING APPLIED\n ";
	return -1;
    }

    double factor;
    int numdata = 1;
    if (OPS_GetDoubleInput(&numdata, &factor) < 0) {
	opserr << "WARNING rayleigh alphaM? betaK? betaK0? betaKc? - could not read betaK? \n";
	return -1;
    }

    Vector modalDampingValues(numEigen);
    for (int i=0; i<numEigen; i++) {
	modalDampingValues(i) = factor;
    }

    Domain* theDomain = OPS_GetDomain();
    if (theDomain != 0) {
	theDomain->setModalDampingFactors(&modalDampingValues, false);
    }

    return 0;
}

int OPS_neesMetaData()
{
    if (cmds == 0) return 0;
    if (OPS_GetNumRemainingInputArgs() < 1) {
	opserr << "WARNING missing args \n";
	return -1;
    }

    SimulationInformation* simulationInfo = cmds->getSimulationInformation();
    if (simulationInfo == 0) return -1;

    while (OPS_GetNumRemainingInputArgs() > 0) {
	const char* flag = OPS_GetString();

	if ((strcmp(flag,"-title") == 0) || (strcmp(flag,"-Title") == 0)
	    || (strcmp(flag,"-TITLE") == 0)) {
	    if (OPS_GetNumRemainingInputArgs() > 0) {
		simulationInfo->setTitle(OPS_GetString());
	    }
	} else if ((strcmp(flag,"-contact") == 0) || (strcmp(flag,"-Contact") == 0)
		   || (strcmp(flag,"-CONTACT") == 0)) {
	    if (OPS_GetNumRemainingInputArgs() > 0) {
		simulationInfo->setContact(OPS_GetString());
	    }
	} else if ((strcmp(flag,"-description") == 0) || (strcmp(flag,"-Description") == 0)
		   || (strcmp(flag,"-DESCRIPTION") == 0)) {
	    if (OPS_GetNumRemainingInputArgs() > 0) {
		simulationInfo->setDescription(OPS_GetString());
	    }
	} else if ((strcmp(flag,"-modelType") == 0) || (strcmp(flag,"-ModelType") == 0)
		   || (strcmp(flag,"-MODELTYPE") == 0)) {
	    if (OPS_GetNumRemainingInputArgs() > 0) {
		simulationInfo->addModelType(OPS_GetString());
	    }
	} else if ((strcmp(flag,"-analysisType") == 0) || (strcmp(flag,"-AnalysisType") == 0)
		   || (strcmp(flag,"-ANALYSISTYPE") == 0)) {
	    if (OPS_GetNumRemainingInputArgs() > 0) {
		simulationInfo->addAnalysisType(OPS_GetString());
	    }
	} else if ((strcmp(flag,"-elementType") == 0) || (strcmp(flag,"-ElementType") == 0)
		   || (strcmp(flag,"-ELEMENTTYPE") == 0)) {
	    if (OPS_GetNumRemainingInputArgs() > 0) {
		simulationInfo->addElementType(OPS_GetString());
	    }
	} else if ((strcmp(flag,"-materialType") == 0) || (strcmp(flag,"-MaterialType") == 0)
		   || (strcmp(flag,"-MATERIALTYPE") == 0)) {
	    if (OPS_GetNumRemainingInputArgs() > 0) {
		simulationInfo->addMaterialType(OPS_GetString());
	    }
	} else if ((strcmp(flag,"-loadingType") == 0) || (strcmp(flag,"-LoadingType") == 0)
		   || (strcmp(flag,"-LOADINGTYPE") == 0)) {
	    if (OPS_GetNumRemainingInputArgs() > 0) {
		simulationInfo->addLoadingType(OPS_GetString());
	    }
	} else {
	    opserr << "WARNING unknown arg type: " << flag << endln;
	    return -1;
	}
    }
    return 0;
}

int OPS_neesUpload()
{
    if (cmds == 0) return 0;
    if (OPS_GetNumRemainingInputArgs() < 2) {
	opserr << "WARNING neesUpload -user isername? -pass passwd? -proj projID? -exp expID? -title title? -description description\n";
	return -1;
    }
    int projID =0;
    int expID =0;
    const char *userName =0;
    const char *userPasswd =0;

    SimulationInformation* simulationInfo = cmds->getSimulationInformation();
    if (simulationInfo == 0) return -1;

    int numdata = 1;
    while (OPS_GetNumRemainingInputArgs() > 1) {
	const char* flag = OPS_GetString();

	if (strcmp(flag,"-user") == 0) {
	    userName = OPS_GetString();

	} else if (strcmp(flag,"-pass") == 0) {
	    userPasswd = OPS_GetString();

	} else if (strcmp(flag,"-projID") == 0) {
	    if (OPS_GetIntInput(&numdata, &projID) < 0) {
		opserr << "WARNING neesUpload -invalid expID\n";
		return -1;
	    }

	} else if (strcmp(flag,"-expID") == 0) {
	    if (OPS_GetIntInput(&numdata, &expID) < 0) {
		opserr << "WARNING neesUpload -invalid expID\n";
		return -1;
	    }

	} else if (strcmp(flag,"-title") == 0) {
	    simulationInfo->setTitle(OPS_GetString());

	} else if (strcmp(flag,"-description") == 0) {
	    simulationInfo->setDescription(OPS_GetString());

	}
    }

    simulationInfo->neesUpload(userName, userPasswd, projID, expID);

    return 0;
}

int OPS_defaultUnits()
{
    if (OPS_GetNumRemainingInputArgs() < 8) {
        opserr << "WARNING defaultUnits - missing a unit type want: defaultUnits -Force type? -Length type? -Time type? -Temperature type?\n";
        return -1;
    }

    const char *force = 0;
    const char *length = 0;
    const char *time = 0;
    const char *temperature = 0;

    while (OPS_GetNumRemainingInputArgs() > 0) {
        const char* unitType = OPS_GetString();

        if ((strcmp(unitType, "-force") == 0) || (strcmp(unitType, "-Force") == 0)
            || (strcmp(unitType, "-FORCE") == 0)) {
            force = OPS_GetString();
        }
        else if ((strcmp(unitType, "-length") == 0) || (strcmp(unitType, "-Length") == 0)
            || (strcmp(unitType, "-LENGTH") == 0)) {
            length = OPS_GetString();
        }
        else if ((strcmp(unitType, "-time") == 0) || (strcmp(unitType, "-Time") == 0)
            || (strcmp(unitType, "-TIME") == 0)) {
            time = OPS_GetString();
        }
        else if ((strcmp(unitType, "-temperature") == 0) || (strcmp(unitType, "-Temperature") == 0)
            || (strcmp(unitType, "-TEMPERATURE") == 0) || (strcmp(unitType, "-temp") == 0)
            || (strcmp(unitType, "-Temp") == 0) || (strcmp(unitType, "-TEMP") == 0)) {
            temperature = OPS_GetString();
        }
        else {
            opserr << "WARNING defaultUnits - unrecognized unit: " << unitType << " want: defaultUnits -Force type? -Length type? -Time type? -Temperature type?\n";
            return -1;
        }
    }

    if (length == 0 || force == 0 || time == 0 || temperature == 0) {
        opserr << "defaultUnits - missing a unit type want: defaultUnits -Force type? -Length type? -Time type? -Temperature type?\n";
        return -1;
    }

    double lb, kip, n, kn, mn, kgf, tonf;
    double in, ft, mm, cm, m;
    double sec, msec;
    double F, C;

    if ((strcmp(force, "lb") == 0) || (strcmp(force, "lbs") == 0)) {
        lb = 1.0;
    }
    else if ((strcmp(force, "kip") == 0) || (strcmp(force, "kips") == 0)) {
        lb = 0.001;
    }
    else if ((strcmp(force, "N") == 0)) {
        lb = 4.4482216152605;
    }
    else if ((strcmp(force, "kN") == 0) || (strcmp(force, "KN") == 0) || (strcmp(force, "kn") == 0)) {
        lb = 0.0044482216152605;
    }
    else if ((strcmp(force, "mN") == 0) || (strcmp(force, "MN") == 0) || (strcmp(force, "mn") == 0)) {
        lb = 0.0000044482216152605;
    }
    else if ((strcmp(force, "kgf") == 0)) {
        lb = 9.80665*4.4482216152605;
    }
    else if ((strcmp(force, "tonf") == 0)) {
        lb = 9.80665 / 1000.0*4.4482216152605;
    }
    else {
        lb = 1.0;
        opserr << "defaultUnits - unknown force type, valid options: lb, kip, N, kN, MN, kgf, tonf\n";
        return -1;
    }

    if ((strcmp(length, "in") == 0) || (strcmp(length, "inch") == 0)) {
        in = 1.0;
    }
    else if ((strcmp(length, "ft") == 0) || (strcmp(length, "feet") == 0)) {
        in = 1.0 / 12.0;
    }
    else if ((strcmp(length, "mm") == 0)) {
        in = 25.4;
    }
    else if ((strcmp(length, "cm") == 0)) {
        in = 2.54;
    }
    else if ((strcmp(length, "m") == 0)) {
        in = 0.0254;
    }
    else {
        in = 1.0;
        opserr << "defaultUnits - unknown length type, valid options: in, ft, mm, cm, m\n";
        return -1;
    }

    if ((strcmp(time, "sec") == 0) || (strcmp(time, "Sec") == 0)) {
        sec = 1.0;
    }
    else if ((strcmp(time, "msec") == 0) || (strcmp(time, "mSec") == 0)) {
        sec = 1000.0;
    }
    else {
        sec = 1.0;
        opserr << "defaultUnits - unknown time type, valid options: sec, msec\n";
        return -1;
    }

    if ((strcmp(temperature, "F") == 0) || (strcmp(temperature, "degF") == 0)) {
        F = 1.0;
    }
    else if ((strcmp(temperature, "C") == 0) || (strcmp(temperature, "degC") == 0)) {
        F = 9.0 / 5.0 + 32.0;
    }
    else {
        F = 1.0;
        opserr << "defaultUnits - unknown temperature type, valid options: F, C\n";
        return -1;
    }

    kip = lb / 0.001;
    n = lb / 4.4482216152605;
    kn = lb / 0.0044482216152605;
    mn = lb / 0.0000044482216152605;
    kgf = lb / (9.80665*4.4482216152605);
    tonf = lb / (9.80665 / 1000.0*4.4482216152605);

    ft = in * 12.0;
    mm = in / 25.44;
    cm = in / 2.54;
    m = in / 0.0254;

    msec = sec * 0.001;

    C = (F - 32.0)*5.0 / 9.0;

    char string[50];
    DL_Interpreter* theInter = cmds->getInterpreter();
    if (theInter == 0) return -1;

    sprintf(string, "lb = %.18e", lb);   theInter->runCommand(string);
    sprintf(string, "lbf = %.18e", lb);   theInter->runCommand(string);
    sprintf(string, "kip = %.18e", kip);   theInter->runCommand(string);
    sprintf(string, "N = %.18e", n);   theInter->runCommand(string);
    sprintf(string, "kN = %.18e", kn);   theInter->runCommand(string);
    sprintf(string, "Newton = %.18e", n);   theInter->runCommand(string);
    sprintf(string, "kNewton = %.18e", kn);   theInter->runCommand(string);
    sprintf(string, "MN = %.18e", mn);   theInter->runCommand(string);
    sprintf(string, "kgf = %.18e", kgf);   theInter->runCommand(string);
    sprintf(string, "tonf = %.18e", tonf);   theInter->runCommand(string);

    //sprintf(string, "in = %.18e", in);   theInter->runCommand(string);  // "in" is a keyword in Python
    sprintf(string, "inch = %.18e", in);   theInter->runCommand(string);
    sprintf(string, "ft = %.18e", ft);   theInter->runCommand(string);
    sprintf(string, "mm = %.18e", mm);   theInter->runCommand(string);
    sprintf(string, "cm = %.18e", cm);   theInter->runCommand(string);
    sprintf(string, "m = %.18e", m);   theInter->runCommand(string);
    sprintf(string, "meter = %.18e", m);   theInter->runCommand(string);

    sprintf(string, "sec = %.18e", sec);   theInter->runCommand(string);
    sprintf(string, "msec = %.18e", msec);   theInter->runCommand(string);

    sprintf(string, "F = %.18e", F);   theInter->runCommand(string);
    sprintf(string, "degF = %.18e", F);   theInter->runCommand(string);
    sprintf(string, "C = %.18e", C);   theInter->runCommand(string);
    sprintf(string, "degC = %.18e", C);   theInter->runCommand(string);

    double g = 32.174049*ft / (sec*sec);
    sprintf(string, "g = %.18e", g);   theInter->runCommand(string);
    sprintf(string, "kg = %.18e", n*sec*sec / m);   theInter->runCommand(string);
    sprintf(string, "Mg = %.18e", 1e3*n*sec*sec / m);   theInter->runCommand(string);
    sprintf(string, "slug = %.18e", lb*sec*sec / ft);   theInter->runCommand(string);
    sprintf(string, "Pa = %.18e", n / (m*m));   theInter->runCommand(string);
    sprintf(string, "kPa = %.18e", 1e3*n / (m*m));   theInter->runCommand(string);
    sprintf(string, "MPa = %.18e", 1e6*n / (m*m));   theInter->runCommand(string);
    sprintf(string, "psi = %.18e", lb / (in*in));   theInter->runCommand(string);
    sprintf(string, "ksi = %.18e", kip / (in*in));   theInter->runCommand(string);
    sprintf(string, "psf = %.18e", lb / (ft*ft));   theInter->runCommand(string);
    sprintf(string, "ksf = %.18e", kip / (ft*ft));   theInter->runCommand(string);
    sprintf(string, "pcf = %.18e", lb / (ft*ft*ft));   theInter->runCommand(string);
    sprintf(string, "in2 = %.18e", in*in);   theInter->runCommand(string);
    sprintf(string, "ft2 = %.18e", ft*ft);   theInter->runCommand(string);
    sprintf(string, "mm2 = %.18e", mm*mm);   theInter->runCommand(string);
    sprintf(string, "cm2 = %.18e", cm*cm);   theInter->runCommand(string);
    sprintf(string, "m2 = %.18e", m*m);   theInter->runCommand(string);
    sprintf(string, "in4 = %.18e", in*in*in*in);   theInter->runCommand(string);
    sprintf(string, "ft4 = %.18e", ft*ft*ft*ft);   theInter->runCommand(string);
    sprintf(string, "mm4 = %.18e", mm*mm*mm*mm);   theInter->runCommand(string);
    sprintf(string, "cm4 = %.18e", cm*cm*cm*cm);   theInter->runCommand(string);
    sprintf(string, "m4 = %.18e", m*m*m*m);   theInter->runCommand(string);
    sprintf(string, "pi = %.18e", 2.0*asin(1.0));   theInter->runCommand(string);
    sprintf(string, "PI = %.18e", 2.0*asin(1.0));   theInter->runCommand(string);

    SimulationInformation* simulationInfo = cmds->getSimulationInformation();
    if (simulationInfo == 0) return -1;

    simulationInfo->setForceUnit(force);
    simulationInfo->setLengthUnit(length);
    simulationInfo->setTimeUnit(time);
    simulationInfo->setTemperatureUnit(temperature);

    return 0;
}

int OPS_totalCPU()
{
    if (cmds == 0) return 0;
    EquiSolnAlgo* theAlgorithm = cmds->getAlgorithm();
    if (theAlgorithm == 0) {
	opserr << "WARNING no algorithm is set\n";
	return -1;
    }

    double value = theAlgorithm->getTotalTimeCPU();
    int numdata = 1;
    if (OPS_SetDoubleOutput(&numdata, &value) < 0) {
	opserr << "WARNING failed to set output\n";
	return -1;
    }

    return 0;
}

int OPS_solveCPU()
{
    if (cmds == 0) return 0;
    EquiSolnAlgo* theAlgorithm = cmds->getAlgorithm();
    if (theAlgorithm == 0) {
	opserr << "WARNING no algorithm is set\n";
	return -1;
    }

    double value = theAlgorithm->getSolveTimeCPU();
    int numdata = 1;
    if (OPS_SetDoubleOutput(&numdata, &value) < 0) {
	opserr << "WARNING failed to set output\n";
	return -1;
    }

    return 0;
}

int OPS_accelCPU()
{
    if (cmds == 0) return 0;
    EquiSolnAlgo* theAlgorithm = cmds->getAlgorithm();
    if (theAlgorithm == 0) {
	opserr << "WARNING no algorithm is set\n";
	return -1;
    }

    double value = theAlgorithm->getAccelTimeCPU();
    int numdata = 1;
    if (OPS_SetDoubleOutput(&numdata, &value) < 0) {
	opserr << "WARNING failed to set output\n";
	return -1;
    }

    return 0;
}

int OPS_numFact()
{
    if (cmds == 0) return 0;
    EquiSolnAlgo* theAlgorithm = cmds->getAlgorithm();
    if (theAlgorithm == 0) {
	opserr << "WARNING no algorithm is set\n";
	return -1;
    }

    double value = theAlgorithm->getNumFactorizations();
    int numdata = 1;
    if (OPS_SetDoubleOutput(&numdata, &value) < 0) {
	opserr << "WARNING failed to set output\n";
	return -1;
    }

    return 0;
}

int OPS_numIter()
{
    if (cmds == 0) return 0;
    EquiSolnAlgo* theAlgorithm = cmds->getAlgorithm();
    if (theAlgorithm == 0) {
	opserr << "WARNING no algorithm is set\n";
	return -1;
    }

    int value = theAlgorithm->getNumIterations();
    int numdata = 1;
    if (OPS_SetIntOutput(&numdata, &value) < 0) {
	opserr << "WARNING failed to set output\n";
	return -1;
    }

    return value;
}

int OPS_systemSize()
{
    if (cmds == 0) return 0;
    LinearSOE* theSOE = cmds->getSOE();
    if (theSOE == 0) {
	opserr << "WARNING no system is set\n";
	return -1;
    }

    double value = theSOE->getNumEqn();
    int numdata = 1;
    if (OPS_SetDoubleOutput(&numdata, &value) < 0) {
	opserr << "WARNING failed to set output\n";
	return -1;
    }

    return 0;
}

// Sensitivity:BEGIN /////////////////////////////////////////////

int OPS_computeGradients()
{
    Integrator* theIntegrator = 0;
    if(cmds->getStaticIntegrator() != 0) {
    	theIntegrator = cmds->getStaticIntegrator();
    } else if(cmds->getTransientIntegrator() != 0) {
    	theIntegrator = cmds->getTransientIntegrator();
    }

    if (theIntegrator == 0) {
    	opserr << "WARNING: No integrator is created\n";
    	return -1;
    }

    if (theIntegrator->computeSensitivities() < 0) {
	opserr << "WARNING: failed to compute sensitivities\n";
	return -1;
    }
    
    return 0;
}

int OPS_sensitivityAlgorithm()
{
    if (cmds == 0) return 0;

    int analysisTypeTag = 1;
    Integrator* theIntegrator = 0;
    if(cmds->getStaticIntegrator() != 0) {
    	theIntegrator = cmds->getStaticIntegrator();
    } else if(cmds->getTransientIntegrator() != 0) {
    	theIntegrator = cmds->getTransientIntegrator();
    }
    
    // 1: compute at each step (default); 2: compute by command; 
    if (OPS_GetNumRemainingInputArgs() < 1) {
    	opserr << "ERROR: Wrong number of parameters to sensitivity algorithm." << "\n";
    	return -1;
    }
    if (theIntegrator == 0) {
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
// Sensitivity:END /////////////////////////////////////////////
