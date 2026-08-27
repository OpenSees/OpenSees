/* ****************************************************************** **
**    OpenSees - Open System for Earthquake Engineering Simulation    **
**          Pacific Earthquake Engineering Research Center            **
**                                                                    **
**                                                                    **
** (C) Copyright 1999, The Regents of the University of California    **
** All Rights Reserved.                                               **
**                                                                    **
** Commercial use of this program without express permission of the   **
** University of California, Berkeley, is strictly prohibited.  See   **
** file 'COPYRIGHT'  in main directory for information on usage and   **
** redistribution,  and for a DISCLAIMER OF ALL WARRANTIES.           **
**                                                                    **
** Developed by:                                                      **
**   Frank McKenna (fmckenna@ce.berkeley.edu)                         **
**   Gregory L. Fenves (fenves@ce.berkeley.edu)                       **
**   Filip C. Filippou (filippou@ce.berkeley.edu)                     **
**                                                                    **
** ****************************************************************** */

// $Revision: 1.3 $
// $Date: 2008-09-16 18:17:57 $
// $Source: /usr/local/cvs/OpenSees/SRC/analysis/algorithm/equiSolnAlgo/AcceleratedNewton.cpp,v $

// Written: MHS
// Created: Oct 2001
//
// Description: This file contains the class definition for 
// AcceleratedNewton.  AcceleratedNewton is a class which uses a Krylov
// subspace accelerator on the modified Newton method.
// The accelerator is described by Carlson and Miller in
// "Design and Application of a 1D GWMFE Code"
// from SIAM Journal of Scientific Computing (Vol. 19, No. 3,
// pp. 728-765, May 1998)

#include <AcceleratedNewton.h>
#include <Accelerator.h>
#include <AnalysisModel.h>
#include <StaticAnalysis.h>
#include <IncrementalIntegrator.h>
#include <LinearSOE.h>
#include <Channel.h>
#include <FEM_ObjectBroker.h>
#include <ConvergenceTest.h>
#include <Matrix.h>
#include <Vector.h>
#include <ID.h>

#include <KrylovAccelerator.h>
#include <RaphsonAccelerator.h>
#include <SecantAccelerator1.h>
#include <SecantAccelerator2.h>
#include <SecantAccelerator3.h>
#include <PeriodicAccelerator.h>

#include <elementAPI.h>
#include <string.h>
#include <stdlib.h>
#include <fstream>

// OPS_* accelerated-Newton parsers (shared Tcl/Python); test set by caller.

static void
parseAcceleratorArgs(int &incrementTangent, int &iterateTangent,
		     int &maxDim, int &factorOnce,
		     int &numTerms, bool &cutOut, double R[2],
		     bool wantMaxDim, bool wantNumTerms, bool wantCutOut)
{
    bool iterateSeen = false;
    bool factorOnceExplicit = false;

    while (OPS_GetNumRemainingInputArgs() > 0) {
	const char *flag = OPS_GetString();

	if (strcmp(flag, "-iterate") == 0 && OPS_GetNumRemainingInputArgs() > 0) {
	    const char *flag2 = OPS_GetString();
	    iterateSeen = true;
	    if (strcmp(flag2, "current") == 0)
		iterateTangent = CURRENT_TANGENT;
	    else if (strcmp(flag2, "initial") == 0)
		iterateTangent = INITIAL_TANGENT;
	    else if (strcmp(flag2, "noTangent") == 0)
		iterateTangent = NO_TANGENT;
	}
	else if (strcmp(flag, "-increment") == 0 && OPS_GetNumRemainingInputArgs() > 0) {
	    const char *flag2 = OPS_GetString();
	    if (strcmp(flag2, "current") == 0)
		incrementTangent = CURRENT_TANGENT;
	    else if (strcmp(flag2, "initial") == 0)
		incrementTangent = INITIAL_TANGENT;
	    else if (strcmp(flag2, "noTangent") == 0)
		incrementTangent = NO_TANGENT;
	}
	else if (wantMaxDim && strcmp(flag, "-maxDim") == 0
		 && OPS_GetNumRemainingInputArgs() > 0) {
	    int numdata = 1;
	    if (OPS_GetIntInput(&numdata, &maxDim) < 0) {
		opserr << "WARNING accelerated-Newton failed to read maxDim\n";
		// leave maxDim at its previous value and keep parsing
	    }
	}
	else if (wantNumTerms && strcmp(flag, "-numTerms") == 0
		 && OPS_GetNumRemainingInputArgs() > 0) {
	    int numdata = 1;
	    if (OPS_GetIntInput(&numdata, &numTerms) < 0) {
		opserr << "WARNING SecantNewton failed to read numTerms\n";
	    }
	}
	else if (wantCutOut
		 && (strcmp(flag, "-cutOut") == 0 || strcmp(flag, "-cutout") == 0)
		 && OPS_GetNumRemainingInputArgs() > 1) {
	    int numdata = 2;
	    if (OPS_GetDoubleInput(&numdata, R) < 0) {
		opserr << "WARNING SecantNewton failed to read cutOut values R1 and R2\n";
	    } else {
		cutOut = true;
	    }
	}
	else if (strcmp(flag, "-factorIncrementOnce") == 0
		 || strcmp(flag, "-FactorIncrementOnce") == 0
		 || strcmp(flag, "-factorOnce") == 0
		 || strcmp(flag, "-factoronce") == 0
		 || strcmp(flag, "-FactorOnce") == 0) {
	    factorOnce = 1;
	    factorOnceExplicit = true;
	}
    }

    if (!iterateSeen && factorOnceExplicit && factorOnce != 0) {
	iterateTangent = NO_TANGENT;
	opserr << "WARNING accelerated-Newton: -factorOnce without -iterate: "
		  "using -iterate noTangent\n";
    }

    if (factorOnce != 0 && iterateTangent != NO_TANGENT
	&& !(incrementTangent == INITIAL_TANGENT
	     && iterateTangent == INITIAL_TANGENT)) {
	opserr << "WARNING accelerated-Newton: -factorOnce / -increment initial "
		  "is disabled when -iterate is not noTangent. For a fixed "
		  "stiffness matrix with factor-once on the increment, use "
		  "-iterate noTangent.\n";
	factorOnce = 0;
    }

    if (incrementTangent == INITIAL_TANGENT
	&& (iterateTangent == NO_TANGENT || iterateTangent == INITIAL_TANGENT))
	factorOnce = 1;
}

void *
OPS_KrylovNewton()
{
    int incrementTangent = CURRENT_TANGENT;
    int iterateTangent = CURRENT_TANGENT;
    int maxDim = 3;
    int factorOnce = 0;
    int numTerms = 0;
    bool cutOut = false;
    double R[2];

    parseAcceleratorArgs(incrementTangent, iterateTangent, maxDim, factorOnce,
			 numTerms, cutOut, R,
			 /*wantMaxDim=*/true,
			 /*wantNumTerms=*/false,
			 /*wantCutOut=*/false);

    Accelerator *theAccel = new KrylovAccelerator(maxDim, iterateTangent);
    return new AcceleratedNewton(theAccel, incrementTangent, factorOnce);
}

void *
OPS_RaphsonNewton()
{
    int incrementTangent = CURRENT_TANGENT;
    int iterateTangent = CURRENT_TANGENT;
    int maxDim = 0;
    int factorOnce = 0;
    int numTerms = 0;
    bool cutOut = false;
    double R[2];

    parseAcceleratorArgs(incrementTangent, iterateTangent, maxDim, factorOnce,
			 numTerms, cutOut, R,
			 /*wantMaxDim=*/false,
			 /*wantNumTerms=*/false,
			 /*wantCutOut=*/false);

    Accelerator *theAccel = new RaphsonAccelerator(iterateTangent);
    return new AcceleratedNewton(theAccel, incrementTangent, factorOnce);
}

void *
OPS_SecantNewton()
{
    int incrementTangent = CURRENT_TANGENT;
    int iterateTangent = CURRENT_TANGENT;
    int maxDim = 3;
    int factorOnce = 0;
    int numTerms = 2;
    bool cutOut = false;
    double R[2];

    parseAcceleratorArgs(incrementTangent, iterateTangent, maxDim, factorOnce,
			 numTerms, cutOut, R,
			 /*wantMaxDim=*/true,
			 /*wantNumTerms=*/true,
			 /*wantCutOut=*/true);

    Accelerator *theAccel = 0;
    if (numTerms <= 1) {
	if (cutOut)
	    theAccel = new SecantAccelerator1(maxDim, iterateTangent, R[0], R[1]);
	else
	    theAccel = new SecantAccelerator1(maxDim, iterateTangent);
    } else if (numTerms == 2) {
	if (cutOut)
	    theAccel = new SecantAccelerator2(maxDim, iterateTangent, R[0], R[1]);
	else
	    theAccel = new SecantAccelerator2(maxDim, iterateTangent);
    } else {
	if (cutOut)
	    theAccel = new SecantAccelerator3(maxDim, iterateTangent, R[0], R[1]);
	else
	    theAccel = new SecantAccelerator3(maxDim, iterateTangent);
    }

    return new AcceleratedNewton(theAccel, incrementTangent, factorOnce);
}

void *
OPS_PeriodicNewton()
{
    int incrementTangent = CURRENT_TANGENT;
    int iterateTangent = CURRENT_TANGENT;
    int maxDim = 3;
    int factorOnce = 0;
    int numTerms = 0;
    bool cutOut = false;
    double R[2];

    parseAcceleratorArgs(incrementTangent, iterateTangent, maxDim, factorOnce,
			 numTerms, cutOut, R,
			 /*wantMaxDim=*/true,
			 /*wantNumTerms=*/false,
			 /*wantCutOut=*/false);

    Accelerator *theAccel = new PeriodicAccelerator(maxDim, iterateTangent);
    return new AcceleratedNewton(theAccel, incrementTangent, factorOnce);
}

void *
OPS_MillerNewton()
{
    int incrementTangent = CURRENT_TANGENT;
    int iterateTangent = CURRENT_TANGENT;
    int maxDim = 3;
    int factorOnce = 0;
    int numTerms = 0;
    bool cutOut = false;
    double R[2];

    parseAcceleratorArgs(incrementTangent, iterateTangent, maxDim, factorOnce,
			 numTerms, cutOut, R,
			 /*wantMaxDim=*/true,
			 /*wantNumTerms=*/false,
			 /*wantCutOut=*/false);

    opserr << "WARNING MillerNewton: Miller acceleration is deactivated; "
	   << "using AcceleratedNewton without MillerAccelerator.\n";

    Accelerator *theAccel = 0;
    return new AcceleratedNewton(theAccel, incrementTangent, factorOnce);
}

// Constructor
AcceleratedNewton::AcceleratedNewton(int theTangentToUse, int factOnce)
  :EquiSolnAlgo(EquiALGORITHM_TAGS_AcceleratedNewton),
   theTest(0), tangent(theTangentToUse), factorOnce(factOnce),
   theAccelerator(0), vAccel(0), 
   numFactorizations(0), numIterations(0)
//   totalTimer(), totalTimeReal(0.0), totalTimeCPU(0.0),
//   solveTimer(), solveTimeReal(0.0), solveTimeCPU(0.0),
//   accelTimer(), accelTimeReal(0.0), accelTimeCPU(0.0)
{

}

// Deferred ConvergenceTest (OPS factories); owns theAccel.
AcceleratedNewton::AcceleratedNewton(Accelerator *theAccel,
				     int theTangentToUse,
				     int factOnce)
  :EquiSolnAlgo(EquiALGORITHM_TAGS_AcceleratedNewton),
   theTest(0), tangent(theTangentToUse), factorOnce(factOnce),
   theAccelerator(theAccel), vAccel(0),
   numFactorizations(0), numIterations(0)
{

}

AcceleratedNewton::AcceleratedNewton(ConvergenceTest &theT,
				     Accelerator *theAccel,
				     int theTangentToUse,
				     int factOnce)
  :EquiSolnAlgo(EquiALGORITHM_TAGS_AcceleratedNewton),
   theTest(&theT), tangent(theTangentToUse), factorOnce(factOnce),
   theAccelerator(theAccel), vAccel(0), 
   numFactorizations(0), numIterations(0)
//   totalTimer(), totalTimeReal(0.0), totalTimeCPU(0.0),
//   solveTimer(), solveTimeReal(0.0), solveTimeCPU(0.0),
//   accelTimer(), accelTimeReal(0.0), accelTimeCPU(0.0)
{
 
}

// Destructor
AcceleratedNewton::~AcceleratedNewton()
{
  if (theAccelerator != 0)
    delete theAccelerator;

  if (vAccel != 0)
    delete vAccel;

  //opserr << "AcceleratedNewton::~AcceleratedNewton " << numFactorizations << endln;
}

int
AcceleratedNewton::setConvergenceTest(ConvergenceTest *newTest)
{
  theTest = newTest;
  return 0;
}

int
AcceleratedNewton::domainChanged(void)
{
  // Cached factorization invalid after domain change / setSize - reform increment tangent next solve.
  if (factorOnce == 2)
    factorOnce = 1;
  return 0;
}

int 
AcceleratedNewton::solveCurrentStep(void)
{
  // set up some pointers and check they are valid
  // NOTE this could be taken away if we set Ptrs as protecetd in superclass
  AnalysisModel *theAnaModel = this->getAnalysisModelPtr();
  IncrementalIntegrator *theIntegrator = this->getIncrementalIntegratorPtr();
  LinearSOE *theSOE = this->getLinearSOEptr();
  
  if ((theAnaModel == 0) || (theIntegrator == 0) || (theSOE == 0)
      || (theTest == 0)){
    opserr << "WARNING AcceleratedNewton::solveCurrentStep() - setLinks() has";
    opserr << " not been called - or no ConvergenceTest has been set\n";
    return -5;
  }	

  if (theAccelerator != 0)
    theAccelerator->newStep(*theSOE);

  int numEqns = theSOE->getNumEqn();

  if (vAccel == 0)
    vAccel = new Vector(numEqns);

  if (vAccel->Size() != numEqns) {
    delete vAccel;
    vAccel = new Vector(numEqns);
  }

  if (vAccel == 0) {
    opserr << "WARNING AcceleratedNewton::solveCurrentStep() - ";
    opserr << " could not allocate correction vector vAccel\n";
    return -6;
  }

  //totalTimer.start();

  // Evaluate system residual R(y_0)
  if (theIntegrator->formUnbalance() < 0) {
    opserr << "WARNING AcceleratedNewton::solveCurrentStep() -";
    opserr << "the Integrator failed in formUnbalance()\n";	
    return -2;
  }

  // Increment-side formTangent only; iterate-side uses accelerator/updateTangent.
  if (factorOnce != 2) {
    if (theIntegrator->formTangent(tangent) < 0){
      opserr << "WARNING AcceleratedNewton::solveCurrentStep() -";
      opserr << "the Integrator failed in formTangent()\n";
      return -1;
    }
    if (factorOnce == 1)
      factorOnce = 2;

    // Count factorization of the first tangent
    numFactorizations++;
  }
  
  // set itself as the ConvergenceTest objects EquiSolnAlgo
  theTest->setEquiSolnAlgo(*this);
  if (theTest->start() < 0) {
    opserr << "AcceleratedNewton::solveCurrentStep() -";
    opserr << "the ConvergenceTest object failed in start()\n";
    return -3;
  }
  
  // Loop counter
  int k = 1;

  int result = -1;

  do {

    //solveTimer.start();
    // Solve for displacement increment
    if (theSOE->solve() < 0) {
      opserr << "WARNING AcceleratedNewton::solveCurrentStep() -";
      opserr << "the LinearSysOfEqn failed in solve()\n";	
      return -3;
    }
//    solveTimer.pause();
//    solveTimeReal += solveTimer.getReal();
//    solveTimeCPU  += solveTimer.getCPU();

    // Get the modified Newton increment
    *vAccel = theSOE->getX();

    // Accelerate the displacement increment
    if (theAccelerator != 0) {

//      accelTimer.start();
      if (theAccelerator->accelerate(*vAccel, *theSOE, *theIntegrator) < 0) {
	opserr << "WARNING AcceleratedNewton::solveCurrentStep() -";
	opserr << "the Accelerator failed in accelerate()\n";
	return -1;
      }
//      accelTimer.pause();
//      accelTimeReal += accelTimer.getReal();
//      accelTimeCPU  += accelTimer.getCPU();
    }

    // Update system with accelerated displacement increment v_{k+1}
    if (theIntegrator->update(*vAccel) < 0) {
      opserr << "WARNING AcceleratedNewton::solveCurrentStep() -";
      opserr << "the Integrator failed in update()\n";	
      return -4;
    }	

    // Evaluate residual
    if (theIntegrator->formUnbalance() < 0) {
      opserr << "WARNING AcceleratedNewton::solveCurrentStep() -";
      opserr << "the Integrator failed in formUnbalance()\n";	
      return -2;
    }

    numIterations++;

    // Check convergence criteria
    result = theTest->test();

    if (result == -1) {
      // Let the accelerator update the tangent if needed
      if (theAccelerator != 0) {
	int ret = theAccelerator->updateTangent(*theIntegrator);
	if (ret < 0) {
	  opserr << "WARNING AcceleratedNewton::solveCurrentStep() -";
	  opserr << "the Accelerator failed in updateTangent()\n";
	  return -1;
	}
	if (ret > 0)
	  numFactorizations++;
      }
    }
    //opserr << "ACCEL: " << numFactorizations << endln;
    this->record(k++);

  } while (result == -1);

//  totalTimer.pause();
//  totalTimeReal += totalTimer.getReal();
//  totalTimeCPU  += totalTimer.getCPU();
  
  if (result == -2) {
    opserr << "AcceleratedNewton::solveCurrentStep() -";
    opserr << "The ConvergenceTest object failed in test()\n";
    return -3;
  }
  
  // note - if positive result we are returning what the convergence
  // test returned which should be the number of iterations
  return result;
}

ConvergenceTest *
AcceleratedNewton::getTest(void)
{
  return theTest;
}

int
AcceleratedNewton::sendSelf(int cTag, Channel &theChannel)
{
  static ID data(3);
  data(0) = tangent;
  if (theAccelerator != 0)
    data(1) = theAccelerator->getClassTag();
  else
    data(1) = -1;
  data(2) = factorOnce;

  int res = theChannel.sendID(0, cTag, data);
  if (res < 0) {
    opserr << "AcceleratedNewton::recvSelf() - failed to send data\n";
    return -1;
  }

  if (theAccelerator != 0) {  
    res = theAccelerator->sendSelf(cTag, theChannel);
    if (res < 0) {
      opserr << "AcceleratedNewton::recvSelf() - accelerator to send\n";
      return -1;
    }
  }

  return 0;
}

int
AcceleratedNewton::recvSelf(int cTag, Channel &theChannel, 
			    FEM_ObjectBroker &theBroker)
{
  static ID data(3);
  int res = theChannel.recvID(0, cTag, data);

  if (res < 0) {
    opserr << "AcceleratedNewton::recvSelf() - failed to recv data\n";
    return -1;
  }

  tangent = data(0);
  factorOnce = (data.Size() > 2) ? data(2) : 0;

  if (data(1) != -1) {

    if (theAccelerator != 0)
      delete theAccelerator;

    theAccelerator = theBroker.getAccelerator(data(1));
    if (theAccelerator == 0) {
      opserr << "AcceleratedNewton::recvSelf() - no acccelerator of classTag " << data(1) << " exists\n";
      return -1;
    }
    
    if (res == 0) {
      res = theAccelerator->recvSelf(cTag, theChannel, theBroker);
      if (res < 0) {
	opserr << "AcceleratedNewton::recvSelf() - accelerator failed to recvSelf\n";
	return -1;
      }
    }
  }

  return 0;
}

void
AcceleratedNewton::Print(OPS_Stream &s, int flag)
{
  s << "AcceleratedNewton" << endln;
  LinearSOE *theSOE = this->getLinearSOEptr();
  s << "\tNumber of equations: " << theSOE->getNumEqn() << endln;

  if (theAccelerator != 0)
    theAccelerator->Print(s,flag);
  else
    s << "\tNo accelerator --> Modified Newton" << endln;
}
