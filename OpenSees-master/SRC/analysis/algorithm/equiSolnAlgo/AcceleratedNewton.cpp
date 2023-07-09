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

#include <fstream>

// Constructor
AcceleratedNewton::AcceleratedNewton(int theTangentToUse)
  :EquiSolnAlgo(EquiALGORITHM_TAGS_AcceleratedNewton),
   theTest(0), tangent(theTangentToUse),
   theAccelerator(0), vAccel(0), 
   numFactorizations(0), numIterations(0)
//   totalTimer(), totalTimeReal(0.0), totalTimeCPU(0.0),
//   solveTimer(), solveTimeReal(0.0), solveTimeCPU(0.0),
//   accelTimer(), accelTimeReal(0.0), accelTimeCPU(0.0)
{

}

AcceleratedNewton::AcceleratedNewton(ConvergenceTest &theT,
				     Accelerator *theAccel,
				     int theTangentToUse)
  :EquiSolnAlgo(EquiALGORITHM_TAGS_AcceleratedNewton),
   theTest(&theT), tangent(theTangentToUse),
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

  // Evaluate system Jacobian J = R'(y)|y_0
  if (theIntegrator->formTangent(tangent) < 0){
    opserr << "WARNING AcceleratedNewton::solveCurrentStep() -";
    opserr << "the Integrator failed in formTangent()\n";
    return -1;
  }
  
  // Count factorization of the first tangent
  numFactorizations++;
  
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
  static ID data(2);
  data(0) = tangent;
  if (theAccelerator != 0)
    data(1) = theAccelerator->getClassTag();
  else
    data(1) = -1;

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
  static ID data(2);
  int res = theChannel.recvID(0, cTag, data);

  if (res < 0) {
    opserr << "AcceleratedNewton::recvSelf() - failed to recv data\n";
    return -1;
  }

  tangent = data(0) = tangent;

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
