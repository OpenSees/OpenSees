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

// $Revision: 1.4 $
// $Date: 2007-04-02 23:41:13 $
// $Source: /usr/local/cvs/OpenSees/SRC/analysis/algorithm/equiSolnAlgo/PeriodicNewton.cpp,v $

// Written: MHS
// Created: Oct 2002
//
// Description: This file contains the class definition for 
// PeriodicNewton. PeriodicNewton is a class which uses the
// Newton-Raphson solution algorithm
// to solve the equations. No member functions are declared as virtual as 
// it is not expected that this class will be subclassed.

#include <PeriodicNewton.h>
#include <AnalysisModel.h>
#include <StaticAnalysis.h>
#include <IncrementalIntegrator.h>
#include <LinearSOE.h>
#include <ID.h>
#include <Channel.h>
#include <FEM_ObjectBroker.h>
#include <ConvergenceTest.h>

// Constructor
PeriodicNewton::PeriodicNewton(int theTangentToUse, int mc)
:EquiSolnAlgo(EquiALGORITHM_TAGS_PeriodicNewton),
 tangent(theTangentToUse), maxCount(mc)
{
  
}


PeriodicNewton::PeriodicNewton(ConvergenceTest &theT, int theTangentToUse, int mc)
:EquiSolnAlgo(EquiALGORITHM_TAGS_PeriodicNewton),
 tangent(theTangentToUse), maxCount(mc)
{

}

// Destructor
PeriodicNewton::~PeriodicNewton()
{

}

int 
PeriodicNewton::solveCurrentStep(void)
{
    // set up some pointers and check they are valid
    // NOTE this could be taken away if we set Ptrs as protecetd in superclass
    AnalysisModel       *theAnalysisModel = this->getAnalysisModelPtr();
    IncrementalIntegrator *theIncIntegratorr = this->getIncrementalIntegratorPtr();
    LinearSOE	        *theSOE = this->getLinearSOEptr();

    if ((theAnalysisModel == 0) || (theIncIntegratorr == 0) || (theSOE == 0)
	|| (theTest == 0)){
	opserr << "WARNING PeriodicNewton::solveCurrentStep() - setLinks() has";
	opserr << " not been called - or no ConvergenceTest has been set\n";
	return -5;
    }	

    // we form the tangent
    
    if (theIncIntegratorr->formUnbalance() < 0) {
	opserr << "WARNING PeriodicNewton::solveCurrentStep() -";
	opserr << "the Integrator failed in formUnbalance()\n";	
	return -2;
    }	

    if (theIncIntegratorr->formTangent(tangent) < 0){
	opserr << "WARNING PeriodicNewton::solveCurrentStep() -";
	opserr << "the Integrator failed in formTangent()\n";
	return -1;
    }		    

    // set itself as the ConvergenceTest objects EquiSolnAlgo
    theTest->setEquiSolnAlgo(*this);
    if (theTest->start() < 0) {
      opserr << "PeriodicNewton::solveCurrentStep() -";
      opserr << "the ConvergenceTest object failed in start()\n";
      return -3;
    }

    // repeat until convergence is obtained or reach max num iterations
    int result = -1;
    int count = 0;
	int iter = 0;
    do {
	if (theSOE->solve() < 0) {
	    opserr << "WARNING PeriodicNewton::solveCurrentStep() -";
	    opserr << "the LinearSysOfEqn failed in solve()\n";	
	    return -3;
	}	    

	if (theIncIntegratorr->update(theSOE->getX()) < 0) {
	    opserr << "WARNING PeriodicNewton::solveCurrentStep() -";
	    opserr << "the Integrator failed in update()\n";	
	    return -4;
	}	        

	if (theIncIntegratorr->formUnbalance() < 0) {
	    opserr << "WARNING PeriodicNewton::solveCurrentStep() -";
	    opserr << "the Integrator failed in formUnbalance()\n";	
	    return -2;
	}	

	this->record(count++);
	result = theTest->test();
	
	iter++;
	if (iter > maxCount) {
		if (theIncIntegratorr->formTangent(tangent) < 0){
		opserr << "WARNING PeriodicNewton::solveCurrentStep() -";
		opserr << "the Integrator failed in formTangent()\n";
		return -1;
		}
		iter = 0;
	}

    } while (result == -1);

    if (result == -2) {
      opserr << "PeriodicNewton::solveCurrentStep() -";
      opserr << "the ConvergenceTest object failed in test()\n";
      return -3;
    }

    return result;
}

int
PeriodicNewton::sendSelf(int cTag, Channel &theChannel)
{
  int result = 0;
  int dataTag = this->getDbTag();
  static ID data(3);
  data(0) = theTest->getClassTag();
  data(1) = theTest->getDbTag();
  data(2) = maxCount;

  result = theChannel.sendID(dataTag, cTag, data);
  if (result != 0) {
    opserr << "PeriodicNewton::sendSelf() - failed to send ID\n";
    return result;
  }

  result = theTest->sendSelf(cTag, theChannel);
  if (result != 0) {
    opserr << "PeriodicNewton::sendSelf() - failed to send CTest object\n";
    return result;
  }
  
  return 0;
}

int
PeriodicNewton::recvSelf(int cTag, 
			Channel &theChannel, 
			FEM_ObjectBroker &theBroker)
{
    static ID data(3);
    int result;
    int dataTag = this->getDbTag();

    result = theChannel.recvID(dataTag, cTag, data);    
    if (result != 0) {
      opserr << "PeriodicNewton::recvSelf() - failed to receive ID\n";
      return result;
    }
    int ctType = data(0);
    int ctDb = data(1);
    maxCount = data(2);

    theTest = theBroker.getNewConvergenceTest(ctType);
    theTest->setDbTag(ctDb);
    result = theTest->recvSelf(cTag, theChannel, theBroker);
    if (result != 0) {
      opserr << "PeriodicNewton::recvSelf() - failed to recv CTest object\n";
      return result;
    }
    
    return 0;
}

void
PeriodicNewton::Print(OPS_Stream &s, int flag)
{
    if (flag == 0) {
	s << "PeriodicNewton" << endln;
	s << "Max count: " << maxCount << endln;
    }
}
