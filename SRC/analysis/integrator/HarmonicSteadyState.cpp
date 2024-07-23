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

// Written: Seweryn Kokot, Opole University of Technology, Poland
// Created: 2021
//
// based on LoadControl.cpp
// written: fmk
//
// Description: This file contains the class definition for HarmonicSteadyState.
// HarmonicSteadyState is an algorithmic class for performing a quasi-static harmonic
// steady-state analysis.

#include "HarmonicSteadyState.h"
#include <AnalysisModel.h>
#include <LinearSOE.h>
#include <Vector.h>
#include <Channel.h>
#include <FE_Element.h>
#include <FE_EleIter.h>
#include <Node.h>
#include <DOF_Group.h>
#include <DOF_GrpIter.h>
#include <LoadPattern.h>
#include <LoadPatternIter.h>
#include <Domain.h>
#include<Parameter.h>
#include<ParameterIter.h>
#include<EquiSolnAlgo.h>
#include <elementAPI.h>
#include <iostream>

void* OPS_HarmonicSteadyState()
{
    if(OPS_GetNumRemainingInputArgs() < 2) {
	opserr<<"insufficient arguments\n";
	return 0;
    }

    double lambda;
    int numData = 1;
    if(OPS_GetDoubleInput(&numData,&lambda) < 0) {
	opserr<<"WARNING failed to read double lambda\n";
	return 0;
    }

	double period = 0;
	numData = 1;
	if(OPS_GetDoubleInput(&numData,&period) < 0) {
	    opserr<<"WARNING failed to read double period\n";
	    return 0;
	}

    int numIter = 1;
    double mLambda[2] = {lambda,lambda};
    if(OPS_GetNumRemainingInputArgs() > 2) {
	if(OPS_GetIntInput(&numData,&numIter) < 0) {
	    opserr<<"WARNING failed to read int numIter\n";
	    return 0;
	}
	numData = 2;
	if(OPS_GetDoubleInput(&numData,&mLambda[0]) < 0) {
	    opserr<<"WARNING failed to read double min and max\n";
	    return 0;
	}
    }

    return new HarmonicSteadyState(lambda,period,numIter,mLambda[0],mLambda[1]);
}

HarmonicSteadyState::HarmonicSteadyState(double dLambda, double period, int numIncr, double min, double max,  int classtag)
    : StaticIntegrator(classtag),
 deltaLambda(dLambda), loadPeriod(period),
 specNumIncrStep(numIncr), numIncrLastStep(numIncr),
	  dLambdaMin(min), dLambdaMax(max), gradNumber(0), sensitivityFlag(0)
{
  // to avoid divide-by-zero error on first update() ensure numIncr != 0
  if (numIncr == 0) {
    opserr << "WARNING HarmonicSteadyState::HarmonicSteadyState() - numIncr set to 0, 1 assumed\n";
    specNumIncrStep = 1.0;
    numIncrLastStep = 1.0;
  }
}


HarmonicSteadyState::~HarmonicSteadyState()
{

}

int
HarmonicSteadyState::newStep(void)
{
    AnalysisModel *theModel = this->getAnalysisModel();
    if (theModel == 0) {
	opserr << "HarmonicSteadyState::newStep() - no associated AnalysisModel\n";
	return -1;
    }

    // determine delta lambda for this step based on dLambda and #iter of last step
    double factor = specNumIncrStep/numIncrLastStep;
    deltaLambda *=factor;

    if (deltaLambda < dLambdaMin)
      deltaLambda = dLambdaMin;
    else if (deltaLambda > dLambdaMax)
      deltaLambda = dLambdaMax;

    double currentLambda = theModel->getCurrentDomainTime();

    currentLambda += deltaLambda;
    theModel->applyLoadDomain(currentLambda);

    numIncrLastStep = 0;

     return 0;
}

int
HarmonicSteadyState::update(const Vector &deltaU)
{
    AnalysisModel *myModel = this->getAnalysisModel();
    LinearSOE *theSOE = this->getLinearSOE();
    if (myModel == 0 || theSOE == 0) {
	opserr << "WARNING HarmonicSteadyState::update() ";
	opserr << "No AnalysisModel or LinearSOE has been set\n";
	return -1;
    }

    myModel->incrDisp(deltaU);
    if (myModel->updateDomain() < 0) {
      opserr << "HarmonicSteadyState::update - model failed to update for new dU\n";
      return -1;
    }

    // Set deltaU for the convergence test
    theSOE->setX(deltaU);

    numIncrLastStep++;

    return 0;
}


int
HarmonicSteadyState::setDeltaLambda(double newValue)
{
  // we set the #incr at last step = #incr so get newValue incr
  numIncrLastStep = specNumIncrStep;
  deltaLambda = newValue;
  return 0;
}


int
HarmonicSteadyState::sendSelf(int cTag,
		      Channel &theChannel)
{
  Vector data(6);
  data(0) = deltaLambda;
  data(1) = loadPeriod;
  data(2) = specNumIncrStep;
  data(3) = numIncrLastStep;
  data(4) = dLambdaMin;
  data(5) = dLambdaMax;
  if (theChannel.sendVector(this->getDbTag(), cTag, data) < 0) {
      opserr << "HarmonicSteadyState::sendSelf() - failed to send the Vector\n";
      return -1;
  }
  return 0;
}


int
HarmonicSteadyState::recvSelf(int cTag,
		      Channel &theChannel, FEM_ObjectBroker &theBroker)
{
  Vector data(6);
  if (theChannel.recvVector(this->getDbTag(), cTag, data) < 0) {
      opserr << "HarmonicSteadyState::sendSelf() - failed to send the Vector\n";
      deltaLambda = 0;
      return -1;
  }
  deltaLambda = data(0);
  loadPeriod = data(1);
  specNumIncrStep = data(2);
  numIncrLastStep = data(3);
  dLambdaMin = data(4);
  dLambdaMax = data(5);
  return 0;
}



void
HarmonicSteadyState::Print(OPS_Stream &s, int flag)
{
     AnalysisModel *theModel = this->getAnalysisModel();
    if (theModel != 0) {
	double currentLambda = theModel->getCurrentDomainTime();
	s << "\t HarmonicSteadyState - currentLambda: " << currentLambda;
	s << "  deltaLambda: " << deltaLambda << endln;
	s << "  Load Period: " << loadPeriod << endln;
    } else
	s << "\t HarmonicSteadyState - no associated AnalysisModel\n";

}

int
HarmonicSteadyState::formEleTangent(FE_Element *theEle)
{
  static const double twoPi = 2*3.1415926535897932;
  double twoPiSquareOverPeriodSquare = twoPi*twoPi/(loadPeriod*loadPeriod);
  if (statusFlag == CURRENT_TANGENT) {
    theEle->zeroTangent();
    theEle->addKtToTang();
	theEle->addMtoTang(-twoPiSquareOverPeriodSquare);
  } else if (statusFlag == INITIAL_TANGENT) {
    theEle->zeroTangent();
    theEle->addKiToTang();
	theEle->addMtoTang(-twoPiSquareOverPeriodSquare);
  } else if (statusFlag == HALL_TANGENT)  {
    theEle->zeroTangent();
    theEle->addKtToTang(cFactor);
    theEle->addKiToTang(iFactor);
	theEle->addMtoTang(-twoPiSquareOverPeriodSquare);
  }

  return 0;
}


int
HarmonicSteadyState::formEleResidual(FE_Element* theEle)
{
    if(sensitivityFlag == 0) {  // no sensitivity
	this->StaticIntegrator::formEleResidual(theEle);
    } else {
	theEle->zeroResidual();
	theEle->addResistingForceSensitivity(gradNumber);
    }
    return 0;
}

int
HarmonicSteadyState::formIndependentSensitivityRHS()
{
    return 0;
}

int
HarmonicSteadyState::formSensitivityRHS(int passedGradNumber)
{
    sensitivityFlag = 1;

    // Set a couple of data members
    gradNumber = passedGradNumber;

    // get model
    AnalysisModel* theAnalysisModel = this->getAnalysisModel();
    LinearSOE* theSOE = this->getLinearSOE();

    // Loop through elements
    FE_Element *elePtr;
    FE_EleIter &theEles = theAnalysisModel->getFEs();
    while((elePtr = theEles()) != 0) {
      theSOE->addB(  elePtr->getResidual(this),  elePtr->getID()  );
    }

    // Loop through the loadPatterns and add the dPext/dh contributions
    static Vector oneDimVectorWithOne(1);
    oneDimVectorWithOne(0) = 1.0;
    static ID oneDimID(1);

    Node *aNode;
    DOF_Group *aDofGroup;
    int nodeNumber, dofNumber, relevantID, i, sizeRandomLoads, numRandomLoads;
    LoadPattern *loadPatternPtr;

    Domain *theDomain = theAnalysisModel->getDomainPtr();
    LoadPatternIter &thePatterns = theDomain->getLoadPatterns();
    while((loadPatternPtr = thePatterns()) != 0) {
	const Vector &randomLoads = loadPatternPtr->getExternalForceSensitivity(gradNumber);
	sizeRandomLoads = randomLoads.Size();
	if (sizeRandomLoads == 1) {
	    // No random loads in this load pattern
	}
	else {
	    // Random loads: add contributions to the 'B' vector
	    numRandomLoads = (int)(sizeRandomLoads/2);
	    for (i=0; i<numRandomLoads*2; i=i+2) {
		nodeNumber = (int)randomLoads(i);
		dofNumber = (int)randomLoads(i+1);
		aNode = theDomain->getNode(nodeNumber);
		aDofGroup = aNode->getDOF_GroupPtr();
		const ID &anID = aDofGroup->getID();
		relevantID = anID(dofNumber-1);
		oneDimID(0) = relevantID;
		theSOE->addB(oneDimVectorWithOne, oneDimID);
	    }
	}
    }

    // reset sensitivity flag
    sensitivityFlag = 0;

    return 0;
}

int
HarmonicSteadyState::saveSensitivity(const Vector &v, int gradNum, int numGrads)
{
    // get model
    AnalysisModel* theAnalysisModel = this->getAnalysisModel();

    DOF_GrpIter &theDOFGrps = theAnalysisModel->getDOFs();
    DOF_Group 	*dofPtr;

    while ( (dofPtr = theDOFGrps() ) != 0)  {
//	dofPtr->saveSensitivity(v,0,0,gradNum,numGrads);
	dofPtr->saveDispSensitivity(v,gradNum,numGrads);

    }

    return 0;
}

int
HarmonicSteadyState::commitSensitivity(int gradNum, int numGrads)
{
      // get model
    AnalysisModel* theAnalysisModel = this->getAnalysisModel();

    // Loop through the FE_Elements and set unconditional sensitivities
    FE_Element *elePtr;
    FE_EleIter &theEles = theAnalysisModel->getFEs();
    while((elePtr = theEles()) != 0) {
	elePtr->commitSensitivity(gradNum, numGrads);
    }

    return 0;
}



// false for LC and true for DC
   bool
HarmonicSteadyState::computeSensitivityAtEachIteration()
{

return false;
}





int
HarmonicSteadyState::computeSensitivities(void)
{
//  opserr<<" computeSensitivity::start"<<endln;
    LinearSOE *theSOE = this->getLinearSOE();

    /*
  if (theAlgorithm == 0) {
    opserr << "ERROR the FE algorithm must be defined before ";
    opserr << "the sensitivity algorithm\n";
    return -1;
  }
*/
/*
   // Get pointer to the system of equations (SOE)
	LinearSOE *theSOE = theAlgorithm->getLinearSOEptr();
	if (theSOE == 0) {
	  opserr << "ERROR the FE linearSOE must be defined before ";
	  opserr << "the sensitivity algorithm\n";
	  return -1;
	}

	// Get pointer to incremental integrator
	IncrementalIntegrator *theIncInt = theAlgorithm->getIncrementalIntegratorPtr();
//	IncrementalIntegrator *theIncIntSens=theAlgorithm->getIncrementalIntegratorPtr();//Abbas
	if (theIncInt == 0 ) {
	  opserr << "ERROR the FE integrator must be defined before ";
	  opserr << "the sensitivity algorithm\n";
	  return -1;
	}

	// Form current tangent at converged state
	// (would be nice with an if-statement here in case
	// the current tangent is already formed)
	if (this->formTangent(CURRENT_TANGENT) < 0){
		opserr << "WARNING SensitivityAlgorithm::computeGradients() -";
		opserr << "the Integrator failed in formTangent()\n";
		return -1;
	}
	*/
	// Zero out the old right-hand side of the SOE
	theSOE->zeroB();


	if (this == 0) {
	  opserr << "ERROR SensitivityAlgorithm::computeSensitivities() -";
	  opserr << "the SensitivityIntegrator is NULL\n";
	  return -1;
	}

	// Form the part of the RHS which are independent of parameter
	this->formIndependentSensitivityRHS();
	AnalysisModel *theModel = this->getAnalysisModel();
	Domain *theDomain=theModel->getDomainPtr();
	ParameterIter &paramIter = theDomain->getParameters();

	Parameter *theParam;
	// De-activate all parameters
	while ((theParam = paramIter()) != 0)
	  theParam->activate(false);

	// Now, compute sensitivity wrt each parameter
	int numGrads = theDomain->getNumParameters();
	paramIter = theDomain->getParameters();

	while ((theParam = paramIter()) != 0) {

	  // Activate this parameter
	  theParam->activate(true);

	  // Zero the RHS vector
	  theSOE->zeroB();

	  // Get the grad index for this parameter
	  int gradIndex = theParam->getGradIndex();

	  // Form the RHS
	  this->formSensitivityRHS(gradIndex);

	  // Solve for displacement sensitivity

	  theSOE->solve();
	  // Save sensitivity to nodes
	  this->saveSensitivity( theSOE->getX(), gradIndex, numGrads );



	  // Commit unconditional history variables (also for elastic problems; strain sens may be needed anyway)
	  this->commitSensitivity(gradIndex, numGrads);

	  // De-activate this parameter for next sensitivity calc
	  theParam->activate(false);
	//  opserr<<"HarmonicSteadyState::..........ComputeSensitivities. end"<<endln;
	}

	return 0;
}
