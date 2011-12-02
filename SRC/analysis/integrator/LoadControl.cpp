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
                                                                        
// $Revision: 1.5 $
// $Date: 2005-10-19 21:53:57 $
// $Source: /usr/local/cvs/OpenSees/SRC/analysis/integrator/LoadControl.cpp,v $
                                                                        
                                                                        
// 
// Written: fmk 
// Created: 07/98
// Revision: A
//
// Description: This file contains the class definition for LoadControl.
// LoadControl is an algorithmic class for perfroming a static analysis
// using a load control integration scheme.
//
// What: "@(#) LoadControl.h, revA"



#include <LoadControl.h>
#include <AnalysisModel.h>
#include <LinearSOE.h>
#include <Vector.h>
#include <Channel.h>


LoadControl::LoadControl(double dLambda, int numIncr, double min, double max)
:StaticIntegrator(INTEGRATOR_TAGS_LoadControl),
 deltaLambda(dLambda), 
 specNumIncrStep(numIncr), numIncrLastStep(numIncr),
 dLambdaMin(min), dLambdaMax(max)
{
  // to avoid divide-by-zero error on first update() ensure numIncr != 0
  if (numIncr == 0) {
    opserr << "WARNING LoadControl::LoadControl() - numIncr set to 0, 1 assumed\n";
    specNumIncrStep = 1.0;
    numIncrLastStep = 1.0;
  }
}


LoadControl::~LoadControl()
{
    
}

int 
LoadControl::newStep(void)
{
    AnalysisModel *theModel = this->getAnalysisModelPtr();    
    if (theModel == 0) {
	opserr << "LoadControl::newStep() - no associated AnalysisModel\n";
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
LoadControl::update(const Vector &deltaU)
{
    AnalysisModel *myModel = this->getAnalysisModelPtr();
    LinearSOE *theSOE = this->getLinearSOEPtr();
    if (myModel == 0 || theSOE == 0) {
	opserr << "WARNING LoadControl::update() ";
	opserr << "No AnalysisModel or LinearSOE has been set\n";
	return -1;
    }

    myModel->incrDisp(deltaU);    
    if (myModel->updateDomain() < 0) {
      opserr << "LoadControl::update - model failed to update for new dU\n";
      return -1;
    }

    // Set deltaU for the convergence test
    theSOE->setX(deltaU);

    numIncrLastStep++;

    return 0;
}


int
LoadControl::setDeltaLambda(double newValue)
{
  // we set the #incr at last step = #incr so get newValue incr
  numIncrLastStep = specNumIncrStep;
  deltaLambda = newValue;
  return 0;
}


int
LoadControl::sendSelf(int cTag,
		      Channel &theChannel)
{
  Vector data(5);
  data(0) = deltaLambda;
  data(1) = specNumIncrStep;
  data(2) = numIncrLastStep;
  data(3) = dLambdaMin;
  data(4) = dLambdaMax;
  if (theChannel.sendVector(this->getDbTag(), cTag, data) < 0) {
      opserr << "LoadControl::sendSelf() - failed to send the Vector\n";
      return -1;
  }
  return 0;
}


int
LoadControl::recvSelf(int cTag,
		      Channel &theChannel, FEM_ObjectBroker &theBroker)
{
  Vector data(5);
  if (theChannel.recvVector(this->getDbTag(), cTag, data) < 0) {
      opserr << "LoadControl::sendSelf() - failed to send the Vector\n";
      deltaLambda = 0;
      return -1;
  }      
  deltaLambda = data(0);
  specNumIncrStep = data(1);
  numIncrLastStep = data(2);
  dLambdaMin = data(3);
  dLambdaMax = data(4);
  return 0;
}



void
LoadControl::Print(OPS_Stream &s, int flag)
{
    AnalysisModel *theModel = this->getAnalysisModelPtr();
    if (theModel != 0) {
	double currentLambda = theModel->getCurrentDomainTime();
	s << "\t LoadControl - currentLambda: " << currentLambda;
	s << "  deltaLambda: " << deltaLambda << endln;
    } else 
	s << "\t LoadControl - no associated AnalysisModel\n";
    
}

