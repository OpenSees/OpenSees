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
                                                                        
// $Revision: 1.8 $
// $Date: 2005-10-19 21:53:57 $
// $Source: /usr/local/cvs/OpenSees/SRC/analysis/integrator/MinUnbalDispNorm.cpp,v $
                                                                        
                                                                        
// File: ~/analysis/integrator/MinUnbalDispNorm.C
// 
// Written: fmk 
// Created: 07/99
// Revision: A
//
// What: "@(#) MinUnbalDispNorm.C, revA"


#include <MinUnbalDispNorm.h>
#include <AnalysisModel.h>
#include <LinearSOE.h>
#include <Vector.h>
#include <Channel.h>
#include <math.h>

MinUnbalDispNorm::MinUnbalDispNorm(double lambda1, int specNumIter,
		     double min, double max, int signFirstStep)
:StaticIntegrator(INTEGRATOR_TAGS_MinUnbalDispNorm),
 dLambda1LastStep(lambda1), 
 specNumIncrStep(specNumIter), numIncrLastStep(specNumIter),
 deltaUhat(0), deltaUbar(0), deltaU(0), deltaUstep(0), 
 phat(0), deltaLambdaStep(0.0), currentLambda(0.0), 
 signLastDeltaLambdaStep(1),
 dLambda1min(min), dLambda1max(max), signLastDeterminant(1), signFirstStepMethod(signFirstStep)
{
  // to avoid divide-by-zero error on first update() ensure numIncr != 0
  if (specNumIncrStep == 0) {
    opserr << "WARNING LoadControl::LoadControl() - numIncr set to 0, 1 assumed\n";
    specNumIncrStep = 1.0;
    numIncrLastStep = 1.0;
  }
}

MinUnbalDispNorm::~MinUnbalDispNorm()
{
    // delete any vector object created
    if (deltaUhat != 0)
	delete deltaUhat;
    if (deltaU != 0)
	delete deltaU;
    if (deltaUstep != 0)
	delete deltaUstep;
    if (deltaUbar != 0)
	delete deltaUbar;
    if (phat != 0)
	delete phat;
}

int
MinUnbalDispNorm::newStep(void)
{
    // get pointers to AnalysisModel and LinearSOE
    AnalysisModel *theModel = this->getAnalysisModelPtr();
    LinearSOE *theLinSOE = this->getLinearSOEPtr();    
    if (theModel == 0 || theLinSOE == 0) {
	opserr << "WARNING MinUnbalDispNorm::newStep() ";
	opserr << "No AnalysisModel or LinearSOE has been set\n";
	return -1;
    }

    // get the current load factor
    currentLambda = theModel->getCurrentDomainTime();


    // determine dUhat
    this->formTangent();
    theLinSOE->setB(*phat);
    if (theLinSOE->solve() < 0) {
      opserr << "MinUnbalanceDispNorm::newStep(void) - failed in solver\n";
      return -1;
    }
    (*deltaUhat) = theLinSOE->getX();
    Vector &dUhat = *deltaUhat;

    // determine delta lambda(1) == dlambda
    double factor = specNumIncrStep/numIncrLastStep;
    double dLambda = dLambda1LastStep*factor;

    // check aaint min and max values specified in constructor
    if (dLambda < dLambda1min)
      dLambda = dLambda1min;
    else if (dLambda > dLambda1max)
      dLambda = dLambda1max;

    dLambda1LastStep = dLambda;


    if (signFirstStepMethod == SIGN_LAST_STEP) {
      if (deltaLambdaStep < 0)
	signLastDeltaLambdaStep = -1;
      else
	signLastDeltaLambdaStep = +1;
      
      dLambda *= signLastDeltaLambdaStep; // base sign of load change
                                        // on what was happening last step
    } else {
    
      double det = theLinSOE->getDeterminant();
      int signDeterminant = 1;
      if (det < 0)
	signDeterminant = -1;
    
      dLambda *= signDeterminant * signLastDeterminant;
    
      signLastDeterminant = signDeterminant;
    }

    /*
    double work = (*phat)^(dUhat);
    int signCurrentWork = 1;
    if (work < 0) signCurrentWork = -1;

    if (signCurrentWork != signLastDeltaStep)
    */

    deltaLambdaStep = dLambda;
    currentLambda += dLambda;
    numIncrLastStep = 0;

    // determine delta U(1) == dU
    (*deltaU) = dUhat;
    (*deltaU) *= dLambda;
    (*deltaUstep) = (*deltaU);

    // update model with delta lambda and delta U
    theModel->incrDisp(*deltaU);    
    theModel->applyLoadDomain(currentLambda);    
    if (theModel->updateDomain() < 0) {
      opserr << "MinUnbalDispNorm::newStep - model failed to update for new dU\n";
      return -1;
    }

    return 0;
}

int
MinUnbalDispNorm::update(const Vector &dU)
{
    AnalysisModel *theModel = this->getAnalysisModelPtr();
    LinearSOE *theLinSOE = this->getLinearSOEPtr();    
    if (theModel == 0 || theLinSOE == 0) {
	opserr << "WARNING MinUnbalDispNorm::update() ";
	opserr << "No AnalysisModel or LinearSOE has been set\n";
	return -1;
    }

    (*deltaUbar) = dU; // have to do this as the SOE is gonna change

    // determine dUhat    
    theLinSOE->setB(*phat);
    theLinSOE->solve();
    (*deltaUhat) = theLinSOE->getX();    

    // determine delta lambda(i)
    double a = (*deltaUhat)^(*deltaUbar);
    double b = (*deltaUhat)^(*deltaUhat);
    if (b == 0) {
      opserr << "MinUnbalDispNorm::update() - zero denominator\n";
      return -1;
    }

    double dLambda = -a/b;
    
    // determine delta U(i)
    (*deltaU) = (*deltaUbar);    
    deltaU->addVector(1.0, *deltaUhat,dLambda);
    
    // update dU and dlambda
    (*deltaUstep) += *deltaU;
    deltaLambdaStep += dLambda;
    currentLambda += dLambda;

    // update the model
    theModel->incrDisp(*deltaU);    
    theModel->applyLoadDomain(currentLambda);    

    if (theModel->updateDomain() < 0) {
      opserr << "MinUnbalDispNorm::update - model failed to update for new dU\n";
      return -1;
    }
    
    // set the X soln in linearSOE to be deltaU for convergence Test
    theLinSOE->setX(*deltaU);

    numIncrLastStep++;
    return 0;
}



int 
MinUnbalDispNorm::domainChanged(void)
{
    // we first create the Vectors needed
    AnalysisModel *theModel = this->getAnalysisModelPtr();
    LinearSOE *theLinSOE = this->getLinearSOEPtr();    
    if (theModel == 0 || theLinSOE == 0) {
	opserr << "WARNING MinUnbalDispNorm::update() ";
	opserr << "No AnalysisModel or LinearSOE has been set\n";
	return -1;
    }    
    int size = theModel->getNumEqn(); // ask model in case N+1 space

    if (deltaUhat == 0 || deltaUhat->Size() != size) { // create new Vector
	if (deltaUhat != 0)
	    delete deltaUhat;   // delete the old
	deltaUhat = new Vector(size);
	if (deltaUhat == 0 || deltaUhat->Size() != size) { // check got it
	    opserr << "FATAL MinUnbalDispNorm::domainChanged() - ran out of memory for";
	    opserr << " deltaUhat Vector of size " << size << endln;
	    exit(-1);
	}
    }

    if (deltaUbar == 0 || deltaUbar->Size() != size) { // create new Vector
	if (deltaUbar != 0)
	    delete deltaUbar;   // delete the old
	deltaUbar = new Vector(size);
	if (deltaUbar == 0 || deltaUbar->Size() != size) { // check got it
	    opserr << "FATAL MinUnbalDispNorm::domainChanged() - ran out of memory for";
	    opserr << " deltaUbar Vector of size " << size << endln;
	    exit(-1);
	}
    }

    
    if (deltaU == 0 || deltaU->Size() != size) { // create new Vector
	if (deltaU != 0)
	    delete deltaU;   // delete the old
	deltaU = new Vector(size);
	if (deltaU == 0 || deltaU->Size() != size) { // check got it
	    opserr << "FATAL MinUnbalDispNorm::domainChanged() - ran out of memory for";
	    opserr << " deltaU Vector of size " << size << endln;
	    exit(-1);
	}
    }

    if (deltaUstep == 0 || deltaUstep->Size() != size) { 
	if (deltaUstep != 0)
	    delete deltaUstep;  
	deltaUstep = new Vector(size);
	if (deltaUstep == 0 || deltaUstep->Size() != size) { 
	    opserr << "FATAL MinUnbalDispNorm::domainChanged() - ran out of memory for";
	    opserr << " deltaUstep Vector of size " << size << endln;
	    exit(-1);
	}
    }

    if (phat == 0 || phat->Size() != size) { 
	if (phat != 0)
	    delete phat;  
	phat = new Vector(size);
	if (phat == 0 || phat->Size() != size) { 
	    opserr << "FATAL MinUnbalDispNorm::domainChanged() - ran out of memory for";
	    opserr << " phat Vector of size " << size << endln;
	    exit(-1);
	}
    }    

    // now we have to determine phat
    // do this by incrementing lambda by 1, applying load
    // and getting phat from unbalance.
    currentLambda = theModel->getCurrentDomainTime();
    currentLambda += 1.0;
    theModel->applyLoadDomain(currentLambda);    
    this->formUnbalance(); // NOTE: this assumes unbalance at last was 0
    (*phat) = theLinSOE->getB();
    currentLambda -= 1.0;
    theModel->setCurrentDomainTime(currentLambda);    


    // check there is a reference load
    int haveLoad = 0;
    for (int i=0; i<size; i++)
      if ( (*phat)(i) != 0.0 ) {
	haveLoad = 1;
	i = size;
      }

    if (haveLoad == 0) {
      opserr << "WARNING ArcLength::domainChanged() - zero reference load";
      return -1;
    }

    
    return 0;
}

int
MinUnbalDispNorm::sendSelf(int cTag,
		    Channel &theChannel)
{
  Vector data(8);
  data(0) = dLambda1LastStep;
  data(1) = specNumIncrStep;
  data(2) = numIncrLastStep;
  data(3) = deltaLambdaStep;
  data(4) = currentLambda;
  if (signLastDeltaLambdaStep == 1)
    data(5)  = 1.0;
  else
    data(5) = 0.0;
  data(6) = dLambda1min;
  data(7) = dLambda1max;

  if (theChannel.sendVector(this->getDbTag(), cTag, data) < 0) {
      opserr << "MinUnbalDispNorm::sendSelf() - failed to send the data\n";
      return -1;
  }
  return 0;
}


int
MinUnbalDispNorm::recvSelf(int cTag,
		    Channel &theChannel, FEM_ObjectBroker &theBroker)
{
  Vector data(8);
  if (theChannel.recvVector(this->getDbTag(), cTag, data) < 0) {
      opserr << "MinUnbalDispNorm::sendSelf() - failed to send the data\n";
      return -1;
  }      

  // set the data

  dLambda1LastStep = data(0);
  specNumIncrStep = data(1);
  numIncrLastStep = data(2);
  deltaLambdaStep = data(3);
  currentLambda = data(4);
  if (data(5)== 1.0)
    signLastDeltaLambdaStep = 1;
  else
    signLastDeltaLambdaStep = -1;
  dLambda1min = data(6);
  dLambda1max = data(7);

  return 0;
}

void
MinUnbalDispNorm::Print(OPS_Stream &s, int flag)
{
    AnalysisModel *theModel = this->getAnalysisModelPtr();
    if (theModel != 0) {
	double cLambda = theModel->getCurrentDomainTime();
	s << "\t MinUnbalDispNorm - currentLambda: " << cLambda;
    } else 
	s << "\t MinUnbalDispNorm - no associated AnalysisModel\n";
}
