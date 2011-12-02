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
// $Date: 2003-02-14 23:00:46 $
// $Source: /usr/local/cvs/OpenSees/SRC/analysis/integrator/ArcLength1.cpp,v $
                                                                        
                                                                        
// File: ~/analysis/integrator/ArcLength1.C
// 
// Written: fmk 
// Created: 07/98
// Revision: A
//
// Description: This file contains the class definition for ArcLength1.
// ArcLength1 is an algorithmic class for perfroming a static analysis
// using the arc length scheme, that is within a load step the follwing
// constraint is enforced: dU^TdU + alpha^2*dLambda^2 = ArcLength1^2
// where dU is change in nodal displacements for step, dLambda is
// change in applied load and ArcLength1 is a control parameter.
//
// What: "@(#) ArcLength1.C, revA"


#include <ArcLength1.h>
#include <AnalysisModel.h>
#include <LinearSOE.h>
#include <Vector.h>
#include <Channel.h>
#include <math.h>

ArcLength1::ArcLength1(double arcLength, double alpha)
:StaticIntegrator(INTEGRATOR_TAGS_ArcLength1),
 arcLength2(arcLength*arcLength), alpha2(alpha*alpha),
 deltaUhat(0), deltaUbar(0), deltaU(0), deltaUstep(0), 
 phat(0), deltaLambdaStep(0.0), currentLambda(0.0), 
 signLastDeltaLambdaStep(1)
{

}

ArcLength1::~ArcLength1()
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
ArcLength1::newStep(void)
{
    // get pointers to AnalysisModel and LinearSOE
    AnalysisModel *theModel = this->getAnalysisModelPtr();
    LinearSOE *theLinSOE = this->getLinearSOEPtr();    
    if (theModel == 0 || theLinSOE == 0) {
	opserr << "WARNING ArcLength1::newStep() ";
	opserr << "No AnalysisModel or LinearSOE has been set\n";
	return -1;
    }

    // get the current load factor
    currentLambda = theModel->getCurrentDomainTime();

    if (deltaLambdaStep < 0)
	signLastDeltaLambdaStep = -1;
    else
	signLastDeltaLambdaStep = +1;

    // determine dUhat
    this->formTangent();
    theLinSOE->setB(*phat);
    theLinSOE->solve();
    (*deltaUhat) = theLinSOE->getX();
    Vector &dUhat = *deltaUhat;
    
    // determine delta lambda(1) == dlambda
    double dLambda = sqrt(arcLength2/((dUhat^dUhat)+alpha2));
    dLambda *= signLastDeltaLambdaStep; // base sign of load change
                                        // on what was happening last step
    deltaLambdaStep = dLambda;
    currentLambda += dLambda;

    // determine delta U(1) == dU
    (*deltaU) = dUhat;
    (*deltaU) *= dLambda;
    (*deltaUstep) = (*deltaU);

    // update model with delta lambda and delta U
    theModel->incrDisp(*deltaU);    
    theModel->applyLoadDomain(currentLambda);    
    theModel->updateDomain();

    return 0;
}

int
ArcLength1::update(const Vector &dU)
{
    AnalysisModel *theModel = this->getAnalysisModelPtr();
    LinearSOE *theLinSOE = this->getLinearSOEPtr();    
    if (theModel == 0 || theLinSOE == 0) {
	opserr << "WARNING ArcLength1::update() ";
	opserr << "No AnalysisModel or LinearSOE has been set\n";
	return -1;
    }

    (*deltaUbar) = dU; // have to do this as the SOE is gonna change

    // determine dUhat    
    theLinSOE->setB(*phat);
    theLinSOE->solve();
    (*deltaUhat) = theLinSOE->getX();    

    // determine delta lambda(i)
    double a = (*deltaUstep)^(*deltaUbar);
    double b = (*deltaUstep)^(*deltaUhat) + alpha2*deltaLambdaStep;
    if (b == 0) {
      opserr << "ArcLength1::update() - zero denominator,";
      opserr << " alpha was set to 0.0 and zero reference load\n";
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
    theModel->updateDomain();
    
    // set the X soln in linearSOE to be deltaU for convergence Test
    theLinSOE->setX(*deltaU);

    return 0;
}



int 
ArcLength1::domainChanged(void)
{
    // we first create the Vectors needed
    AnalysisModel *theModel = this->getAnalysisModelPtr();
    LinearSOE *theLinSOE = this->getLinearSOEPtr();    
    if (theModel == 0 || theLinSOE == 0) {
	opserr << "WARNING ArcLength1::update() ";
	opserr << "No AnalysisModel or LinearSOE has been set\n";
	return -1;
    }    
    int size = theModel->getNumEqn(); // ask model in case N+1 space

    if (deltaUhat == 0 || deltaUhat->Size() != size) { // create new Vector
	if (deltaUhat != 0)
	    delete deltaUhat;   // delete the old
	deltaUhat = new Vector(size);
	if (deltaUhat == 0 || deltaUhat->Size() != size) { // check got it
	    opserr << "FATAL ArcLength1::domainChanged() - ran out of memory for";
	    opserr << " deltaUhat Vector of size " << size << endln;
	    exit(-1);
	}
    }

    if (deltaUbar == 0 || deltaUbar->Size() != size) { // create new Vector
	if (deltaUbar != 0)
	    delete deltaUbar;   // delete the old
	deltaUbar = new Vector(size);
	if (deltaUbar == 0 || deltaUbar->Size() != size) { // check got it
	    opserr << "FATAL ArcLength1::domainChanged() - ran out of memory for";
	    opserr << " deltaUbar Vector of size " << size << endln;
	    exit(-1);
	}
    }

    
    if (deltaU == 0 || deltaU->Size() != size) { // create new Vector
	if (deltaU != 0)
	    delete deltaU;   // delete the old
	deltaU = new Vector(size);
	if (deltaU == 0 || deltaU->Size() != size) { // check got it
	    opserr << "FATAL ArcLength1::domainChanged() - ran out of memory for";
	    opserr << " deltaU Vector of size " << size << endln;
	    exit(-1);
	}
    }

    if (deltaUstep == 0 || deltaUstep->Size() != size) { 
	if (deltaUstep != 0)
	    delete deltaUstep;  
	deltaUstep = new Vector(size);
	if (deltaUstep == 0 || deltaUstep->Size() != size) { 
	    opserr << "FATAL ArcLength1::domainChanged() - ran out of memory for";
	    opserr << " deltaUstep Vector of size " << size << endln;
	    exit(-1);
	}
    }

    if (phat == 0 || phat->Size() != size) { 
	if (phat != 0)
	    delete phat;  
	phat = new Vector(size);
	if (phat == 0 || phat->Size() != size) { 
	    opserr << "FATAL ArcLength1::domainChanged() - ran out of memory for";
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
    
    return 0;
}

int
ArcLength1::sendSelf(int cTag,
		    Channel &theChannel)
{
  Vector data(5);
  data(0) = arcLength2;
  data(1) = alpha2;
  data(2) = deltaLambdaStep;
  data(3) = currentLambda;
  data(4)  = signLastDeltaLambdaStep;

  if (theChannel.sendVector(this->getDbTag(), cTag, data) < 0) {
      opserr << "ArcLength1::sendSelf() - failed to send the data\n";
      return -1;
  }
  return 0;
}


int
ArcLength1::recvSelf(int cTag,
		    Channel &theChannel, FEM_ObjectBroker &theBroker)
{
  Vector data(5);
  if (theChannel.recvVector(this->getDbTag(), cTag, data) < 0) {
      opserr << "ArcLength1::sendSelf() - failed to send the data\n";
      return -1;
  }      

  // set the data
  arcLength2 = data(0);
  alpha2 = data(1);
  deltaLambdaStep = data(2);
  currentLambda = data(3);
  signLastDeltaLambdaStep = data(4);
  return 0;
}

void
ArcLength1::Print(OPS_Stream &s, int flag)
{
    AnalysisModel *theModel = this->getAnalysisModelPtr();
    if (theModel != 0) {
	double cLambda = theModel->getCurrentDomainTime();
	s << "\t ArcLength1 - currentLambda: " << cLambda;
	s << "  ArcLength1: " << sqrt(arcLength2) <<  "  alpha: ";
	s << sqrt(alpha2) << endln;
    } else 
	s << "\t ArcLength1 - no associated AnalysisModel\n";
}








