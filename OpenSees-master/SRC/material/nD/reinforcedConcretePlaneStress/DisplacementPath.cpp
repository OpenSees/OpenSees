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
                                                                        
                                                                
                                                                       
// Written: JZhong 
// Created: 01/04
//
// Description: This file contains the class definition for DisplacementPath.
// DisplacementPath is an algorithmic class for performing a static analysis
// using the user-defined displacement path.
// For Detailed explanation of the model, please refer to the book
// entitled "Unified Theory of Concrete Structures,"
// by Thomas T.C. Hsu and Y.L. Mo, John Wiley & Sons, April 2010.

#include <DisplacementPath.h>
#include <AnalysisModel.h>
#include <LinearSOE.h>
#include <Vector.h>
#include <Channel.h>
#include <math.h>
#include <Domain.h>
#include <Node.h>
#include <DOF_Group.h>
#include <ID.h>
#include <stdlib.h>


DisplacementPath::DisplacementPath(int node, int dof, 
					 Vector &incrementVector, 
					 Domain *domain) 
:StaticIntegrator(INTEGRATOR_TAGS_DisplacementPath),
 theNode(node), theDof(dof), theIncrementVector(0), theDomain(domain),
 theDofID(-1),
 deltaUhat(0), deltaUbar(0), deltaU(0), deltaUstep(0), 
 phat(0), deltaLambdaStep(0.0), currentLambda(0.0),
 currentStep(0), theCurrentIncrement(0.0)
{
	// get theIncrementVector
	theIncrementVector = new Vector(incrementVector);
	if ( theIncrementVector==0 || theIncrementVector->Size()==0 ) {
		opserr << "DisplacementPath::DisplacementPath() - ran out of memory\n";
		exit(-1);
	}

	//opserr << " theIncrementVector = \n";
	//for (int i =0; i<theIncrementVector->Size(); i++) {		
	//	opserr << (*theIncrementVector)(i) << endln;
	//}

}

DisplacementPath::~DisplacementPath()
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

	if (theIncrementVector !=0)
	delete theIncrementVector;
}

int
DisplacementPath::newStep(void)
{
	if (theDofID == -1) {
		opserr << "DisplacementPath::newStep() - domainChanged has not been called\n";
		return -1;
	}

    // get pointers to AnalysisModel and LinearSOE
    AnalysisModel *theModel = this->getAnalysisModel();
    LinearSOE *theLinSOE = this->getLinearSOE();    
    if (theModel == 0 || theLinSOE == 0) {
	opserr << "WARNING DisplacementPath::newStep() ";
	opserr << "No AnalysisModel or LinearSOE has been set\n";
	return -1;
    }

	// check theIncrementVector Vector
	if ( theIncrementVector == 0 ) {
		opserr << "DisplacementPath::newStep() - no theIncrementVector associated with object\n";
		return -2;
	}


    // determine increment for this iteration
    if (currentStep < theIncrementVector->Size()) {
		theCurrentIncrement = (*theIncrementVector)(currentStep);
	}
	else {
		theCurrentIncrement = 0.0;
		opserr << "DisplacementPath::newStep() - reach the end of specified load path\n";
		opserr << " - setting theCurrentIncrement = 0.0\n";
	}


    // get the current load factor
    currentLambda = theModel->getCurrentDomainTime();


    // determine dUhat and dUabar
    this->formTangent();
    this->formUnbalance();

	(*deltaUbar) = theLinSOE->getX();
	double dUabar = (*deltaUbar)(theDofID);


    theLinSOE->setB(*phat);

    if (theLinSOE->solve() < 0) {
      opserr << "DisplacementControl::newStep(void) - failed in solver\n";
      return -1;
    }
    
	
	(*deltaUhat) = theLinSOE->getX();
    Vector &dUhat = *deltaUhat;

    double dUahat = dUhat(theDofID);

	//opserr << " newStep( ) " << endln;
    //opserr << " theDofID = " << theDofID << endln;
	//opserr << "dUahat = " << dUahat << endln;

    if (dUahat == 0.0) {
	opserr << "WARNING DisplacementPath::newStep() ";
	opserr << "dUahat is zero -- zero reference displacement at control node DOF\n";
    
	opserr << "currentStep = " << currentStep << endln; // add by zhong
	opserr << " theCurrentIncrement = " << theCurrentIncrement << endln; // zhong

	return -1;
    }
    

    // determine delta lambda(1) == dlambda    
    double dLambda = (theCurrentIncrement-dUabar)/dUahat;

    deltaLambdaStep = dLambda;
    currentLambda += dLambda;
 //   opserr << "DisplacementPath: " << dUahat  << " " << theDofID << endln;
 //   opserr << "DisplacementPath::newStep() : " << deltaLambdaStep << endln;
    // determine delta U(1) == dU
    (*deltaU) = dUhat;
    (*deltaU) *= dLambda;
    (*deltaUstep) = (*deltaU);

    // update model with delta lambda and delta U
    theModel->incrDisp(*deltaU);    
    theModel->applyLoadDomain(currentLambda);    
    if (theModel->updateDomain() < 0) {
      opserr << "DisplacementPath::newStep - model failed to update for new dU\n";
      return -1;
    }

    currentStep++;

    return 0;
}

int
DisplacementPath::update(const Vector &dU)
{
    // opserr << " update is invoked " << endln;
	if (theDofID == -1) {
		opserr << "DisplacementControl::newStep() - domainChanged has not been called\n";
		return -1;
	}
    AnalysisModel *theModel = this->getAnalysisModel();
    LinearSOE *theLinSOE = this->getLinearSOE();    
    if (theModel == 0 || theLinSOE == 0) {
	opserr << "WARNING DisplacementPath::update() ";
	opserr << "No AnalysisModel or LinearSOE has been set\n";
	return -1;
    }

    (*deltaUbar) = dU; // have to do this as the SOE is gonna change
    double dUabar = (*deltaUbar)(theDofID);
    
    // determine dUhat    
    theLinSOE->setB(*phat);
    theLinSOE->solve();
    (*deltaUhat) = theLinSOE->getX();    

	// add by zhong for check purpose
    //int size = deltaUhat->Size();
	//opserr << "\n size of deltaUhat = " << size << endln;
    //for (int i=0; i<size; i++) {
    //   opserr << " dektaUhat(i) = " << (*deltaUhat)(i) << endln;
    //}
	// finish here

    double dUahat = (*deltaUhat)(theDofID);
	//opserr << " theDofID = " << theDofID << endln;
	
	//opserr << "update( ) " << endln;
	//opserr << "dUahat = " << dUahat << endln;

    if (dUahat == 0.0) {
	opserr << "WARNING DisplacementPath::update() ";
	opserr << "dUahat is zero -- zero reference displacement at control node DOF\n";
	return -1;
    }
    
    // determine delta lambda(1) == dlambda    
    double dLambda = -dUabar/dUahat;

    // add by zhong
	//opserr << "\n dUahat = " << dUahat << endln;
	//opserr << " dUabar = " << dUabar << endln;
	//opserr << " dLambda = " << dLambda << endln;
	// finish
    
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
      opserr << "DisplacementPath::update - model failed to update for new dU\n";
      return -1;
    }
	
    
    // set the X soln in linearSOE to be deltaU for convergence Test
    theLinSOE->setX(*deltaU);


    return 0;
}



int 
DisplacementPath::domainChanged(void)
{
    // we first create the Vectors needed
    AnalysisModel *theModel = this->getAnalysisModel();
    LinearSOE *theLinSOE = this->getLinearSOE();    
    if (theModel == 0 || theLinSOE == 0) {
	opserr << "WARNING DisplacementPath::update() ";
	opserr << "No AnalysisModel or LinearSOE has been set\n";
	return -1;
    }    
    int size = theModel->getNumEqn(); // ask model in case N+1 space

    if (deltaUhat == 0 || deltaUhat->Size() != size) { // create new Vector
	if (deltaUhat != 0)
	    delete deltaUhat;   // delete the old
	deltaUhat = new Vector(size);
	if (deltaUhat == 0 || deltaUhat->Size() != size) { // check got it
	    opserr << "FATAL DisplacementPath::domainChanged() - ran out of memory for";
	    opserr << " deltaUhat Vector of size " << size << endln;
	    exit(-1);
	}
    }

    if (deltaUbar == 0 || deltaUbar->Size() != size) { // create new Vector
	if (deltaUbar != 0)
	    delete deltaUbar;   // delete the old
	deltaUbar = new Vector(size);
	if (deltaUbar == 0 || deltaUbar->Size() != size) { // check got it
	    opserr << "FATAL DisplacementPath::domainChanged() - ran out of memory for";
	    opserr << " deltaUbar Vector of size " << size << endln;
	    exit(-1);
	}
    }
    
    if (deltaU == 0 || deltaU->Size() != size) { // create new Vector
	if (deltaU != 0)
	    delete deltaU;   // delete the old
	deltaU = new Vector(size);
	if (deltaU == 0 || deltaU->Size() != size) { // check got it
	    opserr << "FATAL DisplacementPath::domainChanged() - ran out of memory for";
	    opserr << " deltaU Vector of size " << size << endln;
	    exit(-1);
	}
    }

    if (deltaUstep == 0 || deltaUstep->Size() != size) { 
	if (deltaUstep != 0)
	    delete deltaUstep;  
	deltaUstep = new Vector(size);
	if (deltaUstep == 0 || deltaUstep->Size() != size) { 
	    opserr << "FATAL DisplacementPath::domainChanged() - ran out of memory for";
	    opserr << " deltaUstep Vector of size " << size << endln;
	    exit(-1);
	}
    }

    if (phat == 0 || phat->Size() != size) { 
	if (phat != 0)
	    delete phat;  
	phat = new Vector(size);
	if (phat == 0 || phat->Size() != size) { 
	    opserr << "FATAL DisplacementPath::domainChanged() - ran out of memory for";
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
      opserr << "WARNING DisplacementPath::domainChanged() - zero reference load";
      return -1;
    }

    // lastly we determine the id of the nodal dof
    // EXTRA CODE TO DO SOME ERROR CHECKING REQUIRED
    
    Node *theNodePtr = theDomain->getNode(theNode);
	if (theNodePtr == 0) {
		opserr << "DisplacementControl::domainChanged - no node\n";
		return -1;
	}

    DOF_Group *theGroup = theNodePtr->getDOF_GroupPtr();
	if (theGroup == 0) {
	   return 0;
	}
	const ID &theID = theGroup->getID();
    theDofID = theID(theDof);
    
    return 0;
}

int
DisplacementPath::sendSelf(int cTag,
		    Channel &theChannel)
{
    // TO FINISH
  return 0;
}


int
DisplacementPath::recvSelf(int cTag,
		    Channel &theChannel, FEM_ObjectBroker &theBroker)
{
    // TO FINISH
  return 0;
}

void
DisplacementPath::Print(OPS_Stream &s, int flag)
{
    // TO FINISH    
}
