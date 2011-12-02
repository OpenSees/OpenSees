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
// $Date: 2003-11-18 01:59:04 $
// $Source: /usr/local/cvs/OpenSees/SRC/domain/subdomain/ActorSubdomain.cpp,v $
                                                                        
                                                                        

#include <ActorSubdomain.h>
#include <FEM_ObjectBroker.h>
#include <Element.h>
#include <Node.h>
#include <SP_Constraint.h>
#include <MP_Constraint.h>
#include <ElementalLoad.h>
#include <NodalLoad.h>
#include <LoadPattern.h>
#include <Matrix.h>
#include <Vector.h>
#include <DomainDecompositionAnalysis.h>
#include <PartitionedModelBuilder.h>

#include <EquiSolnAlgo.h>
#include <IncrementalIntegrator.h>
#include <LinearSOE.h>

#include <ArrayOfTaggedObjects.h>
#include <ShadowActorSubdomain.h>

ActorSubdomain::ActorSubdomain(Channel &theChannel,
			       FEM_ObjectBroker &theBroker)
:Subdomain(0), Actor(theChannel,theBroker,0),
 msgData(4),lastResponse(0)
{
  // does nothing
}
    
ActorSubdomain::~ActorSubdomain()
{
    // does nothing
    opserr << "ActorSubdomain DIED: " << endln;

}


int
ActorSubdomain::run(void)
{
    Vector theVect(4);
    bool exitYet = false;
    opserr << "ActorSubdomain running: " << endln;

    while (exitYet == false) {

  	this->recvID(msgData);
	int action = msgData(0);

	//		opserr << "ACTORSUBDOMAIN: " << action << endln;
	
	int theType, theOtherType, tag, dbTag, loadPatternTag;
	Element *theEle;
	Node *theNod;
	SP_Constraint *theSP;
	MP_Constraint *theMP;
	LoadPattern *theLoadPattern;
	NodalLoad *theNodalLoad;
	ElementalLoad *theElementalLoad;
	DomainDecompositionAnalysis *theDDAnalysis;
	const Matrix *theMatrix;
	const Vector *theVector;
	PartitionedModelBuilder *theBuilder;
	IncrementalIntegrator *theIntegrator;
	EquiSolnAlgo *theAlgorithm;
	LinearSOE *theSOE;

	const ID *theID;
	
	switch (action) {

	  case ShadowActorSubdomain_setTag:
	    tag = msgData(1); // subdomain tag
	    this->setTag(tag);

	    break;

	  case ShadowActorSubdomain_newStep:
	    this->recvVector(theVect);
	    this->newStep(theVect(0));

	    break;

	  case ShadowActorSubdomain_buildSubdomain:
	    theType = msgData(1);
	    tag = msgData(3); // subdomain tag
	    this->setTag(tag);
	    tag = msgData(2); // numSubdomains
	    theBuilder = theBroker->getPtrNewPartitionedModelBuilder(*this, 
								     theType);
	    this->recvObject(*theBuilder);
	    this->buildSubdomain(tag, *theBuilder);

	    break;

	case ShadowActorSubdomain_getRemoteData:
	    theID = &(this->getExternalNodes());
	    msgData(0) = theID->Size();
	    msgData(1) = this->getNumDOF();

	    this->sendID(msgData);
	    if (theID->Size() != 0)
	      this->sendID(*theID);
	    break;

	  case ShadowActorSubdomain_getCost:
       	    theVect(0) = this->getCost(); // have to use [] for Sun's CC!
	    this->sendVector(theVect);
	    break;	    

 	  case ShadowActorSubdomain_addElement:
	    theType = msgData(1);
	    dbTag = msgData(2);

	    theEle = theBroker->getNewElement(theType);

	    if (theEle != 0) {
		theEle->setDbTag(dbTag);		
		this->recvObject(*theEle);
		bool result = this->addElement(theEle);
		if (result == true)
		    msgData(0) = 0;
		else
		    msgData(0) = -1;
	    } else
		msgData(0) = -1;

	    /*
	    this->recvID(msgData);	    
	    opserr << "ActorSubdomain::addElement() : " << msgData;
	    
	    msgData(0) = 1;
	    msgData(1) = 2;
	    msgData(2) = 3;
	    msgData(3) = 4;
	    this->sendID(msgData);	    
	    */

	    break;

	    
	    
	  case ShadowActorSubdomain_addNode:
	    theType = msgData(1);
	    dbTag = msgData(2);
	    
	    theNod = theBroker->getNewNode(theType);

	    if (theNod != 0) {
		theNod->setDbTag(dbTag);		
		this->recvObject(*theNod);
		bool result = this->addNode(theNod);
		if (result == true)
		    msgData(0) = 0;
		else
		    msgData(0) = -1;
	    } else
		msgData(0) = -1;

	    break;
	    

	    
	  case ShadowActorSubdomain_addExternalNode:
	    theType = msgData(1);
	    dbTag = msgData(2);
	    theNod = theBroker->getNewNode(theType);

	    if (theNod != 0) {
		theNod->setDbTag(dbTag);
		this->recvObject(*theNod);
		bool result = this->addExternalNode(theNod);
		if (result == true)
		    msgData(0) = 0;
		else
		    msgData(0) = -1;
	    } else
		msgData(0) = -1;

	    break;	    

	    
	  case ShadowActorSubdomain_addSP_Constraint:
	    theType = msgData(1);
	    dbTag = msgData(2);
	    
	    theSP = theBroker->getNewSP(theType);
	    
	    if (theSP != 0) {
		theSP->setDbTag(dbTag);
		this->recvObject(*theSP);
		bool result = this->addSP_Constraint(theSP);
		if (result == true)
		    msgData(0) = 0;
		else
		    msgData(0) = -1;
	    } else
		msgData(0) = -1;

	    break;	    
	    
	  case ShadowActorSubdomain_addMP_Constraint:
	    theType = msgData(1);
	    dbTag = msgData(2);
	    theMP = theBroker->getNewMP(theType);

	    if (theMP != 0) {
		theMP->setDbTag(dbTag);
		this->recvObject(*theMP);
		bool result = this->addMP_Constraint(theMP);
		if (result == true)
		    msgData(0) = 0;
		else
		    msgData(0) = -1;
	    } else
		msgData(0) = -1;

	    break;	    
	    
	    
	  case ShadowActorSubdomain_addLoadPattern:
	    theType = msgData(1);
	    dbTag = msgData(2);
	    
	    theLoadPattern = theBroker->getNewLoadPattern(theType);

	    if (theLoadPattern != 0) {
		theLoadPattern->setDbTag(dbTag);
		this->recvObject(*theLoadPattern);
		bool result = this->addLoadPattern(theLoadPattern);
		if (result == true)
		    msgData(0) = 0;
		else
		    msgData(0) = -1;
	    } else
		msgData(0) = -1;

	    break;	    	    

	  case ShadowActorSubdomain_addNodalLoadToPattern:
 	    theType = msgData(1);
	    dbTag = msgData(2);
	    loadPatternTag = msgData(3);
	    
	    theNodalLoad = theBroker->getNewNodalLoad(theType);

	    if (theNodalLoad != 0) {
		theNodalLoad->setDbTag(dbTag);
		this->recvObject(*theNodalLoad);
		bool result = this->addNodalLoad(theNodalLoad, loadPatternTag);
		if (result == true)
		    msgData(0) = 0;
		else
		    msgData(0) = -1;
	    } else
		msgData(0) = -1;

	    break;	    
	    
	    
	  case ShadowActorSubdomain_addElementalLoadToPattern:
	    theType = msgData(1);
	    dbTag = msgData(2);
	    loadPatternTag = msgData(3);
	    
	    theElementalLoad = theBroker->getNewElementalLoad(theType);

	    if (theElementalLoad != 0) {
		theElementalLoad->setDbTag(dbTag);
		this->recvObject(*theElementalLoad);
		bool result = this->addElementalLoad(theElementalLoad, 
						     loadPatternTag);
		if (result == true)
		    msgData(0) = 0;
		else
		    msgData(0) = -1;
	    } else
		msgData(0) = -1;

	    break;	    	    
	    
	  case ShadowActorSubdomain_addSP_ConstraintToPattern:
	    theType = msgData(1);
	    dbTag = msgData(2);
	    loadPatternTag = msgData(3);
	    
	    theSP = theBroker->getNewSP(theType);

	    if (theSP != 0) {
		theSP->setDbTag(dbTag);
		this->recvObject(*theSP);
		bool result = this->addSP_Constraint(theSP, loadPatternTag);

		if (result == true)
		    msgData(0) = 0;
		else
		    msgData(0) = -1;
	    } else
		msgData(0) = -1;

	    break;	    	    	    

	  case ShadowActorSubdomain_removeElement:
	    tag = msgData(1);

	    theEle = this->removeElement(tag);

	    if (theEle != 0) 
		msgData(0) = theEle->getClassTag();
	    else
		msgData(0) = -1;

	    this->sendID(msgData);
	    if (theEle != 0) {
		this->sendObject(*theEle);
		delete theEle;
	    }

	    msgData(0) = 0;

	    break;	    	    	    


	  case ShadowActorSubdomain_removeNode:
	    tag = msgData(1);

	    theNod = this->removeNode(tag);

	    if (theNod != 0) 
		msgData(0) = theNod->getClassTag();
	    else
		msgData(0) = -1;

	    this->sendID(msgData);
	    if (theNod != 0) {
		this->sendObject(*theNod);
		delete theNod;
	    }

	    msgData(0) = 0;

	    break;

	  case ShadowActorSubdomain_removeSP_Constraint:
	    tag = msgData(1);

	    theSP = this->removeSP_Constraint(tag);

	    break;	    
	    
	  case ShadowActorSubdomain_removeMP_Constraint:
	    tag = msgData(1);

	    theMP = this->removeMP_Constraint(tag);

	    break;	    	    

	  case ShadowActorSubdomain_removeLoadPattern:
	    tag = msgData(1);

	    theLoadPattern = this->removeLoadPattern(tag);

	    break;	    	    
	    
	  case ShadowActorSubdomain_removeNodalLoadFromPattern:
	    tag = msgData(1);
	    theType = msgData(2);

	    theNodalLoad = this->removeNodalLoad(tag, theType);

	    break;	    	    	    

	  case ShadowActorSubdomain_removeElementalLoadFromPattern:
	    tag = msgData(1);
	    theType = msgData(2);

	    theElementalLoad = this->removeElementalLoad(tag, theType);

	    break;	    	    	    

	  case ShadowActorSubdomain_removeSP_ConstraintFromPattern:
	    tag = msgData(1);
	    theType = msgData(2);

	    theSP = this->removeSP_Constraint(tag, theType);

	    break;	    	    	    
	    
	    
	    
	  case ShadowActorSubdomain_getElement:
	    tag = msgData(1);

	    theEle = this->getElement(tag);

	    if (theEle != 0) 
		msgData(0) = theEle->getClassTag();
	    else
		msgData(0) = -1;

	    this->sendID(msgData);
	    if (theEle != 0) {
		this->sendObject(*theEle);
	    }

	    msgData(0) = 0;

	    break;	    	    	    	    


	  case ShadowActorSubdomain_Print:
	    this->Print(opserr);
	    this->sendID(msgData);

	    break;	    	    	    	    

	  case ShadowActorSubdomain_applyLoad:
	    this->recvVector(theVect);	    
	    this->applyLoad(theVect(0));
	    break;
	    
	  case ShadowActorSubdomain_setLoadConstant:
	    this->setLoadConstant();
	    break;	    
	    
	  case ShadowActorSubdomain_update:
	    this->update();
	    break;

	  case ShadowActorSubdomain_computeNodalResponse:
	    tag = msgData(1);
	    if (lastResponse == 0)
		lastResponse = new Vector(tag);
	    else if (lastResponse->Size() != tag) {
		delete lastResponse;
		lastResponse = new Vector(tag);
	    }
	    this->recvVector(*lastResponse);
	    this->computeNodalResponse();


	    
	  case ShadowActorSubdomain_commit:
	    this->commit();
	    break;
	    
	  case ShadowActorSubdomain_revertToLastCommit:
	    this->revertToLastCommit();
	    break;	    
	    
	  case ShadowActorSubdomain_revertToStart:
	    this->revertToStart();
	    break;	    	    
	    

	  case ShadowActorSubdomain_setDomainDecompAnalysis:
	    theType = msgData(1);
	    theDDAnalysis = 
		theBroker->getNewDomainDecompAnalysis(theType, *this);

	    if (theDDAnalysis != 0) {
		this->recvObject(*theDDAnalysis);
		this->setDomainDecompAnalysis(*theDDAnalysis);
		msgData(0) = 0;
	    } else
		msgData(0) = -1;
	    
	    break;

	case ShadowActorSubdomain_setAnalysisAlgorithm:
	  theType = msgData(1);
	  theAlgorithm = theBroker->getNewEquiSolnAlgo(theType);
	  if (theAlgorithm != 0) {
	    this->recvObject(*theAlgorithm);
	    this->setAnalysisAlgorithm(*theAlgorithm);
	    msgData(0) = 0;
	  } else
	    msgData(0) = -1;
	    
	  break;
	  
	case ShadowActorSubdomain_setAnalysisIntegrator:
	  theType = msgData(1);
	  theIntegrator = theBroker->getNewIncrementalIntegrator(theType);
	  if (theIntegrator != 0) {
	    this->recvObject(*theAlgorithm);
	    this->setAnalysisIntegrator(*theIntegrator);
	    msgData(0) = 0;
	  } else
	    msgData(0) = -1;
	    
	  break;

	case ShadowActorSubdomain_setAnalysisLinearSOE:
	  theType = msgData(1);
	  theOtherType = msgData(2);
	  theSOE = theBroker->getNewLinearSOE(theType, theOtherType);
	  if (theSOE != 0) {
	    this->recvObject(*theAlgorithm);
	    this->setAnalysisLinearSOE(*theSOE);
	    msgData(0) = 0;
	  } else
	    msgData(0) = -1;
	    
	  break;
	    
	  case ShadowActorSubdomain_domainChange:
	    this->domainChange();

	    tag = this->getNumDOF();
	    if (lastResponse == 0)
		lastResponse = new Vector(tag);
	    else if (lastResponse->Size() != tag) {
		delete lastResponse;
		lastResponse = new Vector(tag);
	    }
	    
	    break;

	  case ShadowActorSubdomain_clearAnalysis:
//	    this->clearAnalysis();
	    break;
	  /*
	  case 50:
	    const Matrix *theMatrix1 = &(this->getStiff());
	    this->sendMatrix(*theMatrix1);
	    break;

	  case 51:
	    const Matrix *theMatrix2 = &(this->getDamp());
	    this->sendMatrix(*theMatrix2);
	    break;
	    
	  case 52:
	    const Matrix *theMatrix3 = &(this->getMass());
	    this->sendMatrix(*theMatrix3);
	    break;	    
	    */
	  case  ShadowActorSubdomain_getTang:
	    theMatrix = &(this->getTang());
	    this->sendMatrix(*theMatrix);
	    break;	    
	    
	  case ShadowActorSubdomain_getResistingForce:
	    theVector = &(this->getResistingForce());
	    this->sendVector(*theVector);
	    break;	    	    

	  case ShadowActorSubdomain_computeTang:
	    tag = msgData(1);
	    this->setTag(tag);
	    this->computeTang();
	    break;


	  case ShadowActorSubdomain_computeResidual:
	    this->computeResidual();
	    break;

	  case ShadowActorSubdomain_DIE:
	    exitYet = true;
	    break;

	  default:
	    opserr << "ActorSubdomain::invalid action " << action << "received\n";
	    msgData(0) = -1;
	}
    }
    return 0;
}


Node *
ActorSubdomain::removeNode(int tag)
{
  TaggedObject *object = internalNodes->removeComponent(tag);
  if (object == 0) {
      object = externalNodes->removeComponent(tag);      
      Node *result = (Node *)object;
      this->domainChange();                
      return result;	  
  }
  else {
      this->domainChange();          
      Node *result = (Node *)object;
      return result;	  
  }
}

bool 
ActorSubdomain::addExternalNode(Node *thePtr, bool check)
{

    if (check == true) {
	// check to see it has not alredy been added
	
	int nodTag = thePtr->getTag();
	TaggedObject *other = externalNodes->getComponentPtr(nodTag);
	if (other != 0)
	    return false;
	
	other = internalNodes->getComponentPtr(nodTag);
	if (other != 0)
	    return false;
    }

    bool result = externalNodes->addComponent(thePtr);
    if (result == true) {
	thePtr->setDomain(this);
	this->domainChange();    
    }
    
    return result;
}

const Vector &
ActorSubdomain::getLastExternalSysResponse(void)
{
    int numDOF = this->getNumDOF();
    numDOF = this->getNumDOF();

    if (lastResponse == 0)
	lastResponse = new Vector(numDOF);
    else if (lastResponse->Size() != numDOF) {
	delete lastResponse;
	lastResponse = new Vector(numDOF);
    }
    
    if (mapBuilt == false)
      this->buildMap();

    ID &theMap = *map;
    Vector &localResponse = *lastResponse;
    int numberDOF = this->getNumDOF();
    for (int i=0; i<numberDOF; i++)
      (*mappedVect)(theMap(i)) = localResponse(i);

    return *mappedVect;

}










