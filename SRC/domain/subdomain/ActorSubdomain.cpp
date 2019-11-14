/* ****************************************************************** **
**    OpenSeess - Open System for Earthquake Engineering Simulation   **
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
                                                                        
// $Revision: 1.19 $
// $Date: 2010-09-16 00:07:11 $
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
#include <ConvergenceTest.h>

#include <EquiSolnAlgo.h>
#include <IncrementalIntegrator.h>
#include <LinearSOE.h>
#include <EigenSOE.h>
#include <Recorder.h>
#include <Parameter.h>
#include <Message.h>

#include <ArrayOfTaggedObjects.h>
#include <ShadowActorSubdomain.h>

// 2 procedurs defined in SP_Constraint.cpp
int SP_Constraint_GetNextTag();
int SP_Constraint_SetNextTag(int);

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
}


int
ActorSubdomain::run(void)
{
    static Vector theVect(4);
	static Vector theVect1(1);
    bool exitYet = false;
    int res = 0;

    while (exitYet == false) {
      int action;
      res = this->recvID(msgData);
      if (res != 0) {
	opserr << "ActorSubdomain::run - error receiving msgData\n";
	exitYet = true;
        action = ShadowActorSubdomain_DIE;
      } else {
	action = msgData(0);
      }

      bool change;
      int theType, theOtherType, tag, dbTag, loadPatternTag;
      int startTag, endTag, axisDirn, numSP, i, numMode, dof;
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
      Matrix *theM;
      Vector *theV;
      ID     *theI, *theNodeTags, *theEleTags;
      PartitionedModelBuilder *theBuilder;
      IncrementalIntegrator *theIntegrator;
      EquiSolnAlgo *theAlgorithm;
      LinearSOE *theSOE;
      EigenSOE *theEigenSOE;
      ConvergenceTest *theTest;
      Recorder *theRecorder;
      bool res, generalized, findSmallest;
      double doubleRes;
      int intRes;
      NodeResponseType nodeResponseType;
      Parameter *theParameter;
      int argc;
      char **argv;
      char *allResponseArgs;
      char *currentLoc;
      int argLength, msgLength;
      Message theMessage;

      const ID *theID;
      
      //     opserr << "ActorSubdomain ACTION: " << action << endln;

      switch (action) {
	  case ShadowActorSubdomain_setTag:
	    tag = msgData(1); // subdomain tag
	    this->setTag(tag);
	    this->Actor::setCommitTag(tag);
	    break;

	  case ShadowActorSubdomain_analysisStep:
	    this->recvVector(theVect);
	    this->analysisStep(theVect(0));
	    break;

	  case ShadowActorSubdomain_eigenAnalysis:
	    numMode = msgData(1);
	    if (msgData(2) == 0)
	      generalized = true;
	    else
	      generalized = false;
	    if (msgData(3) == 0)
	      findSmallest = true;
	    else
	      findSmallest = false;
		
	    this->eigenAnalysis(numMode, generalized, findSmallest);
	    break;
	    /*
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
	    */
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

	    
	  case ShadowActorSubdomain_hasNode:
	    theType = msgData(1);
	    res = this->hasNode(theType);
	    if (res == true)
	      msgData(0) = 0;
	    else
	      msgData(0) = -1;
	    this->sendID(msgData);

	    break;

	  case ShadowActorSubdomain_hasElement:
	    theType = msgData(1);
	    res = this->hasElement(theType);
	    if (res == true)
	      msgData(0) = 0;
	    else
	      msgData(0) = -1;
	    this->sendID(msgData);
	   
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
	    //	    opserr << "ActorSubdomain::add node: " << *theNod;
	    break;


	  case ShadowActorSubdomain_addExternalNode:
	    theType = msgData(1);
	    dbTag = msgData(2);
	    theNod = theBroker->getNewNode(theType);

	    if (theNod != 0) {
		theNod->setDbTag(dbTag);
		this->recvObject(*theNod);
		bool result = this->Subdomain::addExternalNode(theNod);
		delete theNod;
		/*
		Node *dummy = new Node(*theNod);
		delete theNod;
		opserr << *dummy;
		opserr << dummy->getMass();
		*/

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


	  case ShadowActorSubdomain_addSP_ConstraintAXIS:

	    axisDirn = msgData(1);
	    theI = new ID(msgData(2));
	    theV = new Vector(2);
	    SP_Constraint_SetNextTag(msgData(3));

	    endTag = 0;
		
	    this->recvID(*theI);
	    this->recvVector(*theV);

	    msgData(0) = 0;				 
	    numSP = this->addSP_Constraint(axisDirn, (*theV)(0), *theI, (*theV)(1));
	    endTag = SP_Constraint_GetNextTag();

	    msgData(1) = numSP;
	    msgData(2) = endTag;

	    this->domainChange();
	    this->sendID(msgData);
		
	    delete theV;
	    delete theI;
		
	    /* DON'T BOTHER SENDING
	    if (numSP > 0) {
	      theI = new ID(numSP);
	      for (i = 0; i<numSP; i++) {
		theSP = this->getSP_Constraint(i+startTag);
		(*theI)(i) = theSP->getClassTag();
	      }
	      this->sendID(*theI);	
		  opserr << "Actor: sent: " << *theI;
	      for (i = 0; i<numSP; i++) {
		theSP = this->getSP_Constraint(i+startTag);
		if (theSP != 0)
		this->sendObject(*theSP);	
		else
			opserr << "ActorSubdomain::addSP_AXIS :: PROBLEMS\n";
	      }
	      delete theI;
	    }
opserr << "ActorSubdomain::addSP_AXIS :: DONE\n";
        */

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

	  case ShadowActorSubdomain_removeSP_ConstraintNoTag:
	    tag = msgData(1);
	    dof = msgData(2);
	    loadPatternTag = msgData(3);
	    msgData(0) = this->removeSP_Constraint(tag, dof, loadPatternTag);
	    this->sendID(msgData);

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


	  case ShadowActorSubdomain_getNode:
	    tag = msgData(1);

	    theNod = this->getNode(tag);

	    if (theNod != 0) 
		msgData(0) = theNod->getClassTag();
	    else
		msgData(0) = -1;

	    this->sendID(msgData);

	    if (theNod != 0) {
		this->sendObject(*theNod);
	    }

	    msgData(0) = 0;

	    break;	    	    	    	    


	  case ShadowActorSubdomain_Print:
	    this->Print(opserr, msgData(3));
	    this->sendID(msgData);

	    break;	    	    	    	    

	  case ShadowActorSubdomain_PrintNodeAndEle:
	    
	    theNodeTags = 0;
	    theEleTags = 0;

	    if (msgData(1) != 0) {
	      theNodeTags = new ID(msgData(1));
	      this->recvID(*theNodeTags);
	    }
	    if (msgData(2) != 0) {
	      theEleTags = new ID(msgData(2));
	      this->recvID(*theEleTags);
	    }
	      
	    this->Print(opserr, theNodeTags, theEleTags, msgData(3));
	    
	    if (theNodeTags != 0)
	      delete theNodeTags;
	    if (theEleTags != 0)
	      delete theEleTags;

	    this->sendID(msgData);

	    break;	    	    	    	    

	  case ShadowActorSubdomain_applyLoad:
	    this->recvVector(theVect);	    
	    this->applyLoad(theVect(0));
	    break;

	  case ShadowActorSubdomain_setCommittedTime:
	    this->recvVector(theVect);	    
	    this->setCurrentTime(theVect(0));
	    this->setCommittedTime(theVect(0));
	    break;	    
	    
	  case ShadowActorSubdomain_setLoadConstant:
	    this->setLoadConstant();
	    break;	    

	  case ShadowActorSubdomain_update:
	    this->update();
	    break;

	  case ShadowActorSubdomain_updateTimeDt:
	    this->updateTimeDt();
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
            break;

	  case ShadowActorSubdomain_record:
	    this->record();
	    break;
	    
	  case ShadowActorSubdomain_commit:
	    this->commit();
	    break;
	    
	  case ShadowActorSubdomain_revertToLastCommit:
	    this->revertToLastCommit();
	    break;	    
	    
	  case ShadowActorSubdomain_revertToStart:
	    this->revertToStart();
	    this->sendID(msgData);

	    break;	    	    

	  case ShadowActorSubdomain_addRecorder:
	    theType = msgData(1);
	    theRecorder = theBroker->getPtrNewRecorder(theType);
	    if (theRecorder != 0) {
	      this->recvObject(*theRecorder);	      
	      this->addRecorder(*theRecorder);
	    }
	    break;	    	    

	  case ShadowActorSubdomain_removeRecorders:
	    this->removeRecorders();
	    this->barrierCheck(1);
	    break;	    	    

	  case ShadowActorSubdomain_removeRecorder:
	    theType = msgData(1);
	    this->removeRecorder(theType);
	    break;	    	    
	    

	case ShadowActorSubdomain_wipeAnalysis:
	  this->wipeAnalysis();	    
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
	    this->recvObject(*theIntegrator);
	    this->setAnalysisIntegrator(*theIntegrator);
	    msgData(0) = 0;
	  } else
	    msgData(0) = -1;
	  this->sendID(msgData);
	  break;

	case ShadowActorSubdomain_setAnalysisLinearSOE:
	  theType = msgData(1);
	  theSOE = theBroker->getNewLinearSOE(theType);

	  if (theSOE != 0) {
	    this->recvObject(*theSOE);
	    this->setAnalysisLinearSOE(*theSOE);
	    msgData(0) = 0;
	  } else
	    msgData(0) = -1;
	    
	  break;

	case ShadowActorSubdomain_setAnalysisEigenSOE:
	  theType = msgData(1);
	  theEigenSOE = theBroker->getNewEigenSOE(theType);

	  if (theEigenSOE != 0) {
	    this->recvObject(*theEigenSOE);
	    this->setAnalysisEigenSOE(*theEigenSOE);
	    msgData(0) = 0;
	  } else
	    msgData(0) = -1;
	    
	  break;

	case ShadowActorSubdomain_setAnalysisConvergenceTest:
	  theType = msgData(1);
	  theTest = theBroker->getNewConvergenceTest(theType);
	  
	  if (theTest != 0) {
	    this->recvObject(*theTest);
	    this->setAnalysisConvergenceTest(*theTest);
	    msgData(0) = 0;
	  } else
	    msgData(0) = -1;

	  break;
	    
	  case ShadowActorSubdomain_domainChange:
	    this->domainChange();

	    tag = this->getNumDOF();
	    if (tag != 0) {
	      if (lastResponse == 0)
		lastResponse = new Vector(tag);
	      else if (lastResponse->Size() != tag) {
		delete lastResponse;
		lastResponse = new Vector(tag);
	      }
	    }
	    break;

	  case ShadowActorSubdomain_getDomainChangeFlag:
	    change = this->getDomainChangeFlag();
	    if (change == true)
	      msgData(0) = 0;
	    else
	      msgData(0) = 1;
	    this->sendID(msgData);
	    
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

	  case ShadowActorSubdomain_clearAll:
	    this->clearAll();
	    this->sendID(msgData);
	    break;

	  case ShadowActorSubdomain_getNodeDisp:
	    tag = msgData(1);  // nodeTag
	    dbTag = msgData(2); // dof
	    doubleRes = this->getNodeDisp(tag, dbTag, intRes);
	    msgData(0) = intRes;
	    this->sendID(msgData);
	    if (intRes == 0) {
	      theV = new Vector(1);
	      (*theV)(0) = doubleRes;
	      this->sendVector(*theV);
	      delete theV;
	    }
	    break;

	  case ShadowActorSubdomain_setMass:
	    tag = msgData(1);  // nodeTag
	    dbTag = msgData(2); // noRows
	    theOtherType = msgData(3); // noRows
	    theM = new Matrix(dbTag, theOtherType);
	    this->recvMatrix(*theM);
	    intRes = this->setMass(*theM, tag);
	    
	    delete theM;
	    msgData(0) = intRes;
	    this->sendID(msgData);
	    break;


	  case ShadowActorSubdomain_getNodeResponse:
	    tag = msgData(1);  // nodeTag
	    nodeResponseType = (NodeResponseType)msgData(2); 
	    theVector = this->getNodeResponse(tag, nodeResponseType);

	    if (theVector == 0)
	      msgData(0) = 0;
	    else {
	      msgData(0) = 1;
	      msgData(1) = theVector->Size();
	    }
	    this->sendID(msgData);

	    if (theVector != 0)
	      this->sendVector(*theVector);

	    break;

	  case ShadowActorSubdomain_getElementResponse:
	    tag = msgData(1);  // eleTag
	    argc = msgData(2);
	    msgLength = msgData(3);

	    if (msgLength == 0) {
	      opserr << "ElementRecorder::recvSelf() - 0 sized string for responses\n";
	      return -1;
	    }

	    allResponseArgs = new char[msgLength];
	    if (allResponseArgs == 0) {
	      opserr << "ElementRecorder::recvSelf() - out of memory\n";
	      return -1;
	    }

	    theMessage.setData(allResponseArgs, msgLength);
	    if (this->recvMessage(theMessage) < 0) {
	      opserr << "ElementRecorder::recvSelf() - failed to recv message\n";
	      return -1;
	    }

	    //
	    // now break this single array into many
	    // 
	    
	    argv = new char *[argc];
	    if (argv == 0) {
	      opserr << "ElementRecorder::recvSelf() - out of memory\n";
	      return -1;
	    }
	    
	    currentLoc = allResponseArgs;
	    for (int j=0; j<argc; j++) {
	      argv[j] = currentLoc;	      
	      argLength = strlen(currentLoc)+1;
	      currentLoc += argLength;
	    }

	    theVector = this->getElementResponse(tag, (const char**)argv, argc);

	    delete [] argv;
	    delete [] allResponseArgs;

	    if (theVector == 0) 
	      msgData(0) = 0;
	    else {
	      msgData(0) = 1;
	      msgData(1) = theVector->Size();
	    }
	    this->sendID(msgData);

	    if (theVector != 0)
	      this->sendVector(*theVector);
      
	    break;

	  case ShadowActorSubdomain_calculateNodalReactions:
	    if (msgData(0) == 0)
	      this->calculateNodalReactions(true);
	    else
	      this->calculateNodalReactions(false);
	    break;

         case ShadowActorSubdomain_setRayleighDampingFactors:
	   theV = new Vector(4);
	   this->recvVector(*theV);
	   intRes = this->Subdomain::setRayleighDampingFactors((*theV)(0), (*theV)(1), (*theV)(2), (*theV)
(3));
	   delete theV;
	   break;


         case ShadowActorSubdomain_addParameter:
	    theType = msgData(1);
	    dbTag = msgData(2);

	    theParameter = theBroker->getParameter(theType);

	    if (theParameter != 0) {
		theParameter->setDbTag(dbTag);		
		this->recvObject(*theParameter);
		//bool result = true;
		bool result = this->addParameter(theParameter);
		if (result == true)
		    msgData(0) = 0;
		else {
		  opserr << "Actor::addParameter - FAILED\n";
		  msgData(0) = -1;
		}
	    } else
		msgData(0) = -1;

	   break;

       case ShadowActorSubdomain_removeParameter:
	   theType = msgData(1);
	  
	   this->removeParameter(theType);
	 
	    this->sendID(msgData);
	   break;

         case ShadowActorSubdomain_updateParameterINT:
	   theType = msgData(1);  // tag
	   dbTag = msgData(2);    // value

	   msgData(0) = this->Domain::updateParameter(theType, dbTag);
	   this->sendID(msgData);
	   break;

       case ShadowActorSubdomain_updateParameterDOUBLE:
	   theType = msgData(1);  // tag
	 
	   this->recvVector(theVect1);
	
	   msgData(0) = this->Domain::updateParameter(theType, theVect1(0));
	  
	   this->sendID(msgData);
	   break;


	  case ShadowActorSubdomain_DIE:
	    exitYet = true;
	    break;

	  default:
	    opserr << "ActorSubdomain::invalid action " << action << "received\n";
	    msgData(0) = -1;
	    
      }
      //      opserr << "DONE ACTION: " << action << endln;
    }

    //    this->sendID(msgData);
    return 0;
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

int
ActorSubdomain::update(void)
{
  int res = this->Domain::update();

  res = this->barrierCheck(res);

  return res;
}

int
ActorSubdomain::updateTimeDt(void)
{
  static Vector data(2);

  this->recvVector(data);

  double newTime = data(0);
  double dT = data(1);
  int res = this->Domain::update(newTime, dT);
  return this->barrierCheck(res);
}

int
ActorSubdomain::barrierCheck(int myResult)
{
  static ID data(1);
  data(0) = myResult;
  this->sendID(data);
  this->recvID(data);

  return data(0);
}







