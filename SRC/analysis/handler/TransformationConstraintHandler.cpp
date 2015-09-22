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
                                                                        
// $Revision: 1.15 $
// $Date: 2008-11-19 23:42:25 $
// $Source: /usr/local/cvs/OpenSees/SRC/analysis/handler/TransformationConstraintHandler.cpp,v $
                                                                        
                                                                        
// Written: fmk 
// Created: May 1998
// Revision: A
//
// What: "@(#) TransformationConstraintHandler.C, revA"

#include <TransformationConstraintHandler.h>
#include <stdlib.h>

#include <AnalysisModel.h>
#include <Domain.h>
#include <FE_Element.h>
#include <FE_EleIter.h>
#include <DOF_Group.h>
#include <DOF_GrpIter.h>
#include <Node.h>
#include <Element.h>
#include <NodeIter.h>
#include <ElementIter.h>
#include <SP_ConstraintIter.h>
#include <SP_Constraint.h>
#include <MP_ConstraintIter.h>
#include <MP_Constraint.h>
#include <Integrator.h>
#include <ID.h>
#include <Subdomain.h>
#include <Channel.h>
#include <FEM_ObjectBroker.h>
#include <TransformationDOF_Group.h>
#include <TransformationFE.h>

void* OPS_TransformationConstraintHandler()
{
    return new TransformationConstraintHandler;
}

TransformationConstraintHandler::TransformationConstraintHandler()
:ConstraintHandler(HANDLER_TAG_TransformationConstraintHandler),
 theFEs(0), theDOFs(0),numFE(0),numDOF(0),numConstrainedNodes(0)
{

}

TransformationConstraintHandler::~TransformationConstraintHandler()
{
  if (theDOFs != 0) 
    delete [] theDOFs;

  if (theFEs != 0)
    delete [] theFEs;
}

int
TransformationConstraintHandler::handle(const ID *nodesLast)
{
    // first check links exist to a Domain and an AnalysisModel object
    Domain *theDomain = this->getDomainPtr();
    AnalysisModel *theModel = this->getAnalysisModelPtr();
    Integrator *theIntegrator = this->getIntegratorPtr();    
    
    if ((theDomain == 0) || (theModel == 0) || (theIntegrator == 0)) {
	opserr << "WARNING TransformationConstraintHandler::handle() - ";
	opserr << " setLinks() has not been called\n";
	return -1;
    }
    
    // get number ofelements and nodes in the domain 
    // and init the theFEs and theDOFs arrays
    int numMPConstraints = theDomain->getNumMPs();

    //    int numSPConstraints = theDomain->getNumSPs();    
    int numSPConstraints = 0;
    SP_ConstraintIter &theSP1s = theDomain->getDomainAndLoadPatternSPs();
    SP_Constraint *theSP1; 
    while ((theSP1 = theSP1s()) != 0) 
	numSPConstraints++;
    
    numDOF = 0;
    ID transformedNode(0, 64);

    int i;
    
    // create an ID of constrained node tags in MP_Constraints
    ID constrainedNodesMP(0, numMPConstraints);
    MP_Constraint **mps =0;
    if (numMPConstraints != 0) {
	mps = new MP_Constraint *[numMPConstraints];
	if (mps == 0) {
	    opserr << "WARNING TransformationConstraintHandler::handle() - ";
	    opserr << "ran out of memory for MP_Constraints"; 
	    opserr << " array of size " << numMPConstraints << endln;
	    return -3;	    
	}
	MP_ConstraintIter &theMPs = theDomain->getMPs();
	MP_Constraint *theMP; 
	int index = 0;
	while ((theMP = theMPs()) != 0) {
	  int nodeConstrained = theMP->getNodeConstrained();
	  if (transformedNode.getLocation(nodeConstrained) < 0)
	    transformedNode[numDOF++] = nodeConstrained;
	  constrainedNodesMP[index] = nodeConstrained;
	  mps[index] = theMP;
	  index++;
	}	
    }

    // create an ID of constrained node tags in SP_Constraints
    ID constrainedNodesSP(0, numSPConstraints);;
    SP_Constraint **sps =0;
    if (numSPConstraints != 0) {
	sps = new SP_Constraint *[numSPConstraints];
	if (sps == 0) {
	    opserr << "WARNING TransformationConstraintHandler::handle() - ";
	    opserr << "ran out of memory for SP_Constraints"; 
	    opserr << " array of size " << numSPConstraints << endln;
	    if (mps != 0) delete [] mps;
	    if (sps != 0) delete [] sps;
	    return -3;	    
	}
	SP_ConstraintIter &theSPs = theDomain->getDomainAndLoadPatternSPs();
	SP_Constraint *theSP; 
	int index = 0;
	while ((theSP = theSPs()) != 0) {
	  int constrainedNode = theSP->getNodeTag();
	  if (transformedNode.getLocation(constrainedNode) < 0)
	    transformedNode[numDOF++] = constrainedNode;	    
	  constrainedNodesSP[index] = constrainedNode;
	  sps[index] = theSP;
	  index++;
	}	
    }

    // create an array for the DOF_Groups and zero it
    if ((numDOF != 0) && ((theDOFs = new DOF_Group *[numDOF]) == 0)) {
	opserr << "WARNING TransformationConstraintHandler::handle() - ";
        opserr << "ran out of memory for DOF_Groups";
	opserr << " array of size " << numDOF << endln;
	return -3;    
    }    
    for (i=0; i<numDOF; i++) theDOFs[i] = 0;

    //create a DOF_Group for each Node and add it to the AnalysisModel.
    //    :must of course set the initial IDs
    NodeIter &theNod = theDomain->getNodes();
    Node *nodPtr;

    int numDofGrp = 0;
    int count3 = 0;
    int countDOF =0;
    
    numConstrainedNodes = 0;
    numDOF = 0;
    while ((nodPtr = theNod()) != 0) {

        DOF_Group *dofPtr = 0;

	int nodeTag = nodPtr->getTag();
	int numNodalDOF = nodPtr->getNumberDOF();
	int loc = -1;
	int createdDOF = 0;

	loc = constrainedNodesMP.getLocation(nodeTag);
	if (loc >= 0) {

	  TransformationDOF_Group *tDofPtr = 
	    new TransformationDOF_Group(numDofGrp++, nodPtr, mps[loc], this); 

	  createdDOF = 1;
	  dofPtr = tDofPtr;
	  
	  // add any SPs
	  if (numSPConstraints != 0) {
	    loc = constrainedNodesSP.getLocation(nodeTag);
	    if (loc >= 0) {
	      tDofPtr->addSP_Constraint(*(sps[loc]));
	      for (int i = loc+1; i<numSPConstraints; i++) {
		if (constrainedNodesSP(i) == nodeTag)
		  tDofPtr->addSP_Constraint(*(sps[i]));
	      }
	    }
	    // add the DOF to the array	    
	    theDOFs[numDOF++] = dofPtr;	    	    
	    numConstrainedNodes++;
	  }
	}
	
	if (createdDOF == 0) {
	  loc = constrainedNodesSP.getLocation(nodeTag);
	  if (loc >= 0) {
	    TransformationDOF_Group *tDofPtr = 
	      new TransformationDOF_Group(numDofGrp++, nodPtr, this);

	    int numSPs = 1;
	    createdDOF = 1;
	    dofPtr = tDofPtr;
	    tDofPtr->addSP_Constraint(*(sps[loc]));
	
	    // check for more SP_constraints acting on node and add them
	    for (int i = loc+1; i<numSPConstraints; i++) {
	      if (constrainedNodesSP(i) == nodeTag) {
		tDofPtr->addSP_Constraint(*(sps[i]));
		numSPs++;
	      }
	    }
	    // add the DOF to the array
	    theDOFs[numDOF++] = dofPtr;	    	    
	    numConstrainedNodes++;	    
	    countDOF+= numNodalDOF - numSPs;		
	  }
	}

	// create an ordinary DOF_Group object if no dof constrained
	if (createdDOF == 0) {
	    if ((dofPtr = new DOF_Group(numDofGrp++, nodPtr)) == 0) {
		opserr << "WARNING TransformationConstraintHandler::handle() ";
		opserr << "- ran out of memory";
		opserr << " creating DOF_Group " << i << endln;	
		if (mps != 0) delete [] mps;
		if (sps != 0) delete [] sps;
		return -4;    		
	    }
	
	    countDOF+= numNodalDOF;
	}
	
	if (dofPtr == 0) 
	  opserr << "TransformationConstraintHandler::handle() - error in logic\n";
	    
	nodPtr->setDOF_GroupPtr(dofPtr);
	theModel->addDOF_Group(dofPtr);
    }

    // create the FE_Elements for the Elements and add to the AnalysisModel
    ElementIter &theEle = theDomain->getElements();
    Element *elePtr;
    FE_Element *fePtr;

    numFE = 0;
    ID transformedEle(0, 64);

    while ((elePtr = theEle()) != 0) {
      int flag = 0;
      if (elePtr->isSubdomain() == true) {
	Subdomain *theSub = (Subdomain *)elePtr;
	if (theSub->doesIndependentAnalysis() == true) 
	  flag = 1;
      }

      if (flag == 0) {
      
	const ID &nodes = elePtr->getExternalNodes();
	int nodesSize = nodes.Size();
	int isConstrainedNode = 0;
	for (int i=0; i<nodesSize; i++) {
	  int nodeTag = nodes(i);
	  if (numMPConstraints != 0) {
	    int loc = constrainedNodesMP.getLocation(nodeTag);
	    if (loc >= 0) {
	      isConstrainedNode = 1;
	      i = nodesSize;
	    }
	  } 
	  if (numSPConstraints != 0 && isConstrainedNode == 0) {
	    int loc = constrainedNodesSP.getLocation(nodeTag);
	    if (loc >= 0) {
	      isConstrainedNode = 1;		    
	      i = nodesSize;
	    }
	  }
	}
	
	if (isConstrainedNode == 1) {
	  transformedEle[numFE++] = elePtr->getTag();
	}
      }
    }
    
    // create an array for the FE_elements and zero it
    if ((numFE != 0) && ((theFEs  = new FE_Element *[numFE]) == 0)) {
      opserr << "WARNING TransformationConstraintHandler::handle() - ";
      opserr << "ran out of memory for FE_elements"; 
      opserr << " array of size " << numFE << endln;
      return -2;
    }
    
    for (i=0; i<numFE; i++) theFEs[i] = 0;

    ElementIter &theEle1 = theDomain->getElements();
    
    // int numConstraints = numMPConstraints+numSPConstraints;
    int numFeEle = 0;
    int numFE = 0;

    while ((elePtr = theEle1()) != 0) {
      int tag = elePtr->getTag();
      if (elePtr->isSubdomain() == true) {
	Subdomain *theSub = (Subdomain *)elePtr;
	if (theSub->doesIndependentAnalysis() == false) {
	  
	  if (transformedEle.getLocation(tag) < 0) {
	    if ((fePtr = new FE_Element(numFeEle, elePtr)) == 0) {
	      opserr << "WARNING TransformationConstraintHandler::handle()";
	      opserr << " - ran out of memory";
	      opserr << " creating FE_Element " << elePtr->getTag() << endln; 
	      if (mps != 0) delete [] mps;
	      if (sps != 0) delete [] sps;
	      return -5;
	    }	
	  } else {
	    if ((fePtr = new TransformationFE(numFeEle, elePtr)) == 0) {		
	      opserr << "WARNING TransformationConstraintHandler::handle()";
	      opserr << " - ran out of memory";
	      opserr << " creating TransformationFE " << elePtr->getTag() << endln; 
	      if (mps != 0) delete [] mps;
	      if (sps != 0) delete [] sps;
	      return -6;		    
	    }
	    theFEs[numFE++] = fePtr;
	  }

	  numFeEle++;
	  theModel->addFE_Element(fePtr);
	  theSub->setFE_ElementPtr(fePtr);
	}
      } else {
	if (transformedEle.getLocation(tag) < 0) {
	  if ((fePtr = new FE_Element(numFeEle, elePtr)) == 0) {
	    opserr << "WARNING TransformationConstraintHandler::handle()";
	    opserr << " - ran out of memory";
	    opserr << " creating FE_Element " << elePtr->getTag() << endln; 
	    if (mps != 0) delete [] mps;
	    if (sps != 0) delete [] sps;
	    return -5;
	  }	
	} else {
	  if ((fePtr = new TransformationFE(numFeEle, elePtr)) == 0) {		
	    opserr << "WARNING TransformationConstraintHandler::handle()";
	    opserr << " - ran out of memory";
	    opserr << " creating TransformationFE " << elePtr->getTag() << endln; 
	    if (mps != 0) delete [] mps;
	    if (sps != 0) delete [] sps;
	    return -6;		    
	  }
	  theFEs[numFE++] = fePtr;
	}
	
	numFeEle++;
	theModel->addFE_Element(fePtr);
      }
    }

    theModel->setNumEqn(countDOF);
    
    // set the number of eqn in the model
    // now see if we have to set any of the dof's to -3
    //    int numLast = 0;
    if (nodesLast != 0) 
	for (i=0; i<nodesLast->Size(); i++) {
	    int nodeID = (*nodesLast)(i);
	    Node *nodPtr = theDomain->getNode(nodeID);
	    if (nodPtr != 0) {
		DOF_Group *dofPtr = nodPtr->getDOF_GroupPtr();
		
		const ID &id = dofPtr->getID();
		// set all the dof values to -3
		for (int j=0; j < id.Size(); j++) {
		    if (id(j) == -2) {
			dofPtr->setID(j,-3);
			count3++;
		    } else {
			opserr << "WARNING TransformationConstraintHandler::handle() ";
			opserr << " - boundary sp constraint in subdomain";
			opserr << " this should not be - results suspect \n";
			if (mps != 0) delete [] mps;
			if (sps != 0) delete [] sps;
		    }
		}
	    }
	}

    if (mps != 0) delete [] mps;
    if (sps != 0) delete [] sps;

    return count3;
}



void 
TransformationConstraintHandler::clearAll(void)
{
    // delete the arrays
    if (theFEs != 0) delete [] theFEs;
    if (theDOFs != 0) delete [] theDOFs;

    // reset the numbers
    numDOF = 0;
    numFE =  0;
    theFEs = 0;
    theDOFs = 0;

    // for the nodes reset the DOF_Group pointers to 0
    Domain *theDomain = this->getDomainPtr();
    if (theDomain == 0)
	return;

    NodeIter &theNod = theDomain->getNodes();
    Node *nodPtr;
    while ((nodPtr = theNod()) != 0)
	nodPtr->setDOF_GroupPtr(0);
}    

int
TransformationConstraintHandler::sendSelf(int cTag, Channel &theChannel)
{
  return 0;
}

int
TransformationConstraintHandler::recvSelf(int cTag, 
				   Channel &theChannel, 
				   FEM_ObjectBroker &theBroker)  
{
  return 0;
}


int 
TransformationConstraintHandler::applyLoad(void)
{
    return this->enforceSPs();
}

int 
TransformationConstraintHandler::enforceSPs(void)
{
    for (int i=1; i<=numConstrainedNodes; i++) {
	// upward cast - safe as i put it in this location
	TransformationDOF_Group *theDof  =
	    (TransformationDOF_Group *)theDOFs[numDOF-i];
	theDof->enforceSPs(1);
    }
    for (int k=1; k<=numConstrainedNodes; k++) {
	// upward cast - safe as i put it in this location
	TransformationDOF_Group *theDof  =
	    (TransformationDOF_Group *)theDOFs[numDOF-k];
	theDof->enforceSPs(0);
    }

    for (int j=0; j<numFE; j++) {
      FE_Element *theEle = theFEs[j];
      theEle->updateElement();
    }

    return 0;
}

int 
TransformationConstraintHandler::doneNumberingDOF(void)
{
    // iterate through the DOF_Groups telling them that their ID has now been set
    AnalysisModel *theModel1=this->getAnalysisModelPtr();
    DOF_GrpIter &theDOFS = theModel1->getDOFs();
    DOF_Group *dofPtr;
    while ((dofPtr = theDOFS()) != 0) {
       dofPtr->doneID();
    }


    // iterate through the FE_Element getting them to set their IDs
    AnalysisModel *theModel=this->getAnalysisModelPtr();
    FE_EleIter &theEle = theModel->getFEs();
    FE_Element *elePtr;
    while ((elePtr = theEle()) != 0) {
      elePtr->setID();
    }

    return 0;
}
