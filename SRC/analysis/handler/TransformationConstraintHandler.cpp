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
                                                                        
// $Revision: 1.1.1.1 $
// $Date: 2000-09-15 08:23:17 $
// $Source: /usr/local/cvs/OpenSees/SRC/analysis/handler/TransformationConstraintHandler.cpp,v $
                                                                        
                                                                        
// File: ~/analysis/handler/TransformationConstraintHandler.C
// 
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
#include <DOF_Group.h>
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


TransformationConstraintHandler::TransformationConstraintHandler()
:ConstraintHandler(HANDLER_TAG_TransformationConstraintHandler),
 theFEs(0), theDOFs(0),numFE(0),numDOF(0),numConstrainedNodes(0)
{

}

TransformationConstraintHandler::~TransformationConstraintHandler()
{
    // delete the FE_Element and DOF_Group objects
    for (int i=0; i<numFE; i++)
	if (theFEs[i] != 0)
	  delete theFEs[i];

    for (int j=0; j<numDOF; j++)
	if (theDOFs[j] != 0)
	  delete theDOFs[j];
    
    // delete the arrays
    if (theFEs != 0) delete [] theFEs;
    if (theDOFs != 0) delete [] theDOFs;
}

int
TransformationConstraintHandler::handle(const ID *nodesLast)
{
    // first check links exist to a Domain and an AnalysisModel object
    Domain *theDomain = this->getDomainPtr();
    AnalysisModel *theModel = this->getAnalysisModelPtr();
    Integrator *theIntegrator = this->getIntegratorPtr();    
    
    if ((theDomain == 0) || (theModel == 0) || (theIntegrator == 0)) {
	cerr << "WARNING TransformationConstraintHandler::handle() - ";
	cerr << " setLinks() has not been called\n";
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
    
    
    // int numConstraints = numMPConstraints+numSPConstraints;
    numFE = theDomain->getNumElements();
    numDOF = theDomain->getNumNodes();

    // create an array for the FE_elements and zero it
    if ((numFE <= 0) || ((theFEs  = new FE_Element *[numFE]) == 0)) {
	cerr << "WARNING TransformationConstraintHandler::handle() - ";
        cerr << "ran out of memory for FE_elements"; 
	cerr << " array of size " << numFE << endl;
	return -2;
    }
    int i;
    for (i=0; i<numFE; i++) theFEs[i] = 0;

    // create an array for the DOF_Groups and zero it
    if ((numDOF <= 0) || ((theDOFs = new DOF_Group *[numDOF]) == 0)) {
	cerr << "WARNING TransformationConstraintHandler::handle() - ";
        cerr << "ran out of memory for DOF_Groups";
	cerr << " array of size " << numDOF << endl;
	return -3;    
    }    
    for (i=0; i<numDOF; i++) theDOFs[i] = 0;
    
    // create an ID of constrained node tags in MP_Constraints
    ID *constrainedNodesMP =0;
    MP_Constraint **mps =0;
    if (numMPConstraints != 0) {
	constrainedNodesMP = new ID(numMPConstraints);
	mps = new MP_Constraint *[numMPConstraints];
	if (mps == 0) {
	    cerr << "WARNING TransformationConstraintHandler::handle() - ";
	    cerr << "ran out of memory for MP_Constraints"; 
	    cerr << " array of size " << numMPConstraints << endl;
	    return -3;	    
	}
	MP_ConstraintIter &theMPs = theDomain->getMPs();
	MP_Constraint *theMP; 
	int index = 0;
	while ((theMP = theMPs()) != 0) {
	    (*constrainedNodesMP)(index) = theMP->getNodeConstrained();
	    mps[index] = theMP;
	    index++;
	}	
    }

    // create an ID of constrained node tags in SP_Constraints
    ID *constrainedNodesSP =0;
    SP_Constraint **sps =0;
    if (numSPConstraints != 0) {
	constrainedNodesSP = new ID(numSPConstraints);
	sps = new SP_Constraint *[numSPConstraints];
	if (sps == 0) {
	    cerr << "WARNING TransformationConstraintHandler::handle() - ";
	    cerr << "ran out of memory for SP_Constraints"; 
	    cerr << " array of size " << numSPConstraints << endl;
	    return -3;	    
	}
	SP_ConstraintIter &theSPs = theDomain->getDomainAndLoadPatternSPs();
	SP_Constraint *theSP; 
	int index = 0;
	while ((theSP = theSPs()) != 0) {
	    (*constrainedNodesSP)(index) = theSP->getNodeTag();
	    sps[index] = theSP;
	    index++;
	}	
    }
    
    //create a DOF_Group for each Node and add it to the AnalysisModel.
    //    :must of course set the initial IDs
    NodeIter &theNod = theDomain->getNodes();
    Node *nodPtr;

    int numDofGrp = 0;
    int count3 = 0;
    int countDOF =0;
    
    numConstrainedNodes = 0;
    int numUnconstrainedNodes = 0;
    
    while ((nodPtr = theNod()) != 0) {

        DOF_Group *dofPtr = 0;

	int nodeTag = nodPtr->getTag();
	int numNodalDOF = nodPtr->getNumberDOF();
	int loc = -1;
	int createdDOF = 0;

	if (numMPConstraints != 0) {	
	    loc = constrainedNodesMP->getLocation(nodeTag);
	    if (loc >= 0) {

		TransformationDOF_Group *tDofPtr = 
		    new TransformationDOF_Group(numDofGrp++, nodPtr, mps[loc]); 
								 
		createdDOF = 1;
		dofPtr = tDofPtr;
		
		// add any SPs
		if (numSPConstraints != 0) {
		    loc = constrainedNodesSP->getLocation(nodeTag);
		    if (loc >= 0) {
			tDofPtr->addSP_Constraint(*(sps[loc]));
			for (int i = loc+1; i<numSPConstraints; i++) {
			    if ((*constrainedNodesSP)(i) == nodeTag)
				tDofPtr->addSP_Constraint(*(sps[i]));
			}
		    }
		    // add the DOF to the array	    
		    theDOFs[numDOF-numConstrainedNodes-1] = dofPtr;	    	    
		    numConstrainedNodes++;
		}
	    }
	}
	
	if (createdDOF == 0 && numSPConstraints != 0) {
	    loc = constrainedNodesSP->getLocation(nodeTag);
	    if (loc >= 0) {
		TransformationDOF_Group *tDofPtr = 
		    new TransformationDOF_Group(numDofGrp++, nodPtr);
		
		int numSPs = 1;
		createdDOF = 1;
		dofPtr = tDofPtr;
		tDofPtr->addSP_Constraint(*(sps[loc]));

		// check for more SP_constraints acting on node and add them
		for (int i = loc+1; i<numSPConstraints; i++) {
		    if ((*constrainedNodesSP)(i) == nodeTag) {
			tDofPtr->addSP_Constraint(*(sps[i]));
			numSPs++;
		    }
		}
		// add the DOF to the array
		theDOFs[numDOF-numConstrainedNodes-1] = dofPtr;	    	    
		numConstrainedNodes++;	    
		countDOF+= numNodalDOF - numSPs;		
	    }
	}

	// create an ordinary DOF_Group object if no dof constrained
	if (createdDOF == 0) {
	    if ((dofPtr = new DOF_Group(numDofGrp++, nodPtr)) == 0) {
		cerr << "WARNING TransformationConstraintHandler::handle() ";
		cerr << "- ran out of memory";
		cerr << " creating DOF_Group " << i << endl;	
		return -4;    		
	    }
	
	    countDOF+= numNodalDOF;
	    theDOFs[numUnconstrainedNodes++] = dofPtr;	    
	}
	
	if (dofPtr == 0) 
	  cerr << "TransformationConstraintHandler::handle() - error in logic\n";
	    
	nodPtr->setDOF_GroupPtr(dofPtr);
	theModel->addDOF_Group(dofPtr);
    }

    // create the FE_Elements for the Elements and add to the AnalysisModel
    ElementIter &theEle = theDomain->getElements();
    Element *elePtr;

    int numFeEle = 0;
    FE_Element *fePtr;
    while ((elePtr = theEle()) != 0) {

	// we must check to see if the element has a constrained node
	const ID &nodes = elePtr->getExternalNodes();
	int nodesSize = nodes.Size();
	int isConstrainedNode = 0;
	for (int i=0; i<nodesSize; i++) {
	    int nodeTag = nodes(i);
	    if (numMPConstraints != 0) {
		int loc = constrainedNodesMP->getLocation(nodeTag);
		if (loc >= 0) {
		    isConstrainedNode = 1;
		    i = nodesSize;
		}
	    } 
	    if (numSPConstraints != 0 && isConstrainedNode == 0) {
		int loc = constrainedNodesSP->getLocation(nodeTag);
		if (loc >= 0) {
		    isConstrainedNode = 1;		    
		    i = nodesSize;
		}
	    }
	}	    
	    
	// if constrained Transformation otherwise normal FE_Element
	if (isConstrainedNode == 0) {
	    if ((fePtr = new FE_Element(elePtr)) == 0) {
		cerr << "WARNING TransformationConstraintHandler::handle()";
		cerr << " - ran out of memory";
		cerr << " creating FE_Element " << elePtr->getTag() << endl; 
		return -5;
	    }		
	} else {
	    if ((fePtr = new TransformationFE(elePtr, *this)) == 0) {		
		cerr << "WARNING TransformationConstraintHandler::handle()";
		cerr << " - ran out of memory";
		cerr << " creating TransformationFE " << elePtr->getTag() << endl; 
		return -6;		    
	    }
	}
	
	theFEs[numFeEle++] = fePtr;
	theModel->addFE_Element(fePtr);

	if (elePtr->isSubdomain() == true) {
	    Subdomain *theSub = (Subdomain *)elePtr;
	    theSub->setFE_ElementPtr(fePtr);
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
			cerr << "WARNING TransformationConstraintHandler::handle() ";
			cerr << " - boundary sp constraint in subdomain";
			cerr << " this should not be - results suspect \n";
		    }
		}
	    }
	}

    return count3;
}



void 
TransformationConstraintHandler::clearAll(void)
{
    // delete the FE_Element and DOF_Group objects
    for (int i=0; i<numFE; i++)
	if (theFEs[i] != 0)
	    delete theFEs[i];

    for (int j=0; j<numDOF; j++)
	if (theDOFs[j] != 0)
	    delete theDOFs[j];
    
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
TransformationConstraintHandler::enforceSPs(void)
{
    for (int i=1; i<=numConstrainedNodes; i++) {
	// upward cast - safe as i put it in this location
	TransformationDOF_Group *theDof  =
	    (TransformationDOF_Group *)theDOFs[numDOF-i];
	theDof->enforceSPs();
    }
    return 0;
}

int 
TransformationConstraintHandler::doneDOFids(void)
{
    for (int i=1; i<=numConstrainedNodes; i++) {
	// upward cast - safe as i put it in this location
	TransformationDOF_Group *theDof  =
	    (TransformationDOF_Group *)theDOFs[numDOF-i];
	theDof->doneID();
    }
    return 0;
}
