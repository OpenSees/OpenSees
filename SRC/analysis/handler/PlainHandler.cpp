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
// $Source: /usr/local/cvs/OpenSees/SRC/analysis/handler/PlainHandler.cpp,v $
                                                                        
                                                                        
// File: ~/analysis/handler/PlainHandler.C
// 
// Written: fmk 
// Created: 11/96
// Revision: A
//
// Description: This file contains the implementation of PlainHandler.
//
// What: "@(#) PlainHandler.C, revA"

#include <PlainHandler.h>
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

PlainHandler::PlainHandler()
:ConstraintHandler(HANDLER_TAG_PlainHandler),
 theFEs(0), theDOFs(0),numFE(0),numDOF(0)
{

}

PlainHandler::~PlainHandler()
{
  // delete the FE_Element and DOF_Group objects
  if (theFEs != 0) {
    for (int i=0; i<numFE; i++)
      if (theFEs[i] != 0) 
        delete theFEs[i];
    delete [] theFEs;
  }
	    
  if (theDOFs != 0) {
    for (int j=0; j<numDOF; j++)
      if (theDOFs[j] != 0) 
        delete theDOFs[j];
    delete [] theDOFs;
  }
}

int
PlainHandler::handle(const ID *nodesLast)
{
    // first check links exist to a Domain and an AnalysisModel object
    Domain *theDomain = this->getDomainPtr();
    AnalysisModel *theModel = this->getAnalysisModelPtr();
    Integrator *theIntegrator = this->getIntegratorPtr();    
    
    if ((theDomain == 0) || (theModel == 0) || (theIntegrator == 0)) {
	opserr << "WARNING PlainHandler::handle() - ";
	opserr << " setLinks() has not been called\n";
	return -1;
    }

    
    // get number ofelements and nodes in the domain 
    // and init the theFEs and theDOFs arrays
    
    numFE = theDomain->getNumElements();
    numDOF = theDomain->getNumNodes();

    // create an array for the FE_elements and zero it
    if ((numFE <= 0) || ((theFEs  = new FE_Element *[numFE]) == 0)) {
	opserr << "WARNING PlainHandler::handle() - ";
        opserr << "ran out of memory for FE_elements"; 
	opserr << " array of size " << numFE << endln;
	return -2;
    }

    for (int i=0; i<numFE; i++) theFEs[i] = 0;

    // create an array for the DOF_Groups and zero it
    if ((numDOF <= 0) || ((theDOFs = new DOF_Group *[numDOF]) == 0)) {
	opserr << "WARNING PlainHandler::handle() - ";
	opserr << "ran out of memory for DOF_Groups";
	opserr << " array of size " << numDOF << endln;
	return -3;    
    }    
    for (int j=0; j<numDOF; j++) theDOFs[j] = 0;

    // initialse the DOF_Groups and add them to the AnalysisModel.
    //    : must of course set the initial IDs
    NodeIter &theNod = theDomain->getNodes();
    Node *nodPtr;
    SP_Constraint *spPtr;
    DOF_Group *dofPtr;
    
    int numDofGrp = 0;
    int count3 = 0;
    int countDOF =0;
    while ((nodPtr = theNod()) != 0) {
	if ((dofPtr = new DOF_Group(numDofGrp, nodPtr)) == 0) {
	    opserr << "WARNING PlainHandler::handle() - ran out of memory";
	    opserr << " creating DOF_Group " << numDofGrp << endln;	
	    return -4;    		
	}
	// initially set all the ID value to -2
	const ID &id = dofPtr->getID();
	for (int j=0; j < id.Size(); j++) {
	    dofPtr->setID(j,-2);
	    countDOF++;
	}
	// loop through the SP_Constraints to see if any of the
	// DOFs are constrained, if so set initial ID value to -1
	int nodeID = nodPtr->getTag();
	SP_ConstraintIter &theSPs = theDomain->getDomainAndLoadPatternSPs();
	while ((spPtr = theSPs()) != 0)
	    if (spPtr->getNodeTag() == nodeID) {
		if (spPtr->isHomogeneous() == false) {
		    opserr << "WARNING PlainHandler::handle() - ";
		    opserr << " non-homogeneos constraint";
		    opserr << " for node " << spPtr->getNodeTag();
		    opserr << " homo assumed\n";
		}
		const ID &id = dofPtr->getID();
		int dof = spPtr->getDOF_Number();		
		if (id(dof) == -2) {
			dofPtr->setID(spPtr->getDOF_Number(),-1);
			countDOF--;	
		} else {
		    opserr << "WARNING PlainHandler::handle() - ";
		    opserr << " multiple single pointconstraints at DOF " << dof;
		    opserr << " for node " << spPtr->getNodeTag() << endln;
		}
	    }

    	// loop through the MP_Constraints to see if any of the
	// DOFs are constrained, note constraint matrix must be diagonal
	// with 1's on the diagonal
	MP_ConstraintIter &theMPs = theDomain->getMPs();
	MP_Constraint *mpPtr;
	while ((mpPtr = theMPs()) != 0)
	    if (mpPtr->getNodeConstrained() == nodeID) {
		if (mpPtr->isTimeVarying() == true) {
		    opserr << "WARNING PlainHandler::handle() - ";
		    opserr << " time-varying constraint";
		    opserr << " for node " << nodeID;
		    opserr << " non-varyng assumed\n";
		}
		const Matrix &C = mpPtr->getConstraint();
		int numRows = C.noRows();
		int numCols = C.noCols();
		if (numRows != numCols) {
			opserr << "WARNING PlainHandler::handle() - ";
			opserr << " constraint matrix not diagonal, ignoring constraint";
			opserr << " for node " << nodeID << endln;
			opserr << " non-varyng assumed\n";
		} else {
			int ok = 0;
			for (int i=0; i<numRows; i++) {
				if (C(i,i) != 1.0) ok = 1;
				for (int j=0; j<numRows; j++)
					if (i != j)
						if (C(i,j) != 0.0)
							ok = 1;
			}
			if (ok != 0) {
				opserr << "WARNING PlainHandler::handle() - ";
				opserr << " constraint matrix not identity, ignoring constraint";
				opserr << " for node " << nodeID << endln;
				opserr << " non-varyng assumed\n";
			} else {
				const ID &dofs = mpPtr->getConstrainedDOFs();
				const ID &id = dofPtr->getID();				
				for (int i=0; i<dofs.Size(); i++) {
					int dof = dofs(i);	
					if (id(dof) == -2) {
						dofPtr->setID(dof,-4);
						countDOF--;	
					} else {
		    				opserr << "WARNING PlainHandler::handle() - ";
		    				opserr << " constraint at dof " << dof << " already specified for constrained node";
		    				opserr << " in MP_Constraint at node " << nodeID << endln;
					}
					
				}
			}
		}
	}

	nodPtr->setDOF_GroupPtr(dofPtr);
	theDOFs[numDofGrp] = dofPtr;
	numDofGrp++;	
	theModel->addDOF_Group(dofPtr);
    }

    // set the number of eqn in the model
    theModel->setNumEqn(countDOF);

    // now see if we have to set any of the dof's to -3
    //    int numLast = 0;
    if (nodesLast != 0) 
	for (int i=0; i<nodesLast->Size(); i++) {
	    int nodeID = (*nodesLast)(i);
	    Node *nodPtr = theDomain->getNode(nodeID);
	    if (nodPtr != 0) {
		dofPtr = nodPtr->getDOF_GroupPtr();
		
		const ID &id = dofPtr->getID();
		// set all the dof values to -3
		for (int j=0; j < id.Size(); j++) 
		    if (id(j) == -2) {
			dofPtr->setID(j,-3);
			count3++;
		    } else {
			opserr << "WARNING PlainHandler::handle() ";
			opserr << " - boundary sp constraint in subdomain";
			opserr << " this should not be - results suspect \n";
		    }
	    }
	}
    
    // initialise the FE_Elements and add to the AnalysisModel.
    ElementIter &theEle = theDomain->getElements();
    Element *elePtr;

    int numFe = 0;    
    FE_Element *fePtr;
    while ((elePtr = theEle()) != 0) {
	if ((fePtr = new FE_Element(elePtr)) == 0) {
	    opserr << "WARNING PlainHandler::handle() - ran out of memory";
	    opserr << " creating FE_Element " << elePtr->getTag() << endln; 
	    return -5;
	}		

	
	theFEs[numFe] = fePtr;
	numFe++;
	
	theModel->addFE_Element(fePtr);
	if (elePtr->isSubdomain() == true) {
	    Subdomain *theSub = (Subdomain *)elePtr;
	    theSub->setFE_ElementPtr(fePtr);
	}
    }

    return count3;
}


void 
PlainHandler::clearAll(void)
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
PlainHandler::sendSelf(int cTag,
		       Channel &theChannel)
{
    return 0;
}

int
PlainHandler::recvSelf(int cTag, 
		       Channel &theChannel, 
		       FEM_ObjectBroker &theBroker)
{
    return 0;
}
