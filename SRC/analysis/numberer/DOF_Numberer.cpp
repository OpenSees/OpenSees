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
// $Date: 2005-11-28 21:28:17 $
// $Source: /usr/local/cvs/OpenSees/SRC/analysis/numberer/DOF_Numberer.cpp,v $
                                                                        
                                                                        
// Written: fmk 
// Created: 9/96
// Revision: A
//
// Description: This file contains the class implementation for DOF_Numberer.
// DOF_Numberer is an abstract base class, i.e. no objects of it's
// type can be created. 
//
// What: "@(#) DOF_Numberer.h, revA"

#include <DOF_Numberer.h>
#include <AnalysisModel.h>
#include <GraphNumberer.h>
#include <ID.h>
#include <DOF_Group.h>
#include <FE_Element.h>
#include <FE_EleIter.h>
#include <Channel.h>
#include <Channel.h>
#include <FEM_ObjectBroker.h>
#include <Channel.h>
#include <FEM_ObjectBroker.h>

#include <Graph.h>

#include <Domain.h>
#include <MP_Constraint.h>
#include <Node.h>
#include <MP_ConstraintIter.h>
#include <DOF_GrpIter.h>
// Constructor

DOF_Numberer::DOF_Numberer(int clsTag) 
:MovableObject(clsTag),
 theAnalysisModel(0), theGraphNumberer(0)
{

}

DOF_Numberer::DOF_Numberer(GraphNumberer &aGraphNumberer)
:MovableObject(NUMBERER_TAG_DOF_Numberer),
 theAnalysisModel(0), theGraphNumberer(&aGraphNumberer)
{

}    

DOF_Numberer::DOF_Numberer()
:MovableObject(NUMBERER_TAG_DOF_Numberer),
 theAnalysisModel(0), theGraphNumberer(0)
{

}    

// Destructor

DOF_Numberer::~DOF_Numberer() 
{
  if (theGraphNumberer != 0)
    delete theGraphNumberer;
}


// void setLinks(AnalysisModel &theModel)
//	Method to set theAnalysisModel pointer to be address of theModel

void
DOF_Numberer::setLinks(AnalysisModel &theModel)
{
    theAnalysisModel = &theModel;
}



int 
DOF_Numberer::numberDOF(int lastDOF_Group) 
{
    // check we have a model and a numberer
    Domain *theDomain = 0;
    if (theAnalysisModel != 0) theDomain = theAnalysisModel->getDomainPtr();
    if ((theAnalysisModel == 0) || (theDomain == 0)) {
	opserr << "WARNING DOF_Numberer::numberDOF - ";
	opserr << "Pointers are not set\n";
	return -1;
    }
    
    if ((theGraphNumberer == 0)) {
	opserr << "WARNING DOF_Numberer::numberDOF - ";
	opserr << "subclasses must provide own implementation\n";
	return -2;
    }    

    // check we cant do quick return
    
    if (theAnalysisModel->getNumDOF_Groups() == 0)
	return 0;

    // we first number the dofs using the dof group graph

    const ID &orderedRefs = theGraphNumberer->
      number(theAnalysisModel->getDOFGroupGraph(), lastDOF_Group);     

    theAnalysisModel->clearDOFGroupGraph();

    // we now iterate through the DOFs first time setting -2 values

    int eqnNumber = 0;
    
    if (orderedRefs.Size() != theAnalysisModel->getNumDOF_Groups()) {
	opserr << "WARNING DOF_Numberer::numberDOF - ";
	opserr << "Incompatable Sizes\n";
	return -3;
    }
    int result = 0;
    
    int size = orderedRefs.Size();
    for (int i=0; i<size; i++) {
	int dofTag = orderedRefs(i);
	DOF_Group *dofPtr;	
	dofPtr = theAnalysisModel->getDOF_GroupPtr(dofTag);
	if (dofPtr == 0) {
	    opserr << "WARNING DOF_Numberer::numberDOF - ";
	    opserr << "DOF_Group " << dofTag << "not in AnalysisModel!\n";
	    result = -4;
	} else {
	    const ID &theID = dofPtr->getID();
	    int idSize = theID.Size();
	    for (int j=0; j<idSize; j++)
		if (theID(j) == -2) dofPtr->setID(j,eqnNumber++);
	}
    }

    // iterate throgh  the DOFs second time setting -3 values
    for (int k=0; k<size; k++) {
	int dofTag = orderedRefs(k);
	DOF_Group *dofPtr;	
	dofPtr = theAnalysisModel->getDOF_GroupPtr(dofTag);
	if (dofPtr != 0) {
	    const ID &theID = dofPtr->getID();
	    int idSize = theID.Size();
	    for (int j=0; j<idSize; j++)
		if (theID(j) == -3) dofPtr->setID(j,eqnNumber++);
	}
    }


    // iterate through the DOFs one last time setting any -4 values
    // iterate throgh  the DOFs second time setting -3 values
    DOF_GrpIter &tDOFs = theAnalysisModel->getDOFs();
    DOF_Group *dofPtr;
    while ((dofPtr = tDOFs()) != 0) {
    	const ID &theID = dofPtr->getID();
    	int have4s = 0;
	for (int i=0; i<theID.Size(); i++)
	    if (theID(i) == -4) have4s = 1;

	if (have4s == 1) {
		int nodeID = dofPtr->getNodeTag();
		// loop through the MP_Constraints to see if any of the
		// DOFs are constrained, note constraint matrix must be diagonal
		// with 1's on the diagonal
		MP_ConstraintIter &theMPs = theDomain->getMPs();
		MP_Constraint *mpPtr;
		while ((mpPtr = theMPs()) != 0 ) {
			// note keep looping over all in case multiple constraints
			// are used to constrain a node -- can't assume intelli user
	    		if (mpPtr->getNodeConstrained() == nodeID) {
	    			int nodeRetained = mpPtr->getNodeRetained();
	    			Node *nodeRetainedPtr = theDomain->getNode(nodeRetained);
	    			DOF_Group *retainedDOF = nodeRetainedPtr->getDOF_GroupPtr();
	    			const ID&retainedDOFIDs = retainedDOF->getID();
	    			const ID&constrainedDOFs = mpPtr->getConstrainedDOFs();
	    			const ID&retainedDOFs = mpPtr->getRetainedDOFs();
	    			for (int i=0; i<constrainedDOFs.Size(); i++) {
	    				int dofC = constrainedDOFs(i);
	    				int dofR = retainedDOFs(i);
	    				int dofID = retainedDOFIDs(dofR);
	    				dofPtr->setID(dofC, dofID);
	    			}
	    		}
		}		
	}	
    }


    int numEqn = eqnNumber;

    // iterate through the FE_Element getting them to set their IDs
    FE_EleIter &theEle = theAnalysisModel->getFEs();
    FE_Element *elePtr;
    while ((elePtr = theEle()) != 0)
	elePtr->setID();

    // set the numOfEquation in the Model
    theAnalysisModel->setNumEqn(numEqn);

    if (result == 0)
	return numEqn;

    return result;
}






int 
DOF_Numberer::numberDOF(ID &lastDOFs) 
{
    // check we have a model and a numberer
    	Domain *theDomain = 0;
   if (theAnalysisModel != 0) theDomain = theAnalysisModel->getDomainPtr();
   if ((theAnalysisModel == 0) || (theDomain == 0)) {
	opserr << "WARNING DOF_Numberer::numberDOF - ";
	opserr << "Pointers are not set\n";
	return -1;
    }
    
    if ((theGraphNumberer == 0)) {
	opserr << "WARNING DOF_Numberer::numberDOF - ";
	opserr << "subclasses must provide own implementation\n";
	return -2;
    }    

    // check we cant do quick return
    if (theAnalysisModel->getNumDOF_Groups() == 0)
	return 0;

    // we first number the dofs using the dof group graph
	
    const ID &orderedRefs = theGraphNumberer->
	number(theAnalysisModel->getDOFGroupGraph(), lastDOFs);     

    theAnalysisModel->clearDOFGroupGraph();

    // we now iterate through the DOFs first time setting -2 values

    int eqnNumber = 0;
    if (orderedRefs.Size() != theAnalysisModel->getNumDOF_Groups()) {
	opserr << "WARNING DOF_Numberer::numberDOF - ";
	opserr << "Incompatable Sizes\n";
	return -3;
    }

    int result =0;
    int size = orderedRefs.Size();
    for (int i=0; i<size; i++) {
	int dofTag = orderedRefs(i);
	DOF_Group *dofPtr;	
	dofPtr = theAnalysisModel->getDOF_GroupPtr(dofTag);
	if (dofPtr == 0) {
	    opserr << "WARNING DOF_Numberer::numberDOF - ";
	    opserr << "DOF_Group " << dofTag << "not in AnalysisModel!\n";
	    result = -4;
	} else {
	    const ID &theID = dofPtr->getID();
	    int idSize = theID.Size();
	    for (int j=0; j<idSize; j++)
		if (theID(j) == -2) dofPtr->setID(j,eqnNumber++);
	}	
    }

    // iterate throgh  the DOFs first time setting -3 values
    
    for (int k=0; k<size; k++) {
	int dofTag = orderedRefs(k);
	DOF_Group *dofPtr;	
	dofPtr = theAnalysisModel->getDOF_GroupPtr(dofTag);
	if (dofPtr != 0) {
	    const ID &theID = dofPtr->getID();
	    int idSize = theID.Size();
	    for (int j=0; j<idSize; j++)
		if (theID(j) == -3) dofPtr->setID(j,eqnNumber++);
	}
    }

    // iterate through the DOFs one last time setting any -4 values
    // iterate throgh  the DOFs second time setting -3 values
    DOF_GrpIter &tDOFs = theAnalysisModel->getDOFs();
    DOF_Group *dofPtr;
    while ((dofPtr = tDOFs()) != 0) {
    	const ID &theID = dofPtr->getID();
    	int have4s = 0;
	for (int i=0; i<theID.Size(); i++)
	    if (theID(i) == -4) have4s = 1;

	if (have4s == 1) {
		int nodeID = dofPtr->getNodeTag();
		// loop through the MP_Constraints to see if any of the
		// DOFs are constrained, note constraint matrix must be diagonal
		// with 1's on the diagonal
		MP_ConstraintIter &theMPs = theDomain->getMPs();
		MP_Constraint *mpPtr;
		while ((mpPtr = theMPs()) != 0 ) {
			// note keep looping over all in case multiple constraints
			// are used to constrain a node -- can't assume intelli user
	    		if (mpPtr->getNodeConstrained() == nodeID) {
	    			int nodeRetained = mpPtr->getNodeRetained();
	    			Node *nodeRetainedPtr = theDomain->getNode(nodeRetained);
	    			DOF_Group *retainedDOF = nodeRetainedPtr->getDOF_GroupPtr();
	    			const ID&retainedDOFIDs = retainedDOF->getID();
	    			const ID&constrainedDOFs = mpPtr->getConstrainedDOFs();
	    			const ID&retainedDOFs = mpPtr->getRetainedDOFs();
	    			for (int i=0; i<constrainedDOFs.Size(); i++) {
	    				int dofC = constrainedDOFs(i);
	    				int dofR = retainedDOFs(i);
	    				int dofID = retainedDOFIDs(dofR);
	    				dofPtr->setID(dofC, dofID);
	    			}
	    		}
		}		
	}	
    }

    int numEqn = eqnNumber;

    // iterate through the FE_Element getting them to set their IDs
    FE_EleIter &theEle = theAnalysisModel->getFEs();
    FE_Element *elePtr;
    while ((elePtr = theEle()) != 0)
	elePtr->setID();

    // set the numOfEquation in the Model
    theAnalysisModel->setNumEqn(numEqn);
    
    if (result == 0)
	return numEqn;
    
    return result;

}



// AnalysisModel *getAnalysisModelPtr(void)
// 	Method to return a pointer to theAnalysisModel for subclasses.

int
DOF_Numberer::sendSelf(int cTag, Channel &theChannel)
{
    ID data(2);
    int dataTag = this->getDbTag();

    data(0) = -1;
    if (theGraphNumberer != 0) {
	data(0) = theGraphNumberer->getClassTag();
	data(1) = theGraphNumberer->getDbTag();
    }
    theChannel.sendID(dataTag, cTag, data);
    if (theGraphNumberer != 0)
      theGraphNumberer->sendSelf(cTag, theChannel);

    return 0;
}

int
DOF_Numberer::recvSelf(int cTag, Channel &theChannel, 
		       FEM_ObjectBroker &theBroker)
{
  ID data(2);
  int dataTag = this->getDbTag();
  theChannel.recvID(dataTag, cTag, data);    

  // get a graphNumberer
  if (data(0) != -1) {
    theGraphNumberer = theBroker.getPtrNewGraphNumberer(data(0));
    if (theGraphNumberer != 0) {
      theGraphNumberer->setDbTag(data(1));
      theGraphNumberer->recvSelf(cTag, theChannel,theBroker);
    }
    else {
      opserr << "DOF_Numberer::recvSelf() - failed to get GraphNumberer\n";
      return -1;
    }
  }

  return 0;
}


AnalysisModel *
DOF_Numberer::getAnalysisModelPtr(void) const
{
    return theAnalysisModel;
}


GraphNumberer *
DOF_Numberer::getGraphNumbererPtr(void) const
{
    return theGraphNumberer;
}

