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
                                                                        
// $Revision: 1.3 $
// $Date: 2003-02-14 23:00:51 $
// $Source: /usr/local/cvs/OpenSees/SRC/analysis/numberer/PlainNumberer.cpp,v $
                                                                        
                                                                        
// File: ~/analysis/numberer/PlainNumberer.C
// 
// Written: fmk 
// Created: 9/96
// Revision: A
//
// Description: This file contains the class definition for PlainNumberer.
// PlainNumberer is a subclass of DOF_Numberer. The PlainNumberer assigns
// equation numbers to the DOFs on a first come first serve basis; that is 
// it gets the DOF_GrpIter and assigns the equation numbers to the DOFs
// as it iterates through the iter.
//
// What: "@(#) PlainNumberer.C, revA"

#include <PlainNumberer.h>
#include <AnalysisModel.h>

#include <DOF_Group.h>
#include <DOF_GrpIter.h>
#include <FE_Element.h>
#include <FE_EleIter.h>
#include <ID.h>
#include <Channel.h>
#include <FEM_ObjectBroker.h>

#include <Domain.h>
#include <MP_Constraint.h>
#include <Node.h>
#include <MP_ConstraintIter.h>

void* OPS_PlainNumberer()
{
    return new PlainNumberer();
}

PlainNumberer::PlainNumberer() 
:DOF_Numberer(NUMBERER_TAG_PlainNumberer)
{
}

PlainNumberer::~PlainNumberer() 
{
}


// int numberDOF(void)
//	Method to number the unnumbered DOFs in the DOF Groups. It does so
//	on a first come fist serve basis, first assigning all DOFs with a -2
//	a number, then all DOFs with a -3 a number.

int
PlainNumberer::numberDOF(int lastDOF)
{
    int eqnNumber = 0; // start equation number = 0

    // get a pointer to the model & check its not null
    AnalysisModel *theModel = this->getAnalysisModelPtr();
    Domain *theDomain = 0;
    if (theModel != 0) theDomain = theModel->getDomainPtr();

    if (theModel == 0 || theDomain == 0) {
	opserr << "WARNING PlainNumberer::numberDOF(int) -";
	opserr << " - no AnalysisModel - has setLinks() been invoked?\n";
	return -1;
    }
    
    if (lastDOF != -1) {
	opserr << "WARNING PlainNumberer::numberDOF(int lastDOF):";
	opserr << " does not use the lastDOF as requested\n";
    }
    
    // iterate throgh  the DOFs first time setting -2 values
    DOF_GrpIter &theDOFs = theModel->getDOFs();
    DOF_Group *dofPtr;
    
    while ((dofPtr = theDOFs()) != 0) {
	const ID &theID = dofPtr->getID();
	for (int i=0; i<theID.Size(); i++)
	    if (theID(i) == -2) 
	      dofPtr->setID(i,eqnNumber++);
    }

    // iterate throgh  the DOFs second time setting -3 values
    DOF_GrpIter &moreDOFs = theModel->getDOFs();
    
    while ((dofPtr = moreDOFs()) != 0) {
	const ID &theID = dofPtr->getID();
	for (int i=0; i<theID.Size(); i++)
	    if (theID(i) == -3) dofPtr->setID(i,eqnNumber++);
    }

    // iterate through the DOFs one last time setting any -4 values
    DOF_GrpIter &tDOFs = theModel->getDOFs();
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

    eqnNumber--;
    int numEqn = eqnNumber - START_EQN_NUMBER +1;
	
    // iterate through the FE_Element getting them to set their IDs
    FE_EleIter &theEle = theModel->getFEs();
    FE_Element *elePtr;
    while ((elePtr = theEle()) != 0)
	elePtr->setID();

    // set the numOfEquation in the Model
    theModel->setNumEqn(numEqn);

    return numEqn;
}


int
PlainNumberer::sendSelf(int cTag, Channel &theChannel)
{
    return 0;
}

int
PlainNumberer::recvSelf(int cTag, 
			Channel &theChannel, 
			FEM_ObjectBroker &theBroker)
{
    return 0;
}


int
PlainNumberer::numberDOF(ID &lastDOFs)
{
   int eqnNumber = 0; // start equation number = 0
    
    // get a pointer to the model & check its not null
    AnalysisModel *theModel = this->getAnalysisModelPtr();
    Domain *theDomain =0;
    if (theModel != 0) theDomain = theModel->getDomainPtr();

    if (theModel == 0 || theDomain == 0) {
	opserr << "WARNING PlainNumberer::numberDOF(int) -";
	opserr << " - no AnalysisModel - has setLinks() been invoked?\n";
	return -1;
    }
    
    opserr << "WARNING PlainNumberer::numberDOF(ID):";
    opserr << " does not use the lastDOFs as requested\n";
    
    // iterate throgh  the DOFs first time setting -2 values
    DOF_GrpIter &theDOFs = theModel->getDOFs();
    DOF_Group *dofPtr;
    
    while ((dofPtr = theDOFs()) != 0) {
	const ID &theID = dofPtr->getID();
	for (int i=0; i<theID.Size(); i++)
	    if (theID(i) == -2) dofPtr->setID(i,eqnNumber++);
    }

    // iterate throgh  the DOFs first time setting -3 values
    DOF_GrpIter &moreDOFs = theModel->getDOFs();
    
    while ((dofPtr = moreDOFs()) != 0) {
	const ID &theID = dofPtr->getID();
	for (int i=0; i<theID.Size(); i++)
	    if (theID(i) == -3) dofPtr->setID(i,eqnNumber++);
    }
    // iterate through the DOFs one last time setting any -4 values
    // iterate through the DOFs one last time setting any -4 values
    DOF_GrpIter &tDOFs = theModel->getDOFs();
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

    eqnNumber--;
    int numEqn = eqnNumber - START_EQN_NUMBER +1;
	
    // iterate through the FE_Element getting them to set their IDs
    FE_EleIter &theEle = theModel->getFEs();
    FE_Element *elePtr;
    while ((elePtr = theEle()) != 0)
	elePtr->setID();

    // set the numOfEquation in the Model
    theModel->setNumEqn(numEqn);

    return numEqn;
}
