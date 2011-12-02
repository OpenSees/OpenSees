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
// $Date: 2000-09-15 08:23:16 $
// $Source: /usr/local/cvs/OpenSees/SRC/analysis/fe_ele/penalty/PenaltyMP_FE.cpp,v $
                                                                        
                                                                        
// File: ~/analysis/fe_ele/penalty/PenaltyMP_FE.C
//
// Written: fmk 
// Created: 11/96
// Revision: A
//
// Purpose: This file contains the code for implementing the methods
// of the PenaltyMP_FE class interface.
//
// the interface:

#include <PenaltyMP_FE.h>
#include <stdlib.h>

#include <Element.h>
#include <Domain.h>
#include <Node.h>
#include <DOF_Group.h>
#include <Integrator.h>
#include <Subdomain.h>
#include <AnalysisModel.h>
#include <Matrix.h>
#include <Vector.h>
#include <Node.h>
#include <MP_Constraint.h>
#include <DOF_Group.h>

PenaltyMP_FE::PenaltyMP_FE(Domain &theDomain, 
			   MP_Constraint &TheMP, double Alpha)
:FE_Element(2,(TheMP.getConstrainedDOFs()).Size()+
 (TheMP.getRetainedDOFs()).Size()),
 theMP(&TheMP), theConstrainedNode(0) , theRetainedNode(0),
 tang(0), resid(0), C(0), alpha(Alpha)
{
    
    int size;
    const ID &id1 = theMP->getConstrainedDOFs();
    size = id1.Size();
    const ID &id2 = theMP->getRetainedDOFs();    
    size += id2.Size();

    tang = new Matrix(size,size);
    resid = new Vector(size);
    C = new Matrix(id1.Size(),size);

    if (tang == 0 || resid == 0 || C == 0 ||
	tang->noCols() != size || C->noCols() != size || 
	resid->Size() != size) {
	cerr << "FATAL PenaltyMP_FE::PenaltyMP_FE() - out of memory\n";
	exit(-1);
    }
	    
    theRetainedNode = theDomain.getNode(theMP->getNodeRetained());    
    theConstrainedNode = theDomain.getNode(theMP->getNodeConstrained());

    if (theRetainedNode == 0 || theConstrainedNode == 0) {
	cerr << "FATAL PenaltyMP_FE::PenaltyMP_FE() - Constrained or Retained";
	cerr << " Node does not exist in Domain\n";
	exit(-1);
    }	


    // set up the dof groups tags
    DOF_Group *dofGrpPtr = 0;
    dofGrpPtr = theRetainedNode->getDOF_GroupPtr();
    if (dofGrpPtr != 0) 
	myDOF_Groups(0) = dofGrpPtr->getTag();	    
    else 
	cerr << "WARNING PenaltyMP_FE::PenaltyMP_FE() - node no Group yet?\n"; 
    dofGrpPtr = theConstrainedNode->getDOF_GroupPtr();
    if (dofGrpPtr != 0) 
	myDOF_Groups(1) = dofGrpPtr->getTag();	        
    else
	cerr << "WARNING PenaltyMP_FE::PenaltyMP_FE() - node no Group yet?\n"; 
    
    
    if (theMP->isTimeVarying() == false) {
	this->determineTangent();
	// we can free up the space taken by C as it is no longer needed
	if (C != 0)
	    delete C;
	C = 0;
    }
}

PenaltyMP_FE::~PenaltyMP_FE()
{
    if (tang != 0)
	delete tang;
    if (resid != 0)
	delete resid;
    if (C != 0)
	delete C;
}    

// void setID(int index, int value);
//	Method to set the correMPonding index of the ID to value.
int
PenaltyMP_FE::setID(void)
{
    int result = 0;

    // first determine the IDs in myID for those DOFs marked
    // as constrained DOFs, this is obtained from the DOF_Group
    // associated with the constrained node
    DOF_Group *theConstrainedNodesDOFs = theConstrainedNode->getDOF_GroupPtr();
    if (theConstrainedNodesDOFs == 0) {
	cerr << "WARNING PenaltyMP_FE::setID(void)";
	cerr << " - no DOF_Group with Constrained Node\n";
	return -2;
    }    

    const ID &constrainedDOFs = theMP->getConstrainedDOFs();
    const ID &theConstrainedNodesID = theConstrainedNodesDOFs->getID();    
    
    int size1 = constrainedDOFs.Size();
    for (int i=0; i<size1; i++) {
	int constrained = constrainedDOFs(i);
	if (constrained < 0 || 
	    constrained >= theConstrainedNode->getNumberDOF()) {
	    
	    cerr << "WARNING PenaltyMP_FE::setID(void) - unknown DOF ";
	    cerr << constrained << " at Node\n";
	    myID(i) = -1; // modify so nothing will be added to equations
	    result = -3;
	}    	
	else {
	    if (constrained >= theConstrainedNodesID.Size()) {
		cerr << "WARNING PenaltyMP_FE::setID(void) - ";
		cerr << " Nodes DOF_Group too small\n";
		myID(i) = -1; // modify so nothing will be added to equations
		result = -4;
	    }
	    else
		myID(i) = theConstrainedNodesID(constrained);
	}
    }
    
    // now determine the IDs for the retained dof's
    DOF_Group *theRetainedNodesDOFs = theRetainedNode->getDOF_GroupPtr();
    if (theRetainedNodesDOFs == 0) {
	cerr << "WARNING PenaltyMP_FE::setID(void)";
	cerr << " - no DOF_Group with Retained Node\n";
	return -2;
    }    
    
    const ID &RetainedDOFs = theMP->getRetainedDOFs();
    const ID &theRetainedNodesID = theRetainedNodesDOFs->getID();    

    int size2 = RetainedDOFs.Size();
    for (int j=0; j<size2; j++) {
	int retained = RetainedDOFs(j);
	if (retained < 0 || retained >= theRetainedNode->getNumberDOF()) {
	    cerr << "WARNING PenaltyMP_FE::setID(void) - unknown DOF ";
	    cerr << retained << " at Node\n";
	    myID(j+size1) = -1; // modify so nothing will be added
	    result = -3;
	}    	
	else {
	    if (retained >= theRetainedNodesID.Size()) {
		cerr << "WARNING PenaltyMP_FE::setID(void) - ";
		cerr << " Nodes DOF_Group too small\n";
		myID(j+size1) = -1; // modify so nothing will be added 
		result = -4;
	    }
	    else
		myID(j+size1) = theRetainedNodesID(retained);
	}
    }

    myDOF_Groups(0) = theConstrainedNodesDOFs->getTag();
    myDOF_Groups(1) = theRetainedNodesDOFs->getTag();

    return result;
}

const Matrix &
PenaltyMP_FE::getTangent(Integrator *theNewIntegrator)
{
    if (theMP->isTimeVarying() == true)
	this->determineTangent();    
    return *tang;
}

const Vector &
PenaltyMP_FE::getResidual(Integrator *theNewIntegrator)
{
    // zero residual, CD = 0
    return *resid;
}



const Vector &
PenaltyMP_FE::getTangForce(const Vector &disp, double fact)
{
    // does nothing , zero residual for CD = 0
    return *resid;
}

void  
PenaltyMP_FE::determineTangent(void)
{
    // first determine [C] = [-I [Ccr]]
    C->Zero();
    const Matrix &constraint = theMP->getConstraint();
    int noRows = constraint.noRows();
    int noCols = constraint.noCols();
    
    for (int j=0; j<noRows; j++)
	(*C)(j,j) = -1.0;
    
    for (int i=0; i<noRows; i++)
	for (int j=0; j<noCols; j++)
	    (*C)(i,j+noRows) = constraint(i,j);
    

    // now form the tangent: [K] = alpha * [C]^t[C]
    *(tang) = (*C)^(*C);
    *(tang) *= alpha;
}


