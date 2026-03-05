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
                                                                        
// $Revision: 1.0 $
// $Date: 2025-05-09$
// $Source: /usr/local/cvs/OpenSees/SRC/analysis/fe_ele/lagrange/LagrangeEQ_FE.cpp,v $
                                                                        
                                                                        
// Written: Yuli Huang (yulee@berkeley.edu)
// Created: 05/2020
// Revision: A
//
// Purpose: This file contains the code for implementing the methods
// of the LagrangeEQ_FE class interface.

#include <LagrangeEQ_FE.h>
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
#include <EQ_Constraint.h>
#include <DOF_Group.h>

LagrangeEQ_FE::LagrangeEQ_FE(int tag, Domain &theDomain, EQ_Constraint &TheEQ,
			     DOF_Group &theGroup, double Alpha)
:FE_Element(tag, 2+(TheEQ.getRetainedDOFs()).Size(),
 2+(TheEQ.getRetainedDOFs()).Size()),
 alpha(Alpha), theEQ(&TheEQ), 
 theConstrainedNode(0), theRetainedNode(0),
 theDofGroup(&theGroup), tang(0), resid(0)
{
    int size = 2 + (theEQ->getRetainedDOFs()).Size();
    
    tang = new Matrix(size,size);
    resid = new Vector(size);
    if (tang == 0 || resid == 0 || tang->noCols() == 0 || resid->Size() == 0) {
	    opserr << "FATAL LagrangeEQ_FE::LagrangeEQ_FE() - out of memory\n";
	    exit(-1);
    }
    tang->Zero();	
    resid->Zero();

    theConstrainedNode = theDomain.getNode(theEQ->getNodeConstrained());

    if (theConstrainedNode == 0) {
	    opserr << "FATAL LagrangeEQ_FE::LagrangeEQ_FE() - Constrained";
	    opserr << " Node does not exist in Domain\n";
	    opserr << theEQ->getNodeConstrained() << endln;
	    exit(-1);
    }

    DOF_Group *theConstrainedNodeDOFGrpPtr = theConstrainedNode->getDOF_GroupPtr();
    if (theConstrainedNodeDOFGrpPtr != 0) 
        myDOF_Groups(0) = theConstrainedNodeDOFGrpPtr->getTag();	        
    else
        opserr << "WARNING LagrangeEQ_FE::LagrangeEQ_FE() - no DOF_Group with Constrained Node\n"; 

    const ID &nodeRetained = theEQ->getNodeRetained();
    theRetainedNode = new Node*[nodeRetained.Size()];
    for (int i = 0; i < nodeRetained.Size(); ++i) {
        theRetainedNode[i] = theDomain.getNode(nodeRetained(i));
        if (theRetainedNode[i] == 0) {
            opserr << "FATAL LagrangeEQ_FE::LagrangeEQ_FE() - Retained";
            opserr << " Node does not exist in Domain\n";
            opserr << nodeRetained(i) << endln;
            exit(-1);
        }
        DOF_Group *theRetainedNodeDOFGrpPtr = theRetainedNode[i]->getDOF_GroupPtr();
        if (theRetainedNodeDOFGrpPtr != 0) 
            myDOF_Groups(i+1) = theRetainedNodeDOFGrpPtr->getTag();	    
        else 
            opserr << "WARNING LagrangeEQ_FE::LagrangeEQ_FE() - no DOF_Group with Retained Node\n";
    }

    myDOF_Groups(nodeRetained.Size() + 1) = theDofGroup->getTag();

    if (theEQ->isTimeVarying() == false) {
        this->determineTangent();
    }
}

LagrangeEQ_FE::~LagrangeEQ_FE()
{
    if (tang != 0)
	delete tang;
    if (resid != 0)
	delete resid;
}    

// void setID(int index, int value);
//	Method to set the corresponding index of the ID to value.
int
LagrangeEQ_FE::setID(void)
{
    int result = 0;

    // first determine the IDs in myID for those DOFs marked
    // as constrainedDOF DOFs, this is obtained from the DOF_Group
    // associated with the constrainedDOF node
    if (theConstrainedNode == 0) {
        opserr << "WARNING LagrangeEQ_FE::setID(void)";
        opserr << "- no asscoiated Constrained Node\n";
        return -1;
    }
    DOF_Group *theConstrainedNodesDOFs = theConstrainedNode->getDOF_GroupPtr();
    if (theConstrainedNodesDOFs == 0) {
        opserr << "WARNING LagrangeEQ_FE::setID(void)";
        opserr << " - no DOF_Group with Constrained Node\n";
        return -2;
    }    
    const ID &theConstrainedNodesID = theConstrainedNodesDOFs->getID();    

    int constrainedDOF = theEQ->getConstrainedDOFs();
	if (constrainedDOF < 0 || constrainedDOF >= theConstrainedNode->getNumberDOF()) {
	    opserr << "WARNING LagrangeEQ_FE::setID(void) - unknown DOF ";
	    opserr << constrainedDOF << " at Node\n";
	    myID(0) = -1; // modify so nothing will be added to equations
	    result = -3;
	}    	
	else {
	    if (constrainedDOF >= theConstrainedNodesID.Size()) {
		    opserr << "WARNING LagrangeEQ_FE::setID(void) - ";
		    opserr << " Nodes DOF_Group too small\n";
		    myID(0) = -1; // modify so nothing will be added to equations
		    result = -4;
	    }
	    else
		    myID(0) = theConstrainedNodesID(constrainedDOF);
	}
    
    // now determine the IDs for the retained dof's
    if (theRetainedNode == 0) {
        opserr << "WARNING LagrangeEQ_FE::setID(void)";
        opserr << "- no asscoiated Retained Node\n";
        return -1;
    }
    const ID &nodeRetained = theEQ->getNodeRetained();
    const ID &RetainedDOFs = theEQ->getRetainedDOFs();
    for (int i = 0; i < nodeRetained.Size(); ++i) {
        if (theRetainedNode[i] == 0) {
            opserr << "WARNING LagrangeEQ_FE::setID(void)";
            opserr << "- no asscoiated Retained Node\n";
            return -1;
        }
        DOF_Group *theRetainedNodesDOFs = theRetainedNode[i]->getDOF_GroupPtr();
        if (theRetainedNodesDOFs == 0) {
            opserr << "WARNING LagrangeEQ_FE::setID(void)";
            opserr << " - no DOF_Group with Retained Node\n";
            return -2;
        }    
        const ID &theRetainedNodesID = theRetainedNodesDOFs->getID();    
        int retained = RetainedDOFs(i);
        if (retained < 0 || retained >= theRetainedNode[i]->getNumberDOF()) {
            opserr << "WARNING LagrangeEQ_FE::setID(void) - unknown DOF ";
            opserr << retained << " at Node\n";
            myID(i+1) = -1; // modify so nothing will be added
            result = -3;
        }    	
        else {
            if (retained >= theRetainedNodesID.Size()) {
                opserr << "WARNING LagrangeEQ_FE::setID(void) - ";
                opserr << " Nodes DOF_Group too small\n";
                myID(i+1) = -1; // modify so nothing will be added 
                result = -4;
            }
            else
               myID(i+1) = theRetainedNodesID(retained);
        }
    }

    // finally set the ID corresponding to the ID's at the LagrangeDOF_Group
    const ID &theGroupsID = theDofGroup->getID();
	myID(nodeRetained.Size() + 1) = theGroupsID(0);
    
    return result;
}

const Matrix &
LagrangeEQ_FE::getTangent(Integrator *theNewIntegrator)
{
    if (theEQ->isTimeVarying() == true)
	this->determineTangent();

    return *tang;
}

const Vector &
LagrangeEQ_FE::getResidual(Integrator *theNewIntegrator)
{
    // get the solution vector [Uc Ur lambda]
    static Vector UU;
    const ID &RetainedDOFs = theEQ->getRetainedDOFs();
    const ID& id2 = theEQ->getRetainedDOFs();
    int size = 2 + RetainedDOFs.Size();
    UU.resize(size);
    const Vector& Uc = theConstrainedNode->getTrialDisp();
    double Uc0 = theEQ->getConstrainedDOFsInitialDisplacement();
    int cdof = theEQ->getConstrainedDOFs();
    if (cdof < 0 || cdof >= Uc.Size()) {
        opserr << "LagrangeEQ_FE::getResidual FATAL Error: Constrained DOF " << cdof << " out of bounds [0-" << Uc.Size() - 1 << "]\n";
        exit(-1);
    }
    UU(0) = Uc(cdof) - Uc0;

    const Vector& Ur0 = theEQ->getRetainedDOFsInitialDisplacement();
    for (int i = 0; i < RetainedDOFs.Size(); ++i) {
        int rdof = RetainedDOFs(i);
        const Vector& Ur = theRetainedNode[i]->getTrialDisp();
        if (rdof < 0 || rdof >= Ur.Size()) {
            opserr << "LagrangeEQ_FE::getResidual FATAL Error: Retained DOF " << rdof << " out of bounds [0-" << Ur.Size() - 1 << "]\n";
            exit(-1);
        }
        UU(i+1) = Ur(rdof) - Ur0(i);
    }

    const Vector& lambda = theDofGroup->getTrialDisp();
    UU(RetainedDOFs.Size() + 1) = lambda(0);

    // compute residual
    const Matrix& KK = getTangent(theNewIntegrator);
    resid->addMatrixVector(0.0, KK, UU, -1.0);

    // done
    return *resid;
}



const Vector &
LagrangeEQ_FE::getTangForce(const Vector &disp, double fact)
{
 opserr << "WARNING lagrangeEQ_FE::getTangForce() - not yet implemented\n";
 return *resid;
}


const Vector &
LagrangeEQ_FE::getK_Force(const Vector &disp, double fact)
{
 opserr << "WARNING lagrangeEQ_FE::getK_Force() - not yet implemented\n";
 return *resid;
}

const Vector &
LagrangeEQ_FE::getKi_Force(const Vector &disp, double fact)
{
 opserr << "WARNING LagrangeEQ_FE::getKi_Force() - not yet implemented\n";
 return *resid;
}

const Vector &
LagrangeEQ_FE::getC_Force(const Vector &disp, double fact)
{
 opserr << "WARNING lagrangeEQ_FE::getC_Force() - not yet implemented\n";
 return *resid;
}

const Vector &
LagrangeEQ_FE::getM_Force(const Vector &disp, double fact)
{
 opserr << "WARNING lagrangeEQ_FE::getM_Force() - not yet implemented\n";
 return *resid;
}

void  
LagrangeEQ_FE::determineTangent(void)
{
    const Vector &constraint = theEQ->getConstraint();
    int size = constraint.Size();
    
    tang->Zero();    

    (*tang)(size + 1, 0) = -alpha;
    (*tang)(0, size + 1) = -alpha;	
    
    for (int i = 0; i < size; i++) {
        double val = constraint(i) * alpha;
        (*tang)(size + 1, i + 1) = val;
        (*tang)(i + 1, size + 1) = val;
    }
}





