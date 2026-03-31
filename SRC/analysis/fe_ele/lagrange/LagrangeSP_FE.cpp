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
// $Date: 2006-02-08 20:20:00 $
// $Source: /usr/local/cvs/OpenSees/SRC/analysis/fe_ele/lagrange/LagrangeSP_FE.cpp,v $
                                                                        
// Written: fmk 
// Created: 02/99
// Revision: A
//
// Purpose: This file contains the code for iSPlementing the methods
// of the LagrangeSP_FE class interface.
//
// the interface:

#include <LagrangeSP_FE.h>
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
#include <SP_Constraint.h>
#include <DOF_Group.h>

LagrangeSP_FE::LagrangeSP_FE(int tag, Domain &theDomain, SP_Constraint &TheSP,
			     DOF_Group &theGroup, double Alpha)
:FE_Element(tag, 2,2),
 alpha(Alpha), tang(0), resid(0), theSP(&TheSP), theDofGroup(&theGroup)
{
    // create a Matrix and a Vector for the tangent and residual
    tang = new Matrix(2,2);
    resid = new Vector(2);
    if ((tang == 0) || (tang->noCols() == 0) || (resid == 0) ||
	(resid->Size() == 0)) {
	opserr << "WARNING LagrangeSP_FE::LagrangeSP_FE()";
	opserr << "- ran out of memory\n";
	exit(-1);
    }

    // zero the Matrix and Vector
    resid->Zero();
    tang->Zero();

    theNode = theDomain.getNode(theSP->getNodeTag());    
    if (theNode == 0) {
	opserr << "WARNING LagrangeSP_FE::LagrangeSP_FE()";
	opserr << "- no asscoiated Node\n";
	exit(-1);
    }

    // set the tangent
    (*tang)(0,1) = alpha;
    (*tang)(1,0) = alpha;
    
    // set the myDOF_Groups tags indicating the attached id's of the
    // DOF_Group objects
    DOF_Group *theNodesDOFs = theNode->getDOF_GroupPtr();
    if (theNodesDOFs == 0) {
	opserr << "WARNING LagrangeSP_FE::LagrangeSP_FE()";
	opserr << " - no DOF_Group with Constrained Node\n";
	exit(-1);	
    }    

    myDOF_Groups(0) = theNodesDOFs->getTag();
    myDOF_Groups(1) = theDofGroup->getTag();
}

LagrangeSP_FE::~LagrangeSP_FE()
{
    if (tang != 0)
	delete tang;
    if (resid != 0)
	delete resid;
}    

// void setID(int index, int value);
//	Method to set the correSPonding index of the ID to value.
int
LagrangeSP_FE::setID(void)
{
    int result = 0;

    // first determine the IDs in myID for those DOFs marked
    // as constrained DOFs, this is obtained from the DOF_Group
    // associated with the constrained node
    DOF_Group *theNodesDOFs = theNode->getDOF_GroupPtr();
    if (theNodesDOFs == 0) {
	opserr << "WARNING LagrangeSP_FE::setID(void)";
	opserr << " - no DOF_Group with Constrained Node\n";
	return -1;
    }    

    int restrainedDOF = theSP->getDOF_Number();
    const ID &theNodesID = theNodesDOFs->getID();
    
    if (restrainedDOF < 0 || restrainedDOF >= theNodesID.Size()) {
	opserr << "WARNING LagrangeSP_FE::setID(void)";
	opserr << " - restrained DOF invalid\n";
	return -2;
    }    	
    
    myID(0) = theNodesID(restrainedDOF);
    myID(1) = (theDofGroup->getID())(0);
    
    return result;
}

const Matrix &
LagrangeSP_FE::getTangent(Integrator *theIntegrator)
{
    return *tang;
}

const Vector &
LagrangeSP_FE::getResidual(Integrator *theNewIntegrator)
{
    double constraint = theSP->getValue();
    double initialValue = theSP->getInitialValue();
    int constrainedDOF = theSP->getDOF_Number();
    const Vector &nodeDisp = theNode->getTrialDisp();
    const Vector& lambda = theDofGroup->getTrialDisp();

    if (constrainedDOF < 0 || constrainedDOF >= nodeDisp.Size()) {
        opserr << "LagrangeSP_FE::getResidual() -";
        opserr << " constrained DOF " << constrainedDOF << " outside range\n";
        resid->Zero();
        return *resid;
    }
    if (lambda.Size() != 1) {
        opserr << "LagrangeSP_FE::getResidual() -";
        opserr << " Lambda.Size() = " << lambda.Size() << " != 1\n";
        resid->Zero();
        return *resid;
    }
    
    /*
    R = -C*U + G
       .R = generalized residual vector
       .C = constraint matrix
       .U = generalized solution vector (displacement, lagrange multipliers)
       .G = imposed displacement values
    | Ru |    | 0  A | | u |   | 0 |
    |    | = -|      |*|   | + |   |
    | Rl |    | A  0 | | l |   | g |
    */
    (*resid)(0) = alpha * (-lambda(0));
    (*resid)(1) = alpha *(constraint - (nodeDisp(constrainedDOF) - initialValue));
    return *resid;
}




const Vector &
LagrangeSP_FE::getTangForce(const Vector &disp, double fact)
{
    double constraint = theSP->getValue();
    int constrainedID = myID(1);
    if (constrainedID < 0 || constrainedID >= disp.Size()) {
	opserr << "WARNING LagrangeSP_FE::getTangForce() - ";	
	opserr << " constrained DOF " << constrainedID << " outside disp\n";
	(*resid)(1) = constraint*alpha;
	return *resid;
    }
    (*resid)(1) = disp(constrainedID);
    return *resid;    
}

const Vector &
LagrangeSP_FE::getK_Force(const Vector &disp, double fact)
{
 opserr << "WARNING PenaltySP_FE::getK_Force() - not yet implemented\n";
 return *resid;
}

const Vector &
LagrangeSP_FE::getKi_Force(const Vector &disp, double fact)
{
 opserr << "WARNING PenaltySP_FE::getKi_Force() - not yet implemented\n";
 return *resid;
}

const Vector &
LagrangeSP_FE::getC_Force(const Vector &disp, double fact)
{
 opserr << "WARNING PenaltySP_FE::getC_Force() - not yet implemented\n";
 return *resid;
}

const Vector &
LagrangeSP_FE::getM_Force(const Vector &disp, double fact)
{
 opserr << "WARNING PenaltySP_FE::getM_Force() - not yet implemented\n";
 return *resid;
}

