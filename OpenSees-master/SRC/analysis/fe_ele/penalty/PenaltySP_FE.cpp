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
                                                                        
// $Revision: 1.8 $
// $Date: 2009-10-13 21:13:01 $
// $Source: /usr/local/cvs/OpenSees/SRC/analysis/fe_ele/penalty/PenaltySP_FE.cpp,v $
                                                                        
                                                                        
// File: ~/analysis/fe_ele/penalty/PenaltySP_FE.C
//
// Written: fmk 
// Created: 11/96
// Revision: A
//
// Purpose: This file contains the code for implementing the methods
// of the PenaltySP_FE class interface.
//
// the interface:

#include <PenaltySP_FE.h>
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

Matrix PenaltySP_FE::tang(1,1);
Vector PenaltySP_FE::resid(1);

PenaltySP_FE::PenaltySP_FE(int tag, Domain &theDomain, 
			   SP_Constraint &TheSP, double Alpha)
:FE_Element(tag, 1,1), alpha(Alpha),
 theSP(&TheSP), theNode(0)
{
    // get a pointer to the Node
    theNode = theDomain.getNode(theSP->getNodeTag());
    if (theNode == 0) {
	opserr << "FATAL PenaltySP_FE::PenaltySP_FE() - no Node: ";
	opserr << theSP->getNodeTag() << "in domain\n";
	exit(-1);
    }

    // set the DOF_Group tags
    DOF_Group *dofGrpPtr = theNode->getDOF_GroupPtr();
    if (dofGrpPtr != 0) 
	myDOF_Groups(0) = dofGrpPtr->getTag();	    
}


PenaltySP_FE::~PenaltySP_FE()
{

}    

// void setID(int index, int value);
//	Method to set the corresponding index of the ID to value.

int
PenaltySP_FE::setID(void)
{
    DOF_Group *theNodesDOFs = theNode->getDOF_GroupPtr();
    if (theNodesDOFs == 0) {
	opserr << "WARNING PenaltySP_FE::setID(void) - no DOF_Group with Node\n";
	return -2;
    }    
    myDOF_Groups(0) = theNodesDOFs->getTag();
    
    int restrainedDOF = theSP->getDOF_Number();
    if (restrainedDOF < 0 || restrainedDOF >= theNode->getNumberDOF()) {
	opserr << "WARNING PenaltySP_FE::setID(void) - unknown DOF ";
	opserr << restrainedDOF << " at Node\n";
	return -3;
    }    	
    const ID &theNodesID = theNodesDOFs->getID();
    if (restrainedDOF >= theNodesID.Size()) {
	opserr << "WARNING PenaltySP_FE::setID(void) - ";
	opserr << " Nodes DOF_Group too small\n";
	return -4;
    }    		
    
    myID(0) = theNodesID(restrainedDOF);

    return 0;
}


const Matrix &
PenaltySP_FE::getTangent(Integrator *theNewIntegrator)
{
  tang(0,0) = alpha;
  return tang;
}


const Vector &
PenaltySP_FE::getResidual(Integrator *theNewIntegrator)
{
    double constraint = theSP->getValue();
    int constrainedDOF = theSP->getDOF_Number();
    const Vector &nodeDisp = theNode->getTrialDisp();
	
    if (constrainedDOF < 0 || constrainedDOF >= nodeDisp.Size()) {
	opserr << "WARNING PenaltySP_FE::getTangForce() - ";	
	opserr << " constrained DOF " << constrainedDOF << " outside disp\n";
	resid(0) = 0;
    }

    //    (*resid)(0) = alpha * (constraint - nodeDisp(constrainedDOF));    
    // is replace with the following to remove possible problems with
    // subtracting very small numbers

    resid(0) = alpha * (constraint - nodeDisp(constrainedDOF));    

    return resid;
}


const Vector &
PenaltySP_FE::getTangForce(const Vector &disp, double fact)
{
    double constraint = theSP->getValue();
    int constrainedID = myID(0);
    if (constrainedID < 0 || constrainedID >= disp.Size()) {
	opserr << "WARNING PenaltySP_FE::getTangForce() - ";	
	opserr << " constrained DOF " << constrainedID << " outside disp\n";
	resid(0) = 0.0;
	return resid;
    }
    resid(0) = alpha * disp(constrainedID);

    return resid;
}

const Vector &
PenaltySP_FE::getK_Force(const Vector &disp, double fact)
{
  opserr << "WARNING PenaltySP_FE::getK_Force() - not yet implemented\n";
  resid(0) = 0.0;
  return resid;
}

const Vector &
PenaltySP_FE::getKi_Force(const Vector &disp, double fact)
{
  opserr << "WARNING PenaltySP_FE::getKi_Force() - not yet implemented\n";
  resid(0) = 0.0;
  return resid;
}


const Vector &
PenaltySP_FE::getC_Force(const Vector &disp, double fact)
{
  // opserr << "WARNING PenaltySP_FE::getC_Force() - not yet implemented\n";
  resid(0) = 0.0;
  return resid;
}

const Vector &
PenaltySP_FE::getM_Force(const Vector &disp, double fact)
{
  // opserr << "WARNING PenaltySP_FE::getM_Force() - not yet implemented\n";
  resid(0) = 0.0;
  return resid;
}




