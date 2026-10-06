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
// $Source: /usr/local/cvs/OpenSees/SRC/analysis/fe_ele/penalty/PenaltyMP_FE.cpp,v $
                                                                        
                                                                        
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

PenaltyMP_FE::PenaltyMP_FE(int tag, Domain &theDomain, 
			   MP_Constraint &TheMP, double Alpha)
:FE_Element(tag, 2,(TheMP.getConstrainedDOFs()).Size()+
 (TheMP.getRetainedDOFs()).Size()),
 theMP(&TheMP), theConstrainedNode(0) , theRetainedNode(0),
 tang(0), sysTang(0), resid(0), C(0),
 timeVarying(TheMP.isTimeVarying()), alpha(Alpha)
{
    
    int size;
    const ID &id1 = theMP->getConstrainedDOFs();
    size = id1.Size();
    const ID &id2 = theMP->getRetainedDOFs();    
    size += id2.Size();

    tang = new Matrix(size,size);
    sysTang = new Matrix(size,size);
    resid = new Vector(size);
    C = new Matrix(id1.Size(),size);

    if (tang == 0 || sysTang == 0 || resid == 0 || C == 0 ||
	tang->noCols() != size || sysTang->noCols() != size ||
	C->noCols() != size || resid->Size() != size) {
	opserr << "FATAL PenaltyMP_FE::PenaltyMP_FE() - out of memory\n";
	exit(-1);
    }
	    
    theRetainedNode = theDomain.getNode(theMP->getNodeRetained());    
    theConstrainedNode = theDomain.getNode(theMP->getNodeConstrained());

    if (theRetainedNode == 0 || theConstrainedNode == 0) {
	opserr << "FATAL PenaltyMP_FE::PenaltyMP_FE() - Constrained or Retained";
	opserr << " Node does not exist in Domain\n";
	opserr << theMP->getNodeRetained() << " " << theMP->getNodeConstrained() << endln;
	exit(-1);
    }	


    // set up the dof groups tags
    DOF_Group *dofGrpPtr = 0;
    dofGrpPtr = theRetainedNode->getDOF_GroupPtr();
    if (dofGrpPtr != 0) 
	myDOF_Groups(0) = dofGrpPtr->getTag();	    
    else 
	opserr << "WARNING PenaltyMP_FE::PenaltyMP_FE() - node no Group yet?\n"; 
    dofGrpPtr = theConstrainedNode->getDOF_GroupPtr();
    if (dofGrpPtr != 0) 
	myDOF_Groups(1) = dofGrpPtr->getTag();	        
    else
	opserr << "WARNING PenaltyMP_FE::PenaltyMP_FE() - node no Group yet?\n"; 
    
    
    if (!timeVarying) {
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
    if (sysTang != 0)
	delete sysTang;
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
	opserr << "WARNING PenaltyMP_FE::setID(void)";
	opserr << " - no DOF_Group with Constrained Node\n";
	return -2;
    }    

    const ID &constrainedDOFs = theMP->getConstrainedDOFs();
    const ID &theConstrainedNodesID = theConstrainedNodesDOFs->getID();    
    
    int size1 = constrainedDOFs.Size();
    for (int i=0; i<size1; i++) {
	int constrained = constrainedDOFs(i);
	if (constrained < 0 || 
	    constrained >= theConstrainedNode->getNumberDOF()) {
	    
	    opserr << "WARNING PenaltyMP_FE::setID(void) - unknown DOF ";
	    opserr << constrained << " at Node\n";
	    myID(i) = -1; // modify so nothing will be added to equations
	    result = -3;
	}    	
	else {
	    if (constrained >= theConstrainedNodesID.Size()) {
		opserr << "WARNING PenaltyMP_FE::setID(void) - ";
		opserr << " Nodes DOF_Group too small\n";
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
	opserr << "WARNING PenaltyMP_FE::setID(void)";
	opserr << " - no DOF_Group with Retained Node\n";
	return -2;
    }    
    
    const ID &RetainedDOFs = theMP->getRetainedDOFs();
    const ID &theRetainedNodesID = theRetainedNodesDOFs->getID();    

    int size2 = RetainedDOFs.Size();
    for (int j=0; j<size2; j++) {
	int retained = RetainedDOFs(j);
	if (retained < 0 || retained >= theRetainedNode->getNumberDOF()) {
	    opserr << "WARNING PenaltyMP_FE::setID(void) - unknown DOF ";
	    opserr << retained << " at Node\n";
	    myID(j+size1) = -1; // modify so nothing will be added
	    result = -3;
	}    	
	else {
	    if (retained >= theRetainedNodesID.Size()) {
		opserr << "WARNING PenaltyMP_FE::setID(void) - ";
		opserr << " Nodes DOF_Group too small\n";
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
    if (theNewIntegrator != 0)
	theNewIntegrator->formEleTangent(this);
    return *sysTang;
}


void
PenaltyMP_FE::zeroTangent(void)
{
    sysTang->Zero();
}


void
PenaltyMP_FE::addKtToTang(double fact)
{
    // Always accumulate: HALL_TANGENT calls both addKtToTang and
    // addKiToTang, so fact==1.0 must not replace sysTang.
    if (fact == 0.0)
        return;

    const Matrix &Ks = this->getStaticTangent();
    sysTang->addMatrix(1.0, Ks, fact);
}


void
PenaltyMP_FE::addKiToTang(double fact)
{
    this->addKtToTang(fact);
}


void
PenaltyMP_FE::addCtoTang(double fact)
{
    // no damping contribution from penalty constraint
}


void
PenaltyMP_FE::addMtoTang(double fact)
{
    // no mass contribution from penalty constraint
}


const Vector &
PenaltyMP_FE::getResidual(Integrator *theNewIntegrator)
{
    if (theNewIntegrator != 0)
        theNewIntegrator->formEleResidual(this);
    return *resid;
}


void
PenaltyMP_FE::zeroResidual(void)
{
    resid->Zero();
}


void
PenaltyMP_FE::addRtoResidual(double fact)
{
    if (fact == 0.0)
        return;

    // get the solution vector [Uc Ur]
    static Vector UU;
    const ID& id1 = theMP->getConstrainedDOFs();
    const ID& id2 = theMP->getRetainedDOFs();
    int size = id1.Size() + id2.Size();
    UU.resize(size);
    const Vector& Uc = theConstrainedNode->getTrialDisp();
    const Vector& Ur = theRetainedNode->getTrialDisp();
    const Vector& Uc0 = theMP->getConstrainedDOFsInitialDisplacement();
    const Vector& Ur0 = theMP->getRetainedDOFsInitialDisplacement();
    for (int i = 0; i < id1.Size(); ++i) {
        int cdof = id1(i);
        if (cdof < 0 || cdof >= Uc.Size()) {
            opserr << "PenaltyMP_FE::addRtoResidual FATAL Error: Constrained DOF " << cdof << " out of bounds [0-" << Uc.Size() << "]\n";
            exit(-1);
        }
        UU(i) = Uc(cdof) - Uc0(i);
    }
    for (int i = 0; i < id2.Size(); ++i) {
        int rdof = id2(i);
        if (rdof < 0 || rdof >= Ur.Size()) {
            opserr << "PenaltyMP_FE::addRtoResidual FATAL Error: Retained DOF " << rdof << " out of bounds [0-" << Ur.Size() << "]\n";
            exit(-1);
        }
        UU(i+id1.Size()) = Ur(rdof) - Ur0(i);
    }

    // residual contribution = -R with R = K_static*U, same sign as before
    resid->addMatrixVector(1.0, this->getStaticTangent(), UU, -fact);
}


void
PenaltyMP_FE::addRIncInertiaToResidual(double fact)
{
    // no mass/damping on the constraint
    this->addRtoResidual(fact);
}


void
PenaltyMP_FE::addM_Force(const Vector &accel, double fact)
{
    // no-op
}


void
PenaltyMP_FE::addD_Force(const Vector &vel, double fact)
{
    // no-op
}


const Vector &
PenaltyMP_FE::getTangForce(const Vector &disp, double fact)
{
    resid->Zero();

    if (fact == 0.0)
        return *resid;

    // use last integrator's system tangent (includes c1)
    const Matrix &Kt = this->getTangent(this->getLastIntegrator());

    const int size = resid->Size();
    const int dispSize = disp.Size();
    Vector tmp(size);
    for (int i = 0; i < size; i++) {
        int dof = myID(i);
        if (dof >= 0 && dof < dispSize)
            tmp(i) = disp(dof);
    }

    if (resid->addMatrixVector(0.0, Kt, tmp, fact) < 0) {
        opserr << "WARNING PenaltyMP_FE::getTangForce() - ";
        opserr << "- addMatrixVector returned error\n";
    }

    return *resid;
}

const Vector &
PenaltyMP_FE::getK_Force(const Vector &disp, double fact)
{
    resid->Zero();

    if (fact == 0.0)
        return *resid;

    const int size = resid->Size();
    const int dispSize = disp.Size();
    Vector tmp(size);  // Vector(int) zeros on construction
    for (int i = 0; i < size; i++) {
        int dof = myID(i);
        if (dof >= 0 && dof < dispSize)
            tmp(i) = disp(dof);
    }

    if (resid->addMatrixVector(0.0, this->getStaticTangent(), tmp, fact) < 0) {
        opserr << "WARNING PenaltyMP_FE::getK_Force() - ";
        opserr << "- addMatrixVector returned error\n";
    }

    return *resid;
}

const Vector &
PenaltyMP_FE::getKi_Force(const Vector &disp, double fact)
{
    return this->getK_Force(disp, fact);
}

const Vector &
PenaltyMP_FE::getC_Force(const Vector &disp, double fact)
{
    resid->Zero();
    return *resid;
}

const Vector &
PenaltyMP_FE::getM_Force(const Vector &disp, double fact)
{
    resid->Zero();
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
    const Matrix &Cref = *C;
    tang->addMatrixTransposeProduct(0.0, Cref, Cref, alpha);
}
