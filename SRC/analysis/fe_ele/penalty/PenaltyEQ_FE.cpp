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
// $Source: /usr/local/cvs/OpenSees/SRC/analysis/fe_ele/penalty/PenaltyEQ_FE.cpp,v $
                                                                        
                                                                        
// Written: Yuli Huang (yulee@berkeley.edu)
// Created: 05/2020
// Revision: A
//
// Purpose: This file contains the code for implementing the methods
// of the PenaltyEQ_FE class interface.
//
// the interface:

#include <PenaltyEQ_FE.h>
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

PenaltyEQ_FE::PenaltyEQ_FE(int tag, Domain &theDomain, 
			   EQ_Constraint &TheEQ, double Alpha)
:FE_Element(tag, 1+(TheEQ.getRetainedDOFs()).Size(),
 1+(TheEQ.getRetainedDOFs()).Size()),
 theEQ(&TheEQ), theConstrainedNode(0) , theRetainedNode(0),
 tang(0), sysTang(0), resid(0), C(0),
 timeVarying(TheEQ.isTimeVarying()), alpha(Alpha)
{
    int size = 1 + (theEQ->getRetainedDOFs()).Size();

    tang = new Matrix(size,size);
    sysTang = new Matrix(size,size);
    resid = new Vector(size);
    C = new Vector(size);

    if (tang == 0 || sysTang == 0 || resid == 0 || C == 0 ||
            tang->noCols() != size || sysTang->noCols() != size ||
            C->Size() != size || resid->Size() != size) {
        opserr << "FATAL PenaltyEQ_FE::PenaltyEQ_FE() - out of memory\n";
        exit(-1);
    }

    theConstrainedNode = theDomain.getNode(theEQ->getNodeConstrained());

    if (theConstrainedNode == 0) {
        opserr << "FATAL PenaltyEQ_FE::PenaltyEQ_FE() - Constrained";
        opserr << " Node does not exist in Domain\n";
        opserr << theEQ->getNodeConstrained() << endln;
        exit(-1);
    }

    DOF_Group *theConstrainedNodeDOFGrpPtr = theConstrainedNode->getDOF_GroupPtr();
    if (theConstrainedNodeDOFGrpPtr != 0) 
        myDOF_Groups(0) = theConstrainedNodeDOFGrpPtr->getTag();	        
    else
        opserr << "WARNING PenaltyEQ_FE::PenaltyEQ_FE() - node no Group yet?\n"; 

    const ID &nodeRetained = theEQ->getNodeRetained();
    theRetainedNode = new Node*[nodeRetained.Size()];
    for (int i = 0; i < nodeRetained.Size(); ++i) {
        theRetainedNode[i] = theDomain.getNode(nodeRetained(i));
        if (theRetainedNode[i] == 0) {
            opserr << "FATAL PenaltyEQ_FE::PenaltyEQ_FE() - Retained";
            opserr << " Node does not exist in Domain\n";
            opserr << nodeRetained(i) << endln;
            exit(-1);
        }
        DOF_Group *theRetainedNodeDOFGrpPtr = theRetainedNode[i]->getDOF_GroupPtr();
        if (theRetainedNodeDOFGrpPtr != 0) 
            myDOF_Groups(i + 1) = theRetainedNodeDOFGrpPtr->getTag();	    
        else 
            opserr << "WARNING PenaltyEQ_FE::PenaltyEQ_FE() - node no Group yet?\n";
    }
    
    if (!timeVarying) {
        this->determineTangent();
        // we can free up the space taken by C as it is no longer needed
        if (C != 0)
            delete C;
        C = 0;
    }
}

PenaltyEQ_FE::~PenaltyEQ_FE()
{
    if (theRetainedNode != 0)
        delete theRetainedNode;
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
//	Method to set the corresponding index of the ID to value.
int
PenaltyEQ_FE::setID(void)
{
    int result = 0;

    // first determine the IDs in myID for those DOFs marked
    // as constrainedDOF DOFs, this is obtained from the DOF_Group
    // associated with the constrainedDOF node
    if (theConstrainedNode == 0) {
        opserr << "WARNING PenaltyEQ_FE::setID(void)";
        opserr << "- no asscoiated Constrained Node\n";
        return -1;
    }
    DOF_Group *theConstrainedNodesDOFs = theConstrainedNode->getDOF_GroupPtr();
    if (theConstrainedNodesDOFs == 0) {
        opserr << "WARNING PenaltyEQ_FE::setID(void)";
        opserr << " - no DOF_Group with Constrained Node\n";
        return -2;
    }    
    const ID &theConstrainedNodesID = theConstrainedNodesDOFs->getID();    
    
    int constrainedDOF = theEQ->getConstrainedDOFs();
    if (constrainedDOF < 0 || constrainedDOF >= theConstrainedNode->getNumberDOF()) {
        opserr << "WARNING PenaltyEQ_FE::setID(void) - unknown DOF ";
        opserr << constrainedDOF << " at Node\n";
        myID(0) = -1; // modify so nothing will be added to equations
        result = -3;
    }    	
    else {
        if (constrainedDOF >= theConstrainedNodesID.Size()) {
            opserr << "WARNING PenaltyEQ_FE::setID(void) - ";
            opserr << " Nodes DOF_Group too small\n";
            myID(0) = -1; // modify so nothing will be added to equations
            result = -4;
        }
        else
            myID(0) = theConstrainedNodesID(constrainedDOF);
    }
    
    // now determine the IDs for the retained dof's
    if (theRetainedNode == 0) {
        opserr << "WARNING PenaltyEQ_FE::setID(void)";
        opserr << "- no asscoiated Retained Node\n";
        return -1;
    }
    const ID &nodeRetained = theEQ->getNodeRetained();
    const ID &RetainedDOFs = theEQ->getRetainedDOFs();
    for (int i = 0; i < nodeRetained.Size(); ++i) {
        if (theRetainedNode[i] == 0) {
            opserr << "WARNING PenaltyEQ_FE::setID(void)";
            opserr << "- no asscoiated Retained Node\n";
            return -1;
        }
        DOF_Group *theRetainedNodesDOFs = theRetainedNode[i]->getDOF_GroupPtr();
        if (theRetainedNodesDOFs == 0) {
            opserr << "WARNING PenaltyEQ_FE::setID(void)";
            opserr << " - no DOF_Group with Retained Node\n";
            return -2;
        }    
        const ID &theRetainedNodesID = theRetainedNodesDOFs->getID();    
        int retained = RetainedDOFs(i);
        if (retained < 0 || retained >= theRetainedNode[i]->getNumberDOF()) {
            opserr << "WARNING PenaltyEQ_FE::setID(void) - unknown DOF ";
            opserr << retained << " at Node\n";
            myID(i + 1) = -1; // modify so nothing will be added
            result = -3;
        }    	
        else {
            if (retained >= theRetainedNodesID.Size()) {
                opserr << "WARNING PenaltyEQ_FE::setID(void) - ";
                opserr << " Nodes DOF_Group too small\n";
                myID(i + 1) = -1; // modify so nothing will be added 
                result = -4;
            }
            else
               myID(i + 1) = theRetainedNodesID(retained);
        }
    }

    return result;
}

const Matrix &
PenaltyEQ_FE::getTangent(Integrator *theNewIntegrator)
{
    if (theNewIntegrator != 0)
	    theNewIntegrator->formEleTangent(this);
    return *sysTang;
}


void
PenaltyEQ_FE::zeroTangent(void)
{
    sysTang->Zero();
}


void
PenaltyEQ_FE::addKtToTang(double fact)
{
    // Always accumulate: HALL_TANGENT calls both addKtToTang and
    // addKiToTang, so fact==1.0 must not replace sysTang.
    if (fact == 0.0)
        return;

    const Matrix &Ks = this->getStaticTangent();
    sysTang->addMatrix(1.0, Ks, fact);
}


void
PenaltyEQ_FE::addKiToTang(double fact)
{
    this->addKtToTang(fact);
}


void
PenaltyEQ_FE::addCtoTang(double fact)
{
    // no damping contribution from penalty constraint
}


void
PenaltyEQ_FE::addMtoTang(double fact)
{
    // no mass contribution from penalty constraint
}


const Vector &
PenaltyEQ_FE::getResidual(Integrator *theNewIntegrator)
{
    if (theNewIntegrator != 0)
        theNewIntegrator->formEleResidual(this);
    return *resid;
}


void
PenaltyEQ_FE::zeroResidual(void)
{
    resid->Zero();
}


void
PenaltyEQ_FE::addRtoResidual(double fact)
{
    if (fact == 0.0)
        return;

    // get the solution vector [Uc Ur]
    static Vector UU;
    const ID &RetainedDOFs = theEQ->getRetainedDOFs();
    int size = 1 + RetainedDOFs.Size();
    UU.resize(size);
    const Vector& Uc = theConstrainedNode->getTrialDisp();
    double Uc0 = theEQ->getConstrainedDOFsInitialDisplacement();
    int cdof = theEQ->getConstrainedDOFs();
    if (cdof < 0 || cdof >= Uc.Size()) {
        opserr << "PenaltyEQ_FE::addRtoResidual FATAL Error: Constrained DOF " << cdof << " out of bounds [0-" << Uc.Size() - 1 << "]\n";
        exit(-1);
    }
    UU(0) = Uc(cdof) - Uc0;

    const Vector& Ur0 = theEQ->getRetainedDOFsInitialDisplacement();
    for (int i = 0; i < RetainedDOFs.Size(); ++i) {
        int rdof = RetainedDOFs(i);
        const Vector& Ur = theRetainedNode[i]->getTrialDisp();
        if (rdof < 0 || rdof >= Ur.Size()) {
            opserr << "PenaltyEQ_FE::addRtoResidual FATAL Error: Retained DOF " << rdof << " out of bounds [0-" << Ur.Size() - 1 << "]\n";
            exit(-1);
        }
        UU(i + 1) = Ur(rdof) - Ur0(i);
    }

    // residual contribution = -R with R = K_static*U, same sign as before
    resid->addMatrixVector(1.0, this->getStaticTangent(), UU, -fact);
}


void
PenaltyEQ_FE::addRIncInertiaToResidual(double fact)
{
    // no mass/damping on the constraint
    this->addRtoResidual(fact);
}


void
PenaltyEQ_FE::addM_Force(const Vector &accel, double fact)
{
    // no-op
}


void
PenaltyEQ_FE::addD_Force(const Vector &vel, double fact)
{
    // no-op
}


const Vector &
PenaltyEQ_FE::getTangForce(const Vector &disp, double fact)
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
        opserr << "WARNING PenaltyEQ_FE::getTangForce() - ";
        opserr << "- addMatrixVector returned error\n";
    }

    return *resid;
}

const Vector &
PenaltyEQ_FE::getK_Force(const Vector &disp, double fact)
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
        opserr << "WARNING PenaltyEQ_FE::getK_Force() - ";
        opserr << "- addMatrixVector returned error\n";
    }

    return *resid;
}

const Vector &
PenaltyEQ_FE::getKi_Force(const Vector &disp, double fact)
{
    return this->getK_Force(disp, fact);
}

const Vector &
PenaltyEQ_FE::getC_Force(const Vector &disp, double fact)
{
    resid->Zero();
    return *resid;
}

const Vector &
PenaltyEQ_FE::getM_Force(const Vector &disp, double fact)
{
    resid->Zero();
    return *resid;
}

void  
PenaltyEQ_FE::determineTangent(void)
{
    // first determine [C] = [-I [Ccr]]
    C->Zero();
    const Vector &constraint = theEQ->getConstraint();
    int size = constraint.Size();
    
    (*C)(0) = -1.0;
    
    for (int i = 0; i < size; i++)
        (*C)(i + 1) = constraint(i);
    
    // now form the tangent: [K] = alpha * [C]^t[C]
    const Vector &Cref = *C;
    *tang = alpha * (Cref % Cref);
}
