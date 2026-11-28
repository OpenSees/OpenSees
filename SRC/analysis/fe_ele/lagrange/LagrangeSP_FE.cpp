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
// Purpose: This file contains the code for implementing the methods
// of the LagrangeSP_FE class interface.

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
    : FE_Element(tag, 2, 2),
      alpha(Alpha), tang(0), sysTang(0), resid(0), urLoaded(false),
      theSP(&TheSP), theNode(0), theDofGroup(&theGroup)
{
    // create a Matrix and a Vector for the tangent and residual
    tang = new Matrix(2, 2);
    sysTang = new Matrix(2, 2);
    resid = new Vector(2);
    if (tang == 0 || tang->noCols() == 0 ||
        sysTang == 0 || sysTang->noCols() == 0 ||
        resid == 0 || resid->Size() == 0) {
        opserr << "WARNING LagrangeSP_FE::LagrangeSP_FE()";
        opserr << "- ran out of memory\n";
        exit(-1);
    }

    // zero the Matrix and Vector
    resid->Zero();
    tang->Zero();
    sysTang->Zero();

    theNode = theDomain.getNode(theSP->getNodeTag());
    if (theNode == 0) {
        opserr << "WARNING LagrangeSP_FE::LagrangeSP_FE()";
        opserr << "- no associated Node\n";
        exit(-1);
    }

    // set the tangent: static coupling on [u, lambda] is [0, alpha; alpha, 0]
    (*tang)(0, 1) = alpha;
    (*tang)(1, 0) = alpha;

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
    if (sysTang != 0)
        delete sysTang;
    if (resid != 0)
        delete resid;
}

// void setID(int index, int value);
//      Method to set the corresponding index of the ID to value.
int
LagrangeSP_FE::setID(void)
{
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

    return 0;
}

const Matrix &
LagrangeSP_FE::getTangent(Integrator *theIntegrator)
{
    if (theIntegrator != 0) {
        theIntegrator->formEleTangent(this);
        return *sysTang;
    }
    return this->getStaticTangent();
}

void
LagrangeSP_FE::zeroTangent(void)
{
    sysTang->Zero();
    urLoaded = false;
}

void
LagrangeSP_FE::addKtToTang(double fact)
{
    // Copy static coupling, then scale lower-left C by integrator fact (c1).
    // Upper-right C^T stays unscaled. HALL_TANGENT may call addKi after addKt:
    // the second call only accumulates the lower-left block.
    if (fact == 0.0)
        return;

    const Matrix &Ks = this->getStaticTangent();

    if (!urLoaded) {
        *sysTang = Ks;
        if (fact != 1.0)
            (*sysTang)(1, 0) *= fact;
        urLoaded = true;
    } else {
        (*sysTang)(1, 0) += Ks(1, 0) * fact;
    }
}

void
LagrangeSP_FE::addKiToTang(double fact)
{
    this->addKtToTang(fact);
}

void
LagrangeSP_FE::addCtoTang(double fact)
{
    // no damping contribution from lagrange constraint
}

void
LagrangeSP_FE::addMtoTang(double fact)
{
    // no mass contribution from lagrange constraint
}

const Vector &
LagrangeSP_FE::getResidual(Integrator *theNewIntegrator)
{
    if (theNewIntegrator != 0)
        theNewIntegrator->formEleResidual(this);
    return *resid;
}

void
LagrangeSP_FE::zeroResidual(void)
{
    resid->Zero();
}

void
LagrangeSP_FE::addRtoResidual(double fact)
{
    if (fact == 0.0)
        return;

    double constraint = theSP->getValue();
    double initialValue = theSP->getInitialValue();
    int constrainedDOF = theSP->getDOF_Number();
    const Vector &nodeDisp = theNode->getTrialDisp();
    const Vector &lambda = theDofGroup->getTrialDisp();

    if (constrainedDOF < 0 || constrainedDOF >= nodeDisp.Size()) {
        opserr << "LagrangeSP_FE::addRtoResidual() -";
        opserr << " constrained DOF " << constrainedDOF << " outside range\n";
        return;
    }
    if (lambda.Size() != 1) {
        opserr << "LagrangeSP_FE::addRtoResidual() -";
        opserr << " lambda.Size() = " << lambda.Size() << " != 1\n";
        return;
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
    (*resid)(0) += fact * alpha * (-lambda(0));
    (*resid)(1) += fact * alpha * (constraint - (nodeDisp(constrainedDOF) - initialValue));
}

void
LagrangeSP_FE::addRIncInertiaToResidual(double fact)
{
    this->addRtoResidual(fact);
}

void
LagrangeSP_FE::addM_Force(const Vector &accel, double fact)
{
    // no-op
}

void
LagrangeSP_FE::addD_Force(const Vector &vel, double fact)
{
    // no-op
}

const Vector &
LagrangeSP_FE::getTangForce(const Vector &disp, double fact)
{
    resid->Zero();

    if (fact == 0.0)
        return *resid;

    // use last integrator's system tangent (includes c1 on lower-left)
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
        opserr << "WARNING LagrangeSP_FE::getTangForce() - ";
        opserr << "- addMatrixVector returned error\n";
    }

    return *resid;
}

const Vector &
LagrangeSP_FE::getK_Force(const Vector &disp, double fact)
{
    resid->Zero();

    if (fact == 0.0)
        return *resid;

    const int size = resid->Size();
    const int dispSize = disp.Size();
    Vector tmp(size);
    for (int i = 0; i < size; i++) {
        int dof = myID(i);
        if (dof >= 0 && dof < dispSize)
            tmp(i) = disp(dof);
    }

    if (resid->addMatrixVector(0.0, this->getStaticTangent(), tmp, fact) < 0) {
        opserr << "WARNING LagrangeSP_FE::getK_Force() - ";
        opserr << "- addMatrixVector returned error\n";
    }

    return *resid;
}

const Vector &
LagrangeSP_FE::getKi_Force(const Vector &disp, double fact)
{
    return this->getK_Force(disp, fact);
}

const Vector &
LagrangeSP_FE::getC_Force(const Vector &disp, double fact)
{
    resid->Zero();
    return *resid;
}

const Vector &
LagrangeSP_FE::getM_Force(const Vector &disp, double fact)
{
    resid->Zero();
    return *resid;
}
