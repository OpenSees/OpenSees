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
// $Source: /usr/local/cvs/OpenSees/SRC/analysis/fe_ele/lagrange/LagrangeMP_FE.cpp,v $


// Written: fmk
// Created: 02/99
// Revision: A
//
// Purpose: This file contains the code for implementing the methods
// of the LagrangeMP_FE class interface.

#include <LagrangeMP_FE.h>
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

// Note 1:
//    changed by M.Petracca: 2024. should be the constrained dof size.
//    now it's the same because all MP constraints in opensees are one-to-one.
//    but in the future we may implement generic constraints of the form:
//    CDOF = a1*RDOF1 + a2*RDOF2 + ... + an*RDOFn + beta

LagrangeMP_FE::LagrangeMP_FE(int tag, Domain &theDomain, MP_Constraint &TheMP,
                             DOF_Group &theGroup, double Alpha)
    : FE_Element(tag, 3, (tag, TheMP.getConstrainedDOFs()).Size() +
                             (TheMP.getRetainedDOFs()).Size() +
                             (TheMP.getConstrainedDOFs()).Size()),  // *see note 1
      theMP(&TheMP),
      theConstrainedNode(0), theRetainedNode(0),
      theDofGroup(&theGroup), tang(0), sysTang(0), resid(0),
      timeVarying(TheMP.isTimeVarying()), urLoaded(false),
      alpha(Alpha),
      numU(TheMP.getConstrainedDOFs().Size() + TheMP.getRetainedDOFs().Size())
{
    const Matrix &constraint = theMP->getConstraint();
    int noRows = constraint.noRows();
    int noCols = constraint.noCols();
    int size = 2 * noRows + noCols;

    tang = new Matrix(size, size);
    sysTang = new Matrix(size, size);
    resid = new Vector(size);
    if (tang == 0 || sysTang == 0 || resid == 0 ||
        tang->noCols() == 0 || sysTang->noCols() == 0 || resid->Size() == 0) {
        opserr << "FATAL LagrangeMP_FE::LagrangeMP_FE() - out of memory\n";
        exit(-1);
    }
    tang->Zero();
    sysTang->Zero();
    resid->Zero();

    theRetainedNode = theDomain.getNode(theMP->getNodeRetained());
    theConstrainedNode = theDomain.getNode(theMP->getNodeConstrained());

    if (theRetainedNode == 0) {
        opserr << "WARNING LagrangeMP_FE::LagrangeMP_FE()";
        opserr << "- no associated Retained Node\n";
        exit(-1);
    }

    if (theConstrainedNode == 0) {
        opserr << "WARNING LagrangeMP_FE::LagrangeMP_FE()";
        opserr << "- no associated Constrained Node\n";
        exit(-1);
    }

    if (!timeVarying)
        this->determineTangent();

    // set the myDOF_Groups tags indicating the attached id's of the
    // DOF_Group objects
    DOF_Group *theConstrainedNodesDOFs = theConstrainedNode->getDOF_GroupPtr();
    if (theConstrainedNodesDOFs == 0) {
        opserr << "WARNING LagrangeMP_FE::LagrangeMP_FE()";
        opserr << " - no DOF_Group with Constrained Node\n";
        exit(-1);
    }

    DOF_Group *theRetainedNodesDOFs = theRetainedNode->getDOF_GroupPtr();
    if (theRetainedNodesDOFs == 0) {
        opserr << "WARNING LagrangeMP_FE::LagrangeMP_FE()";
        opserr << " - no DOF_Group with Retained Node\n";
        exit(-1);
    }

    myDOF_Groups(0) = theConstrainedNodesDOFs->getTag();
    myDOF_Groups(1) = theRetainedNodesDOFs->getTag();
    myDOF_Groups(2) = theDofGroup->getTag();
}

LagrangeMP_FE::~LagrangeMP_FE()
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
LagrangeMP_FE::setID(void)
{
    int result = 0;

    // first determine the IDs in myID for those DOFs marked
    // as constrained DOFs, this is obtained from the DOF_Group
    // associated with the constrained node
    if (theConstrainedNode == 0) {
        opserr << "WARNING LagrangeMP_FE::setID(void)";
        opserr << "- no associated Constrained Node\n";
        return -1;
    }
    DOF_Group *theConstrainedNodesDOFs = theConstrainedNode->getDOF_GroupPtr();
    if (theConstrainedNodesDOFs == 0) {
        opserr << "WARNING LagrangeMP_FE::setID(void)";
        opserr << " - no DOF_Group with Constrained Node\n";
        return -2;
    }

    const ID &constrainedDOFs = theMP->getConstrainedDOFs();
    const ID &theConstrainedNodesID = theConstrainedNodesDOFs->getID();

    int size1 = constrainedDOFs.Size();
    for (int i = 0; i < size1; i++) {
        int constrained = constrainedDOFs(i);
        if (constrained < 0 ||
            constrained >= theConstrainedNode->getNumberDOF()) {
            opserr << "WARNING LagrangeMP_FE::setID(void) - unknown DOF ";
            opserr << constrained << " at Node\n";
            myID(i) = -1;  // modify so nothing will be added to equations
            result = -3;
        } else {
            if (constrained >= theConstrainedNodesID.Size()) {
                opserr << "WARNING LagrangeMP_FE::setID(void) - ";
                opserr << " Nodes DOF_Group too small\n";
                myID(i) = -1;  // modify so nothing will be added to equations
                result = -4;
            } else {
                myID(i) = theConstrainedNodesID(constrained);
            }
        }
    }

    // now determine the IDs for the retained dof's
    if (theRetainedNode == 0) {
        opserr << "WARNING LagrangeMP_FE::setID(void)";
        opserr << "- no associated Retained Node\n";
        return -1;
    }
    DOF_Group *theRetainedNodesDOFs = theRetainedNode->getDOF_GroupPtr();
    if (theRetainedNodesDOFs == 0) {
        opserr << "WARNING LagrangeMP_FE::setID(void)";
        opserr << " - no DOF_Group with Retained Node\n";
        return -2;
    }

    const ID &RetainedDOFs = theMP->getRetainedDOFs();
    const ID &theRetainedNodesID = theRetainedNodesDOFs->getID();

    int size2 = RetainedDOFs.Size();
    for (int j = 0; j < size2; j++) {
        int retained = RetainedDOFs(j);
        if (retained < 0 || retained >= theRetainedNode->getNumberDOF()) {
            opserr << "WARNING LagrangeMP_FE::setID(void) - unknown DOF ";
            opserr << retained << " at Node\n";
            myID(j + size1) = -1;  // modify so nothing will be added
            result = -3;
        } else {
            if (retained >= theRetainedNodesID.Size()) {
                opserr << "WARNING LagrangeMP_FE::setID(void) - ";
                opserr << " Nodes DOF_Group too small\n";
                myID(j + size1) = -1;  // modify so nothing will be added
                result = -4;
            } else {
                myID(j + size1) = theRetainedNodesID(retained);
            }
        }
    }

    // finally set the ID corresponding to the ID's at the LagrangeDOF_Group
    const ID &theGroupsID = theDofGroup->getID();
    int size3 = theGroupsID.Size();
    for (int k = 0; k < size3; k++)
        myID(k + size1 + size2) = theGroupsID(k);

    return result;
}

const Matrix &
LagrangeMP_FE::getTangent(Integrator *theNewIntegrator)
{
    if (theNewIntegrator != 0) {
        theNewIntegrator->formEleTangent(this);
        return *sysTang;
    }
    return this->getStaticTangent();
}

void
LagrangeMP_FE::zeroTangent(void)
{
    sysTang->Zero();
    urLoaded = false;
}

void
LagrangeMP_FE::addKtToTang(double fact)
{
    // Copy static coupling, then scale lower-left C by integrator fact (c1).
    // Upper-right C^T stays unscaled. HALL_TANGENT may call addKi after addKt:
    // the second call only accumulates the lower-left block.
    if (fact == 0.0)
        return;

    const Matrix &Ks = this->getStaticTangent();
    const int nLambda = Ks.noRows() - numU;

    if (!urLoaded) {
        *sysTang = Ks;
        if (fact != 1.0) {
            for (int i = 0; i < nLambda; i++)
                for (int j = 0; j < numU; j++)
                    (*sysTang)(numU + i, j) *= fact;
        }
        urLoaded = true;
    } else {
        for (int i = 0; i < nLambda; i++)
            for (int j = 0; j < numU; j++)
                (*sysTang)(numU + i, j) += Ks(numU + i, j) * fact;
    }
}

void
LagrangeMP_FE::addKiToTang(double fact)
{
    this->addKtToTang(fact);
}

void
LagrangeMP_FE::addCtoTang(double fact)
{
    // no damping contribution from lagrange constraint
}

void
LagrangeMP_FE::addMtoTang(double fact)
{
    // no mass contribution from lagrange constraint
}

const Vector &
LagrangeMP_FE::getResidual(Integrator *theNewIntegrator)
{
    if (theNewIntegrator != 0)
        theNewIntegrator->formEleResidual(this);
    return *resid;
}

void
LagrangeMP_FE::zeroResidual(void)
{
    resid->Zero();
}

void
LagrangeMP_FE::addRtoResidual(double fact)
{
    if (fact == 0.0)
        return;

    // get the solution vector [Uc Ur lambda]
    static Vector UU;
    const ID &id1 = theMP->getConstrainedDOFs();
    const ID &id2 = theMP->getRetainedDOFs();
    const ID &id3 = theDofGroup->getID();
    int size = id1.Size() + id2.Size() + id3.Size();
    UU.resize(size);

    const Vector &Uc = theConstrainedNode->getTrialDisp();
    const Vector &Ur = theRetainedNode->getTrialDisp();
    const Vector &Uc0 = theMP->getConstrainedDOFsInitialDisplacement();
    const Vector &Ur0 = theMP->getRetainedDOFsInitialDisplacement();
    const Vector &lambda = theDofGroup->getTrialDisp();

    for (int i = 0; i < id1.Size(); ++i) {
        int cdof = id1(i);
        if (cdof < 0 || cdof >= Uc.Size()) {
            opserr << "LagrangeMP_FE::addRtoResidual FATAL Error: Constrained DOF "
                   << cdof << " out of bounds [0-" << Uc.Size() << "]\n";
            exit(-1);
        }
        UU(i) = Uc(cdof) - Uc0(i);
    }
    for (int i = 0; i < id2.Size(); ++i) {
        int rdof = id2(i);
        if (rdof < 0 || rdof >= Ur.Size()) {
            opserr << "LagrangeMP_FE::addRtoResidual FATAL Error: Retained DOF "
                   << rdof << " out of bounds [0-" << Ur.Size() << "]\n";
            exit(-1);
        }
        UU(i + id1.Size()) = Ur(rdof) - Ur0(i);
    }
    for (int i = 0; i < id3.Size(); ++i)
        UU(i + id1.Size() + id2.Size()) = lambda(i);

    /*
    R = -C*U + G
       .R = generalized residual vector
       .C = constraint matrix
       .U = generalized solution vector (displacement, lagrange multipliers)
       .G = constrain values for non-homogeneous MP constraints (not available now)
    | Ru |    | 0  A | | u |   | 0 |
    |    | = -|      |*|   | + |   |
    | Rl |    | A  0 | | l |   | g |
    */

    // compute residual from unscaled static coupling (no c1)
    resid->addMatrixVector(1.0, this->getStaticTangent(), UU, -fact);
}

void
LagrangeMP_FE::addRIncInertiaToResidual(double fact)
{
    this->addRtoResidual(fact);
}

void
LagrangeMP_FE::addM_Force(const Vector &accel, double fact)
{
    // no-op
}

void
LagrangeMP_FE::addD_Force(const Vector &vel, double fact)
{
    // no-op
}

const Vector &
LagrangeMP_FE::getTangForce(const Vector &disp, double fact)
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
        opserr << "WARNING LagrangeMP_FE::getTangForce() - ";
        opserr << "- addMatrixVector returned error\n";
    }

    return *resid;
}

const Vector &
LagrangeMP_FE::getK_Force(const Vector &disp, double fact)
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
        opserr << "WARNING LagrangeMP_FE::getK_Force() - ";
        opserr << "- addMatrixVector returned error\n";
    }

    return *resid;
}

const Vector &
LagrangeMP_FE::getKi_Force(const Vector &disp, double fact)
{
    return this->getK_Force(disp, fact);
}

const Vector &
LagrangeMP_FE::getC_Force(const Vector &disp, double fact)
{
    resid->Zero();
    return *resid;
}

const Vector &
LagrangeMP_FE::getM_Force(const Vector &disp, double fact)
{
    resid->Zero();
    return *resid;
}

void
LagrangeMP_FE::determineTangent(void)
{
    const Matrix &constraint = theMP->getConstraint();
    int noRows = constraint.noRows();
    int noCols = constraint.noCols();
    int n = noRows + noCols;

    tang->Zero();

    for (int j = 0; j < noRows; j++) {
        (*tang)(n + j, j) = -alpha;
        (*tang)(j, n + j) = -alpha;
    }

    for (int i = 0; i < noRows; i++)
        for (int j = 0; j < noCols; j++) {
            double val = constraint(i, j) * alpha;
            (*tang)(n + i, j + noRows) = val;
            (*tang)(noRows + j, n + i) = val;
        }
}
