/* ****************************************************************** **
**    OpenSees - Open System for Earthquake Engineering Simulation    **
**          Pacific Earthquake Engineering Research Center            **
** (C) Copyright 1999, The Regents of the University of California    **
** All Rights Reserved.                                               **
** ****************************************************************** */
// 
// Created: Gianluca Cavallo gi.cav.2586@¶gmail.com 
// derived from SSPbrickUP (C.McGann, UW, 10.2011)
//
// This implementation adds pore-pressure regularization for undrained
// or nearly incompressible response.

#include "SSPTethraUP.h"

#include <elementAPI.h>
#include <Information.h>
#include <ElementResponse.h>
#include <ElementalLoad.h>
#include <ID.h>
#include <Domain.h>
#include <Node.h>
#include <Channel.h>
#include <Message.h>
#include <FEM_ObjectBroker.h>
#include <Renderer.h>
#include <G3Globals.h>
#include <ErrorHandler.h>
#include <NDMaterial.h>
#include <Parameter.h>

#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#define OPS_Export

static int num_SSPTethraUP = 0;

OPS_Export void *
OPS_SSPTethraUP(void)
{
    if (num_SSPTethraUP == 0) {
        num_SSPTethraUP++;
        opserr << "SSPTethraUP element - derived from SSPbrickUP: C.McGann, P.Arduino, P.Mackenzie-Helnwein, U.Washington\n";
    }

    Element *theElement = 0;

    int numRemainingInputArgs = OPS_GetNumRemainingInputArgs();
    if (numRemainingInputArgs < 13) {
        opserr << "Invalid #args, want: element SSPTethraUP eleTag? iNode? jNode? kNode? lNode? matTag? fBulk? fDen? k1? k2? k3? e? alpha? <b1? b2? b3?> <-lumped>\n";
        return 0;
    }

    int iData[6];
    double dData[10];

    dData[7] = 0.0;
    dData[8] = 0.0;
    dData[9] = 0.0;

    int numData = 6;
    if (OPS_GetIntInput(&numData, iData) != 0) {
        opserr << "WARNING invalid integer data: element SSPTethraUP " << iData[0] << endln;
        return 0;
    }

    int matID = iData[5];
    NDMaterial *theMaterial = OPS_getNDMaterial(matID);
    if (theMaterial == 0) {
        opserr << "WARNING element SSPTethraUP " << iData[0] << endln;
        opserr << " Material: " << matID << " not found\n";
        return 0;
    }

    numData = 7;
    if (OPS_GetDoubleInput(&numData, dData) != 0) {
        opserr << "WARNING invalid double data: element SSPTethraUP " << iData[0] << endln;
        return 0;
    }

    int massType = 0;
    int currentParam = 7;

    while (OPS_GetNumRemainingInputArgs() > 0) {
        int numArgs = OPS_GetNumRemainingInputArgs();

        numData = 1;
        double optionalVal;

        // Try to read a number first
        if (OPS_GetDoubleInput(&numData, &optionalVal) == 0) {
            if (currentParam < 10)
                dData[currentParam++] = optionalVal;
            continue;
        }

        // If the failed double read consumed an argument, give it back
        if (OPS_GetNumRemainingInputArgs() < numArgs)
            OPS_ResetCurrentInputArg(-1);

        const char *arg = OPS_GetString();
        if (arg == 0)
            break;

        if (strcmp(arg, "-lumped") == 0) {
            massType = 1;
        } else {
            if (currentParam < 10) {
                dData[currentParam++] = atof(arg);
            }
        }
    }

    theElement = new SSPTethraUP(iData[0], iData[1], iData[2], iData[3], iData[4],
                                 *theMaterial,
                                 dData[0], dData[1], dData[2], dData[3], dData[4],
                                 dData[5], dData[6],
                                 dData[7], dData[8], dData[9], massType);

    if (theElement == 0) {
        opserr << "WARNING could not create element of type SSPTethraUP\n";
        return 0;
    }

    return theElement;
}

SSPTethraUP::SSPTethraUP(int tag, int Nd1, int Nd2, int Nd3, int Nd4,
                         NDMaterial &theMat, double Kf, double Rf,
                         double k1, double k2, double k3,
                         double eVoid, double alpha,
                         double b1, double b2, double b3, int theMassType)
    : Element(tag, ELE_TAG_SSPTethraUP),
      theMaterial(0),
      mExternalNodes(STUP_NUM_NODE),
      mTangentStiffness(STUP_NUM_DOF, STUP_NUM_DOF),
      mInternalForces(STUP_NUM_DOF),
      Q(STUP_NUM_DOF),
      mDamp(STUP_NUM_DOF, STUP_NUM_DOF),
      mMass(STUP_NUM_DOF, STUP_NUM_DOF),
      fBulk(Kf),
      fDens(Rf),
      mAlpha(alpha),
      applyLoad(0),
      mVol(0),
      mPorosity(0),
      massType(theMassType),
      KfEff(Kf),
      alphaEff(alpha),
      hChar(0.0),
      mFluidS(STUP_NUM_NODE, STUP_NUM_NODE),
      Bnot(6, 3*STUP_NUM_NODE),
      mSolidK(3*STUP_NUM_NODE, 3*STUP_NUM_NODE),
      mSolidM(3*STUP_NUM_NODE, 3*STUP_NUM_NODE),
      mPerm(STUP_NUM_NODE, STUP_NUM_NODE),
      mPressStab(STUP_NUM_NODE, STUP_NUM_NODE),
      dNmod(STUP_NUM_NODE, STUP_NUM_DIM)
{
    for (int i = 0; i < STUP_NUM_NODE; i++)
        theNodes[i] = 0;

    mExternalNodes(0) = Nd1;
    mExternalNodes(1) = Nd2;
    mExternalNodes(2) = Nd3;
    mExternalNodes(3) = Nd4;

    fBulk  = Kf;
    fDens  = Rf;
    mAlpha = alpha;

    b[0] = b1;
    b[1] = b2;
    b[2] = b3;

    appliedB[0] = 0.0;
    appliedB[1] = 0.0;
    appliedB[2] = 0.0;

    perm[0] = k1;
    perm[1] = k2;
    perm[2] = k3;

    double e = eVoid;
    if (e <= -1.0)
        e = 0.0;

    mPorosity = e/(1.0 + e);

    if (mPorosity < 0.0)
        mPorosity = 0.0;
    if (mPorosity > 1.0)
        mPorosity = 1.0;

    NDMaterial *theMatCopy = theMat.getCopy("ThreeDimensional");
    if (theMatCopy != 0) {
        theMaterial = (NDMaterial *)theMatCopy;
    } else {
        opserr << "SSPTethraUP::SSPTethraUP - failed to get copy of material model\n";
    }

    if (theMaterial == 0) {
        opserr << "SSPTethraUP::SSPTethraUP - failed to allocate material model pointer\n";
        exit(-1);
    }
}

SSPTethraUP::SSPTethraUP()
    : Element(0, ELE_TAG_SSPTethraUP),
      theMaterial(0),
      mExternalNodes(STUP_NUM_NODE),
      mTangentStiffness(STUP_NUM_DOF, STUP_NUM_DOF),
      mInternalForces(STUP_NUM_DOF),
      Q(STUP_NUM_DOF),
      mDamp(STUP_NUM_DOF, STUP_NUM_DOF),
      mMass(STUP_NUM_DOF, STUP_NUM_DOF),
      fBulk(0),
      fDens(0),
      mAlpha(0),
      applyLoad(0),
      mVol(0),
      mPorosity(0),
      massType(0),
      KfEff(0),
      alphaEff(0),
      hChar(0.0),
      mFluidS(STUP_NUM_NODE, STUP_NUM_NODE),
      Bnot(6, 3*STUP_NUM_NODE),
      mSolidK(3*STUP_NUM_NODE, 3*STUP_NUM_NODE),
      mSolidM(3*STUP_NUM_NODE, 3*STUP_NUM_NODE),
      mPerm(STUP_NUM_NODE, STUP_NUM_NODE),
      mPressStab(STUP_NUM_NODE, STUP_NUM_NODE),
      dNmod(STUP_NUM_NODE, STUP_NUM_DIM)
{
    for (int i = 0; i < STUP_NUM_NODE; i++)
        theNodes[i] = 0;

    b[0] = 0.0;
    b[1] = 0.0;
    b[2] = 0.0;

    perm[0] = 0.0;
    perm[1] = 0.0;
    perm[2] = 0.0;

    appliedB[0] = 0.0;
    appliedB[1] = 0.0;
    appliedB[2] = 0.0;
}

SSPTethraUP::~SSPTethraUP()
{
    if (theMaterial != 0) {
        delete theMaterial;
    }
}

int
SSPTethraUP::getNumExternalNodes(void) const
{
    return STUP_NUM_NODE;
}

const ID &
SSPTethraUP::getExternalNodes(void)
{
    return mExternalNodes;
}

Node **
SSPTethraUP::getNodePtrs(void)
{
    return theNodes;
}

int
SSPTethraUP::getNumDOF(void)
{
    return STUP_NUM_DOF;
}

void
SSPTethraUP::setDomain(Domain *theDomain)
{
    if (theDomain == 0) {
        for (int i = 0; i < STUP_NUM_NODE; i++)
            theNodes[i] = 0;
        this->DomainComponent::setDomain(theDomain);
        return;
    }

    theNodes[0] = theDomain->getNode(mExternalNodes(0));
    theNodes[1] = theDomain->getNode(mExternalNodes(1));
    theNodes[2] = theDomain->getNode(mExternalNodes(2));
    theNodes[3] = theDomain->getNode(mExternalNodes(3));

    for (int i = 0; i < STUP_NUM_NODE; i++) {
        if (theNodes[i] == 0) {
            return;
        }

        int requiredDOF = STUP_NUM_DOF/STUP_NUM_NODE;
        if (theNodes[i]->getNumberDOF() != requiredDOF) {
            opserr << "SSPTethraUP::setDomain - element " << this->getTag()
                   << ": node " << mExternalNodes(i)
                   << " has " << theNodes[i]->getNumberDOF()
                   << " dof, 4 required (use model BasicBuilder -ndm 3 -ndf 4)\n";
            mVol = 0.0;
            return;
        }
    }

    setTetrahedronData();

    if (mVol <= 0.0) {
        opserr << "SSPTethraUP::setDomain - element " << this->getTag()
               << " has zero or invalid volume\n";
        return;
    }

    computeRegularization();
    GetSolidMass();
    GetPermeabilityMatrix();

    this->DomainComponent::setDomain(theDomain);
}

int
SSPTethraUP::commitState(void)
{
    int retVal = 0;

    if ((retVal = this->Element::commitState()) != 0) {
        opserr << "SSPTethraUP::commitState() - failed in base class\n";
    }

    retVal = theMaterial->commitState();
    return retVal;
}

int
SSPTethraUP::revertToLastCommit(void)
{
    return theMaterial->revertToLastCommit();
}

int
SSPTethraUP::revertToStart(void)
{
    return theMaterial->revertToStart();
}

int
SSPTethraUP::update(void)
{
    if (mVol <= 0.0)
        return 0;

    const Vector &mDisp_1 = theNodes[0]->getTrialDisp();
    const Vector &mDisp_2 = theNodes[1]->getTrialDisp();
    const Vector &mDisp_3 = theNodes[2]->getTrialDisp();
    const Vector &mDisp_4 = theNodes[3]->getTrialDisp();

    Vector u(12);

    u(0)  = mDisp_1(0);
    u(1)  = mDisp_1(1);
    u(2)  = mDisp_1(2);

    u(3)  = mDisp_2(0);
    u(4)  = mDisp_2(1);
    u(5)  = mDisp_2(2);

    u(6)  = mDisp_3(0);
    u(7)  = mDisp_3(1);
    u(8)  = mDisp_3(2);

    u(9)  = mDisp_4(0);
    u(10) = mDisp_4(1);
    u(11) = mDisp_4(2);

    Vector strain(6);
    strain = Bnot*u;

    theMaterial->setTrialStrain(strain);

    return 0;
}

const Matrix &
SSPTethraUP::getTangentStiff(void)
{
    mTangentStiffness.Zero();

    if (mVol <= 0.0)
        return mTangentStiffness;

    GetSolidStiffness();

    for (int i = 0; i < 4; i++) {
        int I    = 3*i;
        int Ip1  = 3*i+1;
        int Ip2  = 3*i+2;

        int II   = 4*i;
        int IIp1 = 4*i+1;
        int IIp2 = 4*i+2;

        for (int j = 0; j < 4; j++) {
            int J    = 3*j;
            int Jp1  = 3*j+1;
            int Jp2  = 3*j+2;

            int JJ   = 4*j;
            int JJp1 = 4*j+1;
            int JJp2 = 4*j+2;

            mTangentStiffness(II,JJ)     = mSolidK(I,J);
            mTangentStiffness(IIp1,JJ)   = mSolidK(Ip1,J);
            mTangentStiffness(IIp2,JJ)   = mSolidK(Ip2,J);

            mTangentStiffness(II,JJp1)   = mSolidK(I,Jp1);
            mTangentStiffness(IIp1,JJp1) = mSolidK(Ip1,Jp1);
            mTangentStiffness(IIp2,JJp1) = mSolidK(Ip2,Jp1);

            mTangentStiffness(II,JJp2)   = mSolidK(I,Jp2);
            mTangentStiffness(IIp1,JJp2) = mSolidK(Ip1,Jp2);
            mTangentStiffness(IIp2,JJp2) = mSolidK(Ip2,Jp2);
        }
    }

    return mTangentStiffness;
}

const Matrix &
SSPTethraUP::getInitialStiff(void)
{
    mTangentStiffness.Zero();

    if (mVol <= 0.0)
        return mTangentStiffness;

    Matrix K0(12,12);
    K0.addMatrixTripleProduct(0.0, Bnot, theMaterial->getInitialTangent(), mVol);

    for (int i = 0; i < 4; i++) {
        for (int j = 0; j < 4; j++) {
            for (int a = 0; a < 3; a++) {
                for (int c = 0; c < 3; c++) {
                    mTangentStiffness(4*i+a, 4*j+c) = K0(3*i+a, 3*j+c);
                }
            }
        }
    }

    return mTangentStiffness;
}

const Matrix &
SSPTethraUP::getDamp(void)
{
    mDamp.Zero();

    if (mVol <= 0.0)
        return mDamp;

    Matrix dampC(12,12);
    dampC.Zero();

    GetSolidStiffness();

    if (alphaM != 0.0) {
        dampC.addMatrix(0.0, mSolidM, alphaM);
    }

    if (betaK != 0.0) {
        dampC.addMatrix(1.0, mSolidK, betaK);
    }

    if (betaK0 != 0.0) {
        Matrix K0(12,12);
        K0.addMatrixTripleProduct(0.0, Bnot, theMaterial->getInitialTangent(), mVol);
        dampC.addMatrix(1.0, K0, betaK0);
    }

    if (betaKc != 0.0) {
        dampC.addMatrix(1.0, mSolidK, betaKc);
    }

    // Coupling matrix Q
    Matrix couple(12,4);
    Matrix INp(6,4);

    INp.Zero();

    for (int i = 0; i < 4; i++) {
        INp(0,i) = 0.25;
        INp(1,i) = 0.25;
        INp(2,i) = 0.25;
    }

    couple.Zero();
    couple.addMatrixTransposeProduct(0.0, Bnot, INp, mVol);

    for (int i = 0; i < 4; i++) {
        int I    = 3*i;
        int Ip1  = 3*i+1;
        int Ip2  = 3*i+2;

        int II   = 4*i;
        int IIp1 = 4*i+1;
        int IIp2 = 4*i+2;
        int IIp3 = 4*i+3;

        for (int j = 0; j < 4; j++) {
            int J    = 3*j;
            int Jp1  = 3*j+1;
            int Jp2  = 3*j+2;

            int JJ   = 4*j;
            int JJp1 = 4*j+1;
            int JJp2 = 4*j+2;
            int JJp3 = 4*j+3;

            // Solid damping
            mDamp(II,JJ)     = dampC(I,J);
            mDamp(IIp1,JJ)   = dampC(Ip1,J);
            mDamp(IIp2,JJ)   = dampC(Ip2,J);

            mDamp(II,JJp1)   = dampC(I,Jp1);
            mDamp(IIp1,JJp1) = dampC(Ip1,Jp1);
            mDamp(IIp2,JJp1) = dampC(Ip2,Jp1);

            mDamp(II,JJp2)   = dampC(I,Jp2);
            mDamp(IIp1,JJp2) = dampC(Ip1,Jp2);
            mDamp(IIp2,JJp2) = dampC(Ip2,Jp2);

            // Fluid-solid coupling
            mDamp(JJp3,II)   = -couple(I,j);
            mDamp(JJp3,IIp1) = -couple(Ip1,j);
            mDamp(JJp3,IIp2) = -couple(Ip2,j);

            mDamp(II,JJp3)   = -couple(I,j);
            mDamp(IIp1,JJp3) = -couple(Ip1,j);
            mDamp(IIp2,JJp3) = -couple(Ip2,j);

            // Permeability
            mDamp(IIp3,JJp3) = -mPerm(i,j);
        }
    }

    return mDamp;
}

const Matrix &
SSPTethraUP::getMass(void)
{
    mMass.Zero();

    if (mVol <= 0.0)
        return mMass;

    for (int i = 0; i < 4; i++) {
        int I    = 3*i;
        int Ip1  = 3*i+1;
        int Ip2  = 3*i+2;

        int II   = 4*i;
        int IIp1 = 4*i+1;
        int IIp2 = 4*i+2;
        int IIp3 = 4*i+3;

        for (int j = 0; j < 4; j++) {
            int J    = 3*j;
            int Jp1  = 3*j+1;
            int Jp2  = 3*j+2;

            int JJ   = 4*j;
            int JJp1 = 4*j+1;
            int JJp2 = 4*j+2;
            int JJp3 = 4*j+3;

            // Solid mass
            mMass(II,JJ)     = mSolidM(I,J);
            mMass(IIp1,JJ)   = mSolidM(Ip1,J);
            mMass(IIp2,JJ)   = mSolidM(Ip2,J);

            mMass(II,JJp1)   = mSolidM(I,Jp1);
            mMass(IIp1,JJp1) = mSolidM(Ip1,Jp1);
            mMass(IIp2,JJp1) = mSolidM(Ip2,Jp1);

            mMass(II,JJp2)   = mSolidM(I,Jp2);
            mMass(IIp1,JJp2) = mSolidM(Ip1,Jp2);
            mMass(IIp2,JJp2) = mSolidM(Ip2,Jp2);

            // Pressure block: -(S + Ps)
            if (massType == 1) {
                if (i == j) {
                    double rowSum = 0.0;
                    for (int k = 0; k < 4; k++) {
                        rowSum += mFluidS(i,k) + mPressStab(i,k);
                    }
                    mMass(IIp3,JJp3) = -rowSum;
                } else {
                    mMass(IIp3,JJp3) = 0.0;
                }
            } else {
                mMass(IIp3,JJp3) = -(mFluidS(i,j) + mPressStab(i,j));
            }
        }
    }

    return mMass;
}

void
SSPTethraUP::zeroLoad(void)
{
    applyLoad = 0;

    appliedB[0] = 0.0;
    appliedB[1] = 0.0;
    appliedB[2] = 0.0;

    Q.Zero();
}

int
SSPTethraUP::addLoad(ElementalLoad *theLoad, double loadFactor)
{
    int type;
    const Vector &data = theLoad->getData(type, loadFactor);

    if (type == LOAD_TAG_SelfWeight) {
        applyLoad = 1;

        appliedB[0] += loadFactor*data(0)*b[0];
        appliedB[1] += loadFactor*data(1)*b[1];
        appliedB[2] += loadFactor*data(2)*b[2];

        return 0;
    } else {
        opserr << "SSPTethraUP::addLoad - load type unknown for ele with tag: "
               << this->getTag() << endln;
        return -1;
    }

    return -1;
}

int
SSPTethraUP::addInertiaLoadToUnbalance(const Vector &accel)
{
    double density = theMaterial->getRho();

    if (density == 0.0)
        return 0;

    const Vector &Raccel1 = theNodes[0]->getRV(accel);
    const Vector &Raccel2 = theNodes[1]->getRV(accel);
    const Vector &Raccel3 = theNodes[2]->getRV(accel);
    const Vector &Raccel4 = theNodes[3]->getRV(accel);

    Vector ra(16);

    ra(0)  = Raccel1(0);
    ra(1)  = Raccel1(1);
    ra(2)  = Raccel1(2);
    ra(3)  = 0.0;

    ra(4)  = Raccel2(0);
    ra(5)  = Raccel2(1);
    ra(6)  = Raccel2(2);
    ra(7)  = 0.0;

    ra(8)  = Raccel3(0);
    ra(9)  = Raccel3(1);
    ra(10) = Raccel3(2);
    ra(11) = 0.0;

    ra(12) = Raccel4(0);
    ra(13) = Raccel4(1);
    ra(14) = Raccel4(2);
    ra(15) = 0.0;

    this->getMass();

    Vector MRa(16);
    MRa.addMatrixVector(0.0, mMass, ra, 1.0);

    Q.addVector(1.0, MRa, -1.0);

    return 0;
}

const Vector &
SSPTethraUP::getResistingForce(void)
{
    Vector f1(12);
    Vector f2(4);

    f1.Zero();
    f2.Zero();

    if (mVol <= 0.0) {
        mInternalForces.Zero();
        mInternalForces.addVector(1.0, Q, -1.0);
        return mInternalForces;
    }

    Vector mStress = theMaterial->getStress();

    // Effective stress internal force
    f1.addMatrixTransposeVector(1.0, Bnot, mStress, mVol);

    // Body forces on solid skeleton
    double density = theMaterial->getRho();
    double dV4 = 0.25*mVol;

    if (density != 0.0) {
        if (applyLoad == 0) {
            for (int i = 0; i < 4; i++) {
                f1(3*i)   -= density*b[0]*dV4;
                f1(3*i+1) -= density*b[1]*dV4;
                f1(3*i+2) -= density*b[2]*dV4;
            }
        } else {
            for (int i = 0; i < 4; i++) {
                f1(3*i)   -= density*appliedB[0]*dV4;
                f1(3*i+1) -= density*appliedB[1]*dV4;
                f1(3*i+2) -= density*appliedB[2]*dV4;
            }
        }
    }

    // Fluid body forces
    if (fDens != 0.0) {
        Matrix k(3,3);
        Vector body(3);

        k.Zero();

        k(0,0) = (perm[0] > 0.0) ? perm[0] : 0.0;
        k(1,1) = (perm[1] > 0.0) ? perm[1] : 0.0;
        k(2,2) = (perm[2] > 0.0) ? perm[2] : 0.0;

        if (applyLoad == 0) {
            body(0) = b[0];
            body(1) = b[1];
            body(2) = b[2];
        } else {
            body(0) = appliedB[0];
            body(1) = appliedB[1];
            body(2) = appliedB[2];
        }

        Vector kBody = k*body;
        f2 = dNmod*kBody;
        f2 *= mVol*fDens;
    }

    // Assemble 16-dof internal force vector
    mInternalForces(0)  = f1(0);
    mInternalForces(1)  = f1(1);
    mInternalForces(2)  = f1(2);
    mInternalForces(3)  = f2(0);

    mInternalForces(4)  = f1(3);
    mInternalForces(5)  = f1(4);
    mInternalForces(6)  = f1(5);
    mInternalForces(7)  = f2(1);

    mInternalForces(8)  = f1(6);
    mInternalForces(9)  = f1(7);
    mInternalForces(10) = f1(8);
    mInternalForces(11) = f2(2);

    mInternalForces(12) = f1(9);
    mInternalForces(13) = f1(10);
    mInternalForces(14) = f1(11);
    mInternalForces(15) = f2(3);

    mInternalForces.addVector(1.0, Q, -1.0);

    return mInternalForces;
}

const Vector &
SSPTethraUP::getResistingForceIncInertia()
{
    this->getResistingForce();
    this->getMass();
    this->getDamp();

    const Vector &accel1 = theNodes[0]->getTrialAccel();
    const Vector &accel2 = theNodes[1]->getTrialAccel();
    const Vector &accel3 = theNodes[2]->getTrialAccel();
    const Vector &accel4 = theNodes[3]->getTrialAccel();

    Vector a(16);

    a(0)  = accel1(0);
    a(1)  = accel1(1);
    a(2)  = accel1(2);
    a(3)  = accel1(3);

    a(4)  = accel2(0);
    a(5)  = accel2(1);
    a(6)  = accel2(2);
    a(7)  = accel2(3);

    a(8)  = accel3(0);
    a(9)  = accel3(1);
    a(10) = accel3(2);
    a(11) = accel3(3);

    a(12) = accel4(0);
    a(13) = accel4(1);
    a(14) = accel4(2);
    a(15) = accel4(3);

    mInternalForces.addMatrixVector(1.0, mMass, a, 1.0);

    const Vector &vel1 = theNodes[0]->getTrialVel();
    const Vector &vel2 = theNodes[1]->getTrialVel();
    const Vector &vel3 = theNodes[2]->getTrialVel();
    const Vector &vel4 = theNodes[3]->getTrialVel();

    Vector v(16);

    v(0)  = vel1(0);
    v(1)  = vel1(1);
    v(2)  = vel1(2);
    v(3)  = vel1(3);

    v(4)  = vel2(0);
    v(5)  = vel2(1);
    v(6)  = vel2(2);
    v(7)  = vel2(3);

    v(8)  = vel3(0);
    v(9)  = vel3(1);
    v(10) = vel3(2);
    v(11) = vel3(3);

    v(12) = vel4(0);
    v(13) = vel4(1);
    v(14) = vel4(2);
    v(15) = vel4(3);

    mInternalForces.addMatrixVector(1.0, mDamp, v, 1.0);

    return mInternalForces;
}

int
SSPTethraUP::sendSelf(int commitTag, Channel &theChannel)
{
    int res = 0;
    int dataTag = this->getDbTag();

    static Vector data(14);

    data(0)  = this->getTag();
    data(1)  = fBulk;
    data(2)  = fDens;
    data(3)  = perm[0];
    data(4)  = perm[1];
    data(5)  = perm[2];
    data(6)  = mPorosity;
    data(7)  = mAlpha;
    data(8)  = b[0];
    data(9)  = b[1];
    data(10) = b[2];
    data(11) = theMaterial->getClassTag();

    int matDbTag = theMaterial->getDbTag();

    if (matDbTag == 0) {
        matDbTag = theChannel.getDbTag();
        if (matDbTag != 0)
            theMaterial->setDbTag(matDbTag);
    }

    data(12) = matDbTag;
    data(13) = massType;

    res += theChannel.sendVector(dataTag, commitTag, data);
    if (res < 0) {
        opserr << "WARNING SSPTethraUP::sendSelf() - " << this->getTag()
               << " failed to send Vector\n";
        return res;
    }

    res += theChannel.sendID(dataTag, commitTag, mExternalNodes);
    if (res < 0) {
        opserr << "WARNING SSPTethraUP::sendSelf() - " << this->getTag()
               << " failed to send ID\n";
        return res;
    }

    res = theMaterial->sendSelf(commitTag, theChannel);
    if (res < 0) {
        opserr << "WARNING SSPTethraUP::sendSelf() - " << this->getTag()
               << " failed to send its Material\n";
        return -3;
    }

    return 0;
}

int
SSPTethraUP::recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker)
{
    int res = 0;
    int dataTag = this->getDbTag();

    static Vector data(14);

    res += theChannel.recvVector(dataTag, commitTag, data);
    if (res < 0) {
        opserr << "WARNING SSPTethraUP::recvSelf() - failed to receive Vector\n";
        return res;
    }

    this->setTag((int)data(0));

    fBulk      = data(1);
    fDens      = data(2);
    perm[0]    = data(3);
    perm[1]    = data(4);
    perm[2]    = data(5);
    mPorosity  = data(6);
    mAlpha     = data(7);
    b[0]       = data(8);
    b[1]       = data(9);
    b[2]       = data(10);
    massType   = (int)data(13);

    applyLoad = 0;
    appliedB[0] = 0.0;
    appliedB[1] = 0.0;
    appliedB[2] = 0.0;

    KfEff = fBulk;
    alphaEff = mAlpha;

    res += theChannel.recvID(dataTag, commitTag, mExternalNodes);
    if (res < 0) {
        opserr << "WARNING SSPTethraUP::recvSelf() - " << this->getTag()
               << " failed to receive ID\n";
        return res;
    }

    int matClass = (int)data(11);
    int matDb    = (int)data(12);

    if ((theMaterial == 0) || (theMaterial->getClassTag() != matClass)) {
        if (theMaterial != 0)
            delete theMaterial;

        NDMaterial *theMatCopy = theBroker.getNewNDMaterial(matClass);
        theMaterial = (NDMaterial *)theMatCopy;

        if (theMaterial == 0) {
            opserr << "WARNING SSPTethraUP::recvSelf() - " << this->getTag()
                   << " failed to get a blank Material of type " << matClass << endln;
            return -3;
        }
    }

    theMaterial->setDbTag(matDb);

    res = theMaterial->recvSelf(commitTag, theChannel, theBroker);
    if (res < 0) {
        opserr << "WARNING SSPTethraUP::recvSelf() - " << this->getTag()
               << " failed to receive its Material\n";
        return -3;
    }

    return 0;
}

int
SSPTethraUP::displaySelf(Renderer &theViewer, int displayMode, float fact,
                         const char **modes, int numMode)
{
    static Vector v1(3);
    static Vector v2(3);
    static Vector v3(3);
    static Vector v4(3);

    theNodes[0]->getDisplayCrds(v1, fact, displayMode);
    theNodes[1]->getDisplayCrds(v2, fact, displayMode);
    theNodes[2]->getDisplayCrds(v3, fact, displayMode);
    theNodes[3]->getDisplayCrds(v4, fact, displayMode);

    static Vector values(3);
    values(0) = 0.0;
    values(1) = 0.0;
    values(2) = 0.0;

    int res = 0;
    static Matrix coords(3,3);

    // Face 1: 1-3-2
    for (int i = 0; i < 3; i++) {
        coords(0,i) = v1(i);
        coords(1,i) = v3(i);
        coords(2,i) = v2(i);
    }
    res += theViewer.drawPolygon(coords, values, this->getTag());

    // Face 2: 1-2-4
    for (int i = 0; i < 3; i++) {
        coords(0,i) = v1(i);
        coords(1,i) = v2(i);
        coords(2,i) = v4(i);
    }
    res += theViewer.drawPolygon(coords, values, this->getTag());

    // Face 3: 1-4-3
    for (int i = 0; i < 3; i++) {
        coords(0,i) = v1(i);
        coords(1,i) = v4(i);
        coords(2,i) = v3(i);
    }
    res += theViewer.drawPolygon(coords, values, this->getTag());

    // Face 4: 2-3-4
    for (int i = 0; i < 3; i++) {
        coords(0,i) = v2(i);
        coords(1,i) = v3(i);
        coords(2,i) = v4(i);
    }
    res += theViewer.drawPolygon(coords, values, this->getTag());

    return res;
}

void
SSPTethraUP::Print(OPS_Stream &s, int flag)
{
    if (flag == OPS_PRINT_CURRENTSTATE) {
        s << "SSPTethraUP, element id: " << this->getTag() << endln;
        s << "   Connected external nodes: ";
        for (int i = 0; i < STUP_NUM_NODE; i++) {
            s << mExternalNodes(i) << " ";
        }
        s << endln;
    }

    if (flag == OPS_PRINT_PRINTMODEL_JSON) {
        s << "\t\t\t{";
        s << "\"name\": " << this->getTag() << ", ";
        s << "\"type\": \"SSPTethraUP\", ";
        s << "\"nodes\": [" << mExternalNodes(0) << ", ";
        for (int i = 1; i < 3; i++)
            s << mExternalNodes(i) << ", ";
        s << mExternalNodes(3) << "], ";
        s << "\"bodyForces\": [" << b[0] << ", " << b[1] << ", " << b[2] << "], ";
        s << "\"material\": \"" << theMaterial->getTag() << "\"}";
    }
}

Response*
SSPTethraUP::setResponse(const char **argv, int argc, OPS_Stream &eleInfo)
{
    if (argc < 1)
        return 0;

    if (strcmp(argv[0], "stress3D6") == 0) {
        return new ElementResponse(this, 1, Vector(6));
    }

    if (strcmp(argv[0], "strain3D6") == 0) {
        return new ElementResponse(this, 2, Vector(6));
    }

    if (strcmp(argv[0], "mass") == 0) {
        return new ElementResponse(this, 3, Matrix(STUP_NUM_DOF, STUP_NUM_DOF));
    }

    if ((strcmp(argv[0], "pressure") == 0) ||
        (strcmp(argv[0], "porePressure") == 0)) {
        return new ElementResponse(this, 4, Vector(1));
    }

    return theMaterial->setResponse(argv, argc, eleInfo);
}

int
SSPTethraUP::getResponse(int responseID, Information &eleInfo)
{
    if (responseID == 1) {
        const Vector &mStress = theMaterial->getStress();
        return eleInfo.setVector(mStress);
    }

    if (responseID == 2) {
        const Vector &mStrain = theMaterial->getStrain();
        return eleInfo.setVector(mStrain);
    }

    if (responseID == 3) {
        return eleInfo.setMatrix(this->getMass());
    }

    if (responseID == 4) {
        // Physical pore pressure is velocity of dof 4
        double eP = 0.0;

        if (mVol > 0.0) {
            for (int i = 0; i < 4; i++) {
                const Vector &vel = theNodes[i]->getTrialVel();
                eP += 0.25*vel(3);
            }
        }

        Vector pressure(1);
        pressure(0) = eP;

        return eleInfo.setVector(pressure);
    }

    return theMaterial->getResponse(responseID, eleInfo);
}

int
SSPTethraUP::setParameter(const char **argv, int argc, Parameter &param)
{
    if (argc < 1)
        return -1;

    if (strcmp(argv[0], "xPerm") == 0) {
        return param.addObject(3, this);
    } else if (strcmp(argv[0], "yPerm") == 0) {
        return param.addObject(4, this);
    } else if (strcmp(argv[0], "zPerm") == 0) {
        return param.addObject(6, this);
    } else {
        return theMaterial->setParameter(argv, argc, param);
    }
}

int
SSPTethraUP::updateParameter(int parameterID, Information &info)
{
    if (parameterID == 3) {
        perm[0] = info.theDouble;
        GetPermeabilityMatrix();
        return 0;
    } else if (parameterID == 4) {
        perm[1] = info.theDouble;
        GetPermeabilityMatrix();
        return 0;
    } else if (parameterID == 6) {
        perm[2] = info.theDouble;
        GetPermeabilityMatrix();
        return 0;
    } else {
        return theMaterial->updateParameter(parameterID, info);
    }
}

void
SSPTethraUP::setTetrahedronData(void)
{
    Matrix dNloc(4,3);

    // Parent tetrahedron:
    // node 1: (0,0,0)
    // node 2: (1,0,0)
    // node 3: (0,1,0)
    // node 4: (0,0,1)

    dNloc(0,0) = -1.0;  dNloc(0,1) = -1.0;  dNloc(0,2) = -1.0;
    dNloc(1,0) =  1.0;  dNloc(1,1) =  0.0;  dNloc(1,2) =  0.0;
    dNloc(2,0) =  0.0;  dNloc(2,1) =  1.0;  dNloc(2,2) =  0.0;
    dNloc(3,0) =  0.0;  dNloc(3,1) =  0.0;  dNloc(3,2) =  1.0;

    const Vector &crd1 = theNodes[0]->getCrds();
    const Vector &crd2 = theNodes[1]->getCrds();
    const Vector &crd3 = theNodes[2]->getCrds();
    const Vector &crd4 = theNodes[3]->getCrds();

    Matrix coords(3,4);

    for (int a = 0; a < 3; a++) {
        coords(a,0) = crd1(a);
        coords(a,1) = crd2(a);
        coords(a,2) = crd3(a);
        coords(a,3) = crd4(a);
    }

    Matrix Jmat(3,3);
    Jmat.addMatrixProduct(0.0, coords, dNloc, 1.0);

    double detJ =
        Jmat(0,0)*Jmat(1,1)*Jmat(2,2)
      + Jmat(0,1)*Jmat(1,2)*Jmat(2,0)
      + Jmat(0,2)*Jmat(1,0)*Jmat(2,1)
      - Jmat(0,2)*Jmat(1,1)*Jmat(2,0)
      - Jmat(0,1)*Jmat(1,0)*Jmat(2,2)
      - Jmat(0,0)*Jmat(1,2)*Jmat(2,1);

    if (fabs(detJ) < 1.0e-20) {
        mVol = 0.0;
        return;
    }

    if (detJ < 0.0) {
        opserr << "SSPTethraUP::setTetrahedronData - element " << this->getTag()
               << " has negative Jacobian. Using |detJ|. Check node ordering.\n";
    }

    mVol = fabs(detJ)/6.0;

    Matrix Jinv(3,3);
    Jmat.Invert(Jinv);

    Matrix dN(4,3);
    dN.addMatrixProduct(0.0, dNloc, Jinv, 1.0);

    for (int i = 0; i < 4; i++) {
        dNmod(i,0) = dN(i,0);
        dNmod(i,1) = dN(i,1);
        dNmod(i,2) = dN(i,2);
    }

    Bnot.Zero();

    for (int i = 0; i < 4; i++) {
        int I = 3*i;

        Bnot(0,I)   = dNmod(i,0);
        Bnot(1,I+1) = dNmod(i,1);
        Bnot(2,I+2) = dNmod(i,2);

        Bnot(3,I)   = dNmod(i,1);
        Bnot(3,I+1) = dNmod(i,0);

        Bnot(4,I+1) = dNmod(i,2);
        Bnot(4,I+2) = dNmod(i,1);

        Bnot(5,I)   = dNmod(i,2);
        Bnot(5,I+2) = dNmod(i,0);
    }
}

double
SSPTethraUP::getReferenceModulus(void)
{
    double Kref = 0.0;

    const Matrix &C0 = theMaterial->getInitialTangent();

    for (int i = 0; i < 3; i++) {
        double c = fabs(C0(i,i));
        if (c > Kref)
            Kref = c;
    }

    for (int i = 3; i < 6; i++) {
        double c = 2.0*fabs(C0(i,i));
        if (c > Kref)
            Kref = c;
    }

    if (Kref < 1.0e-12) {
        const Matrix &C = theMaterial->getTangent();

        for (int i = 0; i < 3; i++) {
            double c = fabs(C(i,i));
            if (c > Kref)
                Kref = c;
        }

        for (int i = 3; i < 6; i++) {
            double c = 2.0*fabs(C(i,i));
            if (c > Kref)
                Kref = c;
        }
    }

    if (Kref < 1.0e-12) {
        if (fabs(fBulk) > 1.0e-12)
            Kref = fabs(fBulk);
        else
            Kref = 1.0;
    }

    return Kref;
}

void
SSPTethraUP::computeRegularization(void)
{
    mFluidS.Zero();

    if (mVol <= 0.0)
        return;

    hChar = pow(6.0*mVol, 1.0/3.0);
    if (hChar <= 0.0)
        hChar = 1.0;

    double Kref = getReferenceModulus();

    // Cap fluid bulk modulus to avoid exactly incompressible fluid.
    // This is a numerical regularization. Tune if necessary.
    double KfCap = 1000.0*Kref;
    if (KfCap < 1.0e6)
        KfCap = 1.0e6;

    KfEff = fBulk;

    if (KfEff <= 0.0 || KfEff > KfCap)
        KfEff = KfCap;

    // Pressure stabilization parameter.
    // If user gives alpha > 0, use it.
    // Otherwise compute an automatic stabilization parameter.
    if (mAlpha > 0.0) {
        alphaEff = mAlpha;
    } else {
        double Kden = KfEff;
        if (Kref > Kden)
            Kden = Kref;

        if (Kden <= 0.0)
            Kden = 1.0;

        alphaEff = hChar*hChar/(10.0*Kden);
    }

    // Ensure a strictly positive stabilization
    if (alphaEff <= 0.0) {
        alphaEff = hChar*hChar/(1.0e8*Kref);
        if (alphaEff <= 0.0)
            alphaEff = 1.0e-12;
    }

    // Fluid compressibility matrix S.
    // Exact consistent matrix for linear tetrahedron:
    // int N_i N_i dV = V/10
    // int N_i N_j dV = V/20
    double n = mPorosity;
    if (n < 0.0)
        n = 0.0;
    if (n > 1.0)
        n = 1.0;

    double coeff = n/KfEff;

    double diag = coeff*mVol/10.0;
    double off  = coeff*mVol/20.0;

    // Small artificial compressibility floor.
    // This prevents the pressure mass block from becoming too small.
    double Sfloor = 1.0e-5*mVol/Kref;

    if (Sfloor < 1.0e-14*mVol)
        Sfloor = 1.0e-14*mVol;

    for (int i = 0; i < 4; i++) {
        for (int j = 0; j < 4; j++) {
            if (i == j)
                mFluidS(i,j) = diag + Sfloor;
            else
                mFluidS(i,j) = off;
        }
    }
}

void
SSPTethraUP::GetSolidStiffness(void)
{
    mSolidK.Zero();

    if (mVol <= 0.0)
        return;

    const Matrix &Cmat = theMaterial->getTangent();

    mSolidK.addMatrixTripleProduct(1.0, Bnot, Cmat, mVol);
}

void
SSPTethraUP::GetSolidMass(void)
{
    mSolidM.Zero();

    if (mVol <= 0.0)
        return;

    double density = theMaterial->getRho();

    if (density == 0.0)
        return;

    // Lumped solid mass for linear tetrahedron:
    // integral N_i dV = V/4
    double massTerm = density*mVol/4.0;

    for (int i = 0; i < 4; i++) {
        int I = 3*i;

        mSolidM(I,I)     += massTerm;
        mSolidM(I+1,I+1) += massTerm;
        mSolidM(I+2,I+2) += massTerm;
    }
}

void
SSPTethraUP::GetPermeabilityMatrix(void)
{
    mPerm.Zero();
    mPressStab.Zero();

    if (mVol <= 0.0)
        return;

    Matrix k(3,3);
    k.Zero();

    k(0,0) = (perm[0] > 0.0) ? perm[0] : 0.0;
    k(1,1) = (perm[1] > 0.0) ? perm[1] : 0.0;
    k(2,2) = (perm[2] > 0.0) ? perm[2] : 0.0;

    Matrix dNT(3,4);

    for (int i = 0; i < 4; i++) {
        dNT(0,i) = dNmod(i,0);
        dNT(1,i) = dNmod(i,1);
        dNT(2,i) = dNmod(i,2);
    }

    // H = int grad(Np)^T k grad(Np) dV
    mPerm.addMatrixTripleProduct(1.0, dNT, k, mVol);

    // Ps = alpha int grad(Np)^T grad(Np) dV
    mPressStab.addMatrixTransposeProduct(1.0, dNT, dNT, mVol*alphaEff);
}