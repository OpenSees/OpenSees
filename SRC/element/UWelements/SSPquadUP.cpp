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
                                                                       
//
// Adapted: Luis Miranda, LNEC, Andre Barbosa, OSU; July 2015 - added boundary tractions (see LM changes)
//
// Created: Chris McGann, UW, 05.2011
//
// Description: This file contains the implementation of the SSPquadUP class
//                Stabilized Single-Point Quad element with a u-p formulation 
//                for plane strain analysis of saturated porous media
//
// References:  Zienkiewicz, O.C. and Shiomi, T. (1984). "Dynamic behavior of 
//                saturated porous media; the generalized Biot formulation and 
//                its numerical solution." International Journal for Numerical 
//                Methods in Geomechanics, 8, 71-96.
//              McGann, C.R., Arduino, P., and Mackenzie-Helnwein, P. (2012) "Stabilized single-point
//                4-node quadrilateral element for dynamic analysis of fluid saturated porous media."
//                Acta Geotechnica, 7(4):297-311

#include "SSPquadUP.h"

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

#define OPS_Export

static int num_SSPquadUP = 0;

OPS_Export void *
OPS_SSPquadUP(void)
{
    if (num_SSPquadUP == 0) {
        num_SSPquadUP++;
        opserr<<"SSPquadUP element - Written: C.McGann, P.Arduino, P.Mackenzie-Helnwein, U.Washington\n";
    }

    // Pointer to an element that will be returned
    Element *theElement = 0;

    int numRemainingInputArgs = OPS_GetNumRemainingInputArgs();
    // LM change	
    if (numRemainingInputArgs < 13) {
        opserr << "Invalid #args, want: element SSPquadUP eleTag? iNode? jNode? kNode? lNode? matTag? t? fBulk? fDen? k1? k2? e? alpha? <b1? b2?> <Pup? Plow? Pleft? Pright?>?\n";
        return 0;
    }
        
    int iData[6];
    double dData[13];
    dData[7]  = 0.0;
    dData[8]  = 0.0;
	dData[9]  = 0.0;
	dData[10] = 0.0;
	dData[11] = 0.0;
	dData[12] = 0.0;
    // LM change

    int numData = 6;
    if (OPS_GetIntInput(&numData, iData) != 0) {
        opserr << "WARNING invalid integer data: element SSPquadUP " << iData[0] << endln;
        return 0;
    }

    numData = 7;
    if (OPS_GetDoubleInput(&numData, dData) != 0) {
        opserr << "WARNING invalid double data: element SSPquadUP " << iData[0] << endln;
        return 0;
    }

    int matID = iData[5];
    NDMaterial *theMaterial = OPS_getNDMaterial(matID);
    if (theMaterial == 0) {
        opserr << "WARNING element SSPquadUP " << iData[0] << endln;
        opserr << " Material: " << matID << "not found\n";
        return 0;
    }
         
    // LM change
    if (numRemainingInputArgs == 15) {
    	numData = 2;
    	if (OPS_GetDoubleInput(&numData, &dData[7]) != 0) {
      		opserr << "WARNING invalid optional data: element SSPquadUP " << iData[0] << endln;
	  		return 0;
    	}
  	} else if (numRemainingInputArgs == 19) {
        numData = 6;
        if (OPS_GetDoubleInput(&numData, &dData[7]) != 0) {
            opserr << "WARNING invalid optional data: element SSPquadUP " << iData[0] << endln;
            return 0;
        }
    }
		
    // parsing was successful, allocate the element
    theElement = new SSPquadUP(iData[0], iData[1], iData[2], iData[3],  iData[4], *theMaterial, 
                               dData[0], dData[1], dData[2], dData[3],  dData[4],  dData[5], dData[6], 
                               dData[7], dData[8], dData[9], dData[10], dData[11], dData[12]);
    // LM change

    if (theElement == 0) {
        opserr << "WARNING could not create element of type SSPquadUP\n";
        return 0;
    }

    return theElement;
}

// full constructor
SSPquadUP::SSPquadUP(int tag, int Nd1, int Nd2, int Nd3, int Nd4, NDMaterial &theMat, 
                              double thick, double Kf, double Rf, double k1, double k2,
                              double eVoid, double alpha, double b1, double b2,
							  double Pup, double Plow, double Pleft, double Pright)
  :Element(tag,ELE_TAG_SSPquadUP),
    theMaterial(0),
    mExternalNodes(SQUP_NUM_NODE),
    mTangentStiffness(SQUP_NUM_DOF,SQUP_NUM_DOF),
    mInternalForces(SQUP_NUM_DOF),
    Q(SQUP_NUM_DOF),
    mMass(SQUP_NUM_DOF,SQUP_NUM_DOF),
    mDamp(SQUP_NUM_DOF,SQUP_NUM_DOF),
    mNodeCrd(2,4),
    dN(4,2),
    Mmem(3,8),
    Kstab(8,8),
    mSolidK(8,8),
    mSolidM(8,8),
    mPerm(4,4),
    mThickness(thick),
    fBulk(Kf),
    fDens(Rf),
    mAlpha(alpha),
    mPorosity(0),
    applyLoad(0),
	pressureLoad(12),
	pressureUpperSide(Pup),
	pressureLowerSide(Plow),
	pressureLeftSide(Pleft),
	pressureRightSide(Pright)
{
    mExternalNodes(0) = Nd1;
    mExternalNodes(1) = Nd2;
    mExternalNodes(2) = Nd3;
    mExternalNodes(3) = Nd4;
                
    mThickness = thick;
    fBulk      = Kf;
    fDens      = Rf;
    mAlpha     = alpha;

    b[0] = b1;
    b[1] = b2;
        
    appliedB[0] = 0.0;
    appliedB[1] = 0.0;

    perm[0] = k1;
    perm[1] = k2;

	// LM change
	P[0] = Pup;
	P[1] = Plow;
    P[2] = Pleft;
	P[3] = Pright;
	// LM change

    mPorosity = eVoid/(1.0 + eVoid);

    const char *type = "PlaneStrain";

    // get copy of the material object
    NDMaterial *theMatCopy = theMat.getCopy(type);
    if (theMatCopy != 0) {
        theMaterial = (NDMaterial *)theMatCopy;
    } else {
        opserr << "SSPquadUP::SSPquadUP - failed to get copy of material model\n";;
    }

    // check material
    if (theMaterial == 0) {
        opserr << "SSPquadUP::SSPquadUP - failed to allocate material model pointer\n";
        exit(-1);
    }
}

// null constructor
SSPquadUP::SSPquadUP()
  :Element(0,ELE_TAG_SSPquadUP),
    theMaterial(0),
    mExternalNodes(SQUP_NUM_NODE),
    mTangentStiffness(SQUP_NUM_DOF,SQUP_NUM_DOF),
    mInternalForces(SQUP_NUM_DOF),
    Q(SQUP_NUM_DOF),
    mMass(SQUP_NUM_DOF,SQUP_NUM_DOF),
    mDamp(SQUP_NUM_DOF,SQUP_NUM_DOF),
    mNodeCrd(2,4),
    dN(4,2),
    Mmem(3,8),
    Kstab(8,8),
    mSolidK(8,8),
    mSolidM(8,8),
    mPerm(4,4),
    mThickness(0),
    fBulk(0),
    fDens(0),
    mAlpha(0),
    mPorosity(0),
    applyLoad(0),
    pressureLoad(12),
	pressureUpperSide(0.0),
	pressureLowerSide(0.0),
	pressureLeftSide(0.0),
	pressureRightSide(0.0)
{
}

// destructor
SSPquadUP::~SSPquadUP()
{
    if (theMaterial != 0) {
        delete theMaterial;
    }
}

int 
SSPquadUP::getNumExternalNodes(void) const
{
    return SQUP_NUM_NODE;
}

const ID &
SSPquadUP::getExternalNodes(void)
{
    return mExternalNodes;
}

Node **
SSPquadUP::getNodePtrs(void)
{
    return theNodes;
}

int
SSPquadUP::getNumDOF(void)
{
    return SQUP_NUM_DOF;
}

void
SSPquadUP::setDomain(Domain *theDomain)
{
    theNodes[0] = theDomain->getNode(mExternalNodes(0));
    theNodes[1] = theDomain->getNode(mExternalNodes(1));
    theNodes[2] = theDomain->getNode(mExternalNodes(2));
    theNodes[3] = theDomain->getNode(mExternalNodes(3));

    for (int i = 0; i < 4; i++) {
        if (theNodes[i] == 0) {
            return;  // don't go any further - otherwise segmentation fault
        }
    }

    // initialize coordinate vectors
    const Vector &mIcrd_1 = theNodes[0]->getCrds();
    const Vector &mIcrd_2 = theNodes[1]->getCrds();
    const Vector &mIcrd_3 = theNodes[2]->getCrds();
    const Vector &mIcrd_4 = theNodes[3]->getCrds();

    // coordinate matrix
    mNodeCrd(0,0) = mIcrd_1(0);
    mNodeCrd(1,0) = mIcrd_1(1);
    mNodeCrd(0,1) = mIcrd_2(0);
    mNodeCrd(1,1) = mIcrd_2(1);
    mNodeCrd(0,2) = mIcrd_3(0);
    mNodeCrd(1,2) = mIcrd_3(1);
    mNodeCrd(0,3) = mIcrd_4(0);
    mNodeCrd(1,3) = mIcrd_4(1);

    // establish jacobian terms 
    J0 = ((mNodeCrd(0,1)-mNodeCrd(0,3))*(mNodeCrd(1,2)-mNodeCrd(1,0))+(mNodeCrd(0,2)-mNodeCrd(0,0))*(mNodeCrd(1,3)-mNodeCrd(1,1)))/8;
    J1 = ((mNodeCrd(0,1)-mNodeCrd(0,0))*(mNodeCrd(1,2)-mNodeCrd(1,3))+(mNodeCrd(0,2)-mNodeCrd(0,3))*(mNodeCrd(1,0)-mNodeCrd(1,1)))/24;
    J2 = ((mNodeCrd(0,0)-mNodeCrd(0,3))*(mNodeCrd(1,2)-mNodeCrd(1,1))+(mNodeCrd(0,2)-mNodeCrd(0,1))*(mNodeCrd(1,3)-mNodeCrd(1,0)))/24;

    // establish stabilization terms (based on initial material tangent, only need to compute once)
    GetStab();

    // establish mass matrix for solid phase (constant, only need to compute once)
    GetSolidMass();

    // establish permeability matrix (constant, only need to compute once)
    GetPermeabilityMatrix();

    //LM change
	// Compute consistent nodal loads due to surface pressure (at any side)
    this->setPressureLoadAtNodes();
	//LM change
		
    // call the base-class method
    this->DomainComponent::setDomain(theDomain);
}

int
SSPquadUP::commitState(void)
{
    int retVal = 0;
    // call element commitState to do any base class stuff
    if ((retVal = this->Element::commitState()) != 0) {
        opserr << "SSPquadUP::commitState() - failed in base class\n";
    }
    retVal = theMaterial->commitState();

    return retVal;
}

int
SSPquadUP::revertToLastCommit(void)
{
    return theMaterial->revertToLastCommit();
}

int
SSPquadUP::revertToStart(void)
{
    return theMaterial->revertToStart();
}

int
SSPquadUP::update(void)
// this function updates variables for an incremental step n to n+1
{
    // get trial displacement
    const Vector &mDisp_1 = theNodes[0]->getTrialDisp();
    const Vector &mDisp_2 = theNodes[1]->getTrialDisp();
    const Vector &mDisp_3 = theNodes[2]->getTrialDisp();
    const Vector &mDisp_4 = theNodes[3]->getTrialDisp();
        
    // assemble displacement vector
    Vector u(8);
    u(0) = mDisp_1(0);
    u(1) = mDisp_1(1);
    u(2) = mDisp_2(0);
    u(3) = mDisp_2(1);
    u(4) = mDisp_3(0);
    u(5) = mDisp_3(1);
    u(6) = mDisp_4(0);
    u(7) = mDisp_4(1);

    Vector strain(3);
    strain = Mmem*u;
    theMaterial->setTrialStrain(strain);

    return 0;
}

const Matrix &
SSPquadUP::getTangentStiff(void)
// this function computes the tangent stiffness matrix for the element
{
    // solid phase stiffness matrix
    GetSolidStiffness();

    // assemble full element stiffness matrix [ K  0 ]
    // comprised of K submatrix               [ 0  0 ]
    mTangentStiffness.Zero();
    for (int i = 0; i < 4; i++) {

        int I    = 2*i;
        int Ip1  = 2*i+1;
        int II   = 3*i;
        int IIp1 = 3*i+1;

        for (int j = 0; j < 4; j++) {

            int J    = 2*j;
            int Jp1  = 2*j+1;
            int JJ   = 3*j;
            int JJp1 = 3*j+1;

            // contribution of solid phase stiffness matrix
            mTangentStiffness(II,JJ)     = mSolidK(I,J);
            mTangentStiffness(IIp1,JJ)   = mSolidK(Ip1,J);
            mTangentStiffness(IIp1,JJp1) = mSolidK(Ip1,Jp1);
            mTangentStiffness(II,JJp1)   = mSolidK(I,Jp1);
        }
    }

    return mTangentStiffness;
}

const Matrix &
SSPquadUP::getInitialStiff(void)
// this function computes the initial tangent stiffness matrix for the element
{
    return getTangentStiff();
}

const Matrix &
SSPquadUP::getDamp(void)
{
    Matrix dampC(8,8);

    // solid phase stiffness matrix
    GetSolidStiffness();

    // contribution of stiffness matrix for Rayleigh damping
    if (betaK != 0.0) {
        dampC.addMatrix(1.0, mSolidK, betaK);
    } if (betaK0 != 0.0) {
        dampC.addMatrix(1.0, mSolidK, betaK0);
    } if (betaKc != 0.0) {
        dampC.addMatrix(1.0, mSolidK, betaKc);
    }

    // contribution of mass matrix for Rayleigh damping
    if (alphaM != 0.0) {
        dampC.addMatrix(1.0, mSolidM, alphaM);
    }

    // assemble full element damping matrix   [  C  -Q ]
    // comprised of C, Q, and H submatrices   [ -Q' -H ]
    mDamp.Zero();
    for (int i = 0; i < 4; i++) {

        int I    = 2*i;
        int Ip1  = 2*i+1;
        int II   = 3*i;
        int IIp1 = 3*i+1;
        int IIp2 = 3*i+2;

        for (int j = 0; j < 4; j++) {

            int J    = 2*j;
            int Jp1  = 2*j+1;
            int JJ   = 3*j;
            int JJp1 = 3*j+1;
            int JJp2 = 3*j+2;

            // contribution of solid phase damping matrix
            mDamp(II,JJ)     = dampC(I,J);
            mDamp(IIp1,JJ)   = dampC(Ip1,J);
            mDamp(IIp1,JJp1) = dampC(Ip1,Jp1);
            mDamp(II,JJp1)   = dampC(I,Jp1);

            // contribution of solid-fluid coupling matrix
            mDamp(JJp2,II)   = -J0*mThickness*Mmem(0,I);
            mDamp(JJp2,IIp1) = -J0*mThickness*Mmem(1,Ip1);
            mDamp(II,JJp2)   = -J0*mThickness*Mmem(0,I);
            mDamp(IIp1,JJp2) = -J0*mThickness*Mmem(1,Ip1);
                        
            // contribution of permeability matrix
            mDamp(IIp2,JJp2) = -mPerm(i,j);
        }
    }

    return mDamp;
}

const Matrix &
SSPquadUP::getMass(void)
{
    mMass.Zero();

    // compute compressibility matrix term
    double oneOverQ = -0.25*J0*mThickness*mPorosity/fBulk;

    // get mass density from the material
    double density = theMaterial->getRho();

    // transpose the shape function derivative array
    Matrix dNp(2,4);
    dNp(0,0) = dN(0,0); dNp(0,1) = dN(1,0); dNp(0,2) = dN(2,0); dNp(0,3) = dN(3,0);
    dNp(1,0) = dN(0,1); dNp(1,1) = dN(1,1); dNp(1,2) = dN(2,1); dNp(1,3) = dN(3,1);

    // compute stabilization matrix for incompressible problems
    Matrix Kp(4,4);
    Kp = -4.0*mAlpha*J0*mThickness*dN*dNp;

    // return zero matrix if density is zero
    if (density == 0.0) {
        return mMass;
    }

    // full mass matrix for the element [ M  0 ]
    //  includes M and S submatrices    [ 0 -S ]
    for (int i = 0; i < 4; i++) {

        int I    = 2*i;
        int Ip1  = 2*i+1;
        int II   = 3*i;
        int IIp1 = 3*i+1;
        int IIp2 = 3*i+2;

        for (int j = 0; j < 4; j++) {

            int J    = 2*j;
            int Jp1  = 2*j+1;
            int JJ   = 3*j;
            int JJp1 = 3*j+1;
            int JJp2 = 3*j+2;

            mMass(II,JJ)     = mSolidM(I,J);
            mMass(IIp1,JJ)   = mSolidM(Ip1,J);
            mMass(IIp1,JJp1) = mSolidM(Ip1,Jp1);
            mMass(II,JJp1)   = mSolidM(I,Jp1);

            // contribution of compressibility matrix
            mMass(IIp2,JJp2) = Kp(i,j) + oneOverQ;
        }
    }

    return mMass;
}

void
SSPquadUP::zeroLoad(void)
{
    applyLoad = 0;
    appliedB[0] = 0.0;
    appliedB[1] = 0.0;
  
    Q.Zero();

    return;
}

int
SSPquadUP::addLoad(ElementalLoad *theLoad, double loadFactor)
{
    // body forces can be applied in a load pattern
    int type;
    const Vector &data = theLoad->getData(type, loadFactor);

    if (type == LOAD_TAG_SelfWeight) {
        applyLoad = 1;
        appliedB[0] += loadFactor*data(0)*b[0];
        appliedB[1] += loadFactor*data(1)*b[1];
        return 0;
    } else {
        opserr << "SSPquadUP::addLoad - load type unknown for ele with tag: " << this->getTag() << endln;
        return -1;
    } 

    return -1;
}

int
SSPquadUP::addInertiaLoadToUnbalance(const Vector &accel)
{
	// get mass density from the material
	double density = theMaterial->getRho();

	// do nothing if density is zero
	if (density == 0.0) {
		return 0;
	}

	// Get R * accel from the nodes
	const Vector &Raccel1 = theNodes[0]->getRV(accel);
	const Vector &Raccel2 = theNodes[1]->getRV(accel);
	const Vector &Raccel3 = theNodes[2]->getRV(accel);
	const Vector &Raccel4 = theNodes[3]->getRV(accel);

	if (3 != Raccel1.Size() || 3 != Raccel2.Size() || 3 != Raccel3.Size() || 3 != Raccel4.Size()) {
    	opserr << "SSPquadUP::addInertiaLoadToUnbalance matrix and vector sizes are incompatible\n";
    	return -1;
	}

	static double ra[12];
	ra[0]  = Raccel1(0);
	ra[1]  = Raccel1(1);
	ra[2]  = 0.0;
	ra[3]  = Raccel2(0);
	ra[4]  = Raccel2(1);
	ra[5]  = 0.0;
	ra[6]  = Raccel3(0);
	ra[7]  = Raccel3(1);
	ra[8]  = 0.0;
	ra[9]  = Raccel4(0);
	ra[10] = Raccel4(1);
	ra[11] = 0.0;

	// compute mass matrix
	this->getMass();

	for (int i = 0; i < 12; i++) {
		Q(i) += -mMass(i,i)*ra[i];
	}
	
	return 0;
}

const Vector &
SSPquadUP::getResistingForce(void)
// this function computes the resisting force vector for the element
{
    Vector f1(8);
    Vector f2(4);
    Vector mStress(3); 
        
    // get stress from the material
    mStress = theMaterial->getStress();

    // get trial displacement
    const Vector &mDisp_1 = theNodes[0]->getTrialDisp();
    const Vector &mDisp_2 = theNodes[1]->getTrialDisp();
    const Vector &mDisp_3 = theNodes[2]->getTrialDisp();
    const Vector &mDisp_4 = theNodes[3]->getTrialDisp();

    Vector d(8);
    d(0) = mDisp_1(0);
    d(1) = mDisp_1(1);
    d(2) = mDisp_2(0);
    d(3) = mDisp_2(1);
    d(4) = mDisp_3(0);
    d(5) = mDisp_3(1);
    d(6) = mDisp_4(0);
    d(7) = mDisp_4(1);

    // add stabilization force to internal force vector 
    f1 = Kstab*d;

    // add internal force from the stress
    f1.addMatrixTransposeVector(1.0, Mmem, mStress, 4.0*mThickness*J0);

    // get mass density from the material
    double density = theMaterial->getRho();

    // subtract body forces from internal force vector 
    double xi[4];
    double eta[4];
    xi[0]  = -1.0; xi[1]  =  1.0; xi[2]  = 1.0; xi[3]  = -1.0;
    eta[0] = -1.0; eta[1] = -1.0; eta[2] = 1.0; eta[3] =  1.0;

    if (applyLoad == 0) {
        for (int i = 0; i < 4; i++) {
            f1(2*i)   -= density*b[0]*mThickness*(J0 + J1*xi[i] + J2*eta[i]);
            f1(2*i+1) -= density*b[1]*mThickness*(J0 + J1*xi[i] + J2*eta[i]);
        }
    } else {
        for (int i = 0; i < 4; i++) {
            f1(2*i)   -= density*appliedB[0]*mThickness*(J0 + J1*xi[i] + J2*eta[i]);
            f1(2*i+1) -= density*appliedB[1]*mThickness*(J0 + J1*xi[i] + J2*eta[i]);
        }
    }

    // account for fluid body forces
    Matrix k(2,2);
    Vector body(2);
    // permeability tensor
    k(0,0) = perm[0];
    k(1,1) = perm[1];
    // body force vector
    if (applyLoad == 0) {
        body(0) = b[0];
        body(1) = b[1];
    } else {
        body(0) = appliedB[0];
        body(1) = appliedB[1];
    }
    f2 = 4.0*J0*mThickness*fDens*dN*k*body;

    // assemble full internal force vector for the element
    mInternalForces(0)  = f1(0);
    mInternalForces(1)  = f1(1);
    mInternalForces(2)  = f2(0);
    mInternalForces(3)  = f1(2);
    mInternalForces(4)  = f1(3);
    mInternalForces(5)  = f2(1);
    mInternalForces(6)  = f1(4);
    mInternalForces(7)  = f1(5);
    mInternalForces(8)  = f2(2);
    mInternalForces(9)  = f1(6);
    mInternalForces(10) = f1(7);
    mInternalForces(11) = f2(3);

    //LM change
    // Subtract pressure loading from internal force vector
    if (pressureUpperSide != 0.0 || pressureLowerSide != 0.0 || pressureLeftSide != 0.0 || pressureRightSide != 0.0) {
        mInternalForces.addVector(1.0, pressureLoad, -1.0);
    }
    //LM change
    
    // inertial unbalance load
    mInternalForces.addVector(1.0, Q, -1.0);

    return mInternalForces;
}

const Vector &
SSPquadUP::getResistingForceIncInertia()
{
	// terms stemming from acceleration
	const Vector &accel1 = theNodes[0]->getTrialAccel();
	const Vector &accel2 = theNodes[1]->getTrialAccel();
	const Vector &accel3 = theNodes[2]->getTrialAccel();
	const Vector &accel4 = theNodes[3]->getTrialAccel();
	
	// compute current resisting force
	this->getResistingForce();

	// compute mass matrix
	this->getMass();

	Vector a(12);
	a(0)  = accel1(0);
	a(1)  = accel1(1);
	a(2)  = accel1(2);
	a(3)  = accel2(0);
	a(4)  = accel2(1);
	a(5)  = accel2(2);
	a(6)  = accel3(0);
	a(7)  = accel3(1);
	a(8)  = accel3(2);
	a(9)  = accel4(0);
	a(10) = accel4(1);
	a(11) = accel4(2);

	mInternalForces.addMatrixVector(1.0, mMass, a, 1.0);

	// terms stemming from velocity
	const Vector &vel1 = theNodes[0]->getTrialVel();
	const Vector &vel2 = theNodes[1]->getTrialVel();
	const Vector &vel3 = theNodes[2]->getTrialVel();
	const Vector &vel4 = theNodes[3]->getTrialVel();
	
	Vector v(12);
	v(0)  = vel1(0);
	v(1)  = vel1(1);
	v(2)  = vel1(2);
	v(3)  = vel2(0);
	v(4)  = vel2(1);
	v(5)  = vel2(2);
	v(6)  = vel3(0);
	v(7)  = vel3(1);
	v(8)  = vel3(2);
	v(9)  = vel4(0);
	v(10) = vel4(1);
	v(11) = vel4(2);

	// compute damping matrix
	this->getDamp();

	mInternalForces.addMatrixVector(1.0, mDamp, v, 1.0);

	return mInternalForces;
}

int
SSPquadUP::sendSelf(int commitTag, Channel &theChannel)
{
    int res = 0;
  
    // note: we don't check for dataTag == 0 for Element
    // objects as that is taken care of in a commit by the Domain
    // object - don't want to have to do the check if sending data
    int dataTag = this->getDbTag();
  
    // SSPquadUP packs its data into a Vector and sends this to theChannel
    // along with its dbTag and the commitTag passed in the arguments
    //LM change
	static Vector data(15);
    data(0) = this->getTag();
    data(1) = mThickness;
    data(2) = fBulk;
    data(3) = fDens;
    data(4) = perm[0];
    data(5) = perm[1];
    data(6) = mPorosity;
    data(7) = b[0];
    data(8) = b[1];
	data(9) = P[0];
	data(10) = P[1];
	data(11) = P[2];
	data(12) = P[3];
    data(13) = theMaterial->getClassTag();         
    //LM change
        
    // Now SSPquadUP sends the ids of its materials
    int matDbTag = theMaterial->getDbTag();
  
    static ID idData(12);
  
    // NOTE: we do have to ensure that the material has a database
    // tag if we are sending to a database channel.
    if (matDbTag == 0) {
        matDbTag = theChannel.getDbTag();
        if (matDbTag != 0) {
            theMaterial->setDbTag(matDbTag);
        }
    }
    data(14) = matDbTag;

    res += theChannel.sendVector(dataTag, commitTag, data);
    if (res < 0) {
        opserr << "WARNING SSPquadUP::sendSelf() - " << this->getTag() << " failed to send Vector\n";
        return res;
    }

    // SSPquadUP then sends the tags of its four nodes
    res += theChannel.sendID(dataTag, commitTag, mExternalNodes);
    if (res < 0) {
        opserr << "WARNING SSPquadUP::sendSelf() - " << this->getTag() << " failed to send ID\n";
        return res;
    }

    // finally, SSPquadUP asks its material object to send itself
    res = theMaterial->sendSelf(commitTag, theChannel);
    if (res < 0) {
        opserr << "WARNING SSPquadUP::sendSelf() - " << this->getTag() << " failed to send its Material\n";
        return -3;
    }

    return 0;
}

int
SSPquadUP::recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker)
{
    int res = 0;
    int dataTag = this->getDbTag();

    // SSPquadUP creates a Vector, receives the Vector and then sets the 
    // internal data with the data in the Vector
    static Vector data(15);
    res += theChannel.recvVector(dataTag, commitTag, data);
    if (res < 0) {
        opserr << "WARNING SSPquadUP::recvSelf() - failed to receive Vector\n";
        return res;
    }
  
    this->setTag((int)data(0));
    mThickness = data(1);
    fBulk = data(2);
    fDens = data(3);
    perm[0] = data(4);
    perm[1] = data(5);
    mPorosity = data(6);
    b[0] = data(7);
    b[1] = data(8);
	//LM change
	P[0] = data(9);
	P[1] = data(10);
	P[2] = data(11);
	P[3] = data(12);
    //LM change

    // SSPquadUP now receives the tags of its four external nodes
    res += theChannel.recvID(dataTag, commitTag, mExternalNodes);
    if (res < 0) {
        opserr << "WARNING SSPquadUP::recvSelf() - " << this->getTag() << " failed to receive ID\n";
        return res;
    }

    // finally, SSPquadUP creates a material object of the correct type, sets its
    // database tag, and asks this new object to receive itself
    //LM change
	int matClass = (int)data(13);
    int matDb    = (int)data(14);
    //LM change

    // check if material object exists and that it is the right type
    if ((theMaterial == 0) || (theMaterial->getClassTag() != matClass)) {

        // if old one, delete it
        if (theMaterial != 0)
            delete theMaterial;

        // create new material object
        NDMaterial *theMatCopy = theBroker.getNewNDMaterial(matClass);
        theMaterial = (NDMaterial *)theMatCopy;

        if (theMaterial == 0) {
            opserr << "WARNING SSPquadUP::recvSelf() - " << this->getTag() 
                   << " failed to get a blank Material of type " << matClass << endln;
            return -3;
        }
    }

    // NOTE: we set the dbTag before we receive the material
    theMaterial->setDbTag(matDb);
    res = theMaterial->recvSelf(commitTag, theChannel, theBroker);
    if (res < 0) {
        opserr << "WARNING SSPquadUP::recvSelf() - " << this->getTag() << " failed to receive its Material\n";
        return -3;
    }

    return 0; 
}

int
SSPquadUP::displaySelf(Renderer &theViewer, int displayMode, float fact, const char **modes, int numMode)
{
    // get the end point display coords
    static Vector v1(3);
    static Vector v2(3);
    static Vector v3(3);
    static Vector v4(3);
    theNodes[0]->getDisplayCrds(v1, fact, displayMode);
    theNodes[1]->getDisplayCrds(v2, fact, displayMode);
    theNodes[2]->getDisplayCrds(v3, fact, displayMode);
    theNodes[3]->getDisplayCrds(v4, fact, displayMode);

    // place values in coords matrix
    static Matrix coords(4, 3);
    for (int i = 0; i < 3; i++) {
        coords(0, i) = v1(i);
        coords(1, i) = v2(i);
        coords(2, i) = v3(i);
        coords(3, i) = v4(i);
    }

    // fill RGB vector
    static Vector values(4);
    for (int i = 0; i < 4; i++)
        values(i) = 1.0;

    // draw the polygon
    return theViewer.drawPolygon(coords, values, this->getTag());
}

void
SSPquadUP::Print(OPS_Stream &s, int flag)
{
    if (flag == OPS_PRINT_CURRENTSTATE) {
        opserr << "SSPquadUP, element id:  " << this->getTag() << endln;
        opserr << "   Connected external nodes:  ";
        for (int i = 0; i < SQUP_NUM_NODE; i++) {
            opserr << mExternalNodes(i) << " ";
        }
    }
    
    if (flag == OPS_PRINT_PRINTMODEL_JSON) {
        s << "\t\t\t{";
        s << "\"name\": " << this->getTag() << ", ";
        s << "\"type\": \"SSPquadUP\", ";
        s << "\"nodes\": [" << mExternalNodes(0) << ", ";
        s << mExternalNodes(1) << ", ";
        s << mExternalNodes(2) << ", ";
        s << mExternalNodes(3) << "], ";
        s << "\"thickness\": " << mThickness << ", ";
        s << "\"bodyForces\": [" << b[0] << ", " << b[1] << "], ";
        s << "\"material\": \"" << theMaterial->getTag() << "\"}";
    }
}

Response*
SSPquadUP::setResponse(const char **argv, int argc, OPS_Stream &eleInfo)
{
    // no special recorders for this element, call the method in the material class
    return theMaterial->setResponse(argv, argc, eleInfo);
}

int
SSPquadUP::getResponse(int responseID, Information &eleInfo)
{
    // no special recorders for this element, call the method in the material class
    return theMaterial->getResponse(responseID, eleInfo);
}

int
SSPquadUP::setParameter(const char **argv, int argc, Parameter &param)
{
    if (argc < 1) {
        return -1;
    }
    int res = -1;

    // check for element parameters first
    if (strcmp(argv[0],"hPerm") == 0) {
        return param.addObject(3, this);
    } else if (strcmp(argv[0],"vPerm") == 0) {
        return param.addObject(4, this);
    //LM change	
    // quad pressure loading
    } else if (strcmp(argv[0],"pressureUpperSide") == 0) {
        return param.addObject(9, this);
    } else if (strcmp(argv[0],"pressureLowerSide") == 0) {
        return param.addObject(10, this);     
    } else if (strcmp(argv[0],"pressureLeftSide") == 0) {
        return param.addObject(11, this);
    } else if (strcmp(argv[0],"pressureRightSide") == 0) {
        return param.addObject(12, this);
    } else if (strcmp(argv[0],"b1") == 0) {
      return param.addObject(13,this);
    } else if (strcmp(argv[0],"b2") == 0) {
      return param.addObject(14,this);
    //LM change
    } else {
        // default is to call setParameter in the material
        int matRes = res;
        matRes =  theMaterial->setParameter(argv, argc, param);
        if (matRes != -1) {
            res = matRes;
        }
    }
  
  return res;
}

int
SSPquadUP::updateParameter(int parameterID, Information &info)
{
    int res = -1;
    int matRes = res;
    
    if (parameterID == res) {
        return -1;
    } else if (parameterID == 3) {
        // update element permeability in direction 1
        perm[0] = info.theDouble;
        GetPermeabilityMatrix();
        return 0;
    } else if (parameterID == 4) {
        // update element permeability in direction 2
        perm[1] = info.theDouble;
        GetPermeabilityMatrix();
        return 0;
    //LM change
	} else if (parameterID == 9) {
        // update consistent nodal loads - pressure positive for outward normal
        pressureUpperSide = info.theDouble;
        this->setPressureLoadAtNodes();
        return 0;
	} else if (parameterID == 10) {
        // update consistent nodal loads - pressure positive for outward normal
        pressureLowerSide = info.theDouble;
        this->setPressureLoadAtNodes();
        return 0;
	} else if (parameterID == 11) {
        // update consistent nodal loads - pressure positive for outward normal
        pressureLeftSide = info.theDouble;
        this->setPressureLoadAtNodes();
        return 0;
	} else if (parameterID == 12) {
        // update consistent nodal loads - pressure positive for outward normal
        pressureRightSide = info.theDouble;
        this->setPressureLoadAtNodes();
        return 0;		
    } else if (parameterID == 13) {
      b[0] = info.theDouble;
      return 0;
    } else if (parameterID == 14) {
      b[1] = info.theDouble;
      return 0;
	//LM change
	} else {
        // update the material parameter
        matRes = theMaterial->updateParameter(parameterID, info);
        if (matRes != -1) {
            res = matRes;
        }
        return res;
    }
}

Matrix 
SSPquadUP::DyadicProd(Vector v1, Vector v2)
// computes dyadic product for two vectors (2x1)
{
	Matrix result(2,2);
	result.Zero();

	for (int i = 0; i < v1.Size(); i++) {
		for (int j = 0; j < v2.Size(); j++) 
			result(i,j) = v1(i) * v2(j);
	}
	
	return result;
}

void
SSPquadUP::GetStab(void)
// this function computes the stabilization stiffness matrix for the element
{
	Vector g1(SQUP_NUM_DIM);
	Vector g2(SQUP_NUM_DIM);
	Matrix I(SQUP_NUM_DIM,SQUP_NUM_DIM);
	Matrix FCF(SQUP_NUM_DIM,SQUP_NUM_DIM);
	Matrix Jmat(SQUP_NUM_DIM,SQUP_NUM_DIM);
	Matrix Jinv(SQUP_NUM_DIM,SQUP_NUM_DIM);
	Matrix dNloc(SQUP_NUM_NODE,SQUP_NUM_DIM);
	Matrix Mben(2,SQUP_NUM_DOF);
	double Hss;
	double Hst;
	double Htt;
	
	// shape function derivatives (local crd) at center
	dNloc(0,0) = -0.25;
	dNloc(1,0) =  0.25;
	dNloc(2,0) =  0.25;
	dNloc(3,0) = -0.25;
	dNloc(0,1) = -0.25;
	dNloc(1,1) = -0.25;
	dNloc(2,1) =  0.25;
	dNloc(3,1) =  0.25;

	// jacobian matrix
	Jmat = mNodeCrd*dNloc;
	// inverse of the jacobian matrix
	Jmat.Invert(Jinv);

	// shape function derivatives (global crd)
	dN = dNloc*Jinv;

	// define hourglass stabilization vector  gamma = 0.25*(h - (h^x)*bx - (h^y)*by);
	double hx = mNodeCrd(0,0) - mNodeCrd(0,1) + mNodeCrd(0,2) - mNodeCrd(0,3);
	double hy = mNodeCrd(1,0) - mNodeCrd(1,1) + mNodeCrd(1,2) - mNodeCrd(1,3);
	double gamma[4];
	gamma[0] = 0.25*( 1.0 - hx*dN(0,0) - hy*dN(0,1));
	gamma[1] = 0.25*(-1.0 - hx*dN(1,0) - hy*dN(1,1));
	gamma[2] = 0.25*( 1.0 - hx*dN(2,0) - hy*dN(2,1));
	gamma[3] = 0.25*(-1.0 - hx*dN(3,0) - hy*dN(3,1));

	// define mapping matrices
	Mmem.Zero();
	Mben.Zero();
	for (int i = 0; i < 4; i++) {
		Mmem(0,2*i)   = dN(i,0);
		Mmem(1,2*i+1) = dN(i,1);
		Mmem(2,2*i)   = dN(i,1);
		Mmem(2,2*i+1) = dN(i,0);

		Mben(0,2*i)   = gamma[i];
		Mben(1,2*i+1) = gamma[i];
	}

	// base vectors
	g1(0) = Jmat(0,0);
	g1(1) = Jmat(1,0);
	g2(0) = Jmat(0,1);
	g2(1) = Jmat(1,1);
	// normalize base vectors
	g1.Normalize();
	g2.Normalize();
	
	// compute second moment of area tensor
	double fourThree = 4.0/3.0;
	I = fourThree*mThickness*J0*(DyadicProd(g1,g1) + DyadicProd(g2,g2));

	// stabilization terms
	Hss = ( I(0,0)*Jinv(1,0)*Jinv(1,0) + I(0,1)*Jinv(0,0)*Jinv(1,0) + I(1,1)*Jinv(0,0)*Jinv(0,0) )*0.25;
	Htt = ( I(0,0)*Jinv(1,1)*Jinv(1,1) + I(0,1)*Jinv(0,1)*Jinv(1,1) + I(1,1)*Jinv(0,1)*Jinv(0,1) )*0.25;
	Hst = ( I(0,0)*Jinv(1,1)*Jinv(1,0) + I(0,1)*(Jinv(1,0)*Jinv(0,1) + Jinv(1,1)*Jinv(0,0)) + I(1,1)*Jinv(0,1)*Jinv(0,0) )*0.25;

	// get material tangent
	const Matrix &CmatI = theMaterial->getInitialTangent();

	// compute stabilization matrix
	FCF(0,0) = (CmatI(0,0) - (CmatI(0,1) + CmatI(1,0)) + CmatI(1,1))*Hss;
	FCF(0,1) = (CmatI(0,1) - (CmatI(0,0) + CmatI(1,1)) + CmatI(1,0))*Hst;
	FCF(1,0) = (CmatI(1,0) - (CmatI(0,0) + CmatI(1,1)) + CmatI(0,1))*Hst;
	FCF(1,1) = (CmatI(1,1) - (CmatI(0,1) + CmatI(1,0)) + CmatI(0,0))*Htt;
	
	// compute stiffness matrix for stabilization terms
	Kstab.Zero();
	Kstab.addMatrixTripleProduct(1.0, Mben, FCF, 1.0);

	return;
}

void
SSPquadUP::GetSolidStiffness(void)
// this function computes the stiffness matrix for the solid phase
{
	// get material tangent
	const Matrix &Cmat = theMaterial->getTangent();

	mSolidK = Kstab;
	mSolidK.addMatrixTripleProduct(1.0, Mmem, Cmat, 4.0*J0*mThickness);

	return;
}

void
SSPquadUP::GetSolidMass(void)
// this function computes the mass matrix for the solid phase
{
	// get mass density from the material
	double density = theMaterial->getRho();

	// local coordinates of nodes
	double xi[4];
	double eta[4];
	xi[0]  = -1.0; xi[1]  =  1.0; xi[2]  = 1.0; xi[3]  = -1.0;
	eta[0] = -1.0; eta[1] = -1.0; eta[2] = 1.0; eta[3] =  1.0;

	double massTerm;
	for (int i = 0; i < 4; i++) {
		massTerm = density*mThickness*(J0 + J1*xi[i] + J2*eta[i]);
		mSolidM(2*i,2*i)     += massTerm;
		mSolidM(2*i+1,2*i+1) += massTerm;
	}

	return;
}

void
SSPquadUP::GetPermeabilityMatrix(void)
// this function computes the permeability matrix for the element
{
	mPerm.Zero();
	Matrix k(2,2);
	Matrix dNp(2,4);

	// permeability tensor
	k(0,0) = perm[0];
	k(1,1) = perm[1];

	// transpose the shape function derivative array
	dNp(0,0) = dN(0,0); dNp(0,1) = dN(1,0); dNp(0,2) = dN(2,0); dNp(0,3) = dN(3,0);
	dNp(1,0) = dN(0,1); dNp(1,1) = dN(1,1); dNp(1,2) = dN(2,1); dNp(1,3) = dN(3,1);

	// compute permeability matrix
	mPerm.addMatrixTripleProduct(1.0, dNp, k, 4.0*J0*mThickness);

	return;
}

// LM change
// nodes numbering
// 4-----3
// |     |
// |     |
// 1-----2
// UpperSide - 34 LowerSide - 12 LeftSide 14 RightSide - 23
void
SSPquadUP::setPressureLoadAtNodes(void)
{
    pressureLoad.Zero();

    if (pressureUpperSide == 0.0 && pressureLowerSide == 0.0 && pressureLeftSide == 0.0 && pressureRightSide == 0.0) {
        return;
    }

    const Vector &node1 = theNodes[0]->getCrds();
    const Vector &node2 = theNodes[1]->getCrds();
    const Vector &node3 = theNodes[2]->getCrds();
    const Vector &node4 = theNodes[3]->getCrds();

    double x1 = node1(0);
    double y1 = node1(1);
    double x2 = node2(0);
    double y2 = node2(1);
    double x3 = node3(0);
    double y3 = node3(1);
    double x4 = node4(0);
    double y4 = node4(1);

    double dx12 = x2-x1;
    double dy12 = y2-y1;
    double dx23 = x3-x2;
    double dy23 = y3-y2;
    double dx34 = x4-x3;
    double dy34 = y4-y3;
    double dx41 = x1-x4;
    double dy41 = y1-y4;

    double pressureOver12 = pressureLowerSide*mThickness/2.0;
    double pressureOver23 = pressureRightSide*mThickness/2.0;
    double pressureOver34 = pressureUpperSide*mThickness/2.0;
    double pressureOver41 = pressureLeftSide*mThickness/2.0;

    // Contribution from side 12
    pressureLoad(0)  += pressureOver12*dy12;  // horizontal node 1 load
    pressureLoad(3)  += pressureOver12*dy12;  // horizontal node 2 load
    pressureLoad(1)  += pressureOver12*-dx12; // vertical node 1 load
    pressureLoad(2)  += pressureOver12*-dx12; // vertical node 2 load

    // Contribution from side 12
    pressureLoad(3)  += pressureOver23*dy23;  // horizontal node 2 load
    pressureLoad(6)  += pressureOver23*dy23;  // horizontal node 3 load
    pressureLoad(4)  += pressureOver23*-dx23; // vertical node 2 load
    pressureLoad(7)  += pressureOver23*-dx23; // vertical node 3 load

    // Contribution from side 12
    pressureLoad(6)  += pressureOver34*dy34;  // horizontal node 3 load
    pressureLoad(9)  += pressureOver34*dy34;  // horizontal node 4 load
    pressureLoad(7)  += pressureOver34*-dx34; // vertical node 3 load
    pressureLoad(10) += pressureOver34*-dx34; // vertical node 4 load

    // Contribution from side 12
    pressureLoad(9)  += pressureOver41*dy41;  // horizontal node 4 load
    pressureLoad(0)  += pressureOver41*dy41;  // horizontal node 1 load
    pressureLoad(10) += pressureOver41*-dx41; // vertical node 4 load
    pressureLoad(1)  += pressureOver41*-dx41; // vertical node 1 load
}
// LM change
