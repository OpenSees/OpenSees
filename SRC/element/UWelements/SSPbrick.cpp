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
                                                                       
// Created: Chris McGann, UW, 10.2011
//
// Description: This file contains the implementation of the SSPbrick class

#include "SSPbrick.h"

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

static int num_SSPbrick = 0;

OPS_Export void *
OPS_SSPbrick(void)
{
  if (num_SSPbrick == 0) {
    num_SSPbrick++;
    opserr << "SSPbrick element - Written: C.McGann, P.Arduino, P.Mackenzie-Helnwein, U.Washington\n";
  }
  
  	// Pointer to an element that will be returned
  	Element *theElement = 0;

  	int numRemainingInputArgs = OPS_GetNumRemainingInputArgs();

  	if (numRemainingInputArgs < 10) {
    	opserr << "Invalid #args, want: element SSPbrick eleTag? iNode? jNode? kNode? lNode? mNode? nNode? pNode? qNode? matTag? <b1? b2? b3?>\n";
		return 0;
  	}

  	int iData[10];
	double dData[3];
	dData[0] = 0.0;
	dData[1] = 0.0;
	dData[2] = 0.0;

  	int numData = 10;
  	if (OPS_GetIntInput(&numData, iData) != 0) {
    	opserr << "WARNING invalid integer data: element SSPbrick " << iData[0] << endln;
		return 0;
  	}

  	int matID = iData[9];
  	NDMaterial *theMaterial = OPS_getNDMaterial(matID);
  	if (theMaterial == 0) {
    	opserr << "WARNING element SSPbrick " << iData[0] << endln;
		opserr << " Material: " << matID << "not found\n";
		return 0;
  	}

	if (numRemainingInputArgs == 13) {
    	numData = 3;
    	if (OPS_GetDoubleInput(&numData, dData) != 0) {
      		opserr << "WARNING invalid optional data: element SSPbrick " << iData[0] << endln;
	  		return 0;
    	}
  	}

  	// parsing was successful, allocate the element
  	theElement = new SSPbrick(iData[0], iData[1], iData[2], iData[3], iData[4], iData[5], iData[6], iData[7], iData[8],
                                        *theMaterial, dData[0], dData[1], dData[2]);

  	if (theElement == 0) {
    	opserr << "WARNING could not create element of type SSPbrick\n";
		return 0;
  	}

  	return theElement;
}

// full constructor
SSPbrick::SSPbrick(int tag, int Nd1, int Nd2, int Nd3, int Nd4, int Nd5, int Nd6, int Nd7, int Nd8,
                      NDMaterial &theMat, double b1, double b2, double b3)
  :Element(tag,ELE_TAG_SSPbrick),
  	theMaterial(0),
	mExternalNodes(SSPB_NUM_NODE),
	mTangentStiffness(SSPB_NUM_DOF,SSPB_NUM_DOF),
	mInternalForces(SSPB_NUM_DOF),
	Q(SSPB_NUM_DOF),
	mMass(SSPB_NUM_DOF,SSPB_NUM_DOF),
	mNodeCrd(SSPB_NUM_DIM,SSPB_NUM_NODE),
	mVol(0),
	Bnot(6,SSPB_NUM_DOF),
	Kstab(SSPB_NUM_DOF,SSPB_NUM_DOF),
	xi(8),
	et(8),
	ze(8),
	hut(8),
	hus(8),
	hst(8),
	hstu(8)
{
	mExternalNodes(0) = Nd1;
	mExternalNodes(1) = Nd2;
	mExternalNodes(2) = Nd3;
	mExternalNodes(3) = Nd4;
	mExternalNodes(4) = Nd5;
	mExternalNodes(5) = Nd6;
	mExternalNodes(6) = Nd7;
	mExternalNodes(7) = Nd8;

	b[0] = b1;
	b[1] = b2;
	b[2] = b3;

	applyLoad = 0;
	
	appliedB[0] = 0.0;
	appliedB[1] = 0.0;
	appliedB[2] = 0.0;

	// get copy of the material object
	NDMaterial *theMatCopy = theMat.getCopy("ThreeDimensional");
	if (theMatCopy != 0) {
		theMaterial = (NDMaterial *)theMatCopy;
	} else {
		opserr << "SSPbrick::SSPbrick - failed to get copy of material model\n";;
	}

	// check material
	if (theMaterial == 0) {
		opserr << "SSPbrick::SSPbrick - failed to allocate material model pointer\n";
		exit(-1);
	}

	mInitialize = true;
}

// null constructor
SSPbrick::SSPbrick()
  :Element(0,ELE_TAG_SSPbrick),
  	theMaterial(0),
	mExternalNodes(SSPB_NUM_NODE),
	mTangentStiffness(SSPB_NUM_DOF,SSPB_NUM_DOF),
	mInternalForces(SSPB_NUM_DOF),
	Q(SSPB_NUM_DOF),
	mMass(SSPB_NUM_DOF,SSPB_NUM_DOF),
	mNodeCrd(SSPB_NUM_DIM,SSPB_NUM_NODE),
	mVol(0),
	Bnot(6,SSPB_NUM_DOF),
	Kstab(SSPB_NUM_DOF,SSPB_NUM_DOF),
	xi(8),
	et(8),
	ze(8),
	hut(8),
	hus(8),
	hst(8),
	hstu(8)
{
	b[0] = 0.0;
	b[1] = 0.0;
	b[2] = 0.0;

	applyLoad = 0;
	
	appliedB[0] = 0.0;
	appliedB[1] = 0.0;
	appliedB[2] = 0.0;

	mInitialize = false;
}

// destructor
SSPbrick::~SSPbrick()
{
}

int 
SSPbrick::getNumExternalNodes(void) const
{
    return SSPB_NUM_NODE;
}

const ID &
SSPbrick::getExternalNodes(void)
{
    return mExternalNodes;
}

Node **
SSPbrick::getNodePtrs(void)
{
    return theNodes;
}

int
SSPbrick::getNumDOF(void)
{
    return SSPB_NUM_DOF;
}

void
SSPbrick::setDomain(Domain *theDomain)
{
	theNodes[0] = theDomain->getNode(mExternalNodes(0));
	theNodes[1] = theDomain->getNode(mExternalNodes(1));
	theNodes[2] = theDomain->getNode(mExternalNodes(2));
	theNodes[3] = theDomain->getNode(mExternalNodes(3));
	theNodes[4] = theDomain->getNode(mExternalNodes(4));
	theNodes[5] = theDomain->getNode(mExternalNodes(5));
	theNodes[6] = theDomain->getNode(mExternalNodes(6));
	theNodes[7] = theDomain->getNode(mExternalNodes(7));

	for (int i = 0; i < 8; i++) {
		if (theNodes[i] == 0) {
			return;  // don't go any further - otherwise segmentation fault
		}
	}
	
	xi(0) = -0.125;
	xi(1) =  0.125;
	xi(2) =  0.125;
	xi(3) = -0.125;
	xi(4) = -0.125;
	xi(5) =  0.125;
	xi(6) =  0.125;
	xi(7) = -0.125;

	et(0) = -0.125;
	et(1) = -0.125;
	et(2) =  0.125;
	et(3) =  0.125;
	et(4) = -0.125;
	et(5) = -0.125;
	et(6) =  0.125;
	et(7) =  0.125;

	ze(0) = -0.125;
	ze(1) = -0.125;
	ze(2) = -0.125;
	ze(3) = -0.125;
	ze(4) =  0.125;
	ze(5) =  0.125;
	ze(6) =  0.125;
	ze(7) =  0.125;

	hst(0) =  0.125;
	hst(1) = -0.125;
	hst(2) =  0.125;
	hst(3) = -0.125;
	hst(4) =  0.125;
	hst(5) = -0.125;
	hst(6) =  0.125;
	hst(7) = -0.125;

	hut(0) =  0.125;
	hut(1) =  0.125;
	hut(2) = -0.125;
	hut(3) = -0.125;
	hut(4) = -0.125;
	hut(5) = -0.125;
	hut(6) =  0.125;
	hut(7) =  0.125;

	hus(0) =  0.125;
	hus(1) = -0.125;
	hus(2) = -0.125;
	hus(3) =  0.125;
	hus(4) = -0.125;
	hus(5) =  0.125;
	hus(6) =  0.125;
	hus(7) = -0.125;

	hstu(0) = -0.125;
	hstu(1) =  0.125;
	hstu(2) = -0.125;
	hstu(3) =  0.125;
	hstu(4) =  0.125;
	hstu(5) = -0.125;
	hstu(6) =  0.125;
	hstu(7) = -0.125;

	if (mInitialize) {
		Vector mIcrd_1(3);
    	Vector mIcrd_2(3);
    	Vector mIcrd_3(3);
    	Vector mIcrd_4(3);
    	Vector mIcrd_5(3);
    	Vector mIcrd_6(3);
    	Vector mIcrd_7(3);
    	Vector mIcrd_8(3);
    
    	// initialize coordinate vectors
    	mIcrd_1 = theNodes[0]->getCrds();
    	mIcrd_2 = theNodes[1]->getCrds();
    	mIcrd_3 = theNodes[2]->getCrds();
    	mIcrd_4 = theNodes[3]->getCrds();
    	mIcrd_5 = theNodes[4]->getCrds();
    	mIcrd_6 = theNodes[5]->getCrds();
    	mIcrd_7 = theNodes[6]->getCrds();
    	mIcrd_8 = theNodes[7]->getCrds();
    
    	// initialize coordinate matrix
    	mNodeCrd(0,0) = mIcrd_1(0);
    	mNodeCrd(1,0) = mIcrd_1(1);
    	mNodeCrd(2,0) = mIcrd_1(2);
    	mNodeCrd(0,1) = mIcrd_2(0);
    	mNodeCrd(1,1) = mIcrd_2(1);
    	mNodeCrd(2,1) = mIcrd_2(2);
    	mNodeCrd(0,2) = mIcrd_3(0);
    	mNodeCrd(1,2) = mIcrd_3(1);
    	mNodeCrd(2,2) = mIcrd_3(2);
    	mNodeCrd(0,3) = mIcrd_4(0);
    	mNodeCrd(1,3) = mIcrd_4(1);
    	mNodeCrd(2,3) = mIcrd_4(2);
    	mNodeCrd(0,4) = mIcrd_5(0);
    	mNodeCrd(1,4) = mIcrd_5(1);
    	mNodeCrd(2,4) = mIcrd_5(2);
    	mNodeCrd(0,5) = mIcrd_6(0);
    	mNodeCrd(1,5) = mIcrd_6(1);
    	mNodeCrd(2,5) = mIcrd_6(2);
    	mNodeCrd(0,6) = mIcrd_7(0);
    	mNodeCrd(1,6) = mIcrd_7(1);
    	mNodeCrd(2,6) = mIcrd_7(2);
    	mNodeCrd(0,7) = mIcrd_8(0);
    	mNodeCrd(1,7) = mIcrd_8(1);
    	mNodeCrd(2,7) = mIcrd_8(2);

	    // establish stabilization terms (based on initial state, only need to compute once)
	    GetStab();
	}

	// call the base-class method
	this->DomainComponent::setDomain(theDomain);

	/*
	this->update();
	const Matrix &t = this->getTangentStiff();
	for (int i=0; i<t.noRows(); i++) {
	  if (t(i,i) < 0.) {
	    opserr << "ERROR SSP_Brick bad tangent from element: " << this->getTag() << " ";
	    i = 100;
	  }
	}
	*/
}

int
SSPbrick::commitState(void)
{
	int retVal = 0;
	// call element commitState to do any base class stuff
	if ((retVal = this->Element::commitState()) != 0) {
		opserr << "SSPbrick::commitState() - failed in base class\n";
	}
	retVal = theMaterial->commitState();

	return retVal;
}

int
SSPbrick::revertToLastCommit(void)
{
	return theMaterial->revertToLastCommit();
}

int
SSPbrick::revertToStart(void)
{
	return theMaterial->revertToStart();
}

int
SSPbrick::update(void)
// this function updates variables for an incremental step n to n+1
{
	// get trial displacement
	const Vector &mDisp_1 = theNodes[0]->getTrialDisp();
	const Vector &mDisp_2 = theNodes[1]->getTrialDisp();
	const Vector &mDisp_3 = theNodes[2]->getTrialDisp();
	const Vector &mDisp_4 = theNodes[3]->getTrialDisp();
	const Vector &mDisp_5 = theNodes[4]->getTrialDisp();
	const Vector &mDisp_6 = theNodes[5]->getTrialDisp();
	const Vector &mDisp_7 = theNodes[6]->getTrialDisp();
	const Vector &mDisp_8 = theNodes[7]->getTrialDisp();
	
	// assemble displacement vector
	Vector u(24);
	u(0) =  mDisp_1(0);
	u(1) =  mDisp_1(1);
	u(2) =  mDisp_1(2);
	u(3) =  mDisp_2(0);
	u(4) =  mDisp_2(1);
	u(5) =  mDisp_2(2);
	u(6) =  mDisp_3(0);
	u(7) =  mDisp_3(1);
	u(8) =  mDisp_3(2);
	u(9) =  mDisp_4(0);
	u(10) = mDisp_4(1);
	u(11) = mDisp_4(2);
	u(12) = mDisp_5(0);
	u(13) = mDisp_5(1);
	u(14) = mDisp_5(2);
	u(15) = mDisp_6(0);
	u(16) = mDisp_6(1);
	u(17) = mDisp_6(2);
	u(18) = mDisp_7(0);
	u(19) = mDisp_7(1);
	u(20) = mDisp_7(2);
	u(21) = mDisp_8(0);
	u(22) = mDisp_8(1);
	u(23) = mDisp_8(2);

	// compute strain and send it to the material
	Vector strain(6);
	strain = Bnot*u;
	theMaterial->setTrialStrain(strain);

	return 0;
}

const Matrix &
SSPbrick::getTangentStiff(void)
// this function computes the tangent stiffness matrix for the element
{
	// get material tangent
	const Matrix &Cmat = theMaterial->getTangent();

	// full element stiffness matrix
	mTangentStiffness = Kstab;
	mTangentStiffness.addMatrixTripleProduct(1.0, Bnot, Cmat, mVol);
	
	return mTangentStiffness;
}

const Matrix &
SSPbrick::getInitialStiff(void)
// this function computes the initial tangent stiffness matrix for the element
{
	return getTangentStiff();
}

const Matrix &
SSPbrick::getMass(void)
{
	mMass.Zero();

	// get mass density from the material
	double density = theMaterial->getRho();

	// return zero matrix if density is zero
	if (density == 0.0) {
		return mMass;
	}
	
	// use jacobian determinant to get nodal mass values
	double massTerm;
	for (int i = 0; i < 8; i++) {
		massTerm = density*J[0]*(1.0 + (J[1]*xi(i) + J[2]*et(i) + J[3]*ze(i) + J[7] + J[8] + J[9])/3.0
                     + (J[4]*hut(i) + J[5]*hus(i) + J[6]*hst(i) + J[10]*ze(i) + J[11]*et(i) + J[12]*xi(i) + J[13]*ze(i) + J[14]*et(i) + J[15]*xi(i))/9.0
					 + (J[16]*hstu(i) + J[17]*hut(i) + J[18]*hus(i) + J[19]*hst(i))/27.0);
		mMass(3*i,3*i)     += massTerm;
		mMass(3*i+1,3*i+1) += massTerm;
		mMass(3*i+2,3*i+2) += massTerm;
	}

	return mMass;
}

void
SSPbrick::zeroLoad(void)
{
  applyLoad = 0;
  appliedB[0] = 0.0;
  appliedB[1] = 0.0;
  appliedB[2] = 0.0;
  
  Q.Zero();
  
  return;
}

int
SSPbrick::addLoad(ElementalLoad *theLoad, double loadFactor)
{
  // body forces can be applied in a load pattern
  int type;
  const Vector &data = theLoad->getData(type, loadFactor);
  
  if (type == LOAD_TAG_SelfWeight) {
    
    if (data.Size() != 3) {
      opserr << "SSPbrick::addLoad - too few SelfWeight data points, need 3 for ele " << this->getTag() << endln;
      return -1;
    }
    
    applyLoad = 1;
    appliedB[0] += loadFactor*data(0)*b[0];
    appliedB[1] += loadFactor*data(1)*b[1];
    appliedB[2] += loadFactor*data(2)*b[2];
    return 0;

  } else {
    opserr << "SSPbrick::addLoad - load type unknown for ele with tag: " << this->getTag() << endln;
    return -1;
  } 
  
  return -1;
}

int
SSPbrick::addInertiaLoadToUnbalance(const Vector &accel)
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
	const Vector &Raccel5 = theNodes[4]->getRV(accel);
	const Vector &Raccel6 = theNodes[5]->getRV(accel);
	const Vector &Raccel7 = theNodes[6]->getRV(accel);
	const Vector &Raccel8 = theNodes[7]->getRV(accel);

	static double ra[24];
	ra[0] =  Raccel1(0);
	ra[1] =  Raccel1(1);
	ra[2] =  Raccel1(2);
	ra[3] =  Raccel2(0);
	ra[4] =  Raccel2(1);
	ra[5] =  Raccel2(2);
	ra[6] =  Raccel3(0);
	ra[7] =  Raccel3(1);
	ra[8] =  Raccel3(2);
	ra[9] =  Raccel4(0);
	ra[10] = Raccel4(1);
	ra[11] = Raccel4(2);
	ra[12] = Raccel5(0);
	ra[13] = Raccel5(1);
	ra[14] = Raccel5(2);
	ra[15] = Raccel6(0);
	ra[16] = Raccel6(1);
	ra[17] = Raccel6(2);
	ra[18] = Raccel7(0);
	ra[19] = Raccel7(1);
	ra[20] = Raccel7(2);
	ra[21] = Raccel8(0);
	ra[22] = Raccel8(1);
	ra[23] = Raccel8(2);

	// compute mass matrix
	this->getMass();

	for (int i = 0; i < 24; i++) {
		Q(i) += -mMass(i,i)*ra[i];
	}
	
	return 0;
}

const Vector &
SSPbrick::getResistingForce(void)
// this function computes the resisting force vector for the element
{
	// get stress from the material
	Vector mStress = theMaterial->getStress();

	// get trial displacement
	const Vector &mDisp_1 = theNodes[0]->getTrialDisp();
	const Vector &mDisp_2 = theNodes[1]->getTrialDisp();
	const Vector &mDisp_3 = theNodes[2]->getTrialDisp();
	const Vector &mDisp_4 = theNodes[3]->getTrialDisp();
	const Vector &mDisp_5 = theNodes[4]->getTrialDisp();
	const Vector &mDisp_6 = theNodes[5]->getTrialDisp();
	const Vector &mDisp_7 = theNodes[6]->getTrialDisp();
	const Vector &mDisp_8 = theNodes[7]->getTrialDisp();
	
	// assemble displacement vector
	Vector d(24);
	d(0) =  mDisp_1(0);
	d(1) =  mDisp_1(1);
	d(2) =  mDisp_1(2);
	d(3) =  mDisp_2(0);
	d(4) =  mDisp_2(1);
	d(5) =  mDisp_2(2);
	d(6) =  mDisp_3(0);
	d(7) =  mDisp_3(1);
	d(8) =  mDisp_3(2);
	d(9) =  mDisp_4(0);
	d(10) = mDisp_4(1);
	d(11) = mDisp_4(2);
	d(12) = mDisp_5(0);
	d(13) = mDisp_5(1);
	d(14) = mDisp_5(2);
	d(15) = mDisp_6(0);
	d(16) = mDisp_6(1);
	d(17) = mDisp_6(2);
	d(18) = mDisp_7(0);
	d(19) = mDisp_7(1);
	d(20) = mDisp_7(2);
	d(21) = mDisp_8(0);
	d(22) = mDisp_8(1);
	d(23) = mDisp_8(2);

	// add stabilization force to internal force vector
	mInternalForces = Kstab*d;

	// add internal force from the stress  ->  fint = Kstab*d + 8*Jo*Bnot'*stress
	mInternalForces.addMatrixTransposeVector(1.0, Bnot, mStress, mVol);

	// subtract body forces from internal force vector
	Vector body(3);
	if (applyLoad == 0) {
		double polyJac = 0.0;
		for (int i = 0; i < 8; i++) {
			polyJac = J[0]*(1.0 + (J[1]*xi(i) + J[2]*et(i) + J[3]*ze(i) + J[7] + J[8] + J[9])/3.0
                     + (J[4]*hut(i) + J[5]*hus(i) + J[6]*hst(i) + J[10]*ze(i) + J[11]*et(i) + J[12]*xi(i) + J[13]*ze(i) + J[14]*et(i) + J[15]*xi(i))/9.0
					 + (J[16]*hstu(i) + J[17]*hut(i) + J[18]*hus(i) + J[19]*hst(i))/27.0);
			mInternalForces(3*i)   -= b[0]*polyJac;
			mInternalForces(3*i+1) -= b[1]*polyJac;
			mInternalForces(3*i+2) -= b[2]*polyJac;
		}
	} else {
		double polyJac = 0.0;
		for (int i = 0; i < 8; i++) {
			polyJac = J[0]*(1.0 + (J[1]*xi(i) + J[2]*et(i) + J[3]*ze(i) + J[7] + J[8] + J[9])/3.0
                     + (J[4]*hut(i) + J[5]*hus(i) + J[6]*hst(i) + J[10]*ze(i) + J[11]*et(i) + J[12]*xi(i) + J[13]*ze(i) + J[14]*et(i) + J[15]*xi(i))/9.0
					 + (J[16]*hstu(i) + J[17]*hut(i) + J[18]*hus(i) + J[19]*hst(i))/27.0);
			mInternalForces(3*i)   -= appliedB[0]*polyJac;
			mInternalForces(3*i+1) -= appliedB[1]*polyJac;
			mInternalForces(3*i+2) -= appliedB[2]*polyJac;
		}
	}

	// inertial unbalance load
	mInternalForces.addVector(1.0, Q, -1.0);

	return mInternalForces;
}

const Vector &
SSPbrick::getResistingForceIncInertia()
{
  // get mass density from the material
  double density = theMaterial->getRho();
  

  // if density is zero only add damping terms
  if (density == 0.0) {
    this->getResistingForce();
    
    // add the damping forces if rayleigh damping
    if (betaK != 0.0 || betaK0 != 0.0 || betaKc != 0.0) {
      mInternalForces += this->getRayleighDampingForces();
    }
    
    return mInternalForces;
  }
  
  const Vector &accel1 = theNodes[0]->getTrialAccel();
  const Vector &accel2 = theNodes[1]->getTrialAccel();
  const Vector &accel3 = theNodes[2]->getTrialAccel();
  const Vector &accel4 = theNodes[3]->getTrialAccel();
  const Vector &accel5 = theNodes[4]->getTrialAccel();
  const Vector &accel6 = theNodes[5]->getTrialAccel();
  const Vector &accel7 = theNodes[6]->getTrialAccel();
  const Vector &accel8 = theNodes[7]->getTrialAccel();
  
  static double a[24];
  a[0] =  accel1(0);
  a[1] =  accel1(1);
  a[2] =  accel1(2);
  a[3] =  accel2(0);
  a[4] =  accel2(1);
  a[5] =  accel2(2);
  a[6] =  accel3(0);
  a[7] =  accel3(1);
  a[8] =  accel3(2);
  a[9] =  accel4(0);
  a[10] = accel4(1);
  a[11] = accel4(2);
  a[12] = accel5(0);
  a[13] = accel5(1);
  a[14] = accel5(2);
  a[15] = accel6(0);
  a[16] = accel6(1);
  a[17] = accel6(2);
  a[18] = accel7(0);
  a[19] = accel7(1);
  a[20] = accel7(2);
  a[21] = accel8(0);
  a[22] = accel8(1);
  a[23] = accel8(2);
  
  // compute current resisting force
  this->getResistingForce();
  
  // compute mass matrix
  this->getMass();
  
  for (int i = 0; i < 24; i++) {
    mInternalForces(i) += mMass(i,i)*a[i];
  }
  
  // add the damping forces if rayleigh damping
  if (alphaM != 0 || betaK != 0.0 || betaK0 != 0.0 || betaKc != 0.0) {
    mInternalForces += this->getRayleighDampingForces(); //FMK - this would add mass contribution again!
  }
  
  return mInternalForces;
}

int
SSPbrick::sendSelf(int commitTag, Channel &theChannel)
{
  int res = 0;
  
  // note: we don't check for dataTag == 0 for Element
  // objects as that is taken care of in a commit by the Domain
  // object - don't want to have to do the check if sending data
  int dataTag = this->getDbTag();
  
  // SSPbrick packs its data into a Vector and sends this to theChannel
  // along with its dbTag and the commitTag passed in the arguments
  static Vector data(751);
  data(0) = this->getTag();
  data(1) = b[0];
  data(2) = b[1];
  data(3) = b[2];
  data(4) = mVol;
  data(5) = theMaterial->getClassTag();

  int matDbTag = theMaterial->getDbTag();
  // NOTE: we have to ensure that the material has a database tag if we are sending to a database channel
  if (matDbTag == 0) {
    matDbTag = theChannel.getDbTag();
    if (matDbTag != 0)
      theMaterial->setDbTag(matDbTag);
  }
  data(6) = matDbTag;

  data(7) = alphaM;
  data(8) = betaK;
  data(9) = betaK0;
  data(10) = betaKc;
 
  int cnt = 11;
  for (int i = 0; i < 20; i++) {
    data(cnt+i) = J[i];
  }
  
  cnt = 31;
  for (int i = 0; i < 6; i++) {
    for (int j = 0; j < 24; j++) {
      data(cnt+j) = Bnot(i,j);
    }
    cnt = cnt+24;
  }
  
  cnt = 175;
  for (int i = 0; i < 24; i++) {
    for (int j = 0; j < 24; j++) {
      data(cnt+j) = Kstab(i,j);
    }
    cnt = cnt+24;
  }
  
  res = theChannel.sendVector(dataTag, commitTag, data);
  if (res < 0) {
    opserr << "WARNING SSPbrick::sendSelf() - " << this->getTag() << " failed to send Vector\n";
    return res;
  }
  
  // SSPbrick then sends the tags of its eight nodes
  res = theChannel.sendID(dataTag, commitTag, mExternalNodes);
  if (res < 0) {
    opserr << "WARNING SSPbrick::sendSelf() - " << this->getTag() << " failed to send ID\n";
    return res;
  }
  
  // finally, SSPbrick asks its material object to send itself
  res = theMaterial->sendSelf(commitTag, theChannel);
  if (res < 0) {
    opserr << "WARNING SSPbrick::sendSelf() - " << this->getTag() << " failed to send its Material\n";
    return -3;
  }
  
  return 0;
}

int
SSPbrick::recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker)
{
  int res = 0;
  int dataTag = this->getDbTag();
  
  // SSPbrick creates a Vector, receives the Vector and then sets the 
  // internal data with the data in the Vector
  static Vector data(751);
  res = theChannel.recvVector(dataTag, commitTag, data);
  if (res < 0) {
    opserr << "WARNING SSPbrick::recvSelf() - failed to receive Vector\n";
    return res;
  }
  
  this->setTag((int)data(0));
  b[0] = data(1);
  b[1] = data(2);
  b[2] = data(3);
  mVol = data(4);
 
  alphaM = data(7);
  betaK = data(8);
  betaK0 = data(9);
  betaKc = data(10);
 
  int cnt = 11;
  for (int i = 0; i < 20; i++) {
    J[i] = data(cnt+i);
  }
  
  cnt = 31;
  for (int i = 0; i < 6; i++) {
    for (int j = 0; j < 24; j++) {
      Bnot(i,j) = data(cnt+j);
    }
    cnt = cnt+24;
  }
  
  cnt = 175;
  for (int i = 0; i < 24; i++) {
    for (int j = 0; j < 24; j++) {
      Kstab(i,j) = data(cnt+j);
    }
    cnt = cnt+24;
  }
  
  // SSPbrick now receives the tags of its four external nodes
  res = theChannel.recvID(dataTag, commitTag, mExternalNodes);
  if (res < 0) {
    opserr << "WARNING SSPbrick::recvSelf() - " << this->getTag() << " failed to receive ID\n";
    return res;
  }
  
  // finally, SSPbrick creates a material object of the correct type, sets its
  // database tag, and asks this new object to receive itself
  int matClass = (int)data(5);
  int matDb    = (int)data(6);
  
  // check if material object exists and that it is the right type
  if ((theMaterial == 0) || (theMaterial->getClassTag() != matClass)) {
    
    // if old one, delete it
    if (theMaterial != 0)
      delete theMaterial;
    
    // create new material object
    NDMaterial *theMatCopy = theBroker.getNewNDMaterial(matClass);
    theMaterial = (NDMaterial *)theMatCopy;
    
    if (theMaterial == 0) {
      opserr << "WARNING SSPbrick::recvSelf() - " << this->getTag() 
	     << " failed to get a blank Material of type " << matClass << endln;
      return -3;
    }
  }
  
  // NOTE: we set the dbTag before we receive the material
  theMaterial->setDbTag(matDb);
  res = theMaterial->recvSelf(commitTag, theChannel, theBroker);
  if (res < 0) {
    opserr << "WARNING SSPbrick::recvSelf() - " << this->getTag() << " failed to receive its Material\n";
    return -3;
  }
  
  return 0; 
}

int
SSPbrick::displaySelf(Renderer &theViewer, int displayMode, float fact, const char **modes, int numMode)
{
	// vertex display coordinate vectors
	static Vector v1(3);
	static Vector v2(3);
	static Vector v3(3);
	static Vector v4(3);
	static Vector v5(3);
	static Vector v6(3);
	static Vector v7(3);
	static Vector v8(3);
	theNodes[0]->getDisplayCrds(v1, fact, displayMode);
	theNodes[1]->getDisplayCrds(v2, fact, displayMode);
	theNodes[2]->getDisplayCrds(v3, fact, displayMode);
	theNodes[3]->getDisplayCrds(v4, fact, displayMode);
	theNodes[4]->getDisplayCrds(v5, fact, displayMode);
	theNodes[5]->getDisplayCrds(v6, fact, displayMode);
	theNodes[6]->getDisplayCrds(v7, fact, displayMode);
	theNodes[7]->getDisplayCrds(v8, fact, displayMode);

	// add to coord matrix
	static Matrix coords(8, 3);
	int i;
	for (i = 0; i < 3; i++) {
		coords(0, i) = v1(i);
		coords(1, i) = v2(i);
		coords(2, i) = v3(i);
		coords(3, i) = v4(i);
		coords(4, i) = v5(i);
		coords(5, i) = v6(i);
		coords(6, i) = v7(i);
		coords(7, i) = v8(i);
	}

	// get color vector
	static Vector values(8);
	for (i = 0; i < 8; i++)
		values(i) = 1.0;

	// draw cube
	return theViewer.drawCube(coords, values, this->getTag());
}

void
SSPbrick::Print(OPS_Stream &s, int flag)
{
    if (flag == OPS_PRINT_CURRENTSTATE) {
        opserr << "SSPbrick, element id:  " << this->getTag() << endln;
        opserr << "   Connected external nodes:  ";
        for (int i = 0; i < SSPB_NUM_NODE; i++) {
            opserr << mExternalNodes(i) << " ";
        }
        opserr << endln;
    }
    
    if (flag == OPS_PRINT_PRINTMODEL_JSON) {
        s << "\t\t\t{";
        s << "\"name\": " << this->getTag() << ", ";
        s << "\"type\": \"SSPbrick\", ";
        s << "\"nodes\": [" << mExternalNodes(0) << ", ";
        for (int i = 1; i < 7; i++)
            s << mExternalNodes(i) << ", ";
        s << mExternalNodes(7) << "], ";
        s << "\"bodyForces\": [" << b[0] << ", " << b[1] << ", " << b[2] << "], ";
        s << "\"material\": \"" << theMaterial->getTag() << "\"}";
    }
}

Response*
SSPbrick::setResponse(const char **argv, int argc, OPS_Stream &eleInfo)
{
	// no special recorders for this element, call the method in the material class
	return theMaterial->setResponse(argv, argc, eleInfo);
}

int
SSPbrick::getResponse(int responseID, Information &eleInfo)
{
	// no special recorders for this element, call the method in the material class
	return theMaterial->getResponse(responseID, eleInfo);
}

int
SSPbrick::setParameter(const char **argv, int argc, Parameter &param)
{
	if (argc < 1) {
		return -1;
	}

    int res = -1;

    // no element parameters, call setParameter in the material
    int matRes;
    matRes = theMaterial->setParameter(argv, argc, param);

    if (matRes != -1) {
		res = matRes;
	}
  
  return res;
}

int
SSPbrick::updateParameter(int parameterID, Information &info)
{
	int res = -1;
	int matRes = res;

    if (parameterID == res) {
        return -1;
    } else {
        matRes = theMaterial->updateParameter(parameterID, info);
		if (matRes != -1) {
			res = matRes;
		}
		return res;
    }    
}

void
SSPbrick::GetStab(void)
// this function computes the stabilization stiffness matrix for the element
{
    Matrix Mben(12,24);
	Matrix FCF(12,12);
	Matrix dNloc(8,3);
	Matrix dNmod(8,3);
	Matrix Jmat(3,3);
	Matrix Jinv(3,3);
	Matrix G(8,8);
	Matrix I8(8,8);
	
	// define local coord and hourglass vectors
	Vector gst(8);
	Vector gut(8);
	Vector gus(8);
	Vector gstu(8);

	// shape function derivatives (local crd) at center
	dNloc(0,0) = -0.125;
	dNloc(1,0) =  0.125;
	dNloc(2,0) =  0.125;
	dNloc(3,0) = -0.125;
	dNloc(4,0) = -0.125;
	dNloc(5,0) =  0.125;
	dNloc(6,0) =  0.125;
	dNloc(7,0) = -0.125;
	dNloc(0,1) = -0.125;
	dNloc(1,1) = -0.125;
	dNloc(2,1) =  0.125;
	dNloc(3,1) =  0.125;
	dNloc(4,1) = -0.125;
	dNloc(5,1) = -0.125;
	dNloc(6,1) =  0.125;
	dNloc(7,1) =  0.125;
	dNloc(0,2) = -0.125;
	dNloc(1,2) = -0.125;
	dNloc(2,2) = -0.125;
	dNloc(3,2) = -0.125;
	dNloc(4,2) =  0.125;
	dNloc(5,2) =  0.125;
	dNloc(6,2) =  0.125;
	dNloc(7,2) =  0.125;

	// jacobian matrix
	Jmat = mNodeCrd*dNloc;
	// inverse of jacobian matrix
	Jmat.Invert(Jinv);

	// nodal coordinate vectors
	Vector x(8);
	Vector y(8);
	Vector z(8);

	x(0) = mNodeCrd(0,0);
	x(1) = mNodeCrd(0,1);
	x(2) = mNodeCrd(0,2);
	x(3) = mNodeCrd(0,3);
	x(4) = mNodeCrd(0,4);
	x(5) = mNodeCrd(0,5);
	x(6) = mNodeCrd(0,6);
	x(7) = mNodeCrd(0,7);

	y(0) = mNodeCrd(1,0);
	y(1) = mNodeCrd(1,1);
	y(2) = mNodeCrd(1,2);
	y(3) = mNodeCrd(1,3);
	y(4) = mNodeCrd(1,4);
	y(5) = mNodeCrd(1,5);
	y(6) = mNodeCrd(1,6);
	y(7) = mNodeCrd(1,7);

	z(0) = mNodeCrd(2,0);
	z(1) = mNodeCrd(2,1);
	z(2) = mNodeCrd(2,2);
	z(3) = mNodeCrd(2,3);
	z(4) = mNodeCrd(2,4);
	z(5) = mNodeCrd(2,5);
	z(6) = mNodeCrd(2,6);
	z(7) = mNodeCrd(2,7);

	// define coefficient terms for jacobian determinant
    double 	a1 = x^xi;
    double 	a2 = x^et;
    double 	a3 = x^ze;
    double 	a4 = x^hut;
    double 	a5 = x^hus;
    double 	a6 = x^hst;
    double 	a7 = x^hstu;

    double 	b1 = y^xi;
    double 	b2 = y^et;
    double 	b3 = y^ze;
    double 	b4 = y^hut;
    double 	b5 = y^hus;
    double 	b6 = y^hst;
    double 	b7 = y^hstu;

    double 	c1 = z^xi;
    double 	c2 = z^et;
    double 	c3 = z^ze;
    double 	c4 = z^hut;
    double 	c5 = z^hus;
    double 	c6 = z^hst;
    double 	c7 = z^hstu;

	// define coefficient vectors for jacobian determinant
	Vector e1(3);
	Vector e2(3);
	Vector e3(3);
	Vector e4(3);
	Vector e5(3);
	Vector e6(3);
	Vector e7(3);

	e1(0) = a1;  e1(1) = b1;  e1(2) = c1;
	e2(0) = a2;  e2(1) = b2;  e2(2) = c2;
	e3(0) = a3;  e3(1) = b3;  e3(2) = c3;
	e4(0) = a4;  e4(1) = b4;  e4(2) = c4;
	e5(0) = a5;  e5(1) = b5;  e5(2) = c5;
	e6(0) = a6;  e6(1) = b6;  e6(2) = c6;
	e7(0) = a7;  e7(1) = b7;  e7(2) = c7;

	// jacobian determinant terms
	J[0] = e1^(CrossProduct(e2,e3));
	J[1] = (e1^(CrossProduct(e2,e5))) + (e1^(CrossProduct(e6,e3)));
    J[2] = (e1^(CrossProduct(e2,e4))) + (e6^(CrossProduct(e2,e3)));
    J[3] = (e5^(CrossProduct(e2,e3))) + (e1^(CrossProduct(e4,e3)));
	J[1] = 0.0;
	J[2] = 0.0;
	J[3] = 0.0;
    J[4] = (e7^(CrossProduct(e2,e3))) + (e4^(CrossProduct(e5,e2))) + (e4^(CrossProduct(e3,e6)));
    J[5] = (e1^(CrossProduct(e7,e3))) + (e4^(CrossProduct(e5,e1))) + (e3^(CrossProduct(e5,e6)));
    J[6] = (e1^(CrossProduct(e2,e7))) + (e4^(CrossProduct(e1,e6))) + (e2^(CrossProduct(e5,e6)));
    J[7] = -1.0*e1^(CrossProduct(e5,e6));
    J[8] = -1.0*e4^(CrossProduct(e2,e6));
    J[9] = -1.0*e4^(CrossProduct(e5,e3));
    J[10] = e2^(CrossProduct(e4,e7));
    J[11] = -1.0*e3^(CrossProduct(e4,e7));
    J[12] = e3^(CrossProduct(e5,e7));
    J[13] = -1.0*e1^(CrossProduct(e5,e7));
    J[14] = e1^(CrossProduct(e6,e7));
    J[15] = -1.0*e2^(CrossProduct(e6,e7));
    J[16] = 2.0*e4^(CrossProduct(e5,e6));
    J[17] = e7^(CrossProduct(e5,e6));
    J[18] = e4^(CrossProduct(e7,e6));
    J[19] = e4^(CrossProduct(e5,e7));

	// combined jacobian terms
    double J0789  = 8.0*(J[0]/3.0 + J[7]/5.0 + J[8]/9.0 + J[9]/9.0);
    double J0879  = 8.0*(J[0]/3.0 + J[8]/5.0 + J[7]/9.0 + J[9]/9.0);
    double J0978  = 8.0*(J[0]/3.0 + J[9]/5.0 + J[7]/9.0 + J[8]/9.0);
    double J417   = 8.0*(J[4]/9.0 + J[17]/27.0);
    double J518   = 8.0*(J[5]/9.0 + J[18]/27.0);
    double J619   = 8.0*(J[6]/9.0 + J[19]/27.0);
    double J11215 = 8.0*(J[1]/9.0 + J[12]/15.0 + J[15]/27.0);
    double J11512 = 8.0*(J[1]/9.0 + J[15]/15.0 + J[12]/27.0);
    double J21114 = 8.0*(J[2]/9.0 + J[11]/15.0 + J[14]/27.0);
    double J21411 = 8.0*(J[2]/9.0 + J[14]/15.0 + J[11]/27.0);
    double J31013 = 8.0*(J[3]/9.0 + J[10]/15.0 + J[13]/27.0);
    double J31310 = 8.0*(J[3]/9.0 + J[13]/15.0 + J[10]/27.0);
	double J789   = 8.0*(J[0]/9.0 + J[7]/15.0 + J[8]/15.0 + J[9]/27.0);
    double J897   = 8.0*(J[0]/9.0 + J[8]/15.0 + J[9]/15.0 + J[7]/27.0);
    double J798   = 8.0*(J[0]/9.0 + J[7]/15.0 + J[9]/15.0 + J[8]/27.0);
    double J174   = 8.0*(J[4]/27.0 + 64.0*J[17]/45.0);
    double J185   = 8.0*(J[5]/27.0 + 64.0*J[18]/45.0);
    double J196   = 8.0*(J[6]/27.0 + 64.0*J[19]/45.0);
	double J16    = 8.0*J[16]/27.0;

	// compute element volume 
	mVol = 8.0*(J[0] + (J[7] + J[8] + J[9])/3.0);

	// kinematic vectors 
	Vector bx(8);
	Vector by(8);
	Vector bz(8);

	bx = (8.0*((b2*c3-c2*b3)*xi + (b3*c1-c3*b1)*et + (b1*c2-c1*b2)*ze) + (8.0/3.0)*((b6*c5-c6*b5)*xi + (b4*c6-c4*b6)*et + (b5*c4-c5*b4)*ze 
	     + (b5*c1-c5*b1 + b2*c4-c2*b4)*hst + (b6*c2-c6*b2 + b3*c5-c3*b5)*hut + (b1*c6-c1*b6 + b4*c3-c4*b3)*hus))/mVol;
	by = (8.0*((c2*a3-a2*c3)*xi + (c3*a1-a3*c1)*et + (c1*a2-a1*c2)*ze) + (8.0/3.0)*((c6*a5-a6*c5)*xi + (c4*a6-a4*c6)*et + (c5*a4-a5*c4)*ze 
	     + (c5*a1-a5*c1 + c2*a4-a2*c4)*hst + (c6*a2-a6*c2 + c3*a5-a3*c5)*hut + (c1*a6-a1*c6 + c4*a3-a4*c3)*hus))/mVol;
	bz = (8.0*((a2*b3-b2*a3)*xi + (a3*b1-b3*a1)*et + (a1*b2-b1*a2)*ze) + (8.0/3.0)*((a6*b5-b6*a5)*xi + (a4*b6-b4*a6)*et + (a5*b4-b5*a4)*ze 
	     + (a5*b1-b5*a1 + a2*b4-b2*a4)*hst + (a6*b2-b6*a2 + a3*b5-b3*a5)*hut + (a1*b6-b1*a6 + a4*b3-b4*a3)*hus))/mVol;

	for (int i = 0; i < 8; i++) {
		dNmod(i,0) = bx(i);
		dNmod(i,1) = by(i);
		dNmod(i,2) = bz(i);
	}

	// compute hourglass transformation matrix G
	I8.Zero();
	I8(0,0) = 1.0;
	I8(1,1) = 1.0;
	I8(2,2) = 1.0;
	I8(3,3) = 1.0;
	I8(4,4) = 1.0;
	I8(5,5) = 1.0;
	I8(6,6) = 1.0;
	I8(7,7) = 1.0;

	G = I8 - dNmod*mNodeCrd;

	// compute gamma vectors
	gst =  G*hst;
	gut =  G*hut;
	gus =  G*hus;
	gstu = G*hstu;

	// define kinematic and mapping matrices
	Bnot.Zero();
	Mben.Zero();
	for (int i = 0; i < 8; i++) {
		Bnot(0,3*i)   = dNmod(i,0);
    	Bnot(1,3*i+1) = dNmod(i,1);
    	Bnot(2,3*i+2) = dNmod(i,2);
    	Bnot(3,3*i)   = dNmod(i,1);
    	Bnot(3,3*i+1) = dNmod(i,0);
    	Bnot(4,3*i+1) = dNmod(i,2);
    	Bnot(4,3*i+2) = dNmod(i,1);
    	Bnot(5,3*i)   = dNmod(i,2);
    	Bnot(5,3*i+2) = dNmod(i,0);

		Mben(0,3*i)   = gst(i);
		Mben(1,3*i+1) = gst(i);
		Mben(2,3*i+2) = gst(i);
		Mben(3,3*i)   = gut(i);
		Mben(4,3*i+1) = gut(i);
		Mben(5,3*i+2) = gut(i);
		Mben(6,3*i)   = gus(i);
		Mben(7,3*i+1) = gus(i);
		Mben(8,3*i+2) = gus(i);
		Mben(9,3*i)   = gstu(i);
		Mben(10,3*i+1) = gstu(i);
		Mben(11,3*i+2) = gstu(i);
	}
	
	// define terms for FCF matrix 
    double HstXX = J0879*Jinv(0,0)*Jinv(0,0) + J0789*Jinv(1,0)*Jinv(1,0) + J619*(Jinv(0,0)*Jinv(1,0) + Jinv(1,0)*Jinv(0,0));
    double HstXY = J0879*Jinv(0,0)*Jinv(0,1) + J0789*Jinv(1,0)*Jinv(1,1) + J619*(Jinv(0,0)*Jinv(1,1) + Jinv(1,0)*Jinv(0,1));
    double HstXZ = J0879*Jinv(0,0)*Jinv(0,2) + J0789*Jinv(1,0)*Jinv(1,2) + J619*(Jinv(0,0)*Jinv(1,2) + Jinv(1,0)*Jinv(0,2));
    double HstYY = J0879*Jinv(0,1)*Jinv(0,1) + J0789*Jinv(1,1)*Jinv(1,1) + J619*(Jinv(0,1)*Jinv(1,1) + Jinv(1,1)*Jinv(0,1));
    double HstYZ = J0879*Jinv(0,1)*Jinv(0,2) + J0789*Jinv(1,1)*Jinv(1,2) + J619*(Jinv(0,1)*Jinv(1,2) + Jinv(1,1)*Jinv(0,2));
    double HstZZ = J0879*Jinv(0,2)*Jinv(0,2) + J0789*Jinv(1,2)*Jinv(1,2) + J619*(Jinv(0,2)*Jinv(1,2) + Jinv(1,2)*Jinv(0,2));

    double HutXX = J0978*Jinv(1,0)*Jinv(1,0) + J0879*Jinv(2,0)*Jinv(2,0) + J417*(Jinv(1,0)*Jinv(2,0) + Jinv(2,0)*Jinv(1,0));
    double HutXY = J0978*Jinv(1,0)*Jinv(1,1) + J0879*Jinv(2,0)*Jinv(2,1) + J417*(Jinv(1,0)*Jinv(2,1) + Jinv(2,0)*Jinv(1,1));
    double HutXZ = J0978*Jinv(1,0)*Jinv(1,2) + J0879*Jinv(2,0)*Jinv(2,2) + J417*(Jinv(1,0)*Jinv(2,2) + Jinv(2,0)*Jinv(1,2));
    double HutYY = J0978*Jinv(1,1)*Jinv(1,1) + J0879*Jinv(2,1)*Jinv(2,1) + J417*(Jinv(1,1)*Jinv(2,1) + Jinv(2,1)*Jinv(1,1));
    double HutYZ = J0978*Jinv(1,1)*Jinv(1,2) + J0879*Jinv(2,1)*Jinv(2,2) + J417*(Jinv(1,1)*Jinv(2,2) + Jinv(2,1)*Jinv(1,2));
    double HutZZ = J0978*Jinv(1,2)*Jinv(1,2) + J0879*Jinv(2,2)*Jinv(2,2) + J417*(Jinv(1,2)*Jinv(2,2) + Jinv(2,2)*Jinv(1,2));

    double HusXX = J0978*Jinv(0,0)*Jinv(0,0) + J0789*Jinv(2,0)*Jinv(2,0) + J518*(Jinv(0,0)*Jinv(2,0) + Jinv(2,0)*Jinv(0,0));
    double HusXY = J0978*Jinv(0,0)*Jinv(0,1) + J0789*Jinv(2,0)*Jinv(2,1) + J518*(Jinv(0,0)*Jinv(2,1) + Jinv(2,0)*Jinv(0,1));
    double HusXZ = J0978*Jinv(0,0)*Jinv(0,2) + J0789*Jinv(2,0)*Jinv(2,2) + J518*(Jinv(0,0)*Jinv(2,2) + Jinv(2,0)*Jinv(0,2));
    double HusYY = J0978*Jinv(0,1)*Jinv(0,1) + J0789*Jinv(2,1)*Jinv(2,1) + J518*(Jinv(0,1)*Jinv(2,1) + Jinv(2,1)*Jinv(0,1));
    double HusYZ = J0978*Jinv(0,1)*Jinv(0,2) + J0789*Jinv(2,1)*Jinv(2,2) + J518*(Jinv(0,1)*Jinv(2,2) + Jinv(2,1)*Jinv(0,2));
    double HusZZ = J0978*Jinv(0,2)*Jinv(0,2) + J0789*Jinv(2,2)*Jinv(2,2) + J518*(Jinv(0,2)*Jinv(2,2) + Jinv(2,2)*Jinv(0,2));

    double HstuXX = J897*Jinv(0,0)*Jinv(0,0) + J798*Jinv(1,0)*Jinv(1,0) + J789*Jinv(2,0)*Jinv(2,0) + J185*(Jinv(0,0)*Jinv(2,0) + Jinv(2,0)*Jinv(0,0)) 
	                + J196*(Jinv(0,0)*Jinv(1,0) + Jinv(1,0)*Jinv(0,0)) + J174*(Jinv(1,0)*Jinv(2,0) + Jinv(2,0)*Jinv(1,0));
    double HstuXY = J897*Jinv(0,0)*Jinv(0,1) + J798*Jinv(1,0)*Jinv(1,1) + J789*Jinv(2,0)*Jinv(2,1) + J185*(Jinv(0,0)*Jinv(2,1) + Jinv(2,0)*Jinv(0,1)) 
	                + J196*(Jinv(0,0)*Jinv(1,1) + Jinv(1,0)*Jinv(0,1)) + J174*(Jinv(1,0)*Jinv(2,1) + Jinv(2,0)*Jinv(1,1));
    double HstuXZ = J897*Jinv(0,0)*Jinv(0,2) + J798*Jinv(1,0)*Jinv(1,2) + J789*Jinv(2,0)*Jinv(2,2) + J185*(Jinv(0,0)*Jinv(2,2) + Jinv(2,0)*Jinv(0,2)) 
	                + J196*(Jinv(0,0)*Jinv(1,2) + Jinv(1,0)*Jinv(0,2)) + J174*(Jinv(1,0)*Jinv(2,2) + Jinv(2,0)*Jinv(1,2));
    double HstuYY = J897*Jinv(0,1)*Jinv(0,1) + J798*Jinv(1,1)*Jinv(1,1) + J789*Jinv(2,1)*Jinv(2,1) + J185*(Jinv(0,1)*Jinv(2,1) + Jinv(2,1)*Jinv(0,1)) 
	                + J196*(Jinv(0,1)*Jinv(1,1) + Jinv(1,1)*Jinv(0,1)) + J174*(Jinv(1,1)*Jinv(2,1) + Jinv(2,1)*Jinv(1,1));
    double HstuYZ = J897*Jinv(0,1)*Jinv(0,2) + J798*Jinv(1,1)*Jinv(1,2) + J789*Jinv(2,1)*Jinv(2,2) + J185*(Jinv(0,1)*Jinv(2,2) + Jinv(2,1)*Jinv(0,2)) 
	                + J196*(Jinv(0,1)*Jinv(1,2) + Jinv(1,1)*Jinv(0,2)) + J174*(Jinv(1,1)*Jinv(2,2) + Jinv(2,1)*Jinv(1,2));
    double HstuZZ = J897*Jinv(0,2)*Jinv(0,2) + J798*Jinv(1,2)*Jinv(1,2) + J789*Jinv(2,2)*Jinv(2,2) + J185*(Jinv(0,2)*Jinv(2,2) + Jinv(2,2)*Jinv(0,2)) 
	                + J196*(Jinv(0,2)*Jinv(1,2) + Jinv(1,2)*Jinv(0,2)) + J174*(Jinv(1,2)*Jinv(2,2) + Jinv(2,2)*Jinv(1,2));

    double IttXX = J0879*Jinv(0,0)*Jinv(2,0) + J417*Jinv(0,0)*Jinv(1,0) + J518*Jinv(1,0)*Jinv(1,0) + J619*Jinv(1,0)*Jinv(2,0);
    double IttXY = J0879*Jinv(0,0)*Jinv(2,1) + J417*Jinv(0,0)*Jinv(1,1) + J518*Jinv(1,0)*Jinv(1,1) + J619*Jinv(1,0)*Jinv(2,1);
    double IttXZ = J0879*Jinv(0,0)*Jinv(2,2) + J417*Jinv(0,0)*Jinv(1,2) + J518*Jinv(1,0)*Jinv(1,2) + J619*Jinv(1,0)*Jinv(2,2);
    double IttYX = J0879*Jinv(0,1)*Jinv(2,0) + J417*Jinv(0,1)*Jinv(1,0) + J518*Jinv(1,1)*Jinv(1,0) + J619*Jinv(1,1)*Jinv(2,0);
    double IttYY = J0879*Jinv(0,1)*Jinv(2,1) + J417*Jinv(0,1)*Jinv(1,1) + J518*Jinv(1,1)*Jinv(1,1) + J619*Jinv(1,1)*Jinv(2,1);
    double IttYZ = J0879*Jinv(0,1)*Jinv(2,2) + J417*Jinv(0,1)*Jinv(1,2) + J518*Jinv(1,1)*Jinv(1,2) + J619*Jinv(1,1)*Jinv(2,2);
    double IttZX = J0879*Jinv(0,2)*Jinv(2,0) + J417*Jinv(0,2)*Jinv(1,0) + J518*Jinv(1,2)*Jinv(1,0) + J619*Jinv(1,2)*Jinv(2,0);
    double IttZY = J0879*Jinv(0,2)*Jinv(2,1) + J417*Jinv(0,2)*Jinv(1,1) + J518*Jinv(1,2)*Jinv(1,1) + J619*Jinv(1,2)*Jinv(2,1);
    double IttZZ = J0879*Jinv(0,2)*Jinv(2,2) + J417*Jinv(0,2)*Jinv(1,2) + J518*Jinv(1,2)*Jinv(1,2) + J619*Jinv(1,2)*Jinv(2,2);

    double IssXX = J0789*Jinv(1,0)*Jinv(2,0) + J417*Jinv(0,0)*Jinv(0,0) + J518*Jinv(1,0)*Jinv(0,0) + J619*Jinv(0,0)*Jinv(2,0);
    double IssXY = J0789*Jinv(1,0)*Jinv(2,1) + J417*Jinv(0,0)*Jinv(0,1) + J518*Jinv(1,0)*Jinv(0,1) + J619*Jinv(0,0)*Jinv(2,1);
    double IssXZ = J0789*Jinv(1,0)*Jinv(2,2) + J417*Jinv(0,0)*Jinv(0,2) + J518*Jinv(1,0)*Jinv(0,2) + J619*Jinv(0,0)*Jinv(2,2);
    double IssYX = J0789*Jinv(1,1)*Jinv(2,0) + J417*Jinv(0,1)*Jinv(0,0) + J518*Jinv(1,1)*Jinv(0,0) + J619*Jinv(0,1)*Jinv(2,0);
    double IssYY = J0789*Jinv(1,1)*Jinv(2,1) + J417*Jinv(0,1)*Jinv(0,1) + J518*Jinv(1,1)*Jinv(0,1) + J619*Jinv(0,1)*Jinv(2,1);
    double IssYZ = J0789*Jinv(1,1)*Jinv(2,2) + J417*Jinv(0,1)*Jinv(0,2) + J518*Jinv(1,1)*Jinv(0,2) + J619*Jinv(0,1)*Jinv(2,2);
    double IssZX = J0789*Jinv(1,2)*Jinv(2,0) + J417*Jinv(0,2)*Jinv(0,0) + J518*Jinv(1,2)*Jinv(0,0) + J619*Jinv(0,2)*Jinv(2,0);
    double IssZY = J0789*Jinv(1,2)*Jinv(2,1) + J417*Jinv(0,2)*Jinv(0,1) + J518*Jinv(1,2)*Jinv(0,1) + J619*Jinv(0,2)*Jinv(2,1);
    double IssZZ = J0789*Jinv(1,2)*Jinv(2,2) + J417*Jinv(0,2)*Jinv(0,2) + J518*Jinv(1,2)*Jinv(0,2) + J619*Jinv(0,2)*Jinv(2,2);

    double IuuXX = J0978*Jinv(1,0)*Jinv(0,0) + J417*Jinv(2,0)*Jinv(0,0) + J518*Jinv(1,0)*Jinv(2,0) + J619*Jinv(2,0)*Jinv(2,0);
    double IuuXY = J0978*Jinv(1,0)*Jinv(0,1) + J417*Jinv(2,0)*Jinv(0,1) + J518*Jinv(1,0)*Jinv(2,1) + J619*Jinv(2,0)*Jinv(2,1);
    double IuuXZ = J0978*Jinv(1,0)*Jinv(0,2) + J417*Jinv(2,0)*Jinv(0,2) + J518*Jinv(1,0)*Jinv(2,2) + J619*Jinv(2,0)*Jinv(2,2);
    double IuuYX = J0978*Jinv(1,1)*Jinv(0,0) + J417*Jinv(2,1)*Jinv(0,0) + J518*Jinv(1,1)*Jinv(2,0) + J619*Jinv(2,1)*Jinv(2,0);
    double IuuYY = J0978*Jinv(1,1)*Jinv(0,1) + J417*Jinv(2,1)*Jinv(0,1) + J518*Jinv(1,1)*Jinv(2,1) + J619*Jinv(2,1)*Jinv(2,1);
    double IuuYZ = J0978*Jinv(1,1)*Jinv(0,2) + J417*Jinv(2,1)*Jinv(0,2) + J518*Jinv(1,1)*Jinv(2,2) + J619*Jinv(2,1)*Jinv(2,2);
    double IuuZX = J0978*Jinv(1,2)*Jinv(0,0) + J417*Jinv(2,2)*Jinv(0,0) + J518*Jinv(1,2)*Jinv(2,0) + J619*Jinv(2,2)*Jinv(2,0);
    double IuuZY = J0978*Jinv(1,2)*Jinv(0,1) + J417*Jinv(2,2)*Jinv(0,1) + J518*Jinv(1,2)*Jinv(2,1) + J619*Jinv(2,2)*Jinv(2,1);
    double IuuZZ = J0978*Jinv(1,2)*Jinv(0,2) + J417*Jinv(2,2)*Jinv(0,2) + J518*Jinv(1,2)*Jinv(2,2) + J619*Jinv(2,2)*Jinv(2,2);
	
	double IstXX = J31013*Jinv(0,0)*Jinv(0,0) + J31310*Jinv(1,0)*Jinv(1,0) + J21411*Jinv(2,0)*Jinv(1,0) + J11512*Jinv(2,0)*Jinv(0,0) + J16*(Jinv(0,0)*Jinv(1,0) + Jinv(1,0)*Jinv(0,0));
    double IstXY = J31013*Jinv(0,0)*Jinv(0,1) + J31310*Jinv(1,0)*Jinv(1,1) + J21411*Jinv(2,0)*Jinv(1,1) + J11512*Jinv(2,0)*Jinv(0,1) + J16*(Jinv(0,0)*Jinv(1,1) + Jinv(1,0)*Jinv(0,1));
    double IstXZ = J31013*Jinv(0,0)*Jinv(0,2) + J31310*Jinv(1,0)*Jinv(1,2) + J21411*Jinv(2,0)*Jinv(1,2) + J11512*Jinv(2,0)*Jinv(0,2) + J16*(Jinv(0,0)*Jinv(1,2) + Jinv(1,0)*Jinv(0,2));
    double IstYX = J31013*Jinv(0,1)*Jinv(0,0) + J31310*Jinv(1,1)*Jinv(1,0) + J21411*Jinv(2,1)*Jinv(1,0) + J11512*Jinv(2,1)*Jinv(0,0) + J16*(Jinv(0,1)*Jinv(1,0) + Jinv(1,1)*Jinv(0,0));
    double IstYY = J31013*Jinv(0,1)*Jinv(0,1) + J31310*Jinv(1,1)*Jinv(1,1) + J21411*Jinv(2,1)*Jinv(1,1) + J11512*Jinv(2,1)*Jinv(0,1) + J16*(Jinv(0,1)*Jinv(1,1) + Jinv(1,1)*Jinv(0,1));
    double IstYZ = J31013*Jinv(0,1)*Jinv(0,2) + J31310*Jinv(1,1)*Jinv(1,2) + J21411*Jinv(2,1)*Jinv(1,2) + J11512*Jinv(2,1)*Jinv(0,2) + J16*(Jinv(0,1)*Jinv(1,2) + Jinv(1,1)*Jinv(0,2));
    double IstZX = J31013*Jinv(0,2)*Jinv(0,0) + J31310*Jinv(1,2)*Jinv(1,0) + J21411*Jinv(2,2)*Jinv(1,0) + J11512*Jinv(2,2)*Jinv(0,0) + J16*(Jinv(0,2)*Jinv(1,0) + Jinv(1,2)*Jinv(0,0));
    double IstZY = J31013*Jinv(0,2)*Jinv(0,1) + J31310*Jinv(1,2)*Jinv(1,1) + J21411*Jinv(2,2)*Jinv(1,1) + J11512*Jinv(2,2)*Jinv(0,1) + J16*(Jinv(0,2)*Jinv(1,1) + Jinv(1,2)*Jinv(0,1));
    double IstZZ = J31013*Jinv(0,2)*Jinv(0,2) + J31310*Jinv(1,2)*Jinv(1,2) + J21411*Jinv(2,2)*Jinv(1,2) + J11512*Jinv(2,2)*Jinv(0,2) + J16*(Jinv(0,2)*Jinv(1,2) + Jinv(1,2)*Jinv(0,2));

    double IutXX = J21114*Jinv(0,0)*Jinv(1,0) + J31013*Jinv(0,0)*Jinv(2,0) + J11215*Jinv(1,0)*Jinv(1,0) + J11512*Jinv(2,0)*Jinv(2,0) + J16*(Jinv(1,0)*Jinv(2,0) + Jinv(2,0)*Jinv(1,0));
    double IutXY = J21114*Jinv(0,0)*Jinv(1,1) + J31013*Jinv(0,0)*Jinv(2,1) + J11215*Jinv(1,0)*Jinv(1,1) + J11512*Jinv(2,0)*Jinv(2,1) + J16*(Jinv(1,0)*Jinv(2,1) + Jinv(2,0)*Jinv(1,1));
    double IutXZ = J21114*Jinv(0,0)*Jinv(1,2) + J31013*Jinv(0,0)*Jinv(2,2) + J11215*Jinv(1,0)*Jinv(1,2) + J11512*Jinv(2,0)*Jinv(2,2) + J16*(Jinv(1,0)*Jinv(2,2) + Jinv(2,0)*Jinv(1,2));
    double IutYX = J21114*Jinv(0,1)*Jinv(1,0) + J31013*Jinv(0,1)*Jinv(2,0) + J11215*Jinv(1,1)*Jinv(1,0) + J11512*Jinv(2,1)*Jinv(2,0) + J16*(Jinv(1,1)*Jinv(2,0) + Jinv(2,1)*Jinv(1,0));
    double IutYY = J21114*Jinv(0,1)*Jinv(1,1) + J31013*Jinv(0,1)*Jinv(2,1) + J11215*Jinv(1,1)*Jinv(1,1) + J11512*Jinv(2,1)*Jinv(2,1) + J16*(Jinv(1,1)*Jinv(2,1) + Jinv(2,1)*Jinv(1,1));
    double IutYZ = J21114*Jinv(0,1)*Jinv(1,2) + J31013*Jinv(0,1)*Jinv(2,2) + J11215*Jinv(1,1)*Jinv(1,2) + J11512*Jinv(2,1)*Jinv(2,2) + J16*(Jinv(1,1)*Jinv(2,2) + Jinv(2,1)*Jinv(1,2));
    double IutZX = J21114*Jinv(0,2)*Jinv(1,0) + J31013*Jinv(0,2)*Jinv(2,0) + J11215*Jinv(1,2)*Jinv(1,0) + J11512*Jinv(2,2)*Jinv(2,0) + J16*(Jinv(1,2)*Jinv(2,0) + Jinv(2,2)*Jinv(1,0));
    double IutZY = J21114*Jinv(0,2)*Jinv(1,1) + J31013*Jinv(0,2)*Jinv(2,1) + J11215*Jinv(1,2)*Jinv(1,1) + J11512*Jinv(2,2)*Jinv(2,1) + J16*(Jinv(1,2)*Jinv(2,1) + Jinv(2,2)*Jinv(1,1));
    double IutZZ = J21114*Jinv(0,2)*Jinv(1,2) + J31013*Jinv(0,2)*Jinv(2,2) + J11215*Jinv(1,2)*Jinv(1,2) + J11512*Jinv(2,2)*Jinv(2,2) + J16*(Jinv(1,2)*Jinv(2,2) + Jinv(2,2)*Jinv(1,2));

    double IusXX = J21114*Jinv(0,0)*Jinv(0,0) + J11215*Jinv(1,0)*Jinv(0,0) + J31310*Jinv(1,0)*Jinv(2,0) + J21411*Jinv(2,0)*Jinv(2,0) + J16*(Jinv(0,0)*Jinv(2,0) + Jinv(2,0)*Jinv(0,0));
    double IusXY = J21114*Jinv(0,0)*Jinv(0,1) + J11215*Jinv(1,0)*Jinv(0,1) + J31310*Jinv(1,0)*Jinv(2,1) + J21411*Jinv(2,0)*Jinv(2,1) + J16*(Jinv(0,0)*Jinv(2,1) + Jinv(2,0)*Jinv(0,1));
    double IusXZ = J21114*Jinv(0,0)*Jinv(0,2) + J11215*Jinv(1,0)*Jinv(0,2) + J31310*Jinv(1,0)*Jinv(2,2) + J21411*Jinv(2,0)*Jinv(2,2) + J16*(Jinv(0,0)*Jinv(2,2) + Jinv(2,0)*Jinv(0,2));
    double IusYX = J21114*Jinv(0,1)*Jinv(0,0) + J11215*Jinv(1,1)*Jinv(0,0) + J31310*Jinv(1,1)*Jinv(2,0) + J21411*Jinv(2,1)*Jinv(2,0) + J16*(Jinv(0,1)*Jinv(2,0) + Jinv(2,1)*Jinv(0,0));
    double IusYY = J21114*Jinv(0,1)*Jinv(0,1) + J11215*Jinv(1,1)*Jinv(0,1) + J31310*Jinv(1,1)*Jinv(2,1) + J21411*Jinv(2,1)*Jinv(2,1) + J16*(Jinv(0,1)*Jinv(2,1) + Jinv(2,1)*Jinv(0,1));
    double IusYZ = J21114*Jinv(0,1)*Jinv(0,2) + J11215*Jinv(1,1)*Jinv(0,2) + J31310*Jinv(1,1)*Jinv(2,2) + J21411*Jinv(2,1)*Jinv(2,2) + J16*(Jinv(0,1)*Jinv(2,2) + Jinv(2,1)*Jinv(0,2));
    double IusZX = J21114*Jinv(0,2)*Jinv(0,0) + J11215*Jinv(1,2)*Jinv(0,0) + J31310*Jinv(1,2)*Jinv(2,0) + J21411*Jinv(2,2)*Jinv(2,0) + J16*(Jinv(0,2)*Jinv(2,0) + Jinv(2,2)*Jinv(0,0));
    double IusZY = J21114*Jinv(0,2)*Jinv(0,1) + J11215*Jinv(1,2)*Jinv(0,1) + J31310*Jinv(1,2)*Jinv(2,1) + J21411*Jinv(2,2)*Jinv(2,1) + J16*(Jinv(0,2)*Jinv(2,1) + Jinv(2,2)*Jinv(0,1));
    double IusZZ = J21114*Jinv(0,2)*Jinv(0,2) + J11215*Jinv(1,2)*Jinv(0,2) + J31310*Jinv(1,2)*Jinv(2,2) + J21411*Jinv(2,2)*Jinv(2,2) + J16*(Jinv(0,2)*Jinv(2,2) + Jinv(2,2)*Jinv(0,2));

	// constitutive constants from material
	const Matrix &CmatI = theMaterial->getInitialTangent();
	double C1 = CmatI(0,0);
	double C2 = CmatI(0,1);
	double C3 = CmatI(3,3);
	double C4 = C2 + C3;

	// define the 12x12 matrix FCF in 3x3 blocks
    // block11
    FCF(0,0) = C1*HstXX + C3*(HstYY + HstZZ);
    FCF(0,1) = C4*HstXY;
    FCF(0,2) = C4*HstXZ;
    FCF(1,0) = FCF(0,1);
    FCF(1,1) = C1*HstYY + C3*(HstXX + HstZZ);
    FCF(1,2) = C4*HstYZ;
    FCF(2,0) = FCF(0,2);
    FCF(2,1) = FCF(1,2);
    FCF(2,2) = C1*HstZZ + C3*(HstYY + HstXX);

    // block12
    FCF(0,3) = C1*IttXX + C3*(IttYY + IttZZ);
    FCF(0,4) = C2*IttXY + C3*IttYX;
    FCF(0,5) = C2*IttXZ + C3*IttZX;
    FCF(1,3) = C2*IttYX + C3*IttXY;
    FCF(1,4) = C1*IttYY + C3*(IttXX + IttZZ);
    FCF(1,5) = C2*IttYZ + C3*IttZY;
    FCF(2,3) = C2*IttZX + C3*IttXZ;
    FCF(2,4) = C2*IttZY + C3*IttYZ;
    FCF(2,5) = C1*IttZZ + C3*(IttYY + IttXX);

    // block13
    FCF(0,6) = C1*IssXX + C3*(IssYY + IssZZ);
    FCF(0,7) = C2*IssXY + C3*IssYX;
    FCF(0,8) = C2*IssXZ + C3*IssZX;
    FCF(1,6) = C2*IssYX + C3*IssXY;
    FCF(1,7) = C1*IssYY + C3*(IssXX + IssZZ);
    FCF(1,8) = C2*IssYZ + C3*IssZY;
    FCF(2,6) = C2*IssZX + C3*IssXZ;
    FCF(2,7) = C2*IssZY + C3*IssYZ;
    FCF(2,8) = C1*IssZZ + C3*(IssYY + IssXX);

    /*// block14
    FCF(0,9)  = C1*IstXX + C3*(IstYY + IstZZ);
    FCF(0,10) = C2*IstYX + C3*IstXY;
    FCF(0,11) = C2*IstZX + C3*IstXZ;
    FCF(1,9)  = C2*IstXY + C3*IstYX;
    FCF(1,10) = C1*IstYY + C3*(IstXX + IstZZ);
    FCF(1,11) = C2*IstZY + C3*IstYZ;
    FCF(2,9)  = C2*IstXZ + C3*IstZX;
    FCF(2,10) = C2*IstYZ + C3*IstZY;
    FCF(2,11) = C1*IstZZ + C3*(IstYY + IstXX);*/

	// block14
    FCF(0,9)  = C3*(IstYY + IstZZ);
    FCF(0,10) = C3*IstXY;
    FCF(0,11) = C3*IstXZ;
    FCF(1,9)  = C3*IstYX;
    FCF(1,10) = C3*(IstXX + IstZZ);
    FCF(1,11) = C3*IstYZ;
    FCF(2,9)  = C3*IstZX;
    FCF(2,10) = C3*IstZY;
    FCF(2,11) = C3*(IstYY + IstXX);

    // block21
    FCF(3,0) = C1*IttXX + C3*(IttYY + IttZZ);
    FCF(3,1) = C2*IttYX + C3*IttXY;
    FCF(3,2) = C2*IttZX + C3*IttXZ;
    FCF(4,0) = C2*IttXY + C3*IttYX;
    FCF(4,1) = C1*IttYY + C3*(IttXX + IttZZ);
    FCF(4,2) = C2*IttZY + C3*IttYZ;
    FCF(5,0) = C2*IttXZ + C3*IttZX;
    FCF(5,1) = C2*IttYZ + C3*IttZY;
    FCF(5,2) = C1*IttZZ + C3*(IttYY + IttXX);

    // block22
    FCF(3,3) = C1*HutXX + C3*(HutYY + HutZZ);
    FCF(3,4) = C4*HutXY;
    FCF(3,5) = C4*HutXZ;
    FCF(4,3) = FCF(3,4);
    FCF(4,4) = C1*HutYY + C3*(HutXX + HutZZ);
    FCF(4,5) = C4*HutYZ;
    FCF(5,3) = FCF(3,5);
    FCF(5,4) = FCF(4,5);
    FCF(5,5) = C1*HutZZ + C3*(HutYY + HutXX);

    // block23
    FCF(3,6) = C1*IuuXX + C3*(IuuYY + IuuZZ);
    FCF(3,7) = C2*IuuXY + C3*IuuYX;
    FCF(3,8) = C2*IuuXZ + C3*IuuZX;
    FCF(4,6) = C2*IuuYX + C3*IuuXY;
    FCF(4,7) = C1*IuuYY + C3*(IuuXX + IuuZZ);
    FCF(4,8) = C2*IuuYZ + C3*IuuZY;
    FCF(5,6) = C2*IuuZX + C3*IuuXZ;
    FCF(5,7) = C2*IuuZY + C3*IuuYZ;
    FCF(5,8) = C1*IuuZZ + C3*(IuuYY + IuuXX);

    /*// block24
    FCF(3,9)  = C1*IutXX + C3*(IutYY + IutZZ);
    FCF(3,10) = C2*IutYX + C3*IutXY;
    FCF(3,11) = C2*IutZX + C3*IutXZ;
    FCF(4,9)  = C2*IutXY + C3*IutYX;
    FCF(4,10) = C1*IutYY + C3*(IutXX + IutZZ);
    FCF(4,11) = C2*IutZY + C3*IutYZ;
    FCF(5,9)  = C2*IutXZ + C3*IutZX;
    FCF(5,10) = C2*IutYZ + C3*IutZY;
    FCF(5,11) = C1*IutZZ + C3*(IutYY + IutXX);*/

	// block24
    FCF(3,9)  = C3*(IutYY + IutZZ);
    FCF(3,10) = C3*IutXY;
    FCF(3,11) = C3*IutXZ;
    FCF(4,9)  = C3*IutYX;
    FCF(4,10) = C3*(IutXX + IutZZ);
    FCF(4,11) = C3*IutYZ;
    FCF(5,9)  = C3*IutZX;
    FCF(5,10) = C3*IutZY;
    FCF(5,11) = C3*(IutYY + IutXX);

    // block31
    FCF(6,0) = C1*IssXX + C3*(IssYY + IssZZ);
    FCF(6,1) = C2*IssYX + C3*IssXY;
    FCF(6,2) = C2*IssZX + C3*IssXZ;
    FCF(7,0) = C2*IssXY + C3*IssYX;
    FCF(7,1) = C1*IssYY + C3*(IssXX + IssZZ);
    FCF(7,2) = C2*IssZY + C3*IssYZ;
    FCF(8,0) = C2*IssXZ + C3*IssZX;
    FCF(8,1) = C2*IssYZ + C3*IssZY;
    FCF(8,2) = C1*IssZZ + C3*(IssYY + IssXX);

    // block32
    FCF(6,3) = C1*IuuXX + C3*(IuuYY + IuuZZ);
    FCF(6,4) = C2*IuuYX + C3*IuuXY;
    FCF(6,5) = C2*IuuZX + C3*IuuXZ;
    FCF(7,3) = C2*IuuXY + C3*IuuYX;
    FCF(7,4) = C1*IuuYY + C3*(IuuXX + IuuZZ);
    FCF(7,5) = C2*IuuZY + C3*IuuYZ;
    FCF(8,3) = C2*IuuXZ + C3*IuuZX;
    FCF(8,4) = C2*IuuYZ + C3*IuuZY;
    FCF(8,5) = C1*IuuZZ + C3*(IuuYY + IuuXX);

    // block33
    FCF(6,6) = C1*HusXX + C3*(HusYY + HusZZ);
    FCF(6,7) = C4*HusXY;
    FCF(6,8) = C4*HusXZ;
    FCF(7,6) = FCF(6,7);
    FCF(7,7) = C1*HusYY + C3*(HusXX + HusZZ);
    FCF(7,8) = C4*HusYZ;
    FCF(8,6) = FCF(6,8);
    FCF(8,7) = FCF(7,8);
    FCF(8,8) = C1*HusZZ + C3*(HusYY + HusXX);

    /*// block34
    FCF(6,9)  = C1*IusXX + C3*(IusYY + IusZZ);
    FCF(6,10) = C2*IusYX + C3*IusXY;
    FCF(6,11) = C2*IusZX + C3*IusXZ;
    FCF(7,9)  = C2*IusXY + C3*IusYX;
    FCF(7,10) = C1*IusYY + C3*(IusXX + IusZZ);
    FCF(7,11) = C2*IusZY + C3*IusYZ;
    FCF(8,9)  = C2*IusXZ + C3*IusZX;
    FCF(8,10) = C2*IusYZ + C3*IusZY;
    FCF(8,11) = C1*IusZZ + C3*(IusYY + IusXX);*/

	// block34
    FCF(6,9)  = C3*(IusYY + IusZZ);
    FCF(6,10) = C3*IusXY;
    FCF(6,11) = C3*IusXZ;
    FCF(7,9)  = C3*IusYX;
    FCF(7,10) = C3*(IusXX + IusZZ);
    FCF(7,11) = C3*IusYZ;
    FCF(8,9)  = C3*IusZX;
    FCF(8,10) = C3*IusZY;
    FCF(8,11) = C3*(IusYY + IusXX);

    /*// block41
    FCF(9,0)  = C1*IstXX + C3*(IstYY + IstZZ);
    FCF(9,1)  = C2*IstXY + C3*IstYX;
    FCF(9,2)  = C2*IstXZ + C3*IstZX;
    FCF(10,0) = C2*IstYX + C3*IstXY;
    FCF(10,1) = C1*IstYY + C3*(IstXX + IstZZ);
    FCF(10,2) = C2*IstYZ + C3*IstZY;
    FCF(11,0) = C2*IstZX + C3*IstXZ;
    FCF(11,1) = C2*IstZY + C3*IstYZ;
    FCF(11,2) = C1*IstZZ + C3*(IstYY + IstXX);

    // block42
    FCF(9,3)  = C1*IutXX + C3*(IutYY + IutZZ);
    FCF(9,4)  = C2*IutXY + C3*IutYX;
    FCF(9,5)  = C2*IutXZ + C3*IutZX;
    FCF(10,3) = C2*IutYX + C3*IutXY;
    FCF(10,4) = C1*IutYY + C3*(IutXX + IutZZ);
    FCF(10,5) = C2*IutYZ + C3*IutZY;
    FCF(11,3) = C2*IutZX + C3*IutXZ;
    FCF(11,4) = C2*IutZY + C3*IutYZ;
    FCF(11,5) = C1*IutZZ + C3*(IutYY + IutXX);

    // block43
    FCF(9,6)  = C1*IusXX + C3*(IusYY + IusZZ);
    FCF(9,7)  = C2*IusXY + C3*IusYX;
    FCF(9,8)  = C2*IusXZ + C3*IusZX;
    FCF(10,6) = C2*IusYX + C3*IusXY;
    FCF(10,7) = C1*IusYY + C3*(IusXX + IusZZ);
    FCF(10,8) = C2*IusYZ + C3*IusZY;
    FCF(11,6) = C2*IusZX + C3*IusXZ;
    FCF(11,7) = C2*IusZY + C3*IusYZ;
    FCF(11,8) = C1*IusZZ + C3*(IusYY + IusXX);*/

	// block41
    FCF(9,0)  = C3*(IstYY + IstZZ);
    FCF(9,1)  = C3*IstYX;
    FCF(9,2)  = C3*IstZX;
    FCF(10,0) = C3*IstXY;
    FCF(10,1) = C3*(IstXX + IstZZ);
    FCF(10,2) = C3*IstZY;
    FCF(11,0) = C3*IstXZ;
    FCF(11,1) = C3*IstYZ;
    FCF(11,2) = C3*(IstYY + IstXX);

    // block42
    FCF(9,3)  = C3*(IutYY + IutZZ);
    FCF(9,4)  = C3*IutYX;
    FCF(9,5)  = C3*IutZX;
    FCF(10,3) = C3*IutXY;
    FCF(10,4) = C3*(IutXX + IutZZ);
    FCF(10,5) = C3*IutZY;
    FCF(11,3) = C3*IutXZ;
    FCF(11,4) = C3*IutYZ;
    FCF(11,5) = C3*(IutYY + IutXX);

    // block43
    FCF(9,6)  = C3*(IusYY + IusZZ);
    FCF(9,7)  = C3*IusYX;
    FCF(9,8)  = C3*IusZX;
    FCF(10,6) = C3*IusXY;
    FCF(10,7) = C3*(IusXX + IusZZ);
    FCF(10,8) = C3*IusZY;
    FCF(11,6) = C3*IusXZ;
    FCF(11,7) = C3*IusYZ;
    FCF(11,8) = C3*(IusYY + IusXX);

    /*// block44
    FCF(9,9)   = C1*HstuXX + C3*(HstuYY + HstuZZ);
    FCF(9,10)  = C4*HstuXY;
    FCF(9,11)  = C4*HstuXZ;
    FCF(10,9)  = C4*HstuXY;
    FCF(10,10) = C1*HstuYY + C3*(HstuXX + HstuZZ);
    FCF(10,11) = C4*HstuYZ;
    FCF(11,9)  = C4*HstuXZ;
    FCF(11,10) = C4*HstuYZ;
    FCF(11,11) = C1*HstuZZ + C3*(HstuYY + HstuXX);*/

	// block44
    FCF(9,9)   = C3*(HstuYY + HstuZZ);
    FCF(9,10)  = C3*HstuXY;
    FCF(9,11)  = C3*HstuXZ;
    FCF(10,9)  = FCF(9,10);
    FCF(10,10) = C3*(HstuXX + HstuZZ);
    FCF(10,11) = C3*HstuYZ;
    FCF(11,9)  = FCF(9,11);
    FCF(11,10) = FCF(10,11);
    FCF(11,11) = C3*(HstuYY + HstuXX);

	// enhanced strain portion of the stabilization stiffness matrix
	// define the constitutive coefficients
	double CssXX = C1*Jinv(0,0)*Jinv(0,0) + C3*(Jinv(0,1)*Jinv(0,1) + Jinv(0,2)*Jinv(0,2));
	double CssXY = C4*Jinv(0,0)*Jinv(0,1);
	double CssXZ = C4*Jinv(0,0)*Jinv(0,2);
	double CssYY = C1*Jinv(0,1)*Jinv(0,1) + C3*(Jinv(0,0)*Jinv(0,0) + Jinv(0,2)*Jinv(0,2));
	double CssYZ = C4*Jinv(0,1)*Jinv(0,2);
	double CssZZ = C1*Jinv(0,2)*Jinv(0,2) + C3*(Jinv(0,0)*Jinv(0,0) + Jinv(0,1)*Jinv(0,1));

	double CttXX = C1*Jinv(1,0)*Jinv(1,0) + C3*(Jinv(1,1)*Jinv(1,1) + Jinv(1,2)*Jinv(1,2));
	double CttXY = C4*Jinv(1,0)*Jinv(1,1);
	double CttXZ = C4*Jinv(1,0)*Jinv(1,2);
	double CttYY = C1*Jinv(1,1)*Jinv(1,1) + C3*(Jinv(1,0)*Jinv(1,0) + Jinv(1,2)*Jinv(1,2));
	double CttYZ = C4*Jinv(1,1)*Jinv(1,2);
	double CttZZ = C1*Jinv(1,2)*Jinv(1,2) + C3*(Jinv(1,0)*Jinv(1,0) + Jinv(1,1)*Jinv(1,1));

	double CuuXX = C1*Jinv(2,0)*Jinv(2,0) + C3*(Jinv(2,1)*Jinv(2,1) + Jinv(2,2)*Jinv(2,2));
	double CuuXY = C4*Jinv(2,0)*Jinv(2,1);
	double CuuXZ = C4*Jinv(2,0)*Jinv(2,2);
	double CuuYY = C1*Jinv(2,1)*Jinv(2,1) + C3*(Jinv(2,0)*Jinv(2,0) + Jinv(2,2)*Jinv(2,2));
	double CuuYZ = C4*Jinv(2,1)*Jinv(2,2);
	double CuuZZ = C1*Jinv(2,2)*Jinv(2,2) + C3*(Jinv(2,0)*Jinv(2,0) + Jinv(2,1)*Jinv(2,1));

	double CstXX = C1*Jinv(0,0)*Jinv(1,0) + C3*(Jinv(0,1)*Jinv(1,1) + Jinv(0,2)*Jinv(1,2));
	double CstXY = C2*Jinv(0,0)*Jinv(1,1) + C3*Jinv(0,1)*Jinv(1,0);
	double CstXZ = C2*Jinv(0,0)*Jinv(1,2) + C3*Jinv(0,2)*Jinv(1,0);
	double CstYX = C2*Jinv(0,1)*Jinv(1,0) + C3*Jinv(0,0)*Jinv(1,1);
	double CstYY = C1*Jinv(0,1)*Jinv(1,1) + C3*(Jinv(0,0)*Jinv(1,0) + Jinv(0,2)*Jinv(1,2));
	double CstYZ = C2*Jinv(0,1)*Jinv(1,2) + C3*Jinv(0,2)*Jinv(1,1);
	double CstZX = C2*Jinv(0,2)*Jinv(1,0) + C3*Jinv(0,0)*Jinv(1,2);
	double CstZY = C2*Jinv(0,2)*Jinv(1,1) + C3*Jinv(0,1)*Jinv(1,2);
	double CstZZ = C1*Jinv(0,2)*Jinv(1,2) + C3*(Jinv(0,0)*Jinv(1,0) + Jinv(0,1)*Jinv(1,1));

	double CsuXX = C1*Jinv(0,0)*Jinv(2,0) + C3*(Jinv(0,1)*Jinv(2,1) + Jinv(0,2)*Jinv(2,2));
	double CsuXY = C2*Jinv(0,0)*Jinv(2,1) + C3*Jinv(0,1)*Jinv(2,0);
	double CsuXZ = C2*Jinv(0,0)*Jinv(2,2) + C3*Jinv(0,2)*Jinv(2,0);
	double CsuYX = C2*Jinv(0,1)*Jinv(2,0) + C3*Jinv(0,0)*Jinv(2,1);
	double CsuYY = C1*Jinv(0,1)*Jinv(2,1) + C3*(Jinv(0,0)*Jinv(2,0) + Jinv(0,2)*Jinv(2,2));
	double CsuYZ = C2*Jinv(0,1)*Jinv(2,2) + C3*Jinv(0,2)*Jinv(2,1);
	double CsuZX = C2*Jinv(0,2)*Jinv(2,0) + C3*Jinv(0,0)*Jinv(2,2);
	double CsuZY = C2*Jinv(0,2)*Jinv(2,1) + C3*Jinv(0,1)*Jinv(2,2);
	double CsuZZ = C1*Jinv(0,2)*Jinv(2,2) + C3*(Jinv(0,0)*Jinv(2,0) + Jinv(0,1)*Jinv(2,1));

	double CtuXX = C1*Jinv(1,0)*Jinv(2,0) + C3*(Jinv(1,1)*Jinv(2,1) + Jinv(1,2)*Jinv(2,2));
	double CtuXY = C2*Jinv(1,0)*Jinv(2,1) + C3*Jinv(1,1)*Jinv(2,0);
	double CtuXZ = C2*Jinv(1,0)*Jinv(2,2) + C3*Jinv(1,2)*Jinv(2,0);
	double CtuYX = C2*Jinv(1,1)*Jinv(2,0) + C3*Jinv(1,0)*Jinv(2,1);
	double CtuYY = C1*Jinv(1,1)*Jinv(2,1) + C3*(Jinv(1,0)*Jinv(2,0) + Jinv(1,2)*Jinv(2,2));
	double CtuYZ = C2*Jinv(1,1)*Jinv(2,2) + C3*Jinv(1,2)*Jinv(2,1);
	double CtuZX = C2*Jinv(1,2)*Jinv(2,0) + C3*Jinv(1,0)*Jinv(2,2);
	double CtuZY = C2*Jinv(1,2)*Jinv(2,1) + C3*Jinv(1,1)*Jinv(2,2);
	double CtuZZ = C1*Jinv(1,2)*Jinv(2,2) + C3*(Jinv(1,0)*Jinv(2,0) + Jinv(1,1)*Jinv(2,1));

	// define the integrated matrix [FenT][C][Fen]
	Matrix FeCFe(9,9);
	Matrix FeCFeInv(9,9);

	// block11
	FeCFe(0,0) = CssXX*J0789;
	FeCFe(0,1) = CssXY*J0789;
	FeCFe(0,2) = CssXZ*J0789;
	FeCFe(1,0) = FeCFe(0,1);
	FeCFe(1,1) = CssYY*J0789;
	FeCFe(1,2) = CssYZ*J0789;
	FeCFe(2,0) = FeCFe(0,2);
	FeCFe(2,1) = FeCFe(1,2);
	FeCFe(2,2) = CssZZ*J0789;

	// block12
	FeCFe(0,3) = CstXX*J619;
	FeCFe(0,4) = CstXY*J619;
	FeCFe(0,5) = CstXZ*J619;
	FeCFe(1,3) = CstYX*J619;
	FeCFe(1,4) = CstYY*J619;
	FeCFe(1,5) = CstYZ*J619;
	FeCFe(2,3) = CstZX*J619;
	FeCFe(2,4) = CstZY*J619;
	FeCFe(2,5) = CstZZ*J619;
	
	// block13
	FeCFe(0,6) = CsuXX*J518;
	FeCFe(0,7) = CsuXY*J518;
	FeCFe(0,8) = CsuXZ*J518;
	FeCFe(1,6) = CsuYX*J518;
	FeCFe(1,7) = CsuYY*J518;
	FeCFe(1,8) = CsuYZ*J518;
	FeCFe(2,6) = CsuZX*J518;
	FeCFe(2,7) = CsuZY*J518;
	FeCFe(2,8) = CsuZZ*J518;

	// block21
	FeCFe(3,0) = CstXX*J619;
	FeCFe(3,1) = CstYX*J619;
	FeCFe(3,2) = CstZX*J619;
	FeCFe(4,0) = CstXY*J619;
	FeCFe(4,1) = CstYY*J619;
	FeCFe(4,2) = CstZY*J619;
	FeCFe(5,0) = CstXZ*J619;
	FeCFe(5,1) = CstYZ*J619;
	FeCFe(5,2) = CstZZ*J619;

	// block22
	FeCFe(3,3) = CttXX*J0879;
	FeCFe(3,4) = CttXY*J0879;
	FeCFe(3,5) = CttXZ*J0879;
	FeCFe(4,3) = FeCFe(3,4);
	FeCFe(4,4) = CttYY*J0879;
	FeCFe(4,5) = CttYZ*J0879;
	FeCFe(5,3) = FeCFe(3,5);
	FeCFe(5,4) = FeCFe(4,5);
	FeCFe(5,5) = CttZZ*J0879;

	// block23
	FeCFe(3,6) = CtuXX*J417;
	FeCFe(3,7) = CtuXY*J417;
	FeCFe(3,8) = CtuXZ*J417;
	FeCFe(4,6) = CtuYX*J417;
	FeCFe(4,7) = CtuYY*J417;
	FeCFe(4,8) = CtuYZ*J417;
	FeCFe(5,6) = CtuZX*J417;
	FeCFe(5,7) = CtuZY*J417;
	FeCFe(5,8) = CtuZZ*J417;

	// block31
	FeCFe(6,0) = CsuXX*J518;
	FeCFe(6,1) = CsuYX*J518;
	FeCFe(6,2) = CsuZX*J518;
	FeCFe(7,0) = CsuXY*J518;
	FeCFe(7,1) = CsuYY*J518;
	FeCFe(7,2) = CsuZY*J518;
	FeCFe(8,0) = CsuXZ*J518;
	FeCFe(8,1) = CsuYZ*J518;
	FeCFe(8,2) = CsuZZ*J518;

	// block32
	FeCFe(6,3) = CtuXX*J417;
	FeCFe(6,4) = CtuYX*J417;
	FeCFe(6,5) = CtuZX*J417;
	FeCFe(7,3) = CtuXY*J417;
	FeCFe(7,4) = CtuYY*J417;
	FeCFe(7,5) = CtuZY*J417;
	FeCFe(8,3) = CtuXZ*J417;
	FeCFe(8,4) = CtuYZ*J417;
	FeCFe(8,5) = CtuZZ*J417;

	// block33
	FeCFe(6,6) = CuuXX*J0978;
	FeCFe(6,7) = CuuXY*J0978;
	FeCFe(6,8) = CuuXZ*J0978;
	FeCFe(7,6) = FeCFe(6,7);
	FeCFe(7,7) = CuuYY*J0978;
	FeCFe(7,8) = CuuYZ*J0978;
	FeCFe(8,6) = FeCFe(6,8);
	FeCFe(8,7) = FeCFe(7,8);
	FeCFe(8,8) = CuuZZ*J0978;

	// inverse of [FenT][C][Fen]
	FeCFe.Invert(FeCFeInv);

	// define the integrated matrix [FenT][C][Fhg]
	Matrix FeCFhg(9,12);
	FeCFhg.Zero();

	// block11
	FeCFhg(0,0) = CstXX*J0789 + CssXX*J619;
	FeCFhg(0,1) = CstXY*J0789 + CssXY*J619;
	FeCFhg(0,2) = CstXZ*J0789 + CssXZ*J619;
	FeCFhg(1,0) = CstYX*J0789 + CssXY*J619;
	FeCFhg(1,1) = CstYY*J0789 + CssYY*J619;
	FeCFhg(1,2) = CstYZ*J0789 + CssYZ*J619;
	FeCFhg(2,0) = CstZX*J0789 + CssXZ*J619;
	FeCFhg(2,1) = CstZY*J0789 + CssYZ*J619;
	FeCFhg(2,2) = CstZZ*J0789 + CssZZ*J619;

	// block12
	FeCFhg(0,3) = CsuXX*J619 + CstXX*J518;
	FeCFhg(0,4) = CsuXY*J619 + CstXY*J518;
	FeCFhg(0,5) = CsuXZ*J619 + CstXZ*J518;
	FeCFhg(1,3) = CsuYX*J619 + CstYX*J518;
	FeCFhg(1,4) = CsuYY*J619 + CstYY*J518;
	FeCFhg(1,5) = CsuYZ*J619 + CstYZ*J518;
	FeCFhg(2,3) = CsuZX*J619 + CstZX*J518;
	FeCFhg(2,4) = CsuZY*J619 + CstZY*J518;
	FeCFhg(2,5) = CsuZZ*J619 + CstZZ*J518;

	// block13
	FeCFhg(0,6) = CsuXX*J0789 + CssXX*J518;
	FeCFhg(0,7) = CsuXY*J0789 + CssXY*J518;
	FeCFhg(0,8) = CsuXZ*J0789 + CssXZ*J518;
	FeCFhg(1,6) = CsuYX*J0789 + CssXY*J518;
	FeCFhg(1,7) = CsuYY*J0789 + CssYY*J518;
	FeCFhg(1,8) = CsuYZ*J0789 + CssYZ*J518;
	FeCFhg(2,6) = CsuZX*J0789 + CssXZ*J518;
	FeCFhg(2,7) = CsuZY*J0789 + CssYZ*J518;
	FeCFhg(2,8) = CsuZZ*J0789 + CssZZ*J518;

	// block21
	FeCFhg(3,0) = CstXX*J0879 + CttXX*J619;
	FeCFhg(3,1) = CstYX*J0879 + CttXY*J619;
	FeCFhg(3,2) = CstZX*J0879 + CttXZ*J619;
	FeCFhg(4,0) = CstXY*J0879 + CttXY*J619;
	FeCFhg(4,1) = CstYY*J0879 + CttYY*J619;
	FeCFhg(4,2) = CstZY*J0879 + CttYZ*J619;
	FeCFhg(5,0) = CstXZ*J0879 + CttXZ*J619;
	FeCFhg(5,1) = CstYZ*J0879 + CttYZ*J619;
	FeCFhg(5,2) = CstZZ*J0879 + CttZZ*J619;

	// block22
	FeCFhg(3,3) = CtuXX*J0879 + CttXX*J417;
	FeCFhg(3,4) = CtuXY*J0879 + CttXY*J417;
	FeCFhg(3,5) = CtuXZ*J0879 + CttXZ*J417;
	FeCFhg(4,3) = CtuYX*J0879 + CttXY*J417;
	FeCFhg(4,4) = CtuYY*J0879 + CttYY*J417;
	FeCFhg(4,5) = CtuYZ*J0879 + CttYZ*J417;
	FeCFhg(5,3) = CtuZX*J0879 + CttXZ*J417;
	FeCFhg(5,4) = CtuZY*J0879 + CttYZ*J417;
	FeCFhg(5,5) = CtuZZ*J0879 + CttZZ*J417;

	// block23
	FeCFhg(3,6) = CtuXX*J619 + CstXX*J417;
	FeCFhg(3,7) = CtuXY*J619 + CstYX*J417;
	FeCFhg(3,8) = CtuXZ*J619 + CstZX*J417;
	FeCFhg(4,6) = CtuYX*J619 + CstXY*J417;
	FeCFhg(4,7) = CtuYY*J619 + CstYY*J417;
	FeCFhg(4,8) = CtuYZ*J619 + CstZY*J417;
	FeCFhg(5,6) = CtuZX*J619 + CstXZ*J417;
	FeCFhg(5,7) = CtuZY*J619 + CstYZ*J417;
	FeCFhg(5,8) = CtuZZ*J619 + CstZZ*J417;

	// block31
	FeCFhg(6,0) = CtuXX*J518 + CsuXX*J417;
	FeCFhg(6,1) = CtuYX*J518 + CsuYX*J417;
	FeCFhg(6,2) = CtuZX*J518 + CsuZX*J417;
	FeCFhg(7,0) = CtuXY*J518 + CsuXY*J417;
	FeCFhg(7,1) = CtuYY*J518 + CsuYY*J417;
	FeCFhg(7,2) = CtuZY*J518 + CsuZY*J417;
	FeCFhg(8,0) = CtuXZ*J518 + CsuXZ*J417;
	FeCFhg(8,1) = CtuYZ*J518 + CsuYZ*J417;
	FeCFhg(8,2) = CtuZZ*J518 + CsuZZ*J417;

	// block32
	FeCFhg(6,3) = CtuXX*J0978 + CuuXX*J417;
	FeCFhg(6,4) = CtuYX*J0978 + CuuXY*J417;
	FeCFhg(6,5) = CtuZX*J0978 + CuuXZ*J417;
	FeCFhg(7,3) = CtuXY*J0978 + CuuXY*J417;
	FeCFhg(7,4) = CtuYY*J0978 + CuuYY*J417;
	FeCFhg(7,5) = CtuZY*J0978 + CuuYZ*J417;
	FeCFhg(8,3) = CtuXZ*J0978 + CuuXZ*J417;
	FeCFhg(8,4) = CtuYZ*J0978 + CuuYZ*J417;
	FeCFhg(8,5) = CtuZZ*J0978 + CuuZZ*J417;

	// block33
	FeCFhg(6,6) = CsuXX*J0978 + CuuXX*J518;
	FeCFhg(6,7) = CsuYX*J0978 + CuuXY*J518;
	FeCFhg(6,8) = CsuZX*J0978 + CuuXZ*J518;
	FeCFhg(7,6) = CsuXY*J0978 + CuuXY*J518;
	FeCFhg(7,7) = CsuYY*J0978 + CuuYY*J518;
	FeCFhg(7,8) = CsuZY*J0978 + CuuYZ*J518;
	FeCFhg(8,6) = CsuXZ*J0978 + CuuXZ*J518;
	FeCFhg(8,7) = CsuYZ*J0978 + CuuYZ*J518;
	FeCFhg(8,8) = CsuZZ*J0978 + CuuZZ*J518;

	/*// off-diagonal enhanced strain matrix [Fe]^T[C][F] is same as [Fe]^T[C][Fhg] with exception of last three columns
	Matrix FeCF(9,12);
	FeCF = FeCFhg;

	// block 14
	FeCF(0,9)  = C3*(J21411*(Jinv(2,1)*Jinv(0,1) + Jinv(2,2)*Jinv(0,2)) + J31310*(Jinv(1,1)*Jinv(0,1) + Jinv(1,2)*Jinv(0,2)) + J16*(Jinv(0,1)*Jinv(0,1) + Jinv(0,2)*Jinv(0,2)));
	FeCF(0,10) = C3*(J21411*Jinv(2,0)*Jinv(0,1) + J31310*Jinv(1,0)*Jinv(0,1) + J16*Jinv(0,0)*Jinv(0,1));
	FeCF(0,11) = C3*(J21411*Jinv(2,0)*Jinv(0,2) + J31310*Jinv(1,0)*Jinv(0,2) + J16*Jinv(0,0)*Jinv(0,2));
	FeCF(1,9)  = C3*(J21411*Jinv(2,1)*Jinv(0,0) + J31310*Jinv(1,1)*Jinv(0,0) + J16*Jinv(0,0)*Jinv(0,1));
	FeCF(1,10) = C3*(J21411*(Jinv(2,0)*Jinv(0,0) + Jinv(2,2)*Jinv(0,2)) + J31310*(Jinv(1,0)*Jinv(0,0) + Jinv(1,2)*Jinv(0,2)) + J16*(Jinv(0,0)*Jinv(0,0) + Jinv(0,2)*Jinv(0,2)));
	FeCF(1,11) = C3*(J21411*Jinv(2,1)*Jinv(0,2) + J31310*Jinv(1,1)*Jinv(0,2) + J16*Jinv(0,1)*Jinv(0,2));
	FeCF(2,9)  = C3*(J21411*Jinv(2,2)*Jinv(0,0) + J31310*Jinv(1,2)*Jinv(0,0) + J16*Jinv(0,0)*Jinv(0,2));
	FeCF(2,10) = C3*(J21411*Jinv(2,2)*Jinv(0,1) + J31310*Jinv(1,2)*Jinv(0,1) + J16*Jinv(0,1)*Jinv(0,2));
	FeCF(2,11) = C3*(J21411*(Jinv(2,0)*Jinv(0,0) + Jinv(2,1)*Jinv(0,1)) + J31310*(Jinv(1,0)*Jinv(0,0) + Jinv(1,1)*Jinv(0,1)) + J16*(Jinv(0,0)*Jinv(0,0) + Jinv(0,1)*Jinv(0,1)));

	// block 24
	FeCF(3,9)  = C3*(J11512*(Jinv(2,1)*Jinv(1,1) + Jinv(2,2)*Jinv(1,2)) + J31013*(Jinv(1,1)*Jinv(0,1) + Jinv(1,2)*Jinv(0,2)) + J16*(Jinv(1,1)*Jinv(1,1) + Jinv(1,2)*Jinv(1,2)));
	FeCF(3,10) = C3*(J11512*Jinv(2,0)*Jinv(1,1) + J31013*Jinv(1,1)*Jinv(0,0) + J16*Jinv(1,0)*Jinv(1,1));
	FeCF(3,11) = C3*(J11512*Jinv(2,0)*Jinv(1,2) + J31013*Jinv(1,2)*Jinv(0,0) + J16*Jinv(1,0)*Jinv(1,2));
	FeCF(4,9)  = C3*(J11512*Jinv(2,1)*Jinv(1,0) + J31013*Jinv(1,0)*Jinv(0,1) + J16*Jinv(1,0)*Jinv(1,1));
	FeCF(4,10) = C3*(J11512*(Jinv(2,0)*Jinv(1,0) + Jinv(2,2)*Jinv(1,2)) + J31013*(Jinv(1,0)*Jinv(0,0) + Jinv(1,2)*Jinv(0,2)) + J16*(Jinv(1,0)*Jinv(1,0) + Jinv(1,2)*Jinv(1,2)));
	FeCF(4,11) = C3*(J11512*Jinv(2,1)*Jinv(1,2) + J31013*Jinv(1,2)*Jinv(0,1) + J16*Jinv(1,1)*Jinv(1,2));
	FeCF(5,9)  = C3*(J11512*Jinv(2,2)*Jinv(1,0) + J31013*Jinv(1,0)*Jinv(0,2) + J16*Jinv(1,0)*Jinv(1,2));
	FeCF(5,10) = C3*(J11512*Jinv(2,2)*Jinv(1,1) + J31013*Jinv(1,1)*Jinv(0,2) + J16*Jinv(1,1)*Jinv(1,2));
	FeCF(5,11) = C3*(J11512*(Jinv(2,0)*Jinv(1,0) + Jinv(2,1)*Jinv(1,1)) + J31013*(Jinv(1,0)*Jinv(0,0) + Jinv(1,1)*Jinv(0,1)) + J16*(Jinv(1,0)*Jinv(1,0) + Jinv(1,1)*Jinv(1,1)));

	// block 34
	FeCF(6,9)  = C3*(J11215*(Jinv(2,1)*Jinv(1,1) + Jinv(2,2)*Jinv(1,2)) + J21114*(Jinv(2,1)*Jinv(0,1) + Jinv(2,2)*Jinv(0,2)) + J16*(Jinv(2,1)*Jinv(2,1) + Jinv(2,2)*Jinv(2,2)));
	FeCF(6,10) = C3*(J11215*Jinv(2,1)*Jinv(1,0) + J21114*Jinv(2,1)*Jinv(0,0) + J16*Jinv(2,0)*Jinv(2,1));
	FeCF(6,11) = C3*(J11215*Jinv(2,2)*Jinv(1,0) + J21114*Jinv(2,2)*Jinv(0,0) + J16*Jinv(2,0)*Jinv(2,2));
	FeCF(7,9)  = C3*(J11215*Jinv(2,0)*Jinv(1,1) + J21114*Jinv(2,0)*Jinv(0,1) + J16*Jinv(2,0)*Jinv(2,1));
	FeCF(7,10) = C3*(J11215*(Jinv(2,0)*Jinv(1,0) + Jinv(2,2)*Jinv(1,2)) + J21114*(Jinv(2,0)*Jinv(0,0) + Jinv(2,2)*Jinv(0,2)) + J16*(Jinv(2,0)*Jinv(2,0) + Jinv(2,2)*Jinv(2,2)));
	FeCF(7,11) = C3*(J11215*Jinv(2,2)*Jinv(1,1) + J21114*Jinv(2,2)*Jinv(0,1) + J16*Jinv(2,1)*Jinv(2,2));
	FeCF(8,9)  = C3*(J11215*Jinv(2,0)*Jinv(1,2) + J21114*Jinv(2,0)*Jinv(0,2) + J16*Jinv(2,0)*Jinv(2,2));
	FeCF(8,10) = C3*(J11215*Jinv(2,1)*Jinv(1,2) + J21114*Jinv(2,1)*Jinv(0,2) + J16*Jinv(2,1)*Jinv(2,2));
	FeCF(8,11) = C3*(J11215*(Jinv(2,0)*Jinv(1,0) + Jinv(2,1)*Jinv(1,1)) + J21114*(Jinv(2,0)*Jinv(0,0) + Jinv(2,1)*Jinv(0,1)) + J16*(Jinv(2,0)*Jinv(2,0) + Jinv(2,1)*Jinv(2,1)));

	// transpose the FeCF matrix to get the FCFe matrix
	Matrix FCFe(12,9);
	FCFe = Transpose(9,12,FeCF);*/
	
	// transpose the Kwu matrix
	Matrix KuT(12,9);
	KuT = Transpose(9,12,FeCFhg);

	// compute the stabilization stiffness matrix
	Kstab.Zero();

	/*Matrix KKF(12,12);
	Matrix FKK(12,12);
	Matrix KKFKK(12,12);

	KKF = KuT*FeCFeInv*FeCF;
	FKK = FCFe*FeCFeInv*FeCFhg;
	KKFKK = KuT*FeCFeInv*FeCFhg;*/

	Matrix interior(12,12);

	//interior = FCF - KKF - FKK + KKFKK;
	interior = FCF - KuT*FeCFeInv*FeCFhg;

	Kstab.addMatrixTripleProduct(1.0, Mben, interior, 1.0);
	
	return;
}

Vector
SSPbrick::CrossProduct(Vector v1, Vector v2)
// computes the cross product of two 3x1 vectors,  v1 x v2
{
	Vector result(3);

	result(0) = v1(1)*v2(2) - v1(2)*v2(1);
	result(1) = v1(2)*v2(0) - v1(0)*v2(2);
	result(2) = v1(0)*v2(1) - v1(1)*v2(0);

	return result;
}

Matrix  
SSPbrick::Transpose(int d1, int d2, const Matrix &M)
// transpose the input matrix
{
  	Matrix Mtran(d2,d1);

  	for (int i = 0; i < d1; i++) {
    	for (int j = 0; j < d2; j++) {
        	Mtran(j,i) = M(i,j);
		}
  	}

  	return Mtran;
}
