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
                                                                       
// Created: Chris McGann, UW
//          December 2011
//
// Description: This file contains the implementation of the BeamContact2Dp class

#include "BeamContact2Dp.h"

#include <elementAPI.h>
#include <Information.h>
#include <ElementResponse.h>
#include <ID.h>
#include <Domain.h>
#include <Node.h>
#include <Channel.h>
#include <Message.h>
#include <FEM_ObjectBroker.h>
#include <Renderer.h>
#include <G3Globals.h>
#include <ErrorHandler.h>
#include <ContactMaterial2D.h>
#include <Parameter.h>

#include <math.h>
#include <stdlib.h>
#include <stdio.h>

#define OPS_Export

static int num_BeamContact2Dp = 0;

OPS_Export void *
OPS_BeamContact2Dp(void)
{
  if (num_BeamContact2Dp == 0) {
    num_BeamContact2Dp++;
    opserr<<"BeamContact2Dp element - Written: C.McGann, P.Arduino, P.Mackenzie-Helnwein, U.Washington\n";
  }

  // Pointer to an element that will be returned
  Element *theElement = 0;

  int numRemainingInputArgs = OPS_GetNumRemainingInputArgs();

  if (numRemainingInputArgs < 7) {
    opserr << "Invalid #args, want: element BeamContact2Dp eleTag? iNode? jNode? secondaryNode? matTag? width? penalty? <cSwitch>?\n";
	return 0;
  }

  int iData[5];
  double dData[2];
  int icSwitch = 0;

  int numData = 5;
  if (OPS_GetIntInput(&numData, iData) != 0) {
    opserr << "WARNING invalid integer data: element BeamContact2Dp " << iData[0] << endln;
	return 0;
  }

  numData = 2;
  if (OPS_GetDoubleInput(&numData, dData) !=0) {
    opserr << "WARNING invalid data: element BeamContact2Dp " << iData[0] << endln;
	return 0;
  }
  
  int matID = iData[4];
  NDMaterial *theMaterial = OPS_getNDMaterial(matID);
  if (theMaterial == 0) {
    opserr << "WARNING element BeamContact2Dp " << iData[0] << endln;
	opserr << " Material: " << matID << "not found\n";
	return 0;
  }

  numRemainingInputArgs -= 7;
  while (numRemainingInputArgs >= 1) {
    numData = 1;
    if (OPS_GetIntInput(&numData, &icSwitch) != 0) {
      opserr << "WARNING invalid initial contact flag: element BeamContact2Dp " << iData[0] << endln;
	  return 0;
    }
	numRemainingInputArgs -= 1;
  }

  // Parsing was successful, allocate the element
  theElement = new BeamContact2Dp(iData[0], iData[1], iData[2], iData[3], *theMaterial, dData[0], dData[1],  icSwitch);

  if (theElement == 0) {
    opserr << "WARNING could not create element of type BeamContact2Dp\n";
	return 0;
  }

  return theElement;
}

// constructors
BeamContact2Dp::BeamContact2Dp(int tag, int Nd1, int Nd2, int NdS, NDMaterial &theMat, double width, double pen, int cSwitch)
  :Element(tag,ELE_TAG_BeamContact2Dp),
    theMaterial(0),
	mExternalNodes(BC2Dp_NUM_NODE),
	mTangentStiffness(BC2Dp_NUM_DOF, BC2Dp_NUM_DOF),
	mInternalForces(BC2Dp_NUM_DOF),
	mEye1(BC2Dp_NUM_DIM, BC2Dp_NUM_DIM),
	mEyeS(BC2Dp_NUM_DIM, BC2Dp_NUM_DIM),
	mg_xi(BC2Dp_NUM_DIM),
	mNormal(BC2Dp_NUM_DIM),
	mShape(4),
	mDshape(4),
	mBn(BC2Dp_NUM_DOF),
	mBs(BC2Dp_NUM_DOF),
	ma_1(BC2Dp_NUM_DIM),
	mb_1(BC2Dp_NUM_DIM),
	mc_1(BC2Dp_NUM_DIM),
	mIcrd_a(BC2Dp_NUM_DIM),
	mIcrd_b(BC2Dp_NUM_DIM),
	mIcrd_s(BC2Dp_NUM_DIM),
    mDcrd_a(BC2Dp_NUM_DIM),
	mDcrd_b(BC2Dp_NUM_DIM),
	mDcrd_s(BC2Dp_NUM_DIM),
	mDisp_a_n(3),
	mDisp_b_n(3)
{
    mExternalNodes(0) = Nd1;
	mExternalNodes(1) = Nd2;
	mExternalNodes(2) = NdS;
	
	// input parameters
	mRadius     = width/2.0;
	mPenalty    = pen;
	mIniContact = cSwitch;

	// initialize contact state
	if (mIniContact == 0) {
		inContact          = true;
		was_inContact      = true;
		in_bounds          = true;
	} else {
		inContact          = false;
		was_inContact      = false;
		in_bounds          = true;
	}

	// initialize penetration function and contact force
	mGap    = 0.0;
	mLambda = 0.0;

	// get copy of the material object
	NDMaterial *theMatCopy = theMat.getCopy("ContactMaterial2D");
	if (theMatCopy != 0) {
	  theMaterial = (ContactMaterial2D *)theMatCopy;
	} else {
	  opserr << "BeamContact2Dp::BeamContact2Dp - material needs to be ContactMaterial2D for ele: " << this->getTag() << endln;
	}

	// check material
	if (theMaterial == 0) {
	  opserr << "BeamContact2Dp::BeamContact2Dp - failed allocate material model pointer\n";
	  exit(-1);
	}
}

BeamContact2Dp::BeamContact2Dp()
  :Element(0,ELE_TAG_BeamContact2Dp),
    theMaterial(0),
	mExternalNodes(BC2Dp_NUM_NODE),
	mTangentStiffness(BC2Dp_NUM_DOF, BC2Dp_NUM_DOF),
	mInternalForces(BC2Dp_NUM_DOF),
	mEye1(BC2Dp_NUM_DIM, BC2Dp_NUM_DIM),
	mEyeS(BC2Dp_NUM_DIM, BC2Dp_NUM_DIM),
	mg_xi(BC2Dp_NUM_DIM),
	mNormal(BC2Dp_NUM_DIM),
	mShape(4),
	mDshape(4),
	mBn(BC2Dp_NUM_DOF),
	mBs(BC2Dp_NUM_DOF),
	ma_1(BC2Dp_NUM_DIM),
	mb_1(BC2Dp_NUM_DIM),
	mc_1(BC2Dp_NUM_DIM),
	mIcrd_a(BC2Dp_NUM_DIM),
	mIcrd_b(BC2Dp_NUM_DIM),
	mIcrd_s(BC2Dp_NUM_DIM),
    mDcrd_a(BC2Dp_NUM_DIM),
	mDcrd_b(BC2Dp_NUM_DIM),
	mDcrd_s(BC2Dp_NUM_DIM),
	mDisp_a_n(3),
	mDisp_b_n(3)
{
}

// destructor
BeamContact2Dp::~BeamContact2Dp()
{
    if (theMaterial != 0) {
        delete theMaterial;
    }
}

int 
BeamContact2Dp::getNumExternalNodes(void) const
{
    return BC2Dp_NUM_NODE;
}

const ID &
BeamContact2Dp::getExternalNodes(void)
{
    return mExternalNodes;
}

Node **
BeamContact2Dp::getNodePtrs(void)
{
    return theNodes;
}

int BeamContact2Dp::getNumDOF(void)
{
    return BC2Dp_NUM_DOF;
}

void
BeamContact2Dp::setDomain(Domain *theDomain)
{
	Vector x_c(BC2Dp_NUM_DIM);
	
	mEye1.Zero();
	mEye1(0,0) = 1.0;
	mEye1(1,1) = 1.0;

	mEyeS.Zero();
	mEyeS(0,1) = -1.0;
	mEyeS(1,0) = 1.0;

	theNodes[0] = theDomain->getNode(mExternalNodes(0));
	theNodes[1] = theDomain->getNode(mExternalNodes(1));
	theNodes[2] = theDomain->getNode(mExternalNodes(2));

	for (int i = 0; i < 3; i++) {
        if (theNodes[i] == 0)
            return;  // don't go any further - otherwise segmentation fault
    }

	// initialize coordinate vectors
	mIcrd_a = theNodes[0]->getCrds();
	mIcrd_b = theNodes[1]->getCrds();
	mIcrd_s = theNodes[2]->getCrds();
	mDcrd_a = mIcrd_a;
	mDcrd_b = mIcrd_b;
	mDcrd_s = mIcrd_s;
	mDisp_a_n.Zero();
	mDisp_b_n.Zero();

    // length of beam element
	mLength = (mDcrd_b - mDcrd_a).Norm();

	// initialize tangent vectors at beam nodes
	ma_1 = (mDcrd_b - mDcrd_a)/mLength;
	mb_1 = ma_1;

	// perform projection of secondary node to beam centerline
	mXi = ((mDcrd_b - mDcrd_s)^(mDcrd_b - mDcrd_a))/mLength;  // initial assumption
	mXi = Project(mXi);                                       // actual location

	// update contact state based on projection
	in_bounds = ((mXi > 0.000) && (mXi < 1.000));
	inContact = (was_inContact && in_bounds);

	// centerline projection coordinate
	x_c = mDcrd_a*mShape(0) + ma_1*mLength*mShape(1) + mDcrd_b*mShape(2) + mb_1*mLength*mShape(3);

	// update surface tangent vector, g_xi
	UpdateBase(mXi);

	// adjust cohesion force
	theMaterial->ScaleCohesion(mLength);
	theMaterial->ScaleTensileStrength(mLength);

	// compute vectors Bn and Bs
	ComputeB();

	// call the base class method
	this->DomainComponent::setDomain(theDomain);
}

int
BeamContact2Dp::commitState()
{
	// update projection
	mXi = Project(mXi);

	// update surface tangent vector, g_xi
	UpdateBase(mXi);

	// update vectors Bn and Bs for next step
	ComputeB();

	// update contact state 
    double tol = 0.000001*mRadius;
	was_inContact = (mGap < tol);
	in_bounds     = ((mXi > 0.000) && (mXi < 1.000));
	inContact     = (was_inContact && in_bounds);

	int retVal = 0;
	// call element commitState to do any base class stuff
	if ((retVal = this->Element::commitState()) != 0) {
	    opserr << "BeamContact2Dp::commitState() - failed in base class";
	}
	retVal = theMaterial->commitState();

	return retVal;
}

int
BeamContact2Dp::revertToLastCommit()
{
    return theMaterial->revertToLastCommit();
}

int
BeamContact2Dp::revertToStart()
{
	if (mIniContact == 0) {
		inContact          = true;
		was_inContact      = true;
		in_bounds          = true;
	} else {
		inContact          = false;
		was_inContact      = false;
		in_bounds          = true;
	}

	// reset applicable member variables 
	mDcrd_a = mIcrd_a;
	mDcrd_b = mIcrd_b;
	mDcrd_s = mIcrd_s;
	mDisp_a_n.Zero();
	mDisp_b_n.Zero();

	mLength = (mDcrd_b - mDcrd_a).Norm();

	ma_1 = (mDcrd_b - mDcrd_a)/mLength;
	mb_1 = ma_1;

	mXi = ((mDcrd_b - mDcrd_s)^(mDcrd_b - mDcrd_a))/mLength;
	mXi = Project(mXi);

	in_bounds = ((mXi > 0.000) && (mXi < 1.000));
	inContact = (was_inContact && in_bounds);

	UpdateBase(mXi);
	ComputeB();

	return theMaterial->revertToStart();
}

int
BeamContact2Dp::update(void)
// this function updates variables for an incremental step n to n+1
{
    double tensileStrength;
	Vector a1(BC2Dp_NUM_DIM);
    Vector b1(BC2Dp_NUM_DIM);
	Vector a1_n(BC2Dp_NUM_DIM);
    Vector b1_n(BC2Dp_NUM_DIM);
    Vector disp_a(3);
    Vector disp_b(3);
    Vector disp_L(BC2Dp_NUM_DIM);
    double rot_a;
    double rot_b;
    Vector x_c(BC2Dp_NUM_DIM);

	// update secondary node coordinates
	mDcrd_s = mIcrd_s + theNodes[2]->getTrialDisp();

	// update nodal coordinates
	disp_a = theNodes[0]->getTrialDisp();
	disp_b = theNodes[1]->getTrialDisp();

	for (int i = 0; i < 2; i++) {
	    mDcrd_a(i) = mIcrd_a(i) + disp_a(i);
		mDcrd_b(i) = mIcrd_b(i) + disp_b(i);
	}

	// compute incremental rotation from step n to step n+1
	rot_a = disp_a(2) - mDisp_a_n(2);
	rot_b = disp_b(2) - mDisp_b_n(2);

	// get tangent vectors from last converged step
	a1_n = Geta1();
	b1_n = Getb1();

	// linear update of tangent vectors
	a1 = a1_n + rot_a*mEyeS*a1_n;
	b1 = b1_n + rot_b*mEyeS*b1_n;

	// update centerline projection coordinate
	x_c = mDcrd_a*mShape(0) + a1*mLength*mShape(1) + mDcrd_b*mShape(2) + b1*mLength*mShape(3);

	// update penetration function
	mGap = (mNormal^(mDcrd_s - x_c)) - mRadius;
    double tol = 0.000001*mRadius;
	if (mGap < tol && in_bounds) {
		inContact = true;
	} else {
		//mGap = 0.0;
		inContact = false;
	}

	// update normal contact force
    if (inContact) {
        // compute normal contact force (penalty formulation)
		mLambda = mPenalty*mGap;
        // get tensile strength from contact material
	    tensileStrength = theMaterial->getTensileStrength();
        // check that contact force is < surface adhesion (tensile strength)
        if (mLambda > tensileStrength) {
            mLambda = 0.0;
        }
	} else {
		mLambda = 0.0;
	}

    // determine trial strain vector based on contact state
	if (inContact) {
	    Vector strain(3);
		double slip;
		Vector c1n1(2);
		Vector c2n1(2);

        // tangent at the centerline projection in step n+1
		c1n1 = mDshape(0)*mDcrd_a + mDshape(1)*mLength*ma_1 + mDshape(2)*mDcrd_b + mDshape(3)*mLength*mb_1;

		// update vector c2 for step n+1
		c2n1 = (mDcrd_s - x_c)/((mDcrd_s - x_c).Norm());
		
		// update vector c2 for step n+1
		c2n1(0) = -c1n1(1);
		c2n1(1) = c1n1(0);

		// compute the slip
		slip = mg_xi^(mDcrd_s - x_c - mrho*c2n1);

		// set the strain vector
		strain(0) = mGap;
		strain(1) = slip;
		strain(2) = -mLambda;
		
		theMaterial->setTrialStrain(strain);
	} else {
	    Vector strain(3);

        // set the strain vector
		strain(0) = mGap;
		strain(1) = 0.0;
		strain(2) = -mLambda;
		
		theMaterial->setTrialStrain(strain);
	}

	return 0;
}

double
BeamContact2Dp::Project(double xi)
// this function computes the centerline projection for the current step
{
    double xi_p;
	double H1;
	double H2;
	double H3;
	double H4;
    double dH1;
	double dH2;
	double dH3;
	double dH4;
	double R;
	double DR;
	double dxi;
	Vector a1(BC2Dp_NUM_DIM);
    Vector b1(BC2Dp_NUM_DIM);
	Vector x_c_p(BC2Dp_NUM_DIM);
	Vector t_c(BC2Dp_NUM_DIM);
	Vector ddx_c(BC2Dp_NUM_DIM);

	// initialize to previous projection location
	xi_p = xi;

	// update end point tangents
	UpdateEndFrames();

	// set tangent vectors
	a1 = Geta1();
	b1 = Getb1();

	// Hermitian basis functions and first derivatives
	H1 = 1.0 - 3.0*xi_p*xi_p + 2.0*xi_p*xi_p*xi_p;
	H2 = xi_p - 2.0*xi_p*xi_p + xi_p*xi_p*xi_p;
	H3 = 3.0*xi_p*xi_p - 2*xi_p*xi_p*xi_p;
	H4 = -xi_p*xi_p + xi_p*xi_p*xi_p;
	dH1 = -6.0*xi_p + 6.0*xi_p*xi_p;
	dH2 = 1.0 - 4.0*xi_p + 3.0*xi_p*xi_p;
	dH3 = 6.0*xi_p - 6.0*xi_p*xi_p;
	dH4 = -2.0*xi_p + 3.0*xi_p*xi_p;

    // compute current projection coordinate and tangent
	x_c_p = mDcrd_a*H1 + a1*mLength*H2 + mDcrd_b*H3 + b1*mLength*H4;
	t_c   = mDcrd_a*dH1 + a1*mLength*dH2 + mDcrd_b*dH3 + b1*mLength*dH4;
	
	// compute initial value of residual
	R = (mDcrd_s - x_c_p)^t_c;

	// iterate to determine new value of xi
	int Gapcount = 0;
	while (fabs(R/mLength) > 1.0e-10 && Gapcount < 50) {
	
		// compute current curvature vector
		ddx_c = Get_dxc_xixi(xi_p);

		// increment projection location
		DR   = ((mDcrd_s - x_c_p)^ddx_c) - (t_c^t_c);
		dxi  = -R/DR;
		xi_p = xi_p + dxi;

		// Hermitian basis functions and first derivatives
	    H1 = 1.0 - 3.0*xi_p*xi_p + 2.0*xi_p*xi_p*xi_p;
    	H2 = xi_p - 2.0*xi_p*xi_p + xi_p*xi_p*xi_p;
    	H3 = 3.0*xi_p*xi_p - 2*xi_p*xi_p*xi_p;
    	H4 = -xi_p*xi_p + xi_p*xi_p*xi_p;
    	dH1 = -6.0*xi_p + 6.0*xi_p*xi_p;
    	dH2 = 1.0 - 4.0*xi_p + 3.0*xi_p*xi_p;
    	dH3 = 6.0*xi_p - 6.0*xi_p*xi_p;
    	dH4 = -2.0*xi_p + 3.0*xi_p*xi_p;

		// update projection coordinate and tangent
		x_c_p = mDcrd_a*H1 + a1*mLength*H2 + mDcrd_b*H3 + b1*mLength*H4;
	    t_c   = mDcrd_a*dH1 + a1*mLength*dH2 + mDcrd_b*dH3 + b1*mLength*dH4;

		// compute residual
    	R = (mDcrd_s - x_c_p)^t_c;

		Gapcount += 1;
	}

	// update normal vector for current projection
	mNormal = (mDcrd_s - x_c_p)/((mDcrd_s - x_c_p).Norm());

	// update Hermitian basis functions and derivatives
	mShape(0)  = H1;
	mShape(1)  = H2;
	mShape(2)  = H3;
	mShape(3)  = H4;
	mDshape(0) = dH1;
	mDshape(1) = dH2;
	mDshape(2) = dH3;
	mDshape(3) = dH4;

	return xi_p;
}

int
BeamContact2Dp::UpdateBase(double xi)
// this function computes the surface tangent vector g_xi
{
    Vector t_c(BC2Dp_NUM_DIM);
	Vector ddx_c(BC2Dp_NUM_DIM);
	Vector d_c1(BC2Dp_NUM_DIM);
	Vector c_2(BC2Dp_NUM_DIM);
	Vector d_c2(BC2Dp_NUM_DIM);
	
    // compute current projection tangent
	t_c = Get_dxc_xi(xi);

	// set value of unit tangent vector c1
	mc_1 = t_c/t_c.Norm();

	// compute current projection curvature
	ddx_c = Get_dxc_xixi(xi);

	// compute derivative of c1 with respect to xi
	d_c1 = (1/t_c.Norm())*(ddx_c - ((mc_1^ddx_c)*mc_1));

	// determine orthogonal vector c2
	c_2(0) = -mc_1(1);  c_2(1) = mc_1(0);

	// compute local coordinate transformation term
	mrho = mRadius*(mNormal^c_2);

	// compute change in c2 with xi
	d_c2 = -(d_c1^c_2)*mc_1;

    // compute surface tangent vector
	mg_xi = t_c + mrho*d_c2;

	return 0;
}

Vector
BeamContact2Dp::Get_dxc_xi(double xi)
// this function computes the first derivative of x_c wrt xi
{
    double dH1;
	double dH2;
	double dH3;
	double dH4;
	Vector a1(BC2Dp_NUM_DIM);
	Vector b1(BC2Dp_NUM_DIM);
	Vector dx(BC2Dp_NUM_DIM);

	// first derivatives of Hermitian basis functions
    dH1 = -6.0*xi + 6.0*xi*xi;
	dH2 = 1.0 - 4.0*xi + 3.0*xi*xi;
	dH3 = 6.0*xi - 6.0*xi*xi;
	dH4 = -2.0*xi + 3.0*xi*xi;

	// tangent vectors
	a1 = Geta1();
	b1 = Getb1();

    // compute current projection tangent
	dx = mDcrd_a*dH1 + a1*mLength*dH2 + mDcrd_b*dH3 + b1*mLength*dH4;

	return dx;
}

Vector
BeamContact2Dp::Get_dxc_xixi(double xi)
// this function computes the second derivative of x_c wrt xi
{
    double ddH1;
	double ddH2;
	double ddH3;
	double ddH4;
	Vector a1(BC2Dp_NUM_DIM);
	Vector b1(BC2Dp_NUM_DIM);
	Vector ddx(BC2Dp_NUM_DIM);
	
	// second derivatives of Hermitian basis functions
	ddH1 = -6.0 + 12.0*xi;
	ddH2 = -4.0 + 6.0*xi;
	ddH3 =  6.0 - 12.0*xi;
	ddH4 = -2.0 + 6.0*xi;

	// tangent vectors
	a1 = Geta1();
	b1 = Getb1();

	// compute current curvature vector
	ddx = mDcrd_a*ddH1 + a1*mLength*ddH2 + mDcrd_b*ddH3 + b1*mLength*ddH4;

	return ddx;
}

void
BeamContact2Dp::ComputeB(void)
// this function computes the finite element equation vectors Bn and Bs
{
    double Ka1n;
	double Kb1n;
	double Ka1g;
	double Kb1g;
	Vector a1(BC2Dp_NUM_DIM);
	Vector b1(BC2Dp_NUM_DIM);

	// initialize Bn and Bs
	mBn.Zero();
	mBs.Zero();

	// get tangent vectors
	a1 = Geta1();
	b1 = Getb1();

	// compute terms for vector Bn
    Ka1n = (mEyeS*a1)^mNormal;
	Kb1n = (mEyeS*b1)^mNormal;
	
	mBn(0) = -mShape(0)*mNormal(0);
	mBn(1) = -mShape(0)*mNormal(1);
	mBn(2) = -mShape(1)*mLength*Ka1n;
	mBn(3) = -mShape(2)*mNormal(0);
	mBn(4) = -mShape(2)*mNormal(1);
	mBn(5) = -mShape(3)*mLength*Kb1n;
	mBn(6) = mNormal(0);
	mBn(7) = mNormal(1);

	// compute terms for vector Bs
	Ka1g = (mEyeS*a1)^mg_xi;
	Kb1g = (mEyeS*b1)^mg_xi;

	mBs(0) = -mg_xi(0)*(mShape(0) + mRadius*mDshape(0));
	mBs(1) = -mg_xi(1)*(mShape(0) + mRadius*mDshape(0));
	mBs(2) = -Ka1g*mLength*(mShape(1) + mRadius*mDshape(1));
	mBs(3) = -mg_xi(0)*(mShape(2) + mRadius*mDshape(2));
	mBs(4) = -mg_xi(1)*(mShape(2) + mRadius*mDshape(2));
	mBs(5) = -Kb1g*mLength*(mShape(3) + mRadius*mDshape(3));
    mBs(6) = mg_xi(0);
	mBs(7) = mg_xi(1);

	return;
}

void
BeamContact2Dp::UpdateEndFrames(void)
// this function updates the tangent and displacement at nodes a and b
{
	Vector disp_a(3);
	Vector disp_b(3);
	double rot_a;
	double rot_b;

	// recalculate incremental rotations from step n to n+1
	disp_a = theNodes[0]->getTrialDisp();
	disp_b = theNodes[1]->getTrialDisp();
	rot_a  = disp_a(2) - mDisp_a_n(2);
	rot_b  = disp_b(2) - mDisp_b_n(2);

	// perform update for node a
	ma_1 += rot_a*mEyeS*ma_1;

	// perform update for node b
	mb_1 += rot_b*mEyeS*mb_1;

	// update displacement vectors for next step
	mDisp_a_n = disp_a;
	mDisp_b_n = disp_b;

	return;
} 

Vector
BeamContact2Dp::Geta1(void)
// this function returns the tangent vector at node a from last converged step
{
	Vector a1(BC2Dp_NUM_DIM);

	a1 = ma_1;

	return a1;
}

Vector
BeamContact2Dp::Getb1(void)
// this function returns the tangent vector at node a from last converged step
{
	Vector b1(BC2Dp_NUM_DIM);

	b1 = mb_1;

	return b1;
}

const Matrix &
BeamContact2Dp::getTangentStiff(void)
// this function computes the tangent stiffness matrix for the element
{
	mTangentStiffness.Zero();

	if (inContact) {
		Matrix Cmat = theMaterial->getTangent();

		double Css = Cmat(1,1);
		double Csn = Cmat(1,2);

		for (int i = 0; i < BC2Dp_NUM_DOF; i++) {
			for (int j = 0; j < BC2Dp_NUM_DOF; j++) {

				mTangentStiffness(i,j) = mBs(i)*mBs(j)*Css - mPenalty*(Csn*mBs(i)*mBn(j) - mBn(i)*mBn(j));
			}
		}
	}

	return mTangentStiffness;
}

const Matrix &
BeamContact2Dp::getInitialStiff(void)
// this function computes the initial tangent stiffness matrix for the element
{
	return getTangentStiff();
}

void
BeamContact2Dp::zeroLoad(void)
{
	return;
}

int
BeamContact2Dp::addLoad(ElementalLoad *theLoad, double loadFactor)
{
	return 0;
}

int
BeamContact2Dp::addInertiaLoadToUnbalance(const Vector &accel)
{
	return 0;
}

const Vector &
BeamContact2Dp::getResistingForce()
// this function computes the resisting force vector for the element
{
	mInternalForces.Zero();

	// get contact "stress" vector
	Vector stress = theMaterial->getStress();

	if (inContact) {
		mInternalForces = mLambda*mBn + stress(1)*mBs;
	}

	return mInternalForces;
}

const Vector &
BeamContact2Dp::getResistingForceIncInertia()
{
	return getResistingForce();
}

int
BeamContact2Dp::sendSelf(int commitTag, Channel &theChannel)
{
	int res = 0;

	// NOTE: we don't check for dataTag == 0 for Element
    // objects as that is taken care of in a commit by the Domain
    // object - don't want to have to do the check if sending data
	int dataTag = this->getDbTag();

	// BeamContact2Dp packs it's data into a Vector and sends this to theChannel
    // along with it's dbTag and the commitTag passed in the arguments
	static Vector data(6);
	data(0) = this->getTag();
	data(1) = mRadius;
	data(2) = mPenalty;
	data(3) = mIniContact;
	data(4) = theMaterial->getClassTag();

	int matDbTag = theMaterial->getDbTag();

	// NOTE: we have to ensure that the material has a database tag
	// if we are sending to a database channel
	if (matDbTag == 0) {
		matDbTag = theChannel.getDbTag();
		if (matDbTag != 0)
			theMaterial->setDbTag(matDbTag);
	}
	data(5) = matDbTag;

	res = theChannel.sendVector(dataTag, commitTag, data);
	if (res < 0) {
		opserr << "WARNING BeamContact2Dp::sendSelf() - " << this->getTag() << " failed to send Vector\n";
		return -1;
	}

	// BeamContact2Dp then sends the tags of its four nodes
	res = theChannel.sendID(dataTag, commitTag, mExternalNodes);
	if (res < 0) {
		opserr << "WARNING BeamContact2Dp::sendSelf() - " << this->getTag() << " failed to send Vector\n";
		return -2;
	}

	// finally, BeamContact2Dp asks its material object to send itself
	res = theMaterial->sendSelf(commitTag, theChannel);
	if (res < 0) {
		opserr << "WARNING BeamContact2Dp::sendSelf() - " << this->getTag() << " failed to send its Material\n";
		return -3;
	}

	return 0;
}

int
BeamContact2Dp::recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker)
{
	int res = 0;
	int dataTag = this->getDbTag();

	// BeamContact2Dp creates a vector, receives the vector, and then sets the internal
	// data with the data in the vector
	static Vector data(6);
	res = theChannel.recvVector(dataTag, commitTag, data);
	if (res < 0) {
		opserr << "WARNING BeamContact2Dp::recvSelf() - failed to receive Vector\n";
		return -1;
	}

	this->setTag((int)data(0));
	mRadius     = data(1);
	mPenalty    = data(2);
	mIniContact = (int)data(3);

	// BeamContact2Dp now receives the tags of its four external nodes
	res = theChannel.recvID(dataTag, commitTag, mExternalNodes);
	if (res < 0) {
		opserr << "WARNING BeamContact2Dp::recvSelf() - " << this->getTag() << " failed to receive ID\n";
		return -2;
	}

	// finally, BeamContact2Dp creates a material object of the correct type, sets its
	// database tag, and asks this new object to receive itself
	int matClass = (int)data(4);
	int matDb    = (int)data(5);

	// check if material object exists and that it is the right type
	if ((theMaterial == 0) || (theMaterial->getClassTag() != matClass)) {

		// if old one, delete it
		if (theMaterial != 0)
			delete theMaterial;

		// create new material object
		NDMaterial *theMatCopy = theBroker.getNewNDMaterial(matClass);
		theMaterial = (ContactMaterial2D *)theMatCopy;

		if (theMaterial == 0) {
			opserr << "WARNING BeamContact2Dp::recvSelf() - " << this->getTag() 
			  << " failed to get a blank Material of type " << matClass << endln;
			return -3;
		}
	}

	// NOTE: we set the dbTag before we receive the material
	theMaterial->setDbTag(matDb);
	res = theMaterial->recvSelf(commitTag, theChannel, theBroker);
	if (res < 0) {
		opserr << "WARNING BeamContact2Dp::recvSelf() - " << this->getTag() << " failed to receive its Material\n";
		return -3;
	}

	return 0;
}

int
BeamContact2Dp::displaySelf(Renderer &theViewer, int displayMode, float fact, const char **modes, int numMode)
{
	return 0;
}

void
BeamContact2Dp::Print(OPS_Stream &s, int flag)
{
	opserr << "BeamContact2Dp, element id:  " << this->getTag() << endln;
	opserr << "   Connected external nodes:  ";
	for (int i = 0; i<BC2Dp_NUM_NODE; i++)
	{
		opserr << mExternalNodes(i) << " ";
	}
	return;
}

Response*
BeamContact2Dp::setResponse(const char **argv, int argc, OPS_Stream &eleInfo)
{
	if (strcmp(argv[0],"force") == 0 || strcmp(argv[0],"forces") == 0) {
		// forces on secondary node
		return new ElementResponse(this, 1, Vector(2));

	} else if (strcmp(argv[0],"frictionforce") == 0 || strcmp(argv[0],"frictionforces") == 0) {
		// frictional force vector
		return new ElementResponse(this, 2, Vector(2));

	} else if (strcmp(argv[0],"forcescalar") == 0 || strcmp(argv[0],"forcescalars") == 0) {
		// scalar contact forces
		return new ElementResponse(this, 3, Vector(2));

	} else if (strcmp(argv[0],"masterforce") == 0 || strcmp(argv[0],"masterforces") == 0 ||
		   strcmp(argv[0],"primaryforce") == 0 || strcmp(argv[0],"primaryforces") == 0) {
		// reactions (forces and moments) on primary nodes
		return new ElementResponse(this, 4, Vector(6));

    } else if (strcmp(argv[0],"gap") == 0) {
        // magnitude of gap
        return new ElementResponse(this, 5, Vector(1));

	} else {
		// otherwise response quantity is unknown for the BeamContact2Dp class
		opserr << "BeamContact2Dp::setResponse(const char **argv, int argc, OPS_Stream &eleInfo): "
		  << argv[0] << " unknown recorder request" << endln;
	return 0;
	}
}

int
BeamContact2Dp::getResponse(int responseID, Information &eleInfo)
{
	// initialize variables
	Vector force(2);
	Vector frictForce(2);
	Vector secondaryForce(2);
	Vector primaryForce(6);
    Vector theGap(1);

	// get contact "stress" vector
	Vector stress = theMaterial->getStress();

	if (responseID == 1) {

		// forces on secondar node
		secondaryForce(0) = -mInternalForces(6);
		secondaryForce(1) = -mInternalForces(7);
		return eleInfo.setVector(secondaryForce);
	
	} else if (responseID == 2) {

		// frictional force vector
		frictForce = stress(1)*mg_xi;
		return eleInfo.setVector(frictForce);

    } else if (responseID == 3) {
		
		// scalar contact forces
		force(0) = stress(0);
		force(1) = stress(1);
		return eleInfo.setVector(force);

	} else if (responseID == 4) {

		// reactions (forces and moments) on primary nodes
		for (int i = 0;  i < 3; i++) {

			primaryForce(i)   = -mInternalForces(i);
			primaryForce(i+3) = -mInternalForces(i+3);
		}
		return eleInfo.setVector(primaryForce);

    } else if (responseID == 5) {
        theGap(0) = mGap;
        return eleInfo.setVector(theGap);

	} else {
		// otherwise response quantity is unknown for the BeamContact2Dp class
		opserr << "BeamContact2Dp::getResponse(int responseID = " << responseID << ", Information &eleInfo); "
		  << " unknown request" << endln;
		return -1;
	}
}

int
BeamContact2Dp::setParameter(const char **argv, int argc, Parameter &param)
{
	if (argc < 1)
		return -1;

	if (strcmp(argv[0],"friction") == 0) {
		return param.addObject(1, this);
	}
	
	return -1;
}

int
BeamContact2Dp::updateParameter(int parameterID, Information &info)
{
	int res = -1;
	int matRes =  theMaterial->updateParameter(parameterID, info);
	if (matRes != -1) {
		res = matRes;
	}
	return res;
}
