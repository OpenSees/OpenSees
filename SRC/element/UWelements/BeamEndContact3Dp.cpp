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
                                                                       
// Created: Chris McGann, UW, 01.2011
//
// Description: This file contains the implementation of the BeamEndContact3Dp class

#include "BeamEndContact3Dp.h"

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
#include <Parameter.h>

#include <math.h>
#include <stdlib.h>
#include <stdio.h>

#define OPS_Export

static int num_BeamEndContact3Dp = 0;

OPS_Export void *
OPS_BeamEndContact3Dp(void)
{
  	if (num_BeamEndContact3Dp == 0) {
    	num_BeamEndContact3Dp++;
    	opserr<<"BeamEndContact3Dp element - Written: C.McGann, P.Arduino, P.Mackenzie-Helnwein, U.Washington\n";
  	}

  	// Pointer to an element that will be returned
  	Element *theElement = 0;

  	int numRemainingInputArgs = OPS_GetNumRemainingInputArgs();

  	if (numRemainingInputArgs < 6) {
    	opserr << "Invalid #args, want: element BeamEndContact3Dp eleTag? iNode? jNode? sNode? radius? penalty? <cFlag>?\n";
		return 0;
  	}

  	int iData[4];
  	double dData[2];
  	int icSwitch = 0;

  	int numData = 4;
  	if (OPS_GetIntInput(&numData, iData) != 0) {
    	opserr << "WARNING invalid integer data: element BeamEndContact3Dp " << iData[0] << endln;
		return 0;
  	}

  	numData = 2;
  	if (OPS_GetDoubleInput(&numData, dData) !=0) {
    	opserr << "WARNING invalid double data: element BeamEndContact3Dp " << iData[0] << endln;
		return 0;
  	}
  
  	numRemainingInputArgs -= 6;
  	while (numRemainingInputArgs >= 1) {
    	numData = 1;
    	if (OPS_GetIntInput(&numData, &icSwitch) != 0) {
      		opserr << "WARNING invalid initial contact flag: element BeamEndContact3Dp " << iData[0] << endln;
	  		return 0;
    	}
		numRemainingInputArgs -= 1;
  	}

  	// Parsing was successful, allocate the element
  	theElement = new BeamEndContact3Dp(iData[0], iData[1], iData[2], iData[3], dData[0], dData[1], icSwitch);

  	if (theElement == 0) {
    	opserr << "WARNING could not create element of type BeamEndContact3DpElement\n";
		return 0;
  	}

  	return theElement;
}

// constructors
BeamEndContact3Dp::BeamEndContact3Dp(int tag, int Nd1, int Nd2, int NdS, double rad, double pen, int cSwitch)
  :Element(tag,ELE_TAG_BeamEndContact3Dp),
	mExternalNodes(BEC3p_NUM_NODE),
	mTangentStiffness(BEC3p_NUM_DOF, BEC3p_NUM_DOF),
	mInternalForces(BEC3p_NUM_DOF),
	mEye1(BEC3p_NUM_DIM, BEC3p_NUM_DIM),
	mIniNormal(BEC3p_NUM_DIM),
	mNormal(BEC3p_NUM_DIM),
	mIcrd_a(BEC3p_NUM_DIM),
	mIcrd_s(BEC3p_NUM_DIM),
    mDcrd_a(BEC3p_NUM_DIM),
	mDcrd_s(BEC3p_NUM_DIM)
{
    mExternalNodes(0) = Nd1;
	mExternalNodes(1) = NdS;
	mBeamNode = Nd2;

    mRadius = rad;
	mPenalty = pen;
	mIniContact = cSwitch;

	mGap    = 0.0;
	mLambda = 0.0;

	// set the initial contact state
	if (mIniContact == 0) {
		inContact          = true;
		was_inContact      = true;
		in_bounds          = true;
	} else {
		inContact          = false;
		was_inContact      = false;
		in_bounds          = true;
	}
}

BeamEndContact3Dp::BeamEndContact3Dp()
  :Element(0,ELE_TAG_BeamEndContact3Dp),
	mExternalNodes(BEC3p_NUM_NODE),
	mTangentStiffness(BEC3p_NUM_DOF, BEC3p_NUM_DOF),
	mInternalForces(BEC3p_NUM_DOF),
	mEye1(BEC3p_NUM_DIM, BEC3p_NUM_DIM),
	mIniNormal(BEC3p_NUM_DIM),
	mNormal(BEC3p_NUM_DIM),
	mIcrd_a(BEC3p_NUM_DIM),
	mIcrd_s(BEC3p_NUM_DIM),
    mDcrd_a(BEC3p_NUM_DIM),
	mDcrd_s(BEC3p_NUM_DIM)
{
}

// destructor
BeamEndContact3Dp::~BeamEndContact3Dp()
{
}

int 
BeamEndContact3Dp::getNumExternalNodes(void) const
{
    return BEC3p_NUM_NODE;
}

const ID &
BeamEndContact3Dp::getExternalNodes(void)
{
    return mExternalNodes;
}

Node **
BeamEndContact3Dp::getNodePtrs(void)
{
    return theNodes;
}

int 
BeamEndContact3Dp::getNumDOF(void)
{
    return BEC3p_NUM_DOF;
}

void
BeamEndContact3Dp::setDomain(Domain *theDomain)
{
	double r;
	Vector a1(BEC3p_NUM_DIM);
	
	mEye1.Zero();
	mEye1(0,0) = 1.0;
	mEye1(1,1) = 1.0;
	mEye1(2,2) = 1.0;

	theNodes[0] = theDomain->getNode(mExternalNodes(0));
	theNodes[1] = theDomain->getNode(mExternalNodes(1));

	for (int i = 0; i < 2; i++) {
        if (theNodes[i] == 0)
            return;  // don't go any further - otherwise segmentation fault
    }

	// initialize coordinate vectors
	mIcrd_a = theNodes[0]->getCrds();
	mIcrd_s = theNodes[1]->getCrds();
	mDcrd_a = mIcrd_a;
	mDcrd_s = mIcrd_s;
	x_b = theDomain->getNode(mBeamNode);
	mIcrd_b = x_b->getCrds();

	// initialize the normal vector
	mIniNormal = -1*(mIcrd_b - mIcrd_a)/(mIcrd_b - mIcrd_a).Norm();
	mNormal = mIniNormal;

	// determine initial gap
	mGap = (mDcrd_s - mDcrd_a)^mIniNormal;

	// determine initial projection
	mx_p = mDcrd_s - (mGap*mIniNormal);

	// initialize contact state based on projection
	r = (mx_p - mIcrd_a).Norm();
	in_bounds = (r <= mRadius);
	inContact = (was_inContact && in_bounds);

	// call the base class method
	this->DomainComponent::setDomain(theDomain);
}

int
BeamEndContact3Dp::commitState()
{
	//update contact state 
	double tol = 0.000001*mRadius;
	was_inContact =  (mGap < tol);
	in_bounds =      (((mx_p - mDcrd_a).Norm()) <= mRadius);
	inContact =      (was_inContact && in_bounds);

	int retVal = 0;
	// call element commitState to do any base class stuff
	if ((retVal = this->Element::commitState()) != 0) {
	    opserr << "BeamEndContact3Dp::commitState() - failed in base class";
	}
	return 0;
}

int
BeamEndContact3Dp::revertToLastCommit()
{
	return 0;
}

int
BeamEndContact3Dp::revertToStart()
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
	return 0;
}

int
BeamEndContact3Dp::update(void)
// this function updates variables for an incremental step n to n+1
{
	Vector disp_a(6);
	Vector rot_a(BEC3p_NUM_DIM);
    Vector omega(BEC3p_NUM_DIM);
	Matrix eMap(BEC3p_NUM_DIM,BEC3p_NUM_DIM);

	// update beam node coordinates and rotations
	disp_a = theNodes[0]->getTrialDisp();
	for (int i = 0; i < 3; i++) {
		mDcrd_a(i) = mIcrd_a(i) + disp_a(i);
		rot_a(i) = disp_a(i+3);
	}

	// update secondary node coordinates
	mDcrd_s = mIcrd_s + theNodes[1]->getTrialDisp();

	// update the normal vector
	omega = rot_a - ((mNormal^rot_a)*mNormal);
	// compute exponential map of skew(omega)
	eMap = ExpMap(omega);
	// compute the current normal vector
	mNormal = eMap*mIniNormal;

	// update penetration function
	mGap = (mDcrd_s - mDcrd_a)^mNormal;
	double tol = 0.000001*mRadius;
	if (mGap < tol && in_bounds) {
		inContact = true;
	} else {
		mGap = 0.0;
		inContact = false;
	}

	// update normal contact force
	if (was_inContact) {
		mLambda = mPenalty*mGap;
	} else {
		mLambda = 0.0;
	}

	// update projection coordinate
    mx_p = mDcrd_s - (mGap*mNormal);

	return 0;
}

Matrix
BeamEndContact3Dp::ExpMap(Vector th)
// this function computes the exp map of the given vector
{
	double sf1;
    double sf2;
    double sf3;
    double theta;                                  // vector norm
    Vector theta_vec(BEC3p_NUM_DIM);                // input vector
    Matrix sk_theta(BEC3p_NUM_DIM,BEC3p_NUM_DIM);    // skew of vector
    Matrix theta_theta(BEC3p_NUM_DIM,BEC3p_NUM_DIM); // dyadic product of vector
    Matrix Q(BEC3p_NUM_DIM,BEC3p_NUM_DIM);           // Exonential Map Vector  

	// initialize theta variables
    Q.Zero();
    sk_theta.Zero();
    theta_theta.Zero();
       
    theta_vec = th;
    theta = theta_vec.Norm();
    sk_theta = GetSkew(theta_vec);
    int i, j;
    for (int i = 0; i<3; i++) {
		for (int j = 0; j<3; j++) {
			theta_theta(i,j) = theta_vec(i)*theta_vec(j);
		}
    }
	
	// determine local variables
    sf1 = cos(theta);
    if (theta > 5.0e-3) {
    	sf2 = (sin(theta))/theta;
    } else {
    	// small theta approximation
    	sf2 = 1.0 - theta*theta/6.0 + pow(theta,4.0)/120.0;
    }
    if (theta > 0.1) {
    	sf3 = (1.0-cos(theta))/(theta*theta);
    } else {
    	// small theta approximation
    	sf3 = 0.5 - theta*theta/24.0 + pow(theta, 4.0)/720.0
    	      -pow(theta,6.0)/40320.0 + pow(theta, 8.0)/3628800.0;
    }

	// compute exponental map vector
    Q = sf1*mEye1 + sf2*sk_theta + sf3*theta_theta;

    return Q;
}

Matrix
BeamEndContact3Dp::GetSkew(Vector th)
// this function returns the skew symmetric matrix of the given vector
{
	Matrix skew_th(BEC3p_NUM_DIM,BEC3p_NUM_DIM);

    skew_th(0,0) =  0.0;
    skew_th(0,1) = -th(2);
    skew_th(0,2) =  th(1);
    skew_th(1,0) =  th(2);
    skew_th(1,1) =  0.0;
    skew_th(1,2) = -th(0);
    skew_th(2,0) = -th(1);
    skew_th(2,1) =  th(0);
    skew_th(2,2) =  0.0;

    return skew_th;
}

const Matrix &
BeamEndContact3Dp::getTangentStiff(void)
// this function computes the tangent stiffness matrix for the element
{
	mTangentStiffness.Zero();
	double tangentTerm;

	if (inContact) {

		for (int i = 0; i < 3; i++) {
			for (int j = 0; j < 3; j++) {
				tangentTerm = mPenalty*mNormal(i)*mNormal(j);
				
				mTangentStiffness(i,j)     =  tangentTerm;
				mTangentStiffness(i,6+j)   = -tangentTerm;
				mTangentStiffness(6+i,j)   = -tangentTerm;
				mTangentStiffness(6+i,6+j) =  tangentTerm;
			}
		}
	}
		
	return mTangentStiffness;
}

const Matrix &
BeamEndContact3Dp::getInitialStiff(void)
// this function computes the initial tangent stiffness matrix for the element
{
	return getTangentStiff();
}

void
BeamEndContact3Dp::zeroLoad(void)
{
	return;
}

int
BeamEndContact3Dp::addLoad(ElementalLoad *theLoad, double loadFactor)
{
	return 0;
}

int
BeamEndContact3Dp::addInertiaLoadToUnbalance(const Vector &accel)
{
	return 0;
}

const Vector &
BeamEndContact3Dp::getResistingForce()
// this function computes the resisting force vector for the element
{
	mInternalForces.Zero();

	if (inContact) {
		
		for (int i = 0; i < 3; i++) {
			mInternalForces(i)   = -mLambda*mNormal(i);
			mInternalForces(i+6) =  mLambda*mNormal(i);
		}
	}
	
	return mInternalForces;
}

const Vector &
BeamEndContact3Dp::getResistingForceIncInertia()
{
	return getResistingForce();
}

int
BeamEndContact3Dp::sendSelf(int commitTag, Channel &theChannel)
{
	int res;

	// NOTE: we don't check for dataTag == 0 for Element
    // objects as that is taken care of in a commit by the Domain
    // object - don't want to have to do the check if sending data
	int dataTag = this->getDbTag();

	// BeamEndContact3Dp packs it's data into a Vector and sends this to theChannel
    // along with it's dbTag and the commitTag passed in the arguments
	static Vector data(5);
	data(0) = this->getTag();
	data(1) = mRadius;
	data(2) = mPenalty;
	data(3) = mBeamNode;
	data(4) = mIniContact;

	res = theChannel.sendVector(dataTag, commitTag, data);
	if (res < 0) {
		opserr << "WARNING BeamEndContact3Dp::sendSelf() - " << this->getTag() << " failed to send Vector\n";
		return -1;
	}

	// BeamEndContact3Dp then sends the tags of its nodes
	res = theChannel.sendID(dataTag, commitTag, mExternalNodes);
	if (res < 0) {
		opserr << "WARNING BeamEndContact3Dp::sendSelf() - " << this->getTag() << " failed to send Vector\n";
		return -2;
	}

	return 0;
}

int
BeamEndContact3Dp::recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker)
{
	int res;
	int dataTag = this->getDbTag();

	// BeamEndContact3Dp creates a vector, receives the vector, and then sets the internal
	// data with the data in the vector
	static Vector data(5);
	res = theChannel.recvVector(dataTag, commitTag, data);
	if (res < 0) {
		opserr << "WARNING BeamEndContact3Dp::recvSelf() - failed to receive Vector\n";
		return -1;
	}

	this->setTag((int)data(0));
	mRadius     = data(1);
	mPenalty    = data(2);
	mBeamNode   = data(3);
	mIniContact = (int)data(4);

	// BeamEndContact3Dp now receives the tags of its four external nodes
	res = theChannel.recvID(dataTag, commitTag, mExternalNodes);
	if (res < 0) {
		opserr << "WARNING BeamEndContact3Dp::recvSelf() - " << this->getTag() << " failed to receive ID\n";
		return -2;
	}

	return 0;
}

int
BeamEndContact3Dp::displaySelf(Renderer &theViewer, int displayMode, float fact, const char **modes, int numMode)
{
	return 0;
}

void
BeamEndContact3Dp::Print(OPS_Stream &s, int flag)
{
	opserr << "BeamEndContact3Dp, element id:  " << this->getTag() << endln;
	opserr << "   Connected external nodes:  ";
	for (int i = 0; i<BEC3p_NUM_NODE; i++)
	{
		opserr << mExternalNodes(i) << " ";
	}
	return;
}

Response*
BeamEndContact3Dp::setResponse(const char **argv, int argc, OPS_Stream &eleInfo)
{
	if (strcmp(argv[0],"force") == 0 || strcmp(argv[0],"forces") == 0) {
		// forces on secondary node
		return new ElementResponse(this, 1, Vector(3));

	} else if (strcmp(argv[0],"masterforce") == 0 || strcmp(argv[0],"masterforces") == 0 ||
		   strcmp(argv[0],"primaryforce") == 0 || strcmp(argv[0],"primaryforces") == 0) {
		// reactions (forces and moments) on primary node
		return new ElementResponse(this, 2, Vector(6));

	} else {
		// otherwise response quantity is unknown for the BeamEndContact3Dp class
		opserr << "BeamEndContact3Dp::setResponse(const char **argv, int argc, OPS_Stream &eleInfo): "
		  << argv[0] << " unknown request" << endln;
	return 0;
	}
}

int
BeamEndContact3Dp::getResponse(int responseID, Information &eleInfo)
{
	// initialize variables
	Vector secondaryForce(3);
	Vector primaryForce(6);

	if (responseID == 1) {

		// forces on secondary node
		secondaryForce(0) = -mInternalForces(6);
		secondaryForce(1) = -mInternalForces(7);
		secondaryForce(2) = -mInternalForces(8);
		return eleInfo.setVector(secondaryForce);
	
    } else if (responseID == 2) {

		// reactions (forces and moments) on primary node
		for (int i = 0;  i < 3; i++) {

			primaryForce(i)   = -mInternalForces(i);
			primaryForce(i+3) = -mInternalForces(i+3);
		}
		return eleInfo.setVector(primaryForce);

	} else {
		// otherwise response quantity is unknown for the BeamEndContact3Dp class
		opserr << "BeamEndContact3Dp::getResponse(int responseID = " << responseID << ", Information &eleInfo); "
		  << " unknown request" << endln;
		return -1;
	}
}

int
BeamEndContact3Dp::setParameter(const char **argv, int argc, Parameter &param)
{
	return -1;
}

int
BeamEndContact3Dp::updateParameter(int parameterID, Information &info)
{
	return -1;
}
