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
                                                                       
// Created: Chris McGann, UW, 04.2011
//
// Description: This file contains the implementation of the SSPquad class
//
// Reference:   McGann, C.R., Arduino, P., and Mackenzie-Helnwein, P. (2012) "Stabilized single-point
//                4-node quadrilateral element for dynamic analysis of fluid saturated porous media."
//                Acta Geotechnica, 7(4):297-311

#include "SSPquad.h"

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

static int num_SSPquad = 0;

OPS_Export void *
OPS_SSPquad(void)
{
	if (num_SSPquad == 0) {
    	num_SSPquad++;
    	opserr << "SSPquad element - Written: C.McGann, P.Arduino, P.Mackenzie-Helnwein, U.Washington\n";
  	}

  	// Pointer to an element that will be returned
  	Element *theElement = 0;

  	int numRemainingInputArgs = OPS_GetNumRemainingInputArgs();

  	if (numRemainingInputArgs < 8) {
    	opserr << "Invalid #args, want: element SSPquad eleTag? iNode? jNode? kNode? lNode? matTag? type? thickness? <b1? b2?>?\n";
		return 0;
  	}

  	int iData[6];
  	const char *theType;
    double dData[3] = { 1.0,0.0,0.0 };

  	int numData = 6;
  	if (OPS_GetIntInput(&numData, iData) != 0) {
    	opserr << "WARNING invalid integer data: element SSPquad " << iData[0] << endln;
		return 0;
  	}

	theType = OPS_GetString();

	numData = 1;
	if (OPS_GetDoubleInput(&numData, dData) != 0) {
		opserr << "WARNING invalid thickness data: element SSPquad " << iData[0] << endln;
		return 0;
	}

  	int matID = iData[5];
  	NDMaterial *theMaterial = OPS_getNDMaterial(matID);
  	if (theMaterial == 0) {
    	opserr << "WARNING element SSPquad " << iData[0] << endln;
		opserr << " Material: " << matID << "not found\n";
		return 0;
  	}

	if (numRemainingInputArgs == 10) {
    	numData = 2;
    	if (OPS_GetDoubleInput(&numData, &dData[1]) != 0) {
      		opserr << "WARNING invalid optional data: element SSPquad " << iData[0] << endln;
	  		return 0;
    	}
  	}

  	// parsing was successful, allocate the element
  	theElement = new SSPquad(iData[0], iData[1], iData[2], iData[3], iData[4],
                                 *theMaterial, theType, dData[0], dData[1], dData[2]);

  	if (theElement == 0) {
    	opserr << "WARNING could not create element of type SSPquad\n";
		return 0;
  	}

  	return theElement;
}

// full constructor
SSPquad::SSPquad(int tag, int Nd1, int Nd2, int Nd3, int Nd4, NDMaterial &theMat, 
                          const char *type, double thick, double b1, double b2)
  :Element(tag,ELE_TAG_SSPquad),
  	theMaterial(0),
	mExternalNodes(SSPQ_NUM_NODE),
	mTangentStiffness(SSPQ_NUM_DOF,SSPQ_NUM_DOF),
	mInternalForces(SSPQ_NUM_DOF),
	Q(SSPQ_NUM_DOF),
	mMass(SSPQ_NUM_DOF,SSPQ_NUM_DOF),
	mNodeCrd(2,4),
	Mmem(3,SSPQ_NUM_DOF),
	Kstab(SSPQ_NUM_DOF,SSPQ_NUM_DOF),
	mThickness(thick),
	applyLoad(0)
{
	mExternalNodes(0) = Nd1;
	mExternalNodes(1) = Nd2;
	mExternalNodes(2) = Nd3;
	mExternalNodes(3) = Nd4;
		
	mThickness = thick;

	b[0] = b1;
	b[1] = b2;
	
	appliedB[0] = 0.0;
	appliedB[1] = 0.0;

	// get copy of the material object
	NDMaterial *theMatCopy = theMat.getCopy(type);
	if (theMatCopy != 0) {
		theMaterial = (NDMaterial *)theMatCopy;
	} else {
		opserr << "SSPquad::SSPquad - failed to get copy of material model\n";;
	}

	// check material
	if (theMaterial == 0) {
		opserr << "SSPquad::SSPquad - failed to allocate material model pointer\n";
		exit(-1);
	}
		
	// check the type
	if (strcmp(type,"PlaneStrain") != 0 && strcmp(type,"PlaneStress") != 0) {
		opserr << "SSPquad::SSPquad - improper material type: " << type << "for SSPquad\n";
		exit(-1);
	}
}

// null constructor
SSPquad::SSPquad()
  :Element(0,ELE_TAG_SSPquad),
  	theMaterial(0),
	mExternalNodes(SSPQ_NUM_NODE),
	mTangentStiffness(SSPQ_NUM_DOF,SSPQ_NUM_DOF),
	mInternalForces(SSPQ_NUM_DOF),
	Q(SSPQ_NUM_DOF),
	mMass(SSPQ_NUM_DOF,SSPQ_NUM_DOF),
	mNodeCrd(2,4),
	Mmem(3,SSPQ_NUM_DOF),
	Kstab(SSPQ_NUM_DOF,SSPQ_NUM_DOF),
	mThickness(0),
	applyLoad(0)
{
}

// destructor
SSPquad::~SSPquad()
{
	if (theMaterial != 0) {
        delete theMaterial;
    }
}

int 
SSPquad::getNumExternalNodes(void) const
{
    return SSPQ_NUM_NODE;
}

const ID &
SSPquad::getExternalNodes(void)
{
    return mExternalNodes;
}

Node **
SSPquad::getNodePtrs(void)
{
    return theNodes;
}

int
SSPquad::getNumDOF(void)
{
    return SSPQ_NUM_DOF;
}

void
SSPquad::setDomain(Domain *theDomain)
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

	// call the base-class method
	this->DomainComponent::setDomain(theDomain);
}

int
SSPquad::commitState(void)
{
	int retVal = 0;
	// call element commitState to do any base class stuff
	if ((retVal = this->Element::commitState()) != 0) {
		opserr << "SSPquad::commitState() - failed in base class\n";
	}
	retVal = theMaterial->commitState();

	return retVal;
}

int
SSPquad::revertToLastCommit(void)
{
	return theMaterial->revertToLastCommit();
}

int
SSPquad::revertToStart(void)
{
	return theMaterial->revertToStart();
}

int
SSPquad::update(void)
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
SSPquad::getTangentStiff(void)
// this function computes the tangent stiffness matrix for the element
{
	// get material tangent
	const Matrix &Cmat = theMaterial->getTangent();

	// full element stiffness matrix
	mTangentStiffness = Kstab;
	mTangentStiffness.addMatrixTripleProduct(1.0, Mmem, Cmat, 4.0*J0*mThickness);
	
	return mTangentStiffness;
}

const Matrix &
SSPquad::getInitialStiff(void)
// this function computes the initial tangent stiffness matrix for the element
{
	return getTangentStiff();
}

const Matrix &
SSPquad::getMass(void)
{
	mMass.Zero();

	// get mass density from the material
	double density = theMaterial->getRho();

	// return zero matrix if density is zero
	if (density == 0.0) {
		return mMass;
	}
	
	// local coordinates of nodes
	double xi[4];
	double eta[4];
	xi[0]  = -1.0; xi[1]  =  1.0; xi[2]  = 1.0; xi[3]  = -1.0;
	eta[0] = -1.0; eta[1] = -1.0; eta[2] = 1.0; eta[3] =  1.0;

	double massTerm;
	for (int i = 0; i < 4; i++) {
		massTerm = density*mThickness*(J0 + J1*xi[i] + J2*eta[i]);
		mMass(2*i,2*i)     += massTerm;
		mMass(2*i+1,2*i+1) += massTerm;
	}

	return mMass;
}

void
SSPquad::zeroLoad(void)
{
  applyLoad = 0;
  appliedB[0] = 0.0;
  appliedB[1] = 0.0;
  
  Q.Zero();

  return;
}

int
SSPquad::addLoad(ElementalLoad *theLoad, double loadFactor)
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
		opserr << "SSPquad::addLoad - load type unknown for ele with tag: " << this->getTag() << endln;
		return -1;
	} 

	return -1;
}

int
SSPquad::addInertiaLoadToUnbalance(const Vector &accel)
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

	if (2 != Raccel1.Size() || 2 != Raccel2.Size() || 2 != Raccel3.Size() || 2 != Raccel4.Size()) {
    	opserr << "FourNodeQuad::addInertiaLoadToUnbalance matrix and vector sizes are incompatible\n";
    	return -1;
	}

	static double ra[8];
	ra[0] = Raccel1(0);
	ra[1] = Raccel1(1);
	ra[2] = Raccel2(0);
	ra[3] = Raccel2(1);
	ra[4] = Raccel3(0);
	ra[5] = Raccel3(1);
	ra[6] = Raccel4(0);
	ra[7] = Raccel4(1);

	// compute mass matrix
	this->getMass();

	for (int i = 0; i < 8; i++) {
		Q(i) += -mMass(i,i)*ra[i];
	}
	
	return 0;
}

const Vector &
SSPquad::getResistingForce(void)
// this function computes the resisting force vector for the element
{
	// get stress from the material
	Vector mStress = theMaterial->getStress();

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
	mInternalForces = Kstab*d;

	// add internal force from the stress  ->  fint = Kstab*d + 4*t*Jo*Mmem'*stress
	mInternalForces.addMatrixTransposeVector(1.0, Mmem, mStress, 4.0*mThickness*J0);

	// subtract body forces from internal force vector
	double xi[4];
	double eta[4];
	xi[0]  = -1.0; xi[1]  =  1.0; xi[2]  = 1.0; xi[3]  = -1.0;
	eta[0] = -1.0; eta[1] = -1.0; eta[2] = 1.0; eta[3] =  1.0;

	if (applyLoad == 0) {
		for (int i = 0; i < 4; i++) {
			mInternalForces(2*i)   -= b[0]*mThickness*(J0 + J1*xi[i] + J2*eta[i]);
			mInternalForces(2*i+1) -= b[1]*mThickness*(J0 + J1*xi[i] + J2*eta[i]);
		}
	} else {
		for (int i = 0; i < 4; i++) {
			mInternalForces(2*i)   -= appliedB[0]*mThickness*(J0 + J1*xi[i] + J2*eta[i]);
			mInternalForces(2*i+1) -= appliedB[1]*mThickness*(J0 + J1*xi[i] + J2*eta[i]);
		}
	}

	// inertial unbalance load
	mInternalForces.addVector(1.0, Q, -1.0);

	return mInternalForces;
}

const Vector &
SSPquad::getResistingForceIncInertia()
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
	
	static double a[8];
	a[0] = accel1(0);
	a[1] = accel1(1);
	a[2] = accel2(0);
	a[3] = accel2(1);
	a[4] = accel3(0);
	a[5] = accel3(1);
	a[6] = accel4(0);
	a[7] = accel4(1);

	// compute current resisting force
	this->getResistingForce();

	// compute mass matrix
	this->getMass();

	for (int i = 0; i < 8; i++) {
		mInternalForces(i) += mMass(i,i)*a[i];
	}

	// add the damping forces if rayleigh damping
	if (alphaM != 0 || betaK != 0.0 || betaK0 != 0.0 || betaKc != 0.0) {
		mInternalForces += this->getRayleighDampingForces();
	}

	return mInternalForces;
}

int
SSPquad::sendSelf(int commitTag, Channel &theChannel)
{
  int res = 0;
  
  // note: we don't check for dataTag == 0 for Element
  // objects as that is taken care of in a commit by the Domain
  // object - don't want to have to do the check if sending data
  int dataTag = this->getDbTag();
  
  // SSPquad packs its data into a Vector and sends this to theChannel
  // along with its dbTag and the commitTag passed in the arguments
  static Vector data(10);
  data(0) = this->getTag();
  data(1) = mThickness;
  data(2) = b[0];
  data(3) = b[1];
  data(4) = theMaterial->getClassTag();	      
  
  data(6) = alphaM;
  data(7) = betaK;
  data(8) = betaK0;
  data(9) = betaKc;
  // Now quad sends the ids of its materials
  int matDbTag = theMaterial->getDbTag();
  
  // NOTE: we do have to ensure that the material has a database
  // tag if we are sending to a database channel.
  if (matDbTag == 0) {
    matDbTag = theChannel.getDbTag();
    if (matDbTag != 0)
      theMaterial->setDbTag(matDbTag);
  }
  data(5) = matDbTag;
  
  res += theChannel.sendVector(dataTag, commitTag, data);
  if (res < 0) {
    opserr << "WARNING SSPquad::sendSelf() - " << this->getTag() << " failed to send Vector\n";
    return res;
  }
  
  // SSPquad then sends the tags of its four nodes
  res += theChannel.sendID(dataTag, commitTag, mExternalNodes);
  if (res < 0) {
    opserr << "WARNING SSPquad::sendSelf() - " << this->getTag() << " failed to send ID\n";
    return res;
  }
  
  // finally, SSPquad asks its material object to send itself
  res = theMaterial->sendSelf(commitTag, theChannel);
  if (res < 0) {
    opserr << "WARNING SSPquad::sendSelf() - " << this->getTag() << " failed to send its Material\n";
    return -3;
  }
  
  return 0;
}

int
SSPquad::recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker)
{
  int res = 0;
  int dataTag = this->getDbTag();
  
  // SSPquad creates a Vector, receives the Vector and then sets the 
  // internal data with the data in the Vector
  static Vector data(10);
  res += theChannel.recvVector(dataTag, commitTag, data);
  if (res < 0) {
    opserr << "WARNING SSPquad::recvSelf() - failed to receive Vector\n";
    return res;
  }
  
  this->setTag((int)data(0));
  mThickness = data(1);
  b[0] = data(2);
  b[1] = data(3);
  
  // SSPquad now receives the tags of its four external nodes
  res += theChannel.recvID(dataTag, commitTag, mExternalNodes);
  if (res < 0) {
    opserr << "WARNING SSPquad::recvSelf() - " << this->getTag() << " failed to receive ID\n";
    return res;
  }
  
  // finally, SSPquad creates a material object of the correct type, sets its
  // database tag, and asks this new object to receive itself
  int matClass = (int)data(4);
  int matDb    = (int)data(5);
  
  alphaM = data(6);
  betaK = data(7);
  betaK0 = data(8);
  betaKc = data(9);
  
  // check if material object exists and that it is the right type
  if ((theMaterial == 0) || (theMaterial->getClassTag() != matClass)) {
    
    // if old one, delete it
    if (theMaterial != 0)
      delete theMaterial;
    
    // create new material object
    NDMaterial *theMatCopy = theBroker.getNewNDMaterial(matClass);
    theMaterial = (NDMaterial *)theMatCopy;
    
    if (theMaterial == 0) {
      opserr << "WARNING SSPquad::recvSelf() - " << this->getTag() 
	     << " failed to get a blank Material of type " << matClass << endln;
      return -3;
    }
  }
  
  // NOTE: we set the dbTag before we receive the material
  theMaterial->setDbTag(matDb);
  res = theMaterial->recvSelf(commitTag, theChannel, theBroker);
  if (res < 0) {
    opserr << "WARNING SSPquad::recvSelf() - " << this->getTag() << " failed to receive its Material\n";
    return -3;
  }
  
  return 0; 
}

int
SSPquad::displaySelf(Renderer &theViewer, int displayMode, float fact, const char **modes, int numMode)
{
	return 0;
}

void
SSPquad::Print(OPS_Stream &s, int flag)
{
    if (flag == OPS_PRINT_CURRENTSTATE) {
        opserr << "SSPquad, element id:  " << this->getTag() << endln;
        opserr << "   Connected external nodes:  ";
        for (int i = 0; i < SSPQ_NUM_NODE; i++) {
            opserr << mExternalNodes(i) << " ";
        }
    }
    
    if (flag == OPS_PRINT_PRINTMODEL_JSON) {
        s << "\t\t\t{";
        s << "\"name\": " << this->getTag() << ", ";
        s << "\"type\": \"SSPquad\", ";
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
SSPquad::setResponse(const char **argv, int argc, OPS_Stream &eleInfo)
{
	// no special recorders for this element, call the method in the material class
	return theMaterial->setResponse(argv, argc, eleInfo);
}

int
SSPquad::getResponse(int responseID, Information &eleInfo)
{
	// no special recorders for this element, call the method in the material class
	return theMaterial->getResponse(responseID, eleInfo);
}

int
SSPquad::setParameter(const char **argv, int argc, Parameter &param)
{
	if (argc < 1) {
		return -1;
	}

	int res = -1;

    // no element parameters, call setParameter in the material	
    int matRes = res;
	matRes =  theMaterial->setParameter(argv, argc, param);

    if (matRes != -1) {
		res = matRes;
	}
  
  return res;
}

int
SSPquad::updateParameter(int parameterID, Information &info)
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

Matrix 
SSPquad::DyadicProd(Vector v1, Vector v2)
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
SSPquad::GetStab(void)
// this function computes the stabilization stiffness matrix for the element
{
	Vector g1(SSPQ_NUM_DIM);
	Vector g2(SSPQ_NUM_DIM);
	Matrix I(SSPQ_NUM_DIM,SSPQ_NUM_DIM);
	Matrix FCF(SSPQ_NUM_DIM,SSPQ_NUM_DIM);
	Matrix Jmat(SSPQ_NUM_DIM,SSPQ_NUM_DIM);
	Matrix Jinv(SSPQ_NUM_DIM,SSPQ_NUM_DIM);
	Matrix dNloc(SSPQ_NUM_NODE,SSPQ_NUM_DIM);
	Matrix dN(SSPQ_NUM_NODE,SSPQ_NUM_DIM);
	Matrix Mben(2,SSPQ_NUM_DOF);
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
	Hss = (I(0,0)*Jinv(1,0)*Jinv(1,0) + I(0,1)*Jinv(0,0)*Jinv(1,0) + I(1,1)*Jinv(0,0)*Jinv(0,0))*0.25;
	Htt = (I(0,0)*Jinv(1,1)*Jinv(1,1) + I(0,1)*Jinv(0,1)*Jinv(1,1) + I(1,1)*Jinv(0,1)*Jinv(0,1))*0.25;
	Hst = (I(0,0)*Jinv(1,1)*Jinv(1,0) + I(0,1)*(Jinv(1,0)*Jinv(0,1) + Jinv(1,1)*Jinv(0,0)) + I(1,1)*Jinv(0,1)*Jinv(0,0))*0.25;

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

