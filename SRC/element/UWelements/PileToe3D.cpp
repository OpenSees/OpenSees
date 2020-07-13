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

// Written: Pedro Arduino
// Created: 09/08/14
//
// Revisions
//    09/14 created
//
// Description: This file contains the implementation for the PileToe3D class.

#include "PileToe3D.h"
#include <Information.h>
#include <ElementResponse.h>
#include <CrdTransf.h>

#include <ID.h>

#include <Domain.h>
#include <Node.h>
#include <Channel.h>
#include <FEM_ObjectBroker.h>
#include <Renderer.h>
#include <G3Globals.h>
#include <ErrorHandler.h>
#include <Parameter.h>

#include <math.h>
#include <stdlib.h>
#include <stdio.h>


#include <elementAPI.h>
#define OPS_Export

static int num_PileToe3D = 0;

OPS_Export void *
OPS_PileToe3D(void)
{
  if (num_PileToe3D == 0) {
    num_PileToe3D++;
    //OPS_Error("PileToe3D element - Written: P.Arduino, P.Mackenzie-Helnwein, U.Washington\n", 1);
    opserr <<"PileToe3D element - Written: P.Arduino, P.Mackenzie-Helnwein, U.Washington\n";
  }

  // Pointer to a uniaxial material that will be returned
  Element *theElement = 0;

  int numRemainingInputArgs = OPS_GetNumRemainingInputArgs();

  if (numRemainingInputArgs < 7) {
    opserr << "Invalid #args,  want: element PileToe3D eleTag?  iNode? BiNode? BjNode? radius? k? crdTransf?\n";
    return 0;
  }
   
  int    iData[5];
  double dData[2];

  int numData = 4;
  if (OPS_GetIntInput(&numData, iData) != 0) {
    opserr << "WARNING invalid integer data: element PileToe3D" << endln;
    return 0;
  }

  numData = 1;
  if (OPS_GetDoubleInput(&numData, dData) != 0) {
    opserr << "WARNING invalid radius data: element PileToe3D " << iData[0] << endln;
    return 0;  
  }

  numData = 1;
  if (OPS_GetDoubleInput(&numData, &dData[1]) != 0) {
    opserr << "WARNING invalid  k data: element PileToe3D " << iData[0] << endln;
    return 0;  
  }

  numData = 1;
  if (OPS_GetIntInput(&numData, &iData[4]) != 0) {
    opserr << "WARNING invalid integer crdTransf data: element PileToe3D" << iData[0] << endln;
    return 0;
  }

  int transfTag = iData[4];
  CrdTransf *theTransf = OPS_getCrdTransf(transfTag);
  if (theTransf == 0) {
    opserr << "WARNING element PileToe3D " << iData[0] << endln;
    opserr << " coordTransf: " << transfTag << "not found\n";
    return 0;
  }

  // Parsing was successful, allocate the element
  theElement = new PileToe3D(iData[0], iData[1], iData[2], iData[3], dData[0], dData[1], *theTransf);

  if (theElement == 0) {
    opserr << "WARNING could not create element of type PileToe3D\n";
    return 0;
  }

  return theElement;
}


// constructors:
PileToe3D::PileToe3D(int tag, int Nd1, int BNd1, int BNd2, double rad, double k, CrdTransf &coordTransf)
  :Element(tag,ELE_TAG_PileToe3D),    
   crdTransf(0),
   externalNodes(PT3D_NUM_NODE),
   externalBNodes(2),
   mTangentStiffness(PT3D_NUM_DOF, PT3D_NUM_DOF),
   mInternalForces(PT3D_NUM_DOF)
{
        externalNodes(0) = Nd1;
		externalBNodes(0) = BNd1;
		externalBNodes(1) = BNd2;

        mRadius = rad;
        mSubgradeCoeff = k;
		mCC = mRadius;  // Initial value of mCC --> Total area 
       
        // get copy of the transformation & material object  
        crdTransf = coordTransf.getCopy3d();
       
        // check it:         
        if (!crdTransf) {
		opserr << "Error: PileToe3D:PileToe3D: could not create copy of coordinate transformation object" << endln;
            exit(-1);
        }
              
        // element tag for debugging
        MyTag = tag;

        // set initialization to true for setDomain function
	    mInitialize = true;
}

PileToe3D::PileToe3D()
 :Element(0,ELE_TAG_PileToe3D),    
   crdTransf(0),
   externalNodes(PT3D_NUM_NODE),
   externalBNodes(2),
   mTangentStiffness(PT3D_NUM_DOF, PT3D_NUM_DOF),
   mInternalForces(PT3D_NUM_DOF)
{
    mRadius         = 0.0;
    mSubgradeCoeff  = 0.0;
    mCC             = 0.0;
    mInitialize     = false;
}


//  destructor:
PileToe3D::~PileToe3D()
{
}


int
PileToe3D::getNumExternalNodes(void) const
{
    return PT3D_NUM_NODE;
}

int
PileToe3D::getNumExternalBNodes(void) const
{
    return 2;
}

const ID &
PileToe3D::getExternalNodes(void)
{
    return externalNodes;
}

const ID &
PileToe3D::getExternalBNodes(void)
{
    return externalBNodes;
}

Node **
PileToe3D::getNodePtrs(void)
{
        return theNodes;                        
}

Node **
PileToe3D::getBNodePtrs(void)
{
        return theBNodes;                        
}

int
PileToe3D::getNumDOF(void)
{
    return PT3D_NUM_DOF;
}

void
PileToe3D::setDomain(Domain *theDomain)
{
  theNodes[0]  = theDomain->getNode(externalNodes(0));

  theBNodes[0] = theDomain->getNode(externalBNodes(0));
  theBNodes[1] = theDomain->getNode(externalBNodes(1));
               
  for (int i = 0; i < 1; i++) {
    if (theNodes[i] == 0) {
      opserr << "PileToe3D::setDomain() - no node with tag: " << theNodes[i] << endln;
      return;  // don't go any further - otherwise segmentation fault
    } 
  }

  for (int i = 0; i < 2; i++) {
    if (theBNodes[i] == 0) {
      opserr << "PileToe3D::setDomain() - no beam node with tag: " << theNodes[i] << endln;
      return;  // don't go any further - otherwise segmentation fault
    } 
  }

  // only perform these steps during initial creation of element
  if (mInitialize) {

	// initialize coordinate vectors
	//const Vector &mIcrd_1 = theNodes[0]->getCrds();
	// coordinate matrix
	//mNodeCrd(0,0) = mIcrd_1(0);  
  }

  //Initialize Coordinate Transformation
  if (crdTransf->initialize(theBNodes[0], theBNodes[1])) {
	  // Add some error check
  }

  double L = crdTransf->getInitialLength();

  if (L == 0.0) {
	// Add some error check
  }

  // call the base class method
  this->DomainComponent::setDomain(theDomain);
}        


int
PileToe3D::commitState()
{
     int retVal = 0;
     // call element commitState to do any base class stuff
     if ((retVal = this->Element::commitState()) != 0) {
        opserr << "PileToe3D::commitState () - failed in base class";
      }    
      //retVal = theMaterial->commitState();
      retVal += crdTransf->commitState();
      
	  return retVal;
}


int
PileToe3D::revertToLastCommit()
{
	 int retVal = 0;
	 retVal += crdTransf->revertToLastCommit();
	 return retVal;
}

int
PileToe3D::revertToStart()
{
    int retVal = 0;
    retVal += crdTransf->revertToStart();
    return retVal;
}

int
PileToe3D::update(void)
{
  mInitialize = true;
  // Update the transformation
  crdTransf->update();
  return 0;
}

const Matrix &
PileToe3D::getTangentStiff(void)
{
	//double mPi   = 4.0*atan(1);
	const double mPi = 3.1415926535897;
    double mArea     = mPi* mRadius*mRadius;
	double mII       = mPi*mRadius*mRadius*mRadius*mRadius/4.0;
    // initialize Kt
    mTangentStiffness.Zero();

	/*
	mTangentStiffness(0,0) = mSubgradeCoeff*mArea;
	mTangentStiffness(1,1) = mSubgradeCoeff*mArea;
	mTangentStiffness(2,2) = mSubgradeCoeff*mArea;
	mTangentStiffness(3,3) = mSubgradeCoeff*mII;	
	mTangentStiffness(4,4) = mSubgradeCoeff*mII;
	mTangentStiffness(5,5) = mSubgradeCoeff*mII;	
	
	mTangentStiffness(0,0) = mSubgradeCoeff*mArea;
	mTangentStiffness(4,4) = mSubgradeCoeff*mII;
	mTangentStiffness(5,5) = mSubgradeCoeff*mII;	
	*/
    mTangentStiffness(2,2) = mSubgradeCoeff*mArea;
	mTangentStiffness(3,3) = mSubgradeCoeff*mII;
	mTangentStiffness(4,4) = mSubgradeCoeff*mII;	

    return mTangentStiffness;
}

const Matrix &
PileToe3D::getInitialStiff(void)
{
  return getTangentStiff();
}
   
void
PileToe3D::zeroLoad(void)
{
}

int
PileToe3D::addLoad(ElementalLoad *theLoad, double loadFactor)
{
  return 0;
}

int
PileToe3D::addInertiaLoadToUnbalance(const Vector &accel)
{
  return 0;
}

const Vector &
PileToe3D::getResistingForce()
{    
      // initialize F
      mInternalForces.Zero();

      // get trial displacement
      const Vector &mDisp_1 = theNodes[0]->getTrialDisp();
      mInternalForces = mTangentStiffness*mDisp_1;

      return mInternalForces;
}


const Vector &
PileToe3D::getResistingForceIncInertia()
{      
      return getResistingForce();				
}


int
PileToe3D::sendSelf(int commitTag, Channel &theChannel)
{

  int res = 0;

  // note: we don't check for dataTag == 0 for Element
  // objects as that is taken care of in a commit by the Domain
  // object - don't want to have to do the check if sending data
  int dataTag = this->getDbTag();

  // PileToe3D packs it's data into a Vector and sends this to theChannel
  // along with it's dbTag and the commitTag passed in the arguments

  static Vector data(6);
  data(0) = this->getTag();
  data(1) = mRadius;
  data(2) = mSubgradeCoeff;
  data(3) = mCC;
  data(4) = crdTransf->getClassTag();
  int crdDbTag = crdTransf->getDbTag();
  // NOTE: we do have to ensure that the material has a database
  // tag if we are sending to a database channel.
  if (crdDbTag == 0) {
    crdDbTag = theChannel.getDbTag();
    if (crdDbTag != 0)
      crdTransf->setDbTag(crdDbTag);
  }
  data(5) = crdDbTag;

  res += theChannel.sendVector(dataTag, commitTag, data);
  if (res < 0) {
    opserr <<"WARNING PileToe3D::sendSelf() - " << this->getTag() << " failed to send Vector\n";
    return -1;
  }          

  //PileToe3D then sends the tags of its one node
  res = theChannel.sendID(dataTag, commitTag, externalNodes);
  if (res < 0) {
    opserr <<"WARNING PileToe3D::sendSelf() - " << this->getTag() << " failed to send Vector\n";
    return -2;
  }

  //PileToe3D then sends the tags of its two beam nodes
  res = theChannel.sendID(dataTag, commitTag, externalBNodes);
  if (res < 0) {
    opserr <<"WARNING PileToe3D::sendSelf() - " << this->getTag() << " failed to send Vector\n";
    return -2;
  }


  //PileToe3D asks its crdTransf object to send itself
  res = crdTransf->sendSelf(commitTag, theChannel);
  if (res < 0) {
    opserr <<"WARNING PileToe3D::sendSelf() - " << this->getTag() << " failed to send its crdTransf\n";
    return -4;
  }

  return 0;
}

int
PileToe3D::recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker)
{
  int res;
  int dataTag = this->getDbTag();

  // PileToe3D creates a Vector, receives the Vector and then sets the
  // internal data with the data in the Vector
  static Vector data(6);
  res = theChannel.recvVector(dataTag, commitTag, data);
  if (res < 0) {
    opserr <<"WARNING PileToe3D::recvSelf() - failed to receive Vector\n";
    return -1;
  }          

  this->setTag((int)data(0));
  mRadius            = data(1);
  mSubgradeCoeff     = data(2);
  mCC                = data(3);

  MyTag              = (int)data(0);
 
  // PileToe3D now receives the tags of it's one external node
  res = theChannel.recvID(dataTag, commitTag, externalNodes);
  if (res < 0) {
    opserr <<"WARNING PileToe3D::recvSelf() - " << this->getTag() << " failed to receive ID\n";
    return -2;
  }

  // PileToe3D now receives the tags of it's two beam external nodes
  res = theChannel.recvID(dataTag, commitTag, externalBNodes);
  if (res < 0) {
    opserr <<"WARNING PileToe3D::recvSelf() - " << this->getTag() << " failed to receive ID\n";
    return -2;
  }

  int crdClass = (int)data(4);
  int crdDb = (int)data(5);

  // check if we have a material object already & if we do if of right type
  if ((crdTransf == 0) || (crdTransf->getClassTag() != crdClass)) {
    
    // if old one .. delete it
    if (crdTransf != 0)
      delete crdTransf;
    
    // create a new material object
    crdTransf = theBroker.getNewCrdTransf(crdClass);
    
    if (crdTransf == 0) {
      opserr <<"WARNING PileToe3D::recvSelf() - " << this->getTag()
	     << " failed to get a blank CrdTransf of type " << crdClass << endln;
      return -3;
    }
  }

  crdTransf->setDbTag(crdDb); // note: we set the dbTag before we receive the material
  
  res = crdTransf->recvSelf(commitTag, theChannel, theBroker);
  if (res < 0) {
    opserr <<"WARNING PileToe3D::recvSelf() - "<< this->getTag() << "failed to receive its Material\n";
    return -3;    
  }

  return 0;
}


int
PileToe3D::displaySelf(Renderer &theViewer, int displayMode, float fact, const char **modes, int numMode)
{
  return 0;
}


void
PileToe3D::Print(OPS_Stream &s, int flag)
{
        opserr << "PileToe3D, element id:  " << this->getTag() << endln;
        opserr << "   Connected external nodes:  " ;
        for (int i = 0; i<PT3D_NUM_NODE; i++)
        {
                opserr << externalNodes(i)<< " ";
        }
        return;
}


Response*
PileToe3D::setResponse(const char **argv, int argc, OPS_Stream &eleInfo)
{

    if (strcmp(argv[0],"reaction") == 0 || strcmp(argv[0],"reactions") == 0)
      return new ElementResponse(this, 1, Vector(6));
    // otherwise response quantity is unknown for the BeamContact3D class
    else
      opserr << "BeamContact3D::setResponse(const char **argv, int argc, OPS_Stream &eleInfo): " << argv[0] << " unknown request" << endln;
      return 0;
}


int
PileToe3D::getResponse(int responseID, Information &eleInfo)
{
  Vector mReactions(6);
  if (responseID == 1) {

     // full reactions on primary nodes
     for (int ii=0; ii<6; ii++) {
          mReactions(ii)   = -mInternalForces(ii);
     }
     return eleInfo.setVector(mReactions);

  } else

    opserr << "PileToe3D::getResponse(int responseID=" << responseID << ", Information &eleInfo): " << " unknown request" << endln;
    return -1;
}
