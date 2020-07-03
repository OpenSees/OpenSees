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
                                                                        
                                                                        
// Written: J. Abell
// Created: March 2018
// Modified: 

// Description: This file contains the implementation for the LysmerTriangle class.

#include "LysmerTriangle.h"
#include <Information.h>
#include <Parameter.h>
#include <ElementResponse.h>
#include <ID.h> 
#include <Domain.h>
#include <Node.h>
#include <Channel.h>
#include <FEM_ObjectBroker.h>
#include <Renderer.h>
#include <G3Globals.h>
#include <ErrorHandler.h>
#include <NDMaterial.h>
#include <ElementalLoad.h>

#include <math.h>
#include <stdlib.h>
#include <stdio.h> 

double LysmerTriangle :: oneOverRoot3 = 1.0/sqrt(3.0);
double LysmerTriangle :: GsPts[1][1];
Matrix LysmerTriangle :: Bmat(9,3);
Matrix LysmerTriangle :: tangentStiffness(SL_NUM_DOF, SL_NUM_DOF);
Matrix LysmerTriangle :: tangentDamping(SL_NUM_DOF, SL_NUM_DOF);

#include <elementAPI.h>
static int num_LysmerTriangle = 0;

void *
OPS_LysmerTriangle(void)
{
  if (num_LysmerTriangle == 0) {
    num_LysmerTriangle++;
    opserr<<"LysmerTriangle element - Written: J. A. Abell (UANDES). www.joseabell.com\n";
  }

  // Pointer to a uniaxial material that will be returned
  Element *theElement = 0;

  if (OPS_GetNumRemainingInputArgs() < 7) {
    opserr << "Want: element LysmerTriangle eleTag?  iNode? jNode? kNode? rho Vp Vs? <length> <stage> \n";
    return 0;
  }
    
  int    iData[4];
  double dData[3];
  double eleLength = 0;
  int stage = 0;

  int numData = 4;
  if (OPS_GetIntInput(&numData, iData) != 0) {
    opserr << "WARNING invalid integer data: element LysmerTriangleElement" << endln;
    return 0;
  }

  numData = 3;
  if (OPS_GetDoubleInput(&numData, dData) != 0) {
    opserr << "WARNING invalid data: element LysmerTriangle " << iData[0] << endln;
    return 0; 
  }

  int num_args_remaining = OPS_GetNumRemainingInputArgs();

  // Its optional (but desireable) to input the element-length....
  if (num_args_remaining > 0)
  {
    numData = 1;
    OPS_GetDoubleInput(&numData, &eleLength);
  }

// Its optional to set the element stage... will be set to 0 (damping mode) by default
  if (num_args_remaining > 0)
  {
    numData = 1;
    OPS_GetIntInput(&numData, &stage);
  }

  // Parsing was successful, allocate the material
  theElement = new LysmerTriangle(iData[0], iData[1], iData[2], iData[3], dData[0], dData[1], dData[2], eleLength, stage );

  if (theElement == 0) {
    opserr << "WARNING could not create element of type LysmerTriangleElement\n";
    return 0;
  }

  return theElement;
}

// constructors:
LysmerTriangle::LysmerTriangle(int tag, int Nd1, int Nd2, int Nd3, double rho_, double Vp_, double Vs_, double eleLength, int stage)
 :Element(tag,ELE_TAG_LysmerTriangle),     
   myExternalNodes(SL_NUM_NODE),
   internalForces(SL_NUM_DOF),
   springForces(SL_NUM_DOF),
   rho(rho_), Vp(Vp_), Vs(Vs_), element_length(eleLength),
   g1(SL_NUM_NDF), 
   g2(SL_NUM_NDF),
   myNhat(SL_NUM_NDF), 
   myThat(SL_NUM_NDF), 
   myShat(SL_NUM_NDF), 
   myNI(SL_NUM_NODE),
   dcrd1(SL_NUM_NDF),
   dcrd2(SL_NUM_NDF),
   dcrd3(SL_NUM_NDF),
   gnd_velocity(3),
   stage(stage)
{
    myExternalNodes(0) = Nd1;
    myExternalNodes(1) = Nd2;
    myExternalNodes(2) = Nd3;

	MyTag = tag;

	GsPts[0][0] = 0.5;

	mLoadFactor = 1.0;

  // if(stage != 0)
  // {
  //   opserr << "LysmerTriangle::LysmerTriangle - element at tag # " << this->getTag() << " starting at stage = " << stage << " also L = " << element_length <<  endln;
  // }
}

LysmerTriangle::LysmerTriangle()
  :Element(0,ELE_TAG_LysmerTriangle),     
   	myExternalNodes(SL_NUM_NODE),
   	internalForces(SL_NUM_DOF),
    springForces(SL_NUM_DOF),
    rho(0), Vp(0), Vs(0), element_length(0),
    g1(SL_NUM_NDF), 
    g2(SL_NUM_NDF),
    myNhat(SL_NUM_NDF), 
    myThat(SL_NUM_NDF), 
    myShat(SL_NUM_NDF), 
    myNI(SL_NUM_NODE),
   	dcrd1(SL_NUM_NDF),
   	dcrd2(SL_NUM_NDF),
   	dcrd3(SL_NUM_NDF),
    gnd_velocity(3),
    stage(0)
{
}

//  destructor:
LysmerTriangle::~LysmerTriangle()
{
}

int
LysmerTriangle::getNumExternalNodes(void) const
{
    return SL_NUM_NODE;
}

const ID &
LysmerTriangle::getExternalNodes(void) 
{
    return myExternalNodes;
}

Node **
LysmerTriangle::getNodePtrs(void)
{
    return theNodes;                        
}

int
LysmerTriangle::getNumDOF(void) 
{
    return SL_NUM_DOF;
}

void
LysmerTriangle::setDomain(Domain *theDomain)
{
    theNodes[0] = theDomain->getNode(myExternalNodes(0));
    theNodes[1] = theDomain->getNode(myExternalNodes(1));
    theNodes[2] = theDomain->getNode(myExternalNodes(2));

    for (int i = 0; i < 3; i++) {
    	if (theNodes[i] == 0)
        	return;  // don't go any further - otherwise segmentation fault
    }

    dcrd1 = theNodes[0]->getCrds();
    dcrd2 = theNodes[1]->getCrds();
    dcrd3 = theNodes[2]->getCrds();

    // call the base class method
    this->DomainComponent::setDomain(theDomain);

    UpdateBase(GsPts[0][0], GsPts[0][0]);

    Bmat(0,0) = 0.5;
    Bmat(1,1) = 0.5;
    Bmat(2,2) = 0.5;
    Bmat(3,0) = 0.5;
    Bmat(4,1) = 0.5;
    Bmat(5,2) = 0.5;
    Bmat(6,0) = 0.5;
    Bmat(7,1) = 0.5;
    Bmat(8,2) = 0.5;
}

int
LysmerTriangle::commitState()
{
	 int retVal = 0;

    // call element commitState to do any base class stuff
    if ((retVal = this->Element::commitState()) != 0) {
    	opserr << "LysmerTriangle::commitState () - failed in base class";
    }    

    return 0; 
}

int
LysmerTriangle::revertToLastCommit()
{
	return 0;
}

int
LysmerTriangle::revertToStart()
{
	return 0;
}

int
LysmerTriangle::update(void)
{
	return 0;
}

int
LysmerTriangle::UpdateBase(double Xi, double Eta)
// this function calculates g1, g2, NI, and nHat for given Xi and Eta
{
    g1 = (dcrd2 - dcrd1);
    g2 = (dcrd3 - dcrd1);

  	// shape functions
  	myNI(0) = 0.5;
  	myNI(1) = 0.5;
  	myNI(2) = 0.5;


	  // normal vector to primary surface as cross product of g1 and g2
    myNhat(0) = g1(1)*g2(2) - g1(2)*g2(1);
    myNhat(1) = g1(2)*g2(0) - g1(0)*g2(2);
    myNhat(2) = g1(0)*g2(1) - g1(1)*g2(0);

    A = myNhat.Norm()/2;

    myNhat.Normalize();


    // tangent vector 1 is g1
    myThat(0) = g1(0);
    myThat(1) = g1(1);
    myThat(2) = g1(2);

    myThat.Normalize();


  // tangent vector 2 is N x T
    myShat(0) = myNhat(1)*g1(2) - myNhat(2)*g1(1);
    myShat(1) = myNhat(2)*g1(0) - myNhat(0)*g1(2);
    myShat(2) = myNhat(0)*g1(1) - myNhat(1)*g1(0);

    myShat.Normalize();

    return 0;
}

const Matrix &
LysmerTriangle::getTangentStiff(void)
{
    tangentStiffness.Zero();
    // = 0 (pure damping) 
    // = 1 (pure stiffness) 
    // = 2 (damping and stiffness) 
    // = 3 (damping but preserve elastic forces from springs after gravity analysis)
    if (stage == 1 || stage == 2)  
    {
      double L = 0; // Actual length used
      if (element_length != 0)
        L = element_length;
      else
      {
        // Compute an equivalent length based as average side length
        double a = (dcrd1 - dcrd2).Norm();
        double b = (dcrd2 - dcrd3).Norm();
        double c = (dcrd1 - dcrd3).Norm();
        L = (a+b+c)/3;
      }
      double G = rho*Vs*Vs;
      double M = rho*Vp*Vp;
      double E = G*(3*M-4*G)/(M-G);

      static Matrix subStiff(3,3);
      static Matrix T(3,3);
      static Matrix K(3,3);
      subStiff.Zero();
      tangentStiffness.Zero();
      K.Zero();
      T.Zero();
      // K(0,0) = G/L;
      // K(1,1) = G/L;
      K(2,2) = E/L;

      for (int i = 0; i < 3; ++i)
      {
         T(0,i) = myThat(i);
         T(1,i) = myShat(i);
         T(2,i) = myNhat(i);
      }

      subStiff.addMatrixTripleProduct(1., T, K, A);

      // static bool do_once = true;

      // if(do_once)
      // {
      //   opserr << "L = " << L << endln;
      //   opserr << "G = " << G << endln;
      //   opserr << "M = " << M << endln;
      //   opserr << "E = " << E << endln;
      //   opserr << "A = " << A << endln;
      //   opserr << "subStiff = " << subStiff << endln;
      //   opserr << "T = " << T << endln;
      //   opserr << "K = " << K << endln;
      //   do_once = false;
      // }

      tangentStiffness.addMatrixTripleProduct(1, Bmat, subStiff, 1.0);
    }

    return tangentStiffness;
}

const Matrix &
LysmerTriangle::getInitialStiff(void)
{
    return getTangentStiff();
}
   
const Matrix &
LysmerTriangle::getDamp(void)
{

    // = 0 (pure damping) 
    // = 1 (pure stiffness) 
    // = 2 (damping and stiffness) 
    // = 3 (damping but preserve elastic forces from springs after gravity analysis)
    tangentDamping.Zero();
    if(stage == 0 || stage == 2 || stage == 3)
    {
      static Matrix subDamp(3,3);
      static Matrix T(3,3);
      static Matrix C(3,3);
      subDamp.Zero();
      tangentDamping.Zero();
      T.Zero();
      C.Zero();
      C(0,0) = rho*Vs;
      C(1,1) = rho*Vs;
      C(2,2) = rho*Vp;
  
      for (int i = 0; i < 3; ++i)
      {
         T(0,i) = myThat(i);
         T(1,i) = myShat(i);
         T(2,i) = myNhat(i);
      }
  
      subDamp.addMatrixTripleProduct(1., T, C, A);
      tangentDamping.addMatrixTripleProduct(1,Bmat, subDamp, 1.0);
    }

    return tangentDamping;
} 


void 
LysmerTriangle::zeroLoad(void)
{
  gnd_velocity.Zero();
  return;
}

int 
LysmerTriangle::addLoad(ElementalLoad *theLoad, double loadFactor)
{
	int type;
  const Vector &data = theLoad->getData(type, loadFactor);

  if (type == LOAD_TAG_LysmerVelocityLoader) {
    gnd_velocity += data;
    return 0;
  } else {
    opserr << "LysmerTriangle::addLoad() - ele with tag: " << this->getTag() << " does not accept load type: " << type << endln;
    return -1;
  }

  return -1;
}

int 
LysmerTriangle::addInertiaLoadToUnbalance(const Vector &accel)
{
	return 0;
}

const Vector &
LysmerTriangle::getResistingForce()
{
  // = 0 (pure damping) 
  // = 1 (pure stiffness) 
  // = 2 (damping and stiffness) 
  // = 3 (damping but preserve elastic forces from springs after gravity analysis)
  if (stage == 0)
  {
    springForces.Zero();
  }
  else if ( stage == 1  || stage == 2 )
  {
    static Vector displacements(9);
    springForces.Zero();

    tangentStiffness = getTangentStiff();

    int count = 0;
    for (int node = 0; node < 3; ++node)
    {
      const Vector &d = theNodes[node]->getTrialDisp();
      displacements(count++) = d(0);
      displacements(count++) = d(1);
      displacements(count++) = d(2);
    }

    springForces.addMatrixVector(0, tangentStiffness, displacements, 1.0);
    // internalForces += springForces;
  } 
  if (stage == 3) //whatever is in spring forces is directly added as a reaction force
  {
    internalForces -= springForces;
  }


  return internalForces;
}

const Vector &
LysmerTriangle::getResistingForceIncInertia()
{
  // = 0 (pure damping) 
  // = 1 (pure stiffness) 
  // = 2 (damping and stiffness) 
  // = 3 (damping but preserve elastic forces from springs after gravity analysis)
  if (stage == 0 || stage == 2 || stage == 3) 
  {    
    static Vector velocities(9);
    internalForces.Zero();

    tangentDamping = getDamp();

    int count = 0;
    for (int node = 0; node < 3; ++node)
    {
      const Vector &v = theNodes[node]->getVel();
      velocities(count++) = 0*v(0) + gnd_velocity(0);
      velocities(count++) = 0*v(1) + gnd_velocity(1);
      velocities(count++) = 0*v(2) + gnd_velocity(2);
    }

    internalForces.addMatrixVector(0.0, tangentDamping, velocities, 1.0);
  } 
  else if (stage == 1)
  {
    internalForces.Zero();
  }

  // Now add stiffness terms
  return getResistingForce();
}

int
LysmerTriangle::sendSelf(int commitTag, Channel &theChannel)
{
  int res;

  // note: we don't check for dataTag == 0 for Element
  // objects as that is taken care of in a commit by the Domain
  // object - don't want to have to do the check if sending data
  int dataTag = this->getDbTag();

  // LysmerTriangle packs it's data into a Vector and sends this to theChannel
  // along with it's dbTag and the commitTag passed in the arguments

  static Vector data(9);
  data(0) = this->getTag();
  data(1) = SL_NUM_DOF;
  data(2) = rho;
  data(3) = Vs;
  data(4) = Vp;
  data(5) = mLoadFactor;
  data(6) = element_length;
  data(7) = stage;
  data(8) = A;

  res = theChannel.sendVector(dataTag, commitTag, data);
  if (res < 0) {
    opserr <<"WARNING LysmerTriangle::sendSelf() - " << this->getTag() << " failed to send data\n";
    return -1;
  }           

  // LysmerTriangle then sends the tags of it's four nodes
  res = theChannel.sendID(dataTag, commitTag, myExternalNodes);
  if (res < 0) {
    opserr <<"WARNING LysmerTriangle::sendSelf() - " << this->getTag() << " failed to send myExternalNodes\n";
    return -2;
  }

  res = theChannel.sendVector(dataTag, commitTag, internalForces);
  if (res < 0) {
    opserr <<"WARNING LysmerTriangle::sendSelf() - " << this->getTag() << " failed to send internalForces\n";
    return -2;
  }
  res = theChannel.sendVector(dataTag, commitTag, springForces);
  if (res < 0) {
    opserr <<"WARNING LysmerTriangle::sendSelf() - " << this->getTag() << " failed to send internalForces\n";
    return -2;
  }
  res = theChannel.sendVector(dataTag, commitTag, g1);
  if (res < 0) {
    opserr <<"WARNING LysmerTriangle::sendSelf() - " << this->getTag() << " failed to send g1\n";
    return -2;
  }
  res = theChannel.sendVector(dataTag, commitTag, g2);
  if (res < 0) {
    opserr <<"WARNING LysmerTriangle::sendSelf() - " << this->getTag() << " failed to send g2\n";
    return -2;
  }
  res = theChannel.sendVector(dataTag, commitTag, myNhat);
  if (res < 0) {
    opserr <<"WARNING LysmerTriangle::sendSelf() - " << this->getTag() << " failed to send myNhat\n";
    return -2;
  }
  res = theChannel.sendVector(dataTag, commitTag, myNI);
  if (res < 0) {
    opserr <<"WARNING LysmerTriangle::sendSelf() - " << this->getTag() << " failed to send myNI\n";
    return -2;
  }
  res = theChannel.sendVector(dataTag, commitTag, dcrd1);
  if (res < 0) {
    opserr <<"WARNING LysmerTriangle::sendSelf() - " << this->getTag() << " failed to send dcrd1\n";
    return -2;
  }
  res = theChannel.sendVector(dataTag, commitTag, dcrd2);
  if (res < 0) {
    opserr <<"WARNING LysmerTriangle::sendSelf() - " << this->getTag() << " failed to send dcrd2\n";
    return -2;
  }
  res = theChannel.sendVector(dataTag, commitTag, dcrd3);
  if (res < 0) {
    opserr <<"WARNING LysmerTriangle::sendSelf() - " << this->getTag() << " failed to send dcrd3\n";
    return -2;
  }


  return 0;
}

int
LysmerTriangle::recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker)
{
  int res;
  int dataTag = this->getDbTag();

  // LysmerTriangle creates a Vector, receives the Vector and then sets the 
  // internal data with the data in the Vector
  static Vector data(9);
  res = theChannel.recvVector(dataTag, commitTag, data);
  if (res < 0) {
    opserr <<"WARNING LysmerTriangle::recvSelf() - failed to receive Vector\n";
    return -1;
  }           
  
  MyTag = (int)data(0);
  rho = data(2);
  Vs = data(3);
  Vp = data(4);
  mLoadFactor = data(5);
  element_length = data(6);
  stage = (int )data(7);
  A = data(8);


  this->setTag((int)MyTag);

  // if(stage != 0)
  // {
  //   opserr << "LysmerTriangle::recvSelf - (chantag = " << theChannel.getTag() <<" ) element at tag # " << this->getTag() << " starting at stage = " << stage << " also L = " << element_length << endln;
  // }

  // LysmerTriangle now receives the tags of it's four external nodes
  res = theChannel.recvID(dataTag, commitTag, myExternalNodes);
  if (res < 0) {
    opserr <<"WARNING LysmerTriangle::recvSelf() - " << this->getTag() << " failed to receive ID\n";
    return -2;
  }

  res = theChannel.recvVector(dataTag, commitTag, internalForces);
  if (res < 0) {
    opserr <<"WARNING LysmerTriangle::sendSelf() - " << this->getTag() << " failed to receive internalForces\n";
    return -2;
  }
  res = theChannel.recvVector(dataTag, commitTag, springForces);
  if (res < 0) {
    opserr <<"WARNING LysmerTriangle::sendSelf() - " << this->getTag() << " failed to receive internalForces\n";
    return -2;
  }
  res = theChannel.recvVector(dataTag, commitTag, g1);
  if (res < 0) {
    opserr <<"WARNING LysmerTriangle::sendSelf() - " << this->getTag() << " failed to receive g1\n";
    return -2;
  }
  res = theChannel.recvVector(dataTag, commitTag, g2);
  if (res < 0) {
    opserr <<"WARNING LysmerTriangle::sendSelf() - " << this->getTag() << " failed to receive g2\n";
    return -2;
  }
  res = theChannel.recvVector(dataTag, commitTag, myNhat);
  if (res < 0) {
    opserr <<"WARNING LysmerTriangle::sendSelf() - " << this->getTag() << " failed to receive myNhat\n";
    return -2;
  }
  res = theChannel.recvVector(dataTag, commitTag, myNI);
  if (res < 0) {
    opserr <<"WARNING LysmerTriangle::sendSelf() - " << this->getTag() << " failed to receive myNI\n";
    return -2;
  }
  res = theChannel.recvVector(dataTag, commitTag, dcrd1);
  if (res < 0) {
    opserr <<"WARNING LysmerTriangle::sendSelf() - " << this->getTag() << " failed to receive dcrd1\n";
    return -2;
  }
  res = theChannel.recvVector(dataTag, commitTag, dcrd2);
  if (res < 0) {
    opserr <<"WARNING LysmerTriangle::sendSelf() - " << this->getTag() << " failed to receive dcrd2\n";
    return -2;
  }
  res = theChannel.recvVector(dataTag, commitTag, dcrd3);
  if (res < 0) {
    opserr <<"WARNING LysmerTriangle::sendSelf() - " << this->getTag() << " failed to receive dcrd3\n";
    return -2;
  }


  return 0;
}

int
LysmerTriangle::displaySelf(Renderer &theViewer, int displayMode, float fact, const char **modes, int numMode)
{
  return 0;
}

void
LysmerTriangle::Print(OPS_Stream &s, int flag)
{
    opserr << "LysmerTriangle, element id:  " << this->getTag() << endln;
    opserr << "   Connected external nodes:  " ; 
    for (int i = 0; i<SL_NUM_NODE; i++) {
      opserr << myExternalNodes(i)<< " ";
    }
    opserr << endln ; 
    opserr << "   A:  " << A << endln ; 

    opserr << "   g1  : " <<  g1 << endln;
    opserr << "   g2  : " <<  g2 << endln;
    opserr << "   myNhat  : " <<  myNhat << endln;
    opserr << "   myThat  : " <<  myThat << endln;
    opserr << "   myShat  : " <<  myShat << endln;
    opserr << "   myNI  : " <<  myNI << endln;
    opserr << "   dcrd1  : " <<  dcrd1 << endln;
    opserr << "   dcrd2  : " <<  dcrd2 << endln;
    opserr << "   dcrd3  : " <<  dcrd3 << endln;
    opserr << "   gnd_velocity  : " <<  gnd_velocity << endln;

    return;
}

Response*
LysmerTriangle::setResponse(const char **argv, int argc, Information &eleInfo)
{
    return 0;
}

int 
LysmerTriangle::getResponse(int responseID, Information &eleInfo)
{
	return -1;
}



int
LysmerTriangle::setParameter(const char **argv, int argc, Parameter &param)
{
  if (argc < 1)
    return -1;

  int res = -1;

  if (strcmp(argv[0],"stage") == 0) {
    param.setValue(stage);
    return param.addObject(1, this);
  } else if (strcmp(argv[0],"rho") == 0) {
    param.setValue(rho);
    return param.addObject(2, this);
  } else if (strcmp(argv[0],"Vp") == 0) {
    param.setValue(Vp);
    return param.addObject(3, this);
  } else if (strcmp(argv[0],"Vs") == 0) {
    param.setValue(Vs);
    return param.addObject(4, this);
  }

  return res;
}
    
int
LysmerTriangle::updateParameter(int parameterID, Information &info)
{
  switch(parameterID) {
  case 1:
    
    stage = (int) info.theDouble;
    return 0;
  case 2:
    rho = info.theDouble;
    return 0;
  case 3:
    Vp = info.theDouble;
    return 0;
  case 4:
    Vs = info.theDouble;
    return 0;
  default:
    return -1;
  }
}
