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

// Description: This file contains the implementation for the TriSurfaceLoad class.

#include "TriSurfaceLoad.h"
#include <Information.h>
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

double TriSurfaceLoad::oneOverRoot3 = 1.0/sqrt(3.0);
double TriSurfaceLoad::GsPts[1][1];
Matrix TriSurfaceLoad::tangentStiffness(SL_NUM_DOF, SL_NUM_DOF);
Matrix TriSurfaceLoad::mass(SL_NUM_DOF, SL_NUM_DOF);
Matrix TriSurfaceLoad::damp(SL_NUM_DOF, SL_NUM_DOF);
Vector TriSurfaceLoad::internalForces(SL_NUM_DOF);

#include <elementAPI.h>
static int num_TriSurfaceLoad = 0;

void *
OPS_TriSurfaceLoad(void)
{
  if (num_TriSurfaceLoad == 0) {
    num_TriSurfaceLoad++;
    opserr<<"TriSurfaceLoad element - Written: J. A. Abell (UANDES). Inspired by the makers of SurfaceLoad\n";
  }

  // Pointer to a uniaxial material that will be returned
  Element *theElement = 0;

  if (OPS_GetNumRemainingInputArgs() < 5) {
    opserr << "Want: element TriSurfaceLoad eleTag?  iNode? jNode? kNode? pressure? <rhoH?>\n";
    return 0;
  }
    
  int    iData[4];
  double dData[1];

  int numData = 4;
  if (OPS_GetIntInput(&numData, iData) != 0) {
    opserr << "WARNING invalid integer data: element TriSurfaceLoadElement" << endln;
    return 0;
  }

  numData = 1;
  if (OPS_GetDoubleInput(&numData, dData) != 0) {
    opserr << "WARNING invalid data: element TriSurfaceLoad " << iData[0] << endln;
    return 0;	
  }

  double rhoH = 0;
  int num_args_remaining = OPS_GetNumRemainingInputArgs();

  if (num_args_remaining > 0)
  {
    numData = 1;
    OPS_GetDoubleInput(&numData, &rhoH);
  }


  // Parsing was successful, allocate the material
  theElement = new TriSurfaceLoad(iData[0], iData[1], iData[2], iData[3], dData[0], rhoH);

  if (theElement == 0) {
    opserr << "WARNING could not create element of type TriSurfaceLoadElement\n";
    return 0;
  }

  return theElement;
}

// constructors:
TriSurfaceLoad::TriSurfaceLoad(int tag, int Nd1, int Nd2, int Nd3, double pressure,  double rhoH_)
 :Element(tag,ELE_TAG_TriSurfaceLoad),     
   myExternalNodes(SL_NUM_NODE),
   g1(SL_NUM_NDF), 
   g2(SL_NUM_NDF),
   myNhat(SL_NUM_NDF), 
   myNI(SL_NUM_NODE),
   dcrd1(SL_NUM_NDF),
   dcrd2(SL_NUM_NDF),
   dcrd3(SL_NUM_NDF)
{
    myExternalNodes(0) = Nd1;
    myExternalNodes(1) = Nd2;
    myExternalNodes(2) = Nd3;

	GsPts[0][0] = 0.5;

  my_pressure = pressure;
  rhoH = rhoH_;

	mLoadFactor = 1.0;
}

TriSurfaceLoad::TriSurfaceLoad()
  :Element(0,ELE_TAG_TriSurfaceLoad),     
   	myExternalNodes(SL_NUM_NODE),
   	g1(SL_NUM_NDF), 
   	g2(SL_NUM_NDF),
   	myNhat(SL_NUM_NDF), 
   	myNI(SL_NUM_NODE),
   	dcrd1(SL_NUM_NDF),
   	dcrd2(SL_NUM_NDF),
   	dcrd3(SL_NUM_NDF)
{
  GsPts[0][0] = 0;
  my_pressure = 0;
  rhoH = 0;
  mLoadFactor = 0;
}

//  destructor:
TriSurfaceLoad::~TriSurfaceLoad()
{
}

int
TriSurfaceLoad::getNumExternalNodes(void) const
{
    return SL_NUM_NODE;
}

const ID &
TriSurfaceLoad::getExternalNodes(void) 
{
    return myExternalNodes;
}

Node **
TriSurfaceLoad::getNodePtrs(void)
{
    return theNodes;                        
}

int
TriSurfaceLoad::getNumDOF(void) 
{
    return SL_NUM_DOF;
}

void
TriSurfaceLoad::setDomain(Domain *theDomain)
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
}

int
TriSurfaceLoad::commitState()
{
	int retVal = 0;
    // call element commitState to do any base class stuff
    if ((retVal = this->Element::commitState()) != 0) {
    	opserr << "TriSurfaceLoad::commitState () - failed in base class";
    }    

    return 0; 
}

int
TriSurfaceLoad::revertToLastCommit()
{
	return 0;
}

int
TriSurfaceLoad::revertToStart()
{
	return 0;
}

int
TriSurfaceLoad::update(void)
{
	return 0;
}

int
TriSurfaceLoad::UpdateBase(double Xi, double Eta)
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

    // Normalize
    // double norm = myNhat.Norm();
    myNhat = myNhat / 3;

    return 0;
}

const Matrix &
TriSurfaceLoad::getTangentStiff(void)
{
    tangentStiffness.Zero();
    return tangentStiffness;
}

const Matrix &
TriSurfaceLoad::getInitialStiff(void)
{
    return getTangentStiff();
}
    
void 
TriSurfaceLoad::zeroLoad(void)
{
    return;
}

int 
TriSurfaceLoad::addLoad(ElementalLoad *theLoad, double loadFactor)
{
	int type;
	const Vector &data = theLoad->getData(type, loadFactor);

	if (type == LOAD_TAG_SurfaceLoader) {
		mLoadFactor = loadFactor;
    // opserr << "Trisurface - loadFactor = " << loadFactor << endln;
		return 0;
	} else {
		opserr << "TriSurfaceLoad::addLoad() - ele with tag: " << this->getTag() << " does not accept load type: " << type << endln;
		return -1;
	}

	return -1;
}

int 
TriSurfaceLoad::addInertiaLoadToUnbalance(const Vector &accel)
{
  return 0;
}

const Vector &
TriSurfaceLoad::getResistingForce()
{
	internalForces.Zero();

	// loop over Gauss points
	for(int i = 0; i < 1; i++) {
		this->UpdateBase(GsPts[i][0],GsPts[i][0]);

		// loop over nodes
		for(int j = 0; j < 3; j++) {
			// loop over dof
			for(int k = 0; k < 3; k++) {
				internalForces[j*3+k] = internalForces[j*3+k] - mLoadFactor*my_pressure*myNhat(k)*myNI(j);
			}
		}
	}

	return internalForces;
}

const Vector &
TriSurfaceLoad::getResistingForceIncInertia()
{   
    static Vector accel(SL_NUM_DOF);
    accel.Zero();
    
    internalForces = getResistingForce();

    int pos = 0;
    for (int i = 0; i < SL_NUM_NODE; ++i)
    {
      const Vector &a = theNodes[i]->getAccel();
      accel(pos++) = a(i);
    }
    
    mass = getMass();
    internalForces.addMatrixVector(1.0, mass, accel, -1.0);
  	
    return internalForces;
}

int
TriSurfaceLoad::sendSelf(int commitTag, Channel &theChannel)
{
  int res;

  // note: we don't check for dataTag == 0 for Element
  // objects as that is taken care of in a commit by the Domain
  // object - don't want to have to do the check if sending data
  int dataTag = this->getDbTag();

  // TriSurfaceLoad packs its data into a Vector and sends this to theChannel
  // along with its dbTag and the commitTag passed in the arguments

  static Vector data(4 + 6*SL_NUM_NDF + SL_NUM_NODE);
  data(0) = this->getTag();
  data(1) = my_pressure;
  data(2) = mLoadFactor;
  data(3) = rhoH;

  for (int i = 0; i < SL_NUM_NDF; i++) {
    data(4+             i) = g1(i);
    data(4+  SL_NUM_NDF+i) = g2(i);
    data(4+2*SL_NUM_NDF+i) = myNhat(i);
    data(4+3*SL_NUM_NDF+i) = dcrd1(i);
    data(4+4*SL_NUM_NDF+i) = dcrd2(i);
    data(4+5*SL_NUM_NDF+i) = dcrd3(i);
  }
  for (int i = 0; i < SL_NUM_NODE; i++)
    data(4+6*SL_NUM_NDF+i) = myNI(i);
  
  res = theChannel.sendVector(dataTag, commitTag, data);
  if (res < 0) {
    opserr <<"WARNING TriSurfaceLoad::sendSelf() - " << this->getTag() << " failed to send data\n";
    return -1;
  }           

  // TriSurfaceLoad then sends the tags of its four nodes
  res = theChannel.sendID(dataTag, commitTag, myExternalNodes);
  if (res < 0) {
    opserr <<"WARNING TriSurfaceLoad::sendSelf() - " << this->getTag() << " failed to send myExternalNodes\n";
    return -2;
  }

  return 0;
}

int
TriSurfaceLoad::recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker)
{
  int res;
  int dataTag = this->getDbTag();

  // TriSurfaceLoad creates a Vector, receives the Vector and then sets the 
  // internal data with the data in the Vector
  static Vector data(4 + 6*SL_NUM_NDF + SL_NUM_NODE);
  res = theChannel.recvVector(dataTag, commitTag, data);
  if (res < 0) {
    opserr <<"WARNING TriSurfaceLoad::recvSelf() - failed to receive Vector\n";
    return -1;
  }           
  
  this->setTag(int(data(0)));
  my_pressure = data(1);
  mLoadFactor = data(2);
  rhoH = data(3);

  for (int i = 0; i < SL_NUM_NDF; i++) {
    g1(i)     = data(4+             i);
    g2(i)     = data(4+  SL_NUM_NDF+i);
    myNhat(i) = data(4+2*SL_NUM_NDF+i);
    dcrd1(i)  = data(4+3*SL_NUM_NDF+i);
    dcrd2(i)  = data(4+4*SL_NUM_NDF+i);
    dcrd3(i)  = data(4+5*SL_NUM_NDF+i);
  }
  for (int i = 0; i < SL_NUM_NODE; i++)
    myNI(i) = data(4+6*SL_NUM_NDF+i);
  
  // TriSurfaceLoad now receives the tags of its four external nodes
  res = theChannel.recvID(dataTag, commitTag, myExternalNodes);
  if (res < 0) {
    opserr <<"WARNING TriSurfaceLoad::recvSelf() - " << this->getTag() << " failed to receive ID\n";
    return -2;
  }

  return 0;
}

int
TriSurfaceLoad::displaySelf(Renderer &theViewer, int displayMode, float fact, const char **modes, int numMode)
{
  return 0;
}

void
TriSurfaceLoad::Print(OPS_Stream &s, int flag)
{
    opserr << "TriSurfaceLoad, element id:  " << this->getTag() << endln;
    opserr << "   Connected external nodes:  " ; 
    for (int i = 0; i<SL_NUM_NODE; i++) {
    	opserr << myExternalNodes(i)<< " ";
    }
    return;
}

Response*
TriSurfaceLoad::setResponse(const char **argv, int argc, OPS_Stream &output)
{
  return Element::setResponse(argv, argc, output);
}

int 
TriSurfaceLoad::getResponse(int responseID, Information &eleInfo)
{
  return Element::getResponse(responseID, eleInfo);
}


// Compute a mass matrix using zangar's idea:
const Matrix & TriSurfaceLoad::getMass(void)
{
    double Area = myNhat.Norm();

    mass.Zero();
    
    if(rhoH > 0)
    {
      for (int i = 0; i < SL_NUM_DOF; ++i)
      {
        mass(i,i) = rhoH * Area /3;
      }
    }

    return mass;
}
