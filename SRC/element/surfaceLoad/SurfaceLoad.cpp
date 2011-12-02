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
                                                                        
// $Revision$
// $Date$
// $Source$
                                                                        
// Written: crm
// Created: 04/09
//
// Description: This file contains the implementation for the SurfaceLoad class.

#include "SurfaceLoad.h"
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

#include <math.h>
#include <stdlib.h>
#include <stdio.h> 

double SurfaceLoad :: oneOverRoot3 = 1.0/sqrt(3.0);
double SurfaceLoad :: GsPts[4][2];

#include <elementAPI.h>
static int num_SurfaceLoad = 0;

void *
OPS_SurfaceLoad(void)
{
  if (num_SurfaceLoad == 0) {
    num_SurfaceLoad++;
    OPS_Error("SurfaceLoad element - Written by K.Petek, U.Washington\n", 1);
  }

  // Pointer to a uniaxial material that will be returned
  Element *theElement = 0;

  if (OPS_GetNumRemainingInputArgs() != 6) {
    opserr << "Want: element SurfaceLoad eleTag  iNode jNode kNode lNode pressure\n";
    return 0;
  }
    
  int    iData[5];
  double dData[1];

  int numData = 5;
  if (OPS_GetIntInput(&numData, iData) != 0) {
    opserr << "WARNING invalid integer data: element SurfaceLoadElement" << endln;
    return 0;
  }

  numData = 1;
  if (OPS_GetDoubleInput(&numData, dData) != 0) {
    opserr << "WARNING invalid data: element SurfaceLoad " << iData[0] << endln;
    return 0;	
  }

  // Parsing was successful, allocate the material
  theElement = new SurfaceLoad(iData[0], iData[1], iData[2], iData[3], iData[4], dData[0]);

  if (theElement == 0) {
    opserr << "WARNING could not create element of type SurfaceLoadElement\n";
    return 0;
  }

  return theElement;
}



// constructors:
SurfaceLoad::SurfaceLoad(int tag, int Nd1, int Nd2, int Nd3, int Nd4, double pressure)
 :Element(tag,ELE_TAG_SurfaceLoad),     
   myExternalNodes(SL_NUM_NODE),
   tangentStiffness(SL_NUM_DOF, SL_NUM_DOF),
   internalForces(SL_NUM_DOF),
   g1(SL_NUM_NDF), 
   g2(SL_NUM_NDF),
   myNhat(SL_NUM_NDF), 
   myNI(SL_NUM_NODE),
   dcrd1(SL_NUM_NDF),
   dcrd2(SL_NUM_NDF),
   dcrd3(SL_NUM_NDF),
   dcrd4(SL_NUM_NDF)
{
#ifdef DEBUG
        opserr << "SurfaceLoad::SurfaceLoad(): " << MyTag << endln;
#endif
        myExternalNodes(0) = Nd1;
        myExternalNodes(1) = Nd2;
        myExternalNodes(2) = Nd3;
        myExternalNodes(3) = Nd4;

	MyTag = tag;

	GsPts[0][0] = -oneOverRoot3;
	GsPts[0][1] = -oneOverRoot3;
	GsPts[1][0] = oneOverRoot3;
	GsPts[1][1] = -oneOverRoot3;
	GsPts[2][0] = oneOverRoot3;
	GsPts[2][1] = oneOverRoot3;
	GsPts[3][0] = -oneOverRoot3;
	GsPts[3][1] = oneOverRoot3;

        my_pressure = pressure;
}

SurfaceLoad::SurfaceLoad()
 :Element(0,ELE_TAG_SurfaceLoad),     
   myExternalNodes(SL_NUM_NODE),
   tangentStiffness(SL_NUM_DOF, SL_NUM_DOF),
   internalForces(SL_NUM_DOF),
   g1(SL_NUM_NDF), 
   g2(SL_NUM_NDF),
   myNhat(SL_NUM_NDF), 
   myNI(SL_NUM_NODE),
   dcrd1(SL_NUM_NDF),
   dcrd2(SL_NUM_NDF),
   dcrd3(SL_NUM_NDF),
   dcrd4(SL_NUM_NDF)
{
#ifdef DEBUG
        opserr << "SurfaceLoad::SurfaceLoad(): " << MyTag << endln;
#endif
}


//  destructor:
SurfaceLoad::~SurfaceLoad()
{
#ifdef DEBUG
        opserr << "SurfaceLoad::~SurfaceLoad(): " << MyTag << endln;
#endif
}


int
SurfaceLoad::getNumExternalNodes(void) const
{
#ifdef DEBUG
        opserr << "SurfaceLoad::getNumExternalNodes(): " << MyTag << endln;
#endif
    	return SL_NUM_NODE;
}

const ID &
SurfaceLoad::getExternalNodes(void) 
{
#ifdef DEBUG
        opserr << "SurfaceLoad::getExternalNodes(): " << MyTag << endln;
#endif
    	return myExternalNodes;
}


Node **
SurfaceLoad::getNodePtrs(void)
{
#ifdef DEBUG
        opserr << "SurfaceLoad::getNodePtrs(): " << MyTag << endln;
#endif
        return theNodes;                        
}

int
SurfaceLoad::getNumDOF(void) 
{
#ifdef DEBUG
        opserr << "SurfaceLoad::getNumDOF(): " << MyTag << endln;
#endif
    	return SL_NUM_DOF;
}


void
SurfaceLoad::setDomain(Domain *theDomain)
{
#ifdef DEBUG
        opserr << "SurfaceLoad::setDomain(Domain *theDomain): " << MyTag << endln;
#endif

        theNodes[0] = theDomain->getNode(myExternalNodes(0));
        theNodes[1] = theDomain->getNode(myExternalNodes(1));
        theNodes[2] = theDomain->getNode(myExternalNodes(2));
        theNodes[3] = theDomain->getNode(myExternalNodes(3));

        for (int i = 0; i < 4; i++) {
           if (theNodes[i] == 0)
                return;  // don't go any further - otherwise segmentation fault
        }

        dcrd1 = theNodes[0]->getCrds();
        dcrd2 = theNodes[1]->getCrds();
        dcrd3 = theNodes[2]->getCrds();
        dcrd4 = theNodes[3]->getCrds();

       // call the base class method
        this->DomainComponent::setDomain(theDomain);
}

int
SurfaceLoad::commitState()
{
#ifdef DEBUG
        opserr << "SurfaceLoad::commitState(): " << MyTag << endln;
#endif

        int retVal = 0;
        // call element commitState to do any base class stuff
        if ((retVal = this->Element::commitState()) != 0) {
                opserr << "SurfaceLoad::commitState () - failed in base class";
                }    

        return 0; 
}

int
SurfaceLoad::revertToLastCommit()
{
#ifdef DEBUG
        opserr << "SurfaceLoad::revertToLastCommit(): " << MyTag << endln;
#endif
        return 0;
}

int
SurfaceLoad::revertToStart()
{
#ifdef DEBUG
        opserr << "SurfaceLoad::revertToStart(): " << MyTag << endln;
#endif
        return 0;
}

int
SurfaceLoad::update(void)
{
#ifdef DEBUG
        opserr << "SurfaceLoad::update(): " << MyTag << endln;
#endif
        return 0;
}

int
SurfaceLoad::UpdateBase(double Xi, double Eta)
// this function calculates g1, g2, NI, and nHat for given Xi and Eta
{

#ifdef DEBUG
        opserr << "SurfaceLoad::UpdateBase(): " << MyTag << endln;
#endif

        double oneMinusEta = 1 - Eta;
        double onePlusEta  = 1 + Eta;
        double oneMinusXi  = 1 - Xi;
        double onePlusXi   = 1 + Xi;

        // calculate vectors g1 and g2
        // g1 = d(x_Xi)/dXi, g2 = d(x_Xi)/dEta
        g1 = (oneMinusEta * (dcrd2 - dcrd1) + onePlusEta  * (dcrd3 - dcrd4)) * 0.25;
        g2 = (onePlusXi   * (dcrd3 - dcrd2) + oneMinusXi  * (dcrd4 - dcrd1)) * 0.25;

		// shape functions
		myNI(0) = 0.25 * oneMinusXi * oneMinusEta;
		myNI(1) = 0.25 * onePlusXi  * oneMinusEta;
		myNI(2) = 0.25 * onePlusXi  * onePlusEta;
		myNI(3) = 0.25 * oneMinusXi * onePlusEta;

		// normal vector to master surface as cross product of g1 and g2
        myNhat(0) = g1(1)*g2(2) - g1(2)*g2(1);
        myNhat(1) = g1(2)*g2(0) - g1(0)*g2(2);
        myNhat(2) = g1(0)*g2(1) - g1(1)*g2(0);

        return 0;
}

const Matrix &
SurfaceLoad::getTangentStiff(void)
{
#ifdef DEBUG
        opserr << "SurfaceLoad::getTangentStiff(): " << MyTag << endln;
#endif
        // initialize Kt
        tangentStiffness.Zero();
        return tangentStiffness;
}

const Matrix &
SurfaceLoad::getInitialStiff(void)
{
#ifdef DEBUG
        opserr << "SurfaceLoad::getInitialStiff(): " << MyTag << endln;
#endif
        return getTangentStiff();
}
    
void 
SurfaceLoad::zeroLoad(void)
{
#ifdef DEBUG
        opserr << "SurfaceLoad::zeroLoad(): " << MyTag << endln;
#endif
        return;
}

int 
SurfaceLoad::addLoad(ElementalLoad *theLoad, double loadFactor)
{
#ifdef DEBUG
        opserr << "SurfaceLoad::addLoad(ElementalLoad *theLoad, double loadFactor): " << MyTag << endln;
#endif
	// my_pressure += theLoad->Load(MyTag, ???) * loadFactor;
	return 0;
}

int 
SurfaceLoad::addInertiaLoadToUnbalance(const Vector &accel)
{
#ifdef DEBUG
        opserr << "SurfaceLoad::addInertiaLoadToUnbalance(const Vector &accel): " << MyTag << endln;
#endif
	return 0;
}

const Vector &
SurfaceLoad::getResistingForce()
{
#ifdef DEBUG
        opserr << "SurfaceLoad::getResistingForce(): " << MyTag << endln;
#endif
	internalForces.Zero();

	// loop over Gauss points
	for(int i = 0; i < 4; i++) {
		this->UpdateBase(GsPts[i][0],GsPts[i][1]);

		// loop over nodes
		for(int j = 0; j < 4; j++) {
			// loop over dof
			for(int k = 0; k < 3; k++) {
				internalForces[j*3+k] = internalForces[j*3+k] - my_pressure * myNhat(k) * myNI(j);
			}
		}
	}
#ifdef DEBUG
        opserr << "SurfaceLoad::getResistingForce(): " << internalForces << endln;
#endif

	return internalForces;
}

const Vector &
SurfaceLoad::getResistingForceIncInertia()
{       
#ifdef DEBUG
        opserr << "SurfaceLoad::getResistingForceIncInertia(): " << MyTag << endln;
#endif
  		return theVector;
}

int
SurfaceLoad::sendSelf(int commitTag, Channel &theChannel)
{
#ifdef DEBUG
        opserr << "SurfaceLoad::sendSelf(int commitTag, Channel &theChannel)" << endln;
#endif
        
  int res;

  // note: we don't check for dataTag == 0 for Element
  // objects as that is taken care of in a commit by the Domain
  // object - don't want to have to do the check if sending data
  int dataTag = this->getDbTag();

  // SurfaceLoad packs it's data into a Vector and sends this to theChannel
  // along with it's dbTag and the commitTag passed in the arguments

  static Vector data(2);
  data(0) = this->getTag();
  data(1) = SL_NUM_DOF;

  res = theChannel.sendVector(dataTag, commitTag, data);
  if (res < 0) {
    opserr <<"WARNING SurfaceLoad::sendSelf() - " << this->getTag() << " failed to send Vector\n";
    return -1;
  }           

  // SurfaceLoad then sends the tags of it's four nodes

  res = theChannel.sendID(dataTag, commitTag, myExternalNodes);
  if (res < 0) {
    opserr <<"WARNING SurfaceLoad::sendSelf() - " << this->getTag() << " failed to send Vector\n";
    return -2;
  }

  return 0;
}

int
SurfaceLoad::recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker)
{
#ifdef DEBUG
        opserr << "SurfaceLoad::recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker)" << endln;
#endif

  int res;
  int dataTag = this->getDbTag();

  // SurfaceLoad creates a Vector, receives the Vector and then sets the 
  // internal data with the data in the Vector

  static Vector data(2);
  res = theChannel.recvVector(dataTag, commitTag, data);
  if (res < 0) {
    opserr <<"WARNING SurfaceLoad::recvSelf() - failed to receive Vector\n";
    return -1;
  }           

  // SurfaceLoad now receives the tags of it's four external nodes
  res = theChannel.recvID(dataTag, commitTag, myExternalNodes);
  if (res < 0) {
    opserr <<"WARNING SurfaceLoad::recvSelf() - " << this->getTag() << " failed to receive ID\n";
    return -2;
  }
    return 0;
}

int
SurfaceLoad::displaySelf(Renderer &theViewer, int displayMode, float fact)
{
#ifdef DEBUG
        opserr << "SurfaceLoad::displaySelf(Renderer &theViewer, int displayMode, float fact)" << endln;
#endif
  return 0;
}

void
SurfaceLoad::Print(OPS_Stream &s, int flag)
{
#ifdef DEBUG
        opserr << "SurfaceLoad::Print(OPS_Stream &s, int flag)" << endln;
#endif
        opserr << "SurfaceLoad, element id:  " << this->getTag() << endln;
        opserr << "   Connected external nodes:  " ; 
        for (int i = 0; i<SL_NUM_NODE; i++)
        {
                opserr << myExternalNodes(i)<< " ";
        }
        return;
}

Response*
SurfaceLoad::setResponse(const char **argv, int argc, Information &eleInfo)
{
#ifdef DEBUG
        opserr << "SurfaceLoad::setResponse(const char **argv, int argc, Information &eleInfo): " << MyTag << endln;
#endif
      	return 0;
}

int 
SurfaceLoad::getResponse(int responseID, Information &eleInfo)
{
#ifdef DEBUG
        opserr << "SurfaceLoad::getResponse(int responseID, Information &eleInfo): " << MyTag << endln;
#endif
	    return -1;
}

