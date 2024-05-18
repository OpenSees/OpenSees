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
// Revision: 1.0
// Date: 05/2019
// Source: /usr/local/cvs/OpenSees/SRC/material/uniaxial/ViscoelasticGap.cpp
//
// Written: Patrick J. Hughes, University of California - San Diego
// Created: 05/2019
//
// Description: This file contains the class implementation for the
// ViscoelasticGap uniaxialMaterial.
//
// References:
// Goldsmith, W. (1960). "Impact: The Theory and Physical Behavior of Colliding Solids." 
//   E. Arnold: London.
//
// Variables:
// K: linear stiffness
// C: linear damping coefficient
// gap: initial gap (must be a negative value)

#include <stdlib.h>
#include <Vector.h>
#include <Channel.h>
#include <math.h>
#include <float.h>

#include <elementAPI.h>
#include <OPS_Globals.h>
#include <ViscoelasticGap.h>

static int numViscoelasticGap = 0;

void *
OPS_ViscoelasticGap(void)
{
  // kudos
  if (numViscoelasticGap == 0) {
    numViscoelasticGap++;
    opserr << "ViscoelasticGap model written by Patrick J. Hughes, UC San Diego\n";
  }
  
  // Pointer to a uniaxial material that will be returned
  UniaxialMaterial *theMaterial = 0;

  // check for minimum number of arguments
  if (OPS_GetNumRemainingInputArgs() < 4) {
    opserr << "WARNING incorrect num args want: uniaxialMaterial ViscoelasticGap tag? K? C? gap?" << endln;
    return 0;
  }

  int    iData[1]; // 1 integer input (tag)
  double dData[3]; // 4 double inputs (K, C, gap)
  
  // check material tag validity
  int numData = 1;
  if (OPS_GetIntInput(&numData, iData) != 0) {
    opserr << "WARNING invalid uniaxialMaterial ViscoelasticGap tag" << endln;
    return 0;
  }

  // check material property validity
  numData = OPS_GetNumRemainingInputArgs();
  // all 3 double inputs given
  if (numData >= 3) {
	  numData = 3;
	  if (OPS_GetDoubleInput(&numData, dData) != 0) {
		  opserr << "Invalid data for uniaxialMaterial ViscoelasticGap " << iData[0] << endln;
		  return 0;
	  }
  }
	 
  // parsing was successful - allocate material
  theMaterial = new ViscoelasticGap(iData[0], dData[0], dData[1], dData[2]);
  if (theMaterial == 0) {
    opserr << "WARNING could not create uniaxialMaterial ViscoelasticGap\n";
    return 0;
  }

  return theMaterial;
}

// constructor
ViscoelasticGap::ViscoelasticGap(int tag, double k, double c, double gap0)
:UniaxialMaterial(tag,MAT_TAG_ViscoelasticGap),
 K(k), C(c), gap(gap0)
{
	if (gap >= 0) {
	  opserr << "ViscoelasticGap::ViscoelasticGap -- Initial gap size should be negative for compression-only material\n";
	  gap = -gap;
	  opserr << "Setting gap to negative value, " << gap << endln;
	}
	this->revertToStart(); // initialize state variables
	printFlag = 0; // initialize print flag for impact event
}

// constructor for parallel processing
ViscoelasticGap::ViscoelasticGap()
:UniaxialMaterial(0,MAT_TAG_ViscoelasticGap),
 K(0.0), C(0.0), gap(0.0)
{

}

// destructor
ViscoelasticGap::~ViscoelasticGap()
{
  // does nothing
}

// state determintaion
int 
ViscoelasticGap::setTrialStrain(double strain, double strainRate)
{
	  // set the trial strain & trial strain rate
  	trialStrain = strain;
	trialStrainRate = strainRate;
  
	// gap is open
	if (trialStrain >= gap) {
		if (commitTangent != 0.0) {
			printFlag = 0; // reset print flag
		}
		trialStress = 0.0; // zero force when out of contact
		trialTangent = 0.0; // zero stiffness when out of contact
	
	// gap is closed
	} else {
		
		// new impact event
		if ((commitTangent == 0.0) & (trialStrainRate < 0.0)) {
			if (printFlag == 0) {
				opserr << "ViscoelasticGap impact detected: impact velocity = " << -trialStrainRate << "\n"; // display impact message
				printFlag = 1; // update print flag
			}
		}
		
		// same impact event as previous state
		trialStress = -( K*(-(trialStrain-gap)) + C*(-trialStrainRate) ); // contact force
		trialTangent = K; // contact stiffness
	}
    return 0;
}

// retrieve trial strain
double 
ViscoelasticGap::getStrain(void)
{
  return trialStrain;
}

// retrieve trial strain rate
double 
ViscoelasticGap::getStrainRate(void)
{
  return trialStrainRate;
}

// retrieve trial stress
double 
ViscoelasticGap::getStress(void)
{
  return trialStress;
}

// retrieve trial tangent stiffness
double 
ViscoelasticGap::getTangent(void)
{
  return trialTangent;
}

// retrieve initial tangent stiffness
double 
ViscoelasticGap::getInitialTangent(void)
{
  return 0.0; 
}

// commit state
int 
ViscoelasticGap::commitState(void)
{
	commitStrain = trialStrain;
	commitStrainRate = trialStrainRate;
	commitStress = trialStress;
	commitTangent = trialTangent;
	return 0;
}

// revert to last comitted state
int 
ViscoelasticGap::revertToLastCommit(void)
{
	trialStrain = commitStrain;
	trialStrainRate = commitStrainRate;
	trialStress = commitStress;
	trialTangent = commitTangent;
    return 0;
}

// revert to initial state
int 
ViscoelasticGap::revertToStart(void)
{
    commitStrain = 0.0;
    commitStrainRate = 0.0;
	commitStress = 0.0;
	commitTangent = 0.0;
    trialStrain = 0.0;
    trialStrainRate = 0.0;
	trialStress = 0.0;
	trialTangent = 0.0;
    return 0;
}

// copy state
UniaxialMaterial *
ViscoelasticGap::getCopy(void)
{
    ViscoelasticGap *theCopy = new ViscoelasticGap(this->getTag(),K,C,gap);
	theCopy->commitStrain = commitStrain;
	theCopy->commitStrainRate = commitStrainRate;
	theCopy->commitStress = commitStress;
	theCopy->commitTangent = commitTangent;
	theCopy->trialStrain = trialStrain;
	theCopy->trialStrainRate = trialStrainRate;
	theCopy->trialStress = trialStress;
	theCopy->trialTangent = trialTangent;
    return theCopy;
}

// send ViscoelasticGap object
int 
ViscoelasticGap::sendSelf(int cTag, Channel &theChannel)
{
  int res = 0;
  static Vector data(12);
  data(0)  = this->getTag();
  data(1)  = K;
  data(2)  = C;
  data(3)  = gap;
  data(4)  = commitStrain;
  data(5)  = commitStrainRate;
  data(6)  = commitStress; 
  data(7)  = commitTangent;
  data(8)  = trialStrain;
  data(9) = trialStrainRate;
  data(10) = trialStress; 
  data(11) = trialTangent;
  res = theChannel.sendVector(this->getDbTag(), cTag, data);
  if (res < 0) 
    opserr << "ViscoelasticGap::sendSelf() - failed to send data\n";
  return res;
}

// receive ViscoelasticGap object
int 
ViscoelasticGap::recvSelf(int cTag, Channel &theChannel,
			  FEM_ObjectBroker &theBroker)
{
  int res = 0;
  static Vector data(12);
  res = theChannel.recvVector(this->getDbTag(), cTag, data);
  if (res < 0)
    opserr << "ViscoelasticGap::recvSelf() - failed to recv data\n";
  else {
    this->setTag((int)data(0));
    K = data(1);
    C = data(2);
    gap = data(3);
    commitStrain = data(4);
    commitStrainRate = data(5);
    commitStress = data(6); 
    commitTangent = data(7);
    trialStrain = data(8);
    trialStrainRate = data(9);
    trialStress = data(10); 
    trialTangent = data(11);
  }
  return res;
}

// print ViscoelasticGap object properties
void 
ViscoelasticGap::Print(OPS_Stream &s, int flag)
{
    s << "ViscoelasticGap tag: " << this->getTag() << endln;
    s << "  K: " << K << endln;
    s << "  C: " << C << endln;
	s << "  gap: " << gap << endln;
}
