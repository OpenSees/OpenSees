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
// Date: 03/2019
// Source: /usr/local/cvs/OpenSees/SRC/material/uniaxial/Hertzdamp.cpp
//
// Written: Patrick J. Hughes, University of California - San Diego
// Created: 03/2019
//
// Description: This file contains the class implementation for the
// Hertzdamp uniaxialMaterial.
//
// References:
// Muthukumar, S., and DesRoches, R. (2006). "A Hertz Contact Model with Non-linear Damping for Pounding Simulation." 
//   Earthquake Engineering and Structural Dynamics, 35, 811-828.
// Ye, Kun., Li, L., and Zhu, H. (2008) "A Note on the Hertz Contact Model with Nonlinear Damping for Pounding Simulation."
//	 Earthquake Engineering and Structural Dynamics, 38, 1135-1142.
//
// Variables:
// Kh: nonlinear Hertz contact stiffness
// xiNorm: normalized damping coefficient
// gap: initial gap distance (must be a negative value)
// n: displacement exponent (default is 1.5)

#include <stdlib.h>
#include <Vector.h>
#include <Channel.h>
#include <math.h>
#include <float.h>

#include <elementAPI.h>
#include <OPS_Globals.h>
#include <Hertzdamp.h>

static int numHertzdamp = 0;

void *
OPS_Hertzdamp(void)
{
  // kudos
  if (numHertzdamp == 0) {
    numHertzdamp++;
    opserr << "Hertzdamp model written by Patrick J. Hughes, UC San Diego\n";
  }
  
  // Pointer to a uniaxial material that will be returned
  UniaxialMaterial *theMaterial = 0;

  // check for minimum number of arguments
  if (OPS_GetNumRemainingInputArgs() < 4) {
    opserr << "WARNING incorrect num args want: uniaxialMaterial Hertzdamp tag? Kh? xiNorm? gap? <n?>" << endln;
    return 0;
  }

  int    iData[1]; // 1 integer input (tag)
  double dData[4]; // 4 double inputs (Kh, xiNorm, gap, n)
  
  // check material tag validity
  int numData = 1;
  if (OPS_GetIntInput(&numData, iData) != 0) {
    opserr << "WARNING invalid uniaxialMaterial Hertzdamp tag" << endln;
    return 0;
  }

  // check material property validity
  numData = OPS_GetNumRemainingInputArgs();
  // all 4 double inputs given
  if (numData >= 4) {
	  numData = 4;
	  if (OPS_GetDoubleInput(&numData, dData) != 0) {
		  opserr << "Invalid data for uniaxialMaterial Hertzdamp " << iData[0] << endln;
		  return 0;
	  }
  } else {
	  // 3 double inputs given - use default value for 4th input
	  numData = 3;
	  if (OPS_GetDoubleInput(&numData, dData) != 0) {
		  opserr << "Invalid data for uniaxialMaterial Hertzdamp " << iData[0] << endln;
		  return 0;
	  }
	  dData[3] = 1.5; // default value for displacement exponent is 1.5
  }
	 
  // parsing was successful - allocate material
  theMaterial = new Hertzdamp(iData[0], dData[0], dData[1], dData[2], dData[3]);
  if (theMaterial == 0) {
    opserr << "WARNING could not create uniaxialMaterial Hertzdamp\n";
    return 0;
  }

  return theMaterial;
}

// constructor
Hertzdamp::Hertzdamp(int tag, double kh, double xin, double gap0, double dispexp)
:UniaxialMaterial(tag,MAT_TAG_Hertzdamp),
 Kh(kh), xiNorm(xin), gap(gap0), n(dispexp)
{
	if (gap >= 0) {
	  opserr << "Hertzdamp::Hertzdamp -- Initial gap size should be negative for compression-only material\n";
	  gap = -gap;
	  opserr << "Setting gap to negative value, " << gap << endln;
	}
	this->revertToStart(); // initialize state variables
	xi = 0.0; // initialize impact damping coefficient
	printFlag = 0; // initialize print flag for impact event
}

// constructor for parallel processing
Hertzdamp::Hertzdamp()
:UniaxialMaterial(0,MAT_TAG_Hertzdamp),
 Kh(0.0), xiNorm(0.0), gap(0.0), n(0.0)
{

}

// destructor
Hertzdamp::~Hertzdamp()
{
  // does nothing
}

// state determintaion
int 
Hertzdamp::setTrialStrain(double strain, double strainRate)
{
	  // set the trial strain & trial strain rate
  	trialStrain = strain;
	trialStrainRate = strainRate;
  
	// gap is open
	if (trialStrain >= gap) {
		if (commitTangent != 0.0) {
			xi = 0.0; // reset damping coefficient
			printFlag = 0; // reset print flag
		}
		trialStress = 0.0; // zero force when out of contact
		trialTangent = 0.0; // zero stiffness when out of contact
	
	// gap is closed
	} else {
		
		// new impact event
		if ((commitTangent == 0.0) & (trialStrainRate < 0.0)) {
			xi = xiNorm * Kh/(-trialStrainRate); // compute damping coefficient
			if (printFlag == 0) {
				opserr << "Hertzdamp impact detected: impact velocity = " << -trialStrainRate << ", damping coefficient = " << xi << "\n"; // display impact message
				printFlag = 1; // update print flag
			}
		}
		
		// same impact event as previous state
		trialStress = -(Kh + xi*(-trialStrainRate)) * pow(-(trialStrain-gap),n); // contact force
		trialTangent = -n * (Kh + xi*(-trialStrainRate)) * pow(-(trialStrain-gap),n-1.0); // contact stiffness
	}
    return 0;
}

// retrieve trial strain
double 
Hertzdamp::getStrain(void)
{
  return trialStrain;
}

// retrieve trial strain rate
double 
Hertzdamp::getStrainRate(void)
{
  return trialStrainRate;
}

// retrieve trial stress
double 
Hertzdamp::getStress(void)
{
  return trialStress;
}

// retrieve trial tangent stiffness
double 
Hertzdamp::getTangent(void)
{
  return trialTangent;
}

// retrieve initial tangent stiffness
double 
Hertzdamp::getInitialTangent(void)
{
  return 0.0; 
}

// commit state
int 
Hertzdamp::commitState(void)
{
	commitStrain = trialStrain;
	commitStrainRate = trialStrainRate;
	commitStress = trialStress;
	commitTangent = trialTangent;
	return 0;
}

// revert to last comitted state
int 
Hertzdamp::revertToLastCommit(void)
{
	trialStrain = commitStrain;
	trialStrainRate = commitStrainRate;
	trialStress = commitStress;
	trialTangent = commitTangent;
    return 0;
}

// revert to initial state
int 
Hertzdamp::revertToStart(void)
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
Hertzdamp::getCopy(void)
{
    Hertzdamp *theCopy = new Hertzdamp(this->getTag(),Kh,xiNorm,gap,n);
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

// send Hertzdamp object
int 
Hertzdamp::sendSelf(int cTag, Channel &theChannel)
{
  int res = 0;
  static Vector data(13);
  data(0)  = this->getTag();
  data(1)  = Kh;
  data(2)  = xiNorm;
  data(3)  = gap;
  data(4)  = n;
  data(5)  = commitStrain;
  data(6)  = commitStrainRate;
  data(7)  = commitStress; 
  data(8)  = commitTangent;
  data(9)  = trialStrain;
  data(10) = trialStrainRate;
  data(11) = trialStress; 
  data(12) = trialTangent;
  res = theChannel.sendVector(this->getDbTag(), cTag, data);
  if (res < 0) 
    opserr << "Hertzdamp::sendSelf() - failed to send data\n";
  return res;
}

// receive Hertzdamp object
int 
Hertzdamp::recvSelf(int cTag, Channel &theChannel,
			  FEM_ObjectBroker &theBroker)
{
  int res = 0;
  static Vector data(13);
  res = theChannel.recvVector(this->getDbTag(), cTag, data);
  if (res < 0)
    opserr << "Hertzdamp::recvSelf() - failed to recv data\n";
  else {
    this->setTag((int)data(0));
    Kh = data(1);
    xiNorm = data(2);
    gap = data(3);
    n = data(4);
    commitStrain = data(5);
    commitStrainRate = data(6);
    commitStress = data(7); 
    commitTangent = data(8);
    trialStrain = data(9);
    trialStrainRate = data(10);
    trialStress = data(11); 
    trialTangent = data(12);
  }
  return res;
}

// print Hertzdamp object properties
void 
Hertzdamp::Print(OPS_Stream &s, int flag)
{
    s << "Hertzdamp tag: " << this->getTag() << endln;
    s << "  Kh: " << Kh << endln;
    s << "  xiNorm: " << xiNorm << endln;
	s << "  gap: " << gap << endln;
    s << "  n: " << n << endln;
}
