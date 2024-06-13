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
// Source: /usr/local/cvs/OpenSees/SRC/material/uniaxial/JankowskiImpact.cpp
//
// Written: Patrick J. Hughes, University of California - San Diego
// Created: 05/2019
//
// Description: This file contains the class implementation for the
// JankowskiImpact uniaxialMaterial.
//
// References:
//
// Variables:
// Kh: nonlinear Hertz contact stiffness
// xi: impact damping ratio
// mEff: effective mass of the colliding bodies
// gap: initial gap distance (must be a negative value)
// n: displacement exponent (default is 1.5)

#include <stdlib.h>
#include <Vector.h>
#include <Channel.h>
#include <math.h>
#include <float.h>

#include <elementAPI.h>
#include <OPS_Globals.h>
#include <JankowskiImpact.h>

static int numJankowskiImpact = 0;

void *
OPS_JankowskiImpact(void)
{
  // kudos
  if (numJankowskiImpact == 0) {
    numJankowskiImpact++;
    opserr << "JankowskiImpact model written by Patrick J. Hughes, UC San Diego\n";
  }
  
  // Pointer to a uniaxial material that will be returned
  UniaxialMaterial *theMaterial = 0;

  // check for minimum number of arguments
  if (OPS_GetNumRemainingInputArgs() < 5) {
    opserr << "WARNING incorrect num args want: uniaxialMaterial JankowskiImpact tag? Kh? xi? mEff? gap? <n?>" << endln;
    return 0;
  }

  int    iData[1]; // 1 integer input (tag)
  double dData[5]; // 5 double inputs (Kh, xi, mEff, gap, n)
  
  // check material tag validity
  int numData = 1;
  if (OPS_GetIntInput(&numData, iData) != 0) {
    opserr << "WARNING invalid uniaxialMaterial JankowskiImpact tag" << endln;
    return 0;
  }

  // check material property validity
  numData = OPS_GetNumRemainingInputArgs();
  // all 5 double inputs given
  if (numData >= 5) {
	  numData = 5;
	  if (OPS_GetDoubleInput(&numData, dData) != 0) {
		  opserr << "Invalid data for uniaxialMaterial JankowskiImpact " << iData[0] << endln;
		  return 0;
	  }
  } else {
	  // 4 double inputs given - use default value for 4th input
	  numData = 4;
	  if (OPS_GetDoubleInput(&numData, dData) != 0) {
		  opserr << "Invalid data for uniaxialMaterial JankowskiImpact " << iData[0] << endln;
		  return 0;
	  }
	  dData[4] = 1.5; // default value for displacement exponent is 1.5
  }
	 
  // parsing was successful - allocate material
  theMaterial = new JankowskiImpact(iData[0], dData[0], dData[1], dData[2], dData[3], dData[4]);
  if (theMaterial == 0) {
    opserr << "WARNING could not create uniaxialMaterial JankowskiImpact\n";
    return 0;
  }

  return theMaterial;
}

// constructor
JankowskiImpact::JankowskiImpact(int tag, double kh, double xiin, double meff, double gap0, double dispexp)
:UniaxialMaterial(tag,MAT_TAG_JankowskiImpact),
 Kh(kh), xi(xiin), mEff(meff), gap(gap0), n(dispexp)
{
	if (gap >= 0) {
	  opserr << "JankowskiImpact::JankowskiImpact -- Initial gap size should be negative for compression-only material\n";
	  gap = -gap;
	  opserr << "Setting gap to negative value, " << gap << endln;
	}
	this->revertToStart(); // initialize state variables
	printFlag = 0; // initialize print flag for impact event
}

// constructor for parallel processing
JankowskiImpact::JankowskiImpact()
:UniaxialMaterial(0,MAT_TAG_JankowskiImpact),
 Kh(0.0), xi(0.0), mEff(0.0), gap(0.0), n(0.0)
{

}

// destructor
JankowskiImpact::~JankowskiImpact()
{
  // does nothing
}

// state determintaion
int 
JankowskiImpact::setTrialStrain(double strain, double strainRate)
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
				opserr << "JankowskiImpact impact detected: impact velocity = " << -trialStrainRate << "\n"; // display impact message
				printFlag = 1; // update print flag
			}
		}
		
		// same impact event as previous state
		if (trialStrainRate < 0) {
			trialStress = -( Kh*pow(-(trialStrain-gap),n) + 2*xi*pow(Kh*pow(-(trialStrain-gap),n-1.0)*mEff,0.5)*(-trialStrainRate) ); // contact force - approach period
			trialTangent = -( n*Kh*pow(-(trialStrain-gap),n-1.0) + (n-1.0)*xi*pow(Kh*pow(-(trialStrain-gap),n-3.0)*mEff,0.5)*(-trialStrainRate) ); // contact stiffness - approach period
		} else {
			trialStress = -Kh*pow(-(trialStrain-gap),n); // contact force - restitution period
			trialTangent = -n*Kh*pow(-(trialStrain-gap),n-1.0); // contact stiffness - restitution period
		}
	}
    return 0;
}

// retrieve trial strain
double 
JankowskiImpact::getStrain(void)
{
  return trialStrain;
}

// retrieve trial strain rate
double 
JankowskiImpact::getStrainRate(void)
{
  return trialStrainRate;
}

// retrieve trial stress
double 
JankowskiImpact::getStress(void)
{
  return trialStress;
}

// retrieve trial tangent stiffness
double 
JankowskiImpact::getTangent(void)
{
  return trialTangent;
}

// retrieve initial tangent stiffness
double 
JankowskiImpact::getInitialTangent(void)
{
  return 0.0; 
}

// commit state
int 
JankowskiImpact::commitState(void)
{
	commitStrain = trialStrain;
	commitStrainRate = trialStrainRate;
	commitStress = trialStress;
	commitTangent = trialTangent;
	return 0;
}

// revert to last comitted state
int 
JankowskiImpact::revertToLastCommit(void)
{
	trialStrain = commitStrain;
	trialStrainRate = commitStrainRate;
	trialStress = commitStress;
	trialTangent = commitTangent;
    return 0;
}

// revert to initial state
int 
JankowskiImpact::revertToStart(void)
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
JankowskiImpact::getCopy(void)
{
    JankowskiImpact *theCopy = new JankowskiImpact(this->getTag(),Kh,xi,mEff,gap,n);
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

// send JankowskiImpact object
int 
JankowskiImpact::sendSelf(int cTag, Channel &theChannel)
{
  int res = 0;
  static Vector data(14);
  data(0)  = this->getTag();
  data(1)  = Kh;
  data(2)  = xi;
  data(3)  = mEff;
  data(4)  = gap;
  data(5)  = n;
  data(6)  = commitStrain;
  data(7)  = commitStrainRate;
  data(8)  = commitStress; 
  data(9)  = commitTangent;
  data(10)  = trialStrain;
  data(11) = trialStrainRate;
  data(12) = trialStress; 
  data(13) = trialTangent;
  res = theChannel.sendVector(this->getDbTag(), cTag, data);
  if (res < 0) 
    opserr << "JankowskiImpact::sendSelf() - failed to send data\n";
  return res;
}

// receive JankowskiImpact object
int 
JankowskiImpact::recvSelf(int cTag, Channel &theChannel,
			  FEM_ObjectBroker &theBroker)
{
  int res = 0;
  static Vector data(14);
  res = theChannel.recvVector(this->getDbTag(), cTag, data);
  if (res < 0)
    opserr << "JankowskiImpact::recvSelf() - failed to recv data\n";
  else {
    this->setTag((int)data(0));
    Kh = data(1);
    xi = data(2);
    mEff = data(3);
    gap = data(4);
    n = data(5);
    commitStrain = data(6);
    commitStrainRate = data(7);
    commitStress = data(8); 
    commitTangent = data(9);
    trialStrain = data(10);
    trialStrainRate = data(11);
    trialStress = data(12); 
    trialTangent = data(13);
  }
  return res;
}

// print JankowskiImpact object properties
void 
JankowskiImpact::Print(OPS_Stream &s, int flag)
{
    s << "JankowskiImpact tag: " << this->getTag() << endln;
    s << "  Kh: " << Kh << endln;
    s << "  xi: " << xi << endln;
	s << "  mEff: " << mEff << endln;
	s << "  gap: " << gap << endln;
    s << "  n: " << n << endln;
}
