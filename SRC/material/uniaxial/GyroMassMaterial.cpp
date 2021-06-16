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

// $Revision:
// $Date: 2021/6/16 $
// $Source: C:\Users\Edbert Lumbantobing\Documents\GitHub\OpenSees\SRC\GyroMassMaterial.cpp $

// Written: Edbert Rainer Lumbantobing (Tokyo Institute of Technology)
// Created: June 2021
// Revision:
//
// Description: This file contains the class implementation for GyroMassMaterial that can be used to model an inerter.
// Inerter is a mechanical element which output force is proportional to the relative acceleration between its terminals (between two nodes). 
// F = b (a1 - a2)
// b = inertance
// a1 = acceleration at node 1
// a2 = acceleration at node 2

// This code has been partly published by Hessabi (2017)
// Reference: Hessabi R. (2017). Application of real-time hybrid simulation method in experimental identification of gyromass dampers (Doctoral thesis). University of Toronto, Toronto.


#include <elementAPI.h>
#include "GyroMassMaterial.h"

#include <Vector.h>
#include <Channel.h>
#include <math.h>
#include <float.h>
#include <Domain.h>

static int numGyroMassMaterial = 0;

void *
OPS_GyroMassMaterial()
{
  // print out some KUDO's
  if (numGyroMassMaterial == 0) {
    opserr << "GyroMassMaterial uniaxial material - Written by Edbert Lumbantobing at Tokyo Tech - Use at your own risk\n";
    numGyroMassMaterial =1;
  }

  // Pointer to a uniaxial material that will be returned
  UniaxialMaterial *theMaterial = 0;

  //
  // parse the input line for the material parameters
  //

  int    iData[1];
  double dData[1];
  int numData;
  numData = 1;
  if (OPS_GetIntInput(&numData, iData) != 0) {
    opserr << "WARNING invalid uniaxialMaterial GyroMassMaterial tag" << endln;
    return 0;
  }

  numData = 1;
  if (OPS_GetDoubleInput(&numData, dData) != 0) {
    opserr << "WARNING invalid b\n";
    return 0;	
  }

  // 
  // create a new material
  //

  theMaterial = new GyroMassMaterial(iData[0], dData[0]);       

  if (theMaterial == 0) {
    opserr << "WARNING could not create uniaxialMaterial of type GyroMassMaterial\n";
    return 0;
  }

  // return the material
  return theMaterial;
}


GyroMassMaterial::GyroMassMaterial(int tag, double b)
:UniaxialMaterial(tag, MAT_TAG_GyroMassMaterial),
 trialStrain(0.0), trialStress(0.0),
 commitStrain(0.0), commitStress(0.0),
 trialStrainRate(0.0), trialAccel(0.0), b(b), currentTime(0.0), commitStrainRate(0.0), commitAccel(0.0), commitTime(0.0)
{

}

GyroMassMaterial::GyroMassMaterial()
:UniaxialMaterial(0, MAT_TAG_GyroMassMaterial),
 trialStrain(0.0), trialStress(0.0),
 commitStrain(0.0), commitStress(0.0),
 trialStrainRate(0.0), trialAccel(0.0), b(0.0), currentTime(0.0), commitStrainRate(0.0), commitAccel(0.0), commitTime(0.0)
{

}

GyroMassMaterial::~GyroMassMaterial()
{
  // does nothing
}

int 
GyroMassMaterial::setTrialStrain(double strain, double strainRate)
{
    double dt;
    currentTime = ops_TheActiveDomain->getCurrentTime();
    dt = currentTime - commitTime;
    trialStrain = strain;
    newmark(0.5,0.25,dt,trialStrain,commitStrain,trialStrainRate,commitStrainRate,trialAccel,commitAccel);
    trialStress = b * trialAccel;

    return 0;
}

int
GyroMassMaterial::newmark(double gama, double beta, double dt, double trialStrain, double commitStrain, double& trialStrainRate, double commitStrainRate, double& trialAccel, double commitAccel)
{
    double dDis, dVel, dAccel = 0;
    dDis = trialStrain - commitStrain;
    dVel = gama * dDis / (beta * dt) - gama * commitStrainRate / beta + commitAccel * dt * (1 - gama / (2 * beta));
    dAccel = dDis / (beta * pow(dt, 2)) - commitStrainRate / (beta * dt) - commitAccel / (2 * beta);
    trialStrainRate = commitStrainRate + dVel;
    trialAccel = commitAccel + dAccel;
    
    return 0;
}

double 
GyroMassMaterial::getStrain(void)
{
  return trialStrain;
}

double 
GyroMassMaterial::getStress(void)
{ 
  return trialStress;
}


double 
GyroMassMaterial::getTangent(void)
{
  return 0.0;
}

double
GyroMassMaterial::getInitialTangent(void)
{
    return 0.0;
}

double
GyroMassMaterial::getDampTangent(void)
{
    double dt;
    currentTime = ops_TheActiveDomain->getCurrentTime();
    dt = currentTime - commitTime;
    return 2 * b / dt;
}

int 
GyroMassMaterial::commitState(void)
{
    commitStrain = trialStrain;
    commitStress = trialStress;
    commitStrainRate = trialStrainRate;
    commitAccel = trialAccel;
    commitTime = currentTime;

    return 0;
}	


int 
GyroMassMaterial::revertToLastCommit(void)
{
  trialStrain = commitStrain;
  trialStress = commitStress;
  trialStrainRate = commitStrainRate;
  trialAccel = commitAccel;
  currentTime = commitTime;

  return 0;
}


int 
GyroMassMaterial::revertToStart(void)
{
  trialStrain = commitStrain = 0.0;
  trialStress = commitStress = 0.0;
  trialStrainRate = commitStrainRate = 0;
  trialAccel = commitAccel = 0;
  currentTime = commitTime = 0;

  return 0;
}


UniaxialMaterial *
GyroMassMaterial::getCopy(void)
{
  GyroMassMaterial *theCopy =
  new GyroMassMaterial(this->getTag(),b);
  theCopy->trialStrain = this->trialStrain;
  theCopy->trialStress = this->trialStress;
  theCopy->trialStrainRate = this->trialStrainRate;
  theCopy->trialAccel = this->trialAccel;
  theCopy->currentTime = this->currentTime;
  
  return theCopy;
}


int 
GyroMassMaterial::sendSelf(int cTag, Channel &theChannel)
{
  int res = 0;
  static Vector data(7);
  data(0) = this->getTag();
  data(1) = commitStrain;
  data(2) = commitStress;
  data(3) = b;
  data(4) = commitTime;
  data(5) = commitStrainRate;
  data(6) = commitAccel;

  res = theChannel.sendVector(this->getDbTag(), cTag, data);
  if (res < 0) 
    opserr << "GyroMassMaterial::sendSelf() - failed to send data\n";

  return res;
}

int 
GyroMassMaterial::recvSelf(int cTag, Channel &theChannel, 
				 FEM_ObjectBroker &theBroker)
{
  int res = 0;
  static Vector data(7);
  res = theChannel.recvVector(this->getDbTag(), cTag, data);
  if (res < 0) 
    opserr << "GyroMassMaterial::recvSelf() - failed to recv data\n";
  else {
    this->setTag(data(0));
    commitStrain = data(1);
    commitStress = data(2);
    b = data(3);
    commitTime = data(4);
    commitStrainRate = data(5);
    commitAccel = data(6);
    trialStrain = commitStrain;
    trialStress = commitStress;
    trialStrainRate = commitStrainRate;
    trialAccel = commitAccel;
    currentTime = commitTime;
  }

  return res;
}

void 
GyroMassMaterial::Print(OPS_Stream &s, int flag)
{
  s << "GyroMassMaterial tag: " << this->getTag() << endln;
  s << "  b: " << b << endln;
  s << "  Inertial force: " << trialStress << endln;
}


