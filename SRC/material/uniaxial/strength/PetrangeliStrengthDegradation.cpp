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

// Written: MHS
// Created: July 2000
//
// Description: This file contains the implementation of 
// PetrangeliStrengthDegradation, which models hysteretic shear strength
// degradation according to Ranzo and Petrangeli (1998).  Shear strength
// degradation is a function of beam section axial strain.

#include <PetrangeliStrengthDegradation.h>
#include <G3Globals.h>
#include <Vector.h>
#include <Channel.h>

#include <elementAPI.h>

void *
OPS_PetrangeliStrengthDegradation(void)
{
  StrengthDegradation *theDegradation = 0;

  if (OPS_GetNumRemainingInputArgs() < 3) {
    opserr << "Invalid number of args, want: strengthDegradation Petrangeli tag? e1? V2? e2?" << endln;
    return 0;
  }

  int iData[1];
  double dData[3];
  
  int numData = 1;
  if (OPS_GetIntInput(&numData, iData) != 0) {
    opserr << "WARNING invalid tag for strengthDegradation Petrangeli" << endln;
    return 0;
  }

  numData = 3;
  if (OPS_GetDoubleInput(&numData, dData) != 0) {
    opserr << "WARNING invalid data for strengthDegradation Petrangeli" << endln;
    return 0;
  }

  theDegradation = new PetrangeliStrengthDegradation(iData[0], dData[0], dData[1], dData[2]);
  if (theDegradation == 0) {
    opserr << "WARNING could not create PetrangeliStrengthDegradation\n";
    return 0;
  }

  return theDegradation;
}

PetrangeliStrengthDegradation::PetrangeliStrengthDegradation
(int tag, double E1, double v2, double E2):
  StrengthDegradation(tag,DEG_TAG_STRENGTH_Petrangeli),
  V2(v2), e1(E1), e2(E2)
{
  if (e2 <= e1)
    opserr << "PetrangeliStrengthDegradation::PetrangeliStrengthDegradation -- e2 is <= e1" << endln;
  
  this->revertToStart();
  this->revertToLastCommit();
}
 
PetrangeliStrengthDegradation::PetrangeliStrengthDegradation():
  StrengthDegradation(0,DEG_TAG_STRENGTH_Petrangeli),
  V2(0.0), e1(0.0), e2(0.0)
{

}

PetrangeliStrengthDegradation::~PetrangeliStrengthDegradation()
{

}

const char*
PetrangeliStrengthDegradation::getMeasure(void)
{
  return "axialStrain";
}

int
PetrangeliStrengthDegradation::setTrialMeasure(double measure)
{
  Tmeasure = measure;
  
  if (Tmeasure < Cmeasure)
    Tmeasure = Cmeasure;
  
  return 0;
}

double
PetrangeliStrengthDegradation::getValue(void)
{
  if (Tmeasure < e1)
    return 1.0;
  else if (Tmeasure > e2)
    return V2;
  else {
    double m = (V2-1.0)/(e2-e1);
    return 1.0 + m*(Tmeasure-e1);
  }
}

int
PetrangeliStrengthDegradation::commitState(void)
{
  Cmeasure = Tmeasure;
  
  return 0;
}
 
int
PetrangeliStrengthDegradation::revertToLastCommit(void)
{
  Tmeasure = Cmeasure;
  
  return 0;
}

int
PetrangeliStrengthDegradation::revertToStart(void)
{
  Cmeasure = 0.0;
  
  return 0;
}

StrengthDegradation*
PetrangeliStrengthDegradation::getCopy(void)
{
  PetrangeliStrengthDegradation *theCopy =
    new PetrangeliStrengthDegradation (this->getTag(), e1, V2, e2);
  
  theCopy->Cmeasure = Cmeasure;
  
  return theCopy;
}

int
PetrangeliStrengthDegradation::sendSelf(int commitTag, Channel &theChannel)
{
  static Vector data(5);
  
  data(0) = this->getTag();
  data(1) = e1;
  data(2) = V2;
  data(3) = e2;
  data(4) = Cmeasure;
  
  int res = theChannel.sendVector(this->getDbTag(), commitTag, data);
  
  if (res < 0) 
    opserr << "PetrangeliStrengthDegradation::sendSelf() - failed to send data\n";
  
  return res;
}

int
PetrangeliStrengthDegradation::recvSelf(int commitTag, Channel &theChannel, 
					FEM_ObjectBroker &theBroker)
{
  static Vector data(5);
  int res = theChannel.recvVector(this->getDbTag(), commitTag, data);
  
  if (res < 0) {
    opserr << "PetrangeliStrengthDegradation::recvSelf() - failed to receive data\n";
    this->setTag(0);      
  }
  else {
    this->setTag(int(data(0)));
    e1 = data(1);
    V2 = data(2);
    e2 = data(3);
    Cmeasure = data(4);
  }
  
  return res;
}

void
PetrangeliStrengthDegradation::Print(OPS_Stream &s, int flag)
{
  s << "PetrangeliStrengthDegradation, tag: " << this->getTag() << endln;
  s << "e1: " << e1 << endln;
  s << "V2: " << V2 << endln;
  s << "e2: " << e2 << endln;
}
