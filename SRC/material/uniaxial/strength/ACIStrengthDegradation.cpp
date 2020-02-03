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
// ACIStrengthDegradation, which models hysteretic shear strength
// degradation as a function of beam section curvature ductility.

#include <ACIStrengthDegradation.h>
#include <G3Globals.h>
#include <Vector.h>
#include <Channel.h>

#include <math.h>
#include <float.h>

#include <elementAPI.h>

void *
OPS_ACIStrengthDegradation(void)
{
  StrengthDegradation *theDegradation = 0;

  if (OPS_GetNumRemainingInputArgs() < 3) {
    opserr << "Invalid number of args, want: strengthDegradation ACI tag? Ky? D1? v2? D2?" << endln;
    return 0;
  }

  int iData[1];
  double dData[4];
  
  int numData = 1;
  if (OPS_GetIntInput(&numData, iData) != 0) {
    opserr << "WARNING invalid tag for strengthDegradation ACI" << endln;
    return 0;
  }

  numData = 4;
  if (OPS_GetDoubleInput(&numData, dData) != 0) {
    opserr << "WARNING invalid data for strengthDegradation ACI" << endln;
    return 0;
  }

  theDegradation = new ACIStrengthDegradation(iData[0], dData[0], dData[1], dData[2], dData[3]);
  if (theDegradation == 0) {
    opserr << "WARNING could not create ACIStrengthDegradation\n";
    return 0;
  }

  return theDegradation;
}

ACIStrengthDegradation::ACIStrengthDegradation
(int tag, double Ky, double D1, double v2, double D2):
  StrengthDegradation(tag,DEG_TAG_STRENGTH_ACI),
  V2(v2), d1(D1), d2(D2)
{
  // Store Ky here and invert later
  oneOverKy = fabs(Ky);
  
  if (oneOverKy < DBL_EPSILON)
    opserr << "ACIStrengthDegradation::ACIStrengthDegradation -- yield curvature is zero" << endln;
  
  // Invert Ky
  oneOverKy = 1.0/oneOverKy;
  
  if (d2 <= d1)
    opserr << "ACIStrengthDegradation::ACIStrengthDegradation -- d2 is <= d1" << endln;
  
  this->revertToStart();
  this->revertToLastCommit();
}

ACIStrengthDegradation::ACIStrengthDegradation():
  StrengthDegradation(0,DEG_TAG_STRENGTH_ACI),
  V2(0.0), d1(0.0), d2(0.0), oneOverKy(0.0)
{

}
 
ACIStrengthDegradation::~ACIStrengthDegradation()
{

}

const char*
ACIStrengthDegradation::getMeasure(void)
{
  return "curvatureZ";
}

int
ACIStrengthDegradation::setTrialMeasure(double measure)
{
  Tductility = fabs(measure)*oneOverKy;
  
  if (Tductility < Cductility)
    Tductility = Cductility;
  
  return 0;
}

double
ACIStrengthDegradation::getValue(void)
{
  if (Tductility < 1.0)
    return 1.0;
  else if (Tductility > d2)
    return V2;
  else {
    double m = (V2-1.0)/(d2-d1);
    return 1.0 + m*(Tductility-d1);
  }
}

int
ACIStrengthDegradation::commitState(void)
{
  Cductility = Tductility;
  
  return 0;
}
 
int
ACIStrengthDegradation::revertToLastCommit(void)
{
  Tductility = Cductility;

  return 0;
}

int
ACIStrengthDegradation::revertToStart(void)
{
  Cductility = 0.0;
  
  return 0;
}

StrengthDegradation*
ACIStrengthDegradation::getCopy(void)
{
  ACIStrengthDegradation *theCopy =
    new ACIStrengthDegradation (this->getTag(), 1/oneOverKy, d1, V2, d2);
  
  theCopy->Cductility = Cductility;
  
  return theCopy;
}

int
ACIStrengthDegradation::sendSelf(int commitTag, Channel &theChannel)
{
  static Vector data(6);
  
  data(0) = this->getTag();
  data(1) = oneOverKy;
  data(2) = V2;
  data(3) = d1;
  data(4) = d2;
  data(5) = Cductility;
  
  int res = theChannel.sendVector(this->getDbTag(), commitTag, data);
  
  if (res < 0) 
    opserr << "ACIStrengthDegradation::sendSelf() - failed to send data\n";
  
  return res;
}

int
ACIStrengthDegradation::recvSelf(int commitTag, Channel &theChannel, 
				 FEM_ObjectBroker &theBroker)
{
  static Vector data(6);
  int res = theChannel.recvVector(this->getDbTag(), commitTag, data);
  
  if (res < 0) {
    opserr << "ACIStrengthDegradation::recvSelf() - failed to receive data\n";
    this->setTag(0);      
  }
  else {
    this->setTag(int(data(0)));
    oneOverKy = data(1);
    V2 = data(2);
    d1 = data(3);
    d2 = data(4);
    Cductility = data(5);
  }
  
  return res;
}

void
ACIStrengthDegradation::Print(OPS_Stream &s, int flag)
{
  s << "ACIStrengthDegradation, tag: " << this->getTag() << endln;
  s << "ky: " << 1/oneOverKy << endln;
  s << "V2: " << V2 << endln;
  s << "d1: " << d1 << endln;
  s << "d2: " << d2 << endln;
}
