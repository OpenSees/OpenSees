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
// SectionStrengthDegradation, which models hysteretic shear strength
// degradation according to Ranzo and Section (1998).  Shear strength
// degradation is a function of beam section axial strain.

#include <SectionStrengthDegradation.h>
#include <SectionForceDeformation.h>
#include <G3Globals.h>
#include <Vector.h>
#include <Channel.h>
#include <float.h>
#include <math.h>

SectionStrengthDegradation::SectionStrengthDegradation
(int tag, double E1, double v2, double E2, int c):
  StrengthDegradation(tag,DEG_TAG_STRENGTH_Section),
  code(c), V2(v2), e1(E1), e2(E2)
{
  if (e2 <= e1)
    opserr << "SectionStrengthDegradation::SectionStrengthDegradation -- e2 is <= e1" << endln;
  
  slope = (V2-1.0)/(e2-e1);
  
  this->revertToStart();
  this->revertToLastCommit();
}

SectionStrengthDegradation::SectionStrengthDegradation
(int tag, double Ey, double u1, double v2, double u2, int c):
  StrengthDegradation(tag,DEG_TAG_STRENGTH_Section),
  code(c), V2(v2), e1(u1*Ey), e2(u2*Ey)
{
  if (e2 <= e1)
    opserr << "SectionStrengthDegradation::SectionStrengthDegradation -- e2 is <= e1" << endln;
  
  slope = (V2-1.0)/(e2-e1);
  
  this->revertToStart();
  this->revertToLastCommit();
}

SectionStrengthDegradation::SectionStrengthDegradation():
  StrengthDegradation(0,DEG_TAG_STRENGTH_Section),
  code(0), V2(0.0), e1(0.0), e2(0.0)
{
  
}

SectionStrengthDegradation::~SectionStrengthDegradation()
{
  
}

const char*
SectionStrengthDegradation::getMeasure(void)
{
  switch (code) {
  case SECTION_RESPONSE_P:
    return "axialStrain";
  case SECTION_RESPONSE_MZ:
    return "curvatureZ";
  case SECTION_RESPONSE_VY:
    return "shearStrainY";
  case SECTION_RESPONSE_MY:
    return "curvatureY";
  case SECTION_RESPONSE_VZ:
    return "shearStrainZ";
  case SECTION_RESPONSE_T:
    return "rotationX";
  default:
    return "nothing";
  }
}

int
SectionStrengthDegradation::setTrialMeasure(double measure)
{
  Tmeasure = fabs(measure);
  
  if (Tmeasure < Cmeasure)
    Tmeasure = Cmeasure;
  
  return 0;
}

double
SectionStrengthDegradation::getValue(void)
{
  if (Tmeasure < e1)
    return 1.0;
  else if (Tmeasure > e2)
    return V2;
  else {
    return 1.0 + slope*(Tmeasure-e1);
  }
}

int
SectionStrengthDegradation::commitState(void)
{
  Cmeasure = Tmeasure;
  
  return 0;
}
 
int
SectionStrengthDegradation::revertToLastCommit(void)
{
  Tmeasure = Cmeasure;
  
  return 0;
}

int
SectionStrengthDegradation::revertToStart(void)
{
  Cmeasure = 0.0;
  
  return 0;
}

StrengthDegradation*
SectionStrengthDegradation::getCopy(void)
{
  SectionStrengthDegradation *theCopy =
    new SectionStrengthDegradation(this->getTag(), e1, V2, e2, code);
  
  theCopy->Cmeasure = Cmeasure;
  
  return theCopy;
}

int
SectionStrengthDegradation::sendSelf(int commitTag, Channel &theChannel)
{
  static Vector data(6);
  
  data(0) = this->getTag();
  data(1) = e1;
  data(2) = V2;
  data(3) = e2;
  data(4) = code;
  data(5) = Cmeasure;
  
  int res = theChannel.sendVector(this->getDbTag(), commitTag, data);
  
  if (res < 0) 
    opserr << "SectionStrengthDegradation::sendSelf() - failed to send data\n";
  
  return res;
}

int
SectionStrengthDegradation::recvSelf(int commitTag, Channel &theChannel, 
				     FEM_ObjectBroker &theBroker)
{
  static Vector data(6);
  int res = theChannel.recvVector(this->getDbTag(), commitTag, data);
  
  if (res < 0) {
    opserr << "SectionStrengthDegradation::recvSelf() - failed to receive data\n";
    this->setTag(0);      
  }
  else {
    this->setTag(int(data(0)));
    e1 = data(1);
    V2 = data(2);
    e2 = data(3);
    slope = (V2-1.0)/(e2-e1);
    code = int(data(4));
    Cmeasure = data(5);
  }
  
  return res;
}

void
SectionStrengthDegradation::Print(OPS_Stream &s, int flag)
{
  s << "SectionStrengthDegradation, tag: " << this->getTag() << endln;
  s << "\tcode: " << code << endln;
  s << "\te1: " << e1 << endln;
  s << "\tV2: " << V2 << endln;
  s << "\te2: " << e2 << endln;
}
