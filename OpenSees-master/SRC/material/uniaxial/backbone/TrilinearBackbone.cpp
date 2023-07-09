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
                                                                        
// $Revision: 1.2 $
// $Date: 2008-11-09 06:05:48 $
// $Source: /usr/local/cvs/OpenSees/SRC/material/uniaxial/backbone/TrilinearBackbone.cpp,v $

// Written: MHS
// Created: Aug 2000
//
// Description: This file contains the implementation of 
// TrilinearBackbone, which is a trilinear backbone.

#include <TrilinearBackbone.h>
#include <Vector.h>
#include <Channel.h>
#include <string.h>

#include <elementAPI.h>

void *
OPS_TrilinearBackbone(void)
{
  HystereticBackbone *theBackbone = 0;

  if (OPS_GetNumRemainingInputArgs() < 7) {
    opserr << "Invalid number of args, want: hystereticBackbone Trilinear tag? e1? s1? e2? s2? e3? s3?" << endln;
    return 0;
  }

  int iData[1];
  double dData[6];
  
  int numData = 1;
  if (OPS_GetIntInput(&numData, iData) != 0) {
    opserr << "WARNING invalid tag for hystereticBackbone Trilinear" << endln;
    return 0;
  }

  numData = 6;
  if (OPS_GetDoubleInput(&numData, dData) != 0) {
    opserr << "WARNING invalid data for hystereticBackbone Trilinear" << endln;
    return 0;
  }

  theBackbone = new TrilinearBackbone(iData[0], dData[0], dData[1], dData[2], dData[3], dData[4], dData[5]);
  if (theBackbone == 0) {
    opserr << "WARNING could not create TrilinearBackbone\n";
    return 0;
  }

  return theBackbone;
}

void *
OPS_BilinearBackbone(void)
{
  HystereticBackbone *theBackbone = 0;

  if (OPS_GetNumRemainingInputArgs() < 5) {
    opserr << "Invalid number of args, want: hystereticBackbone Bilinear tag? e1? s1? e2? s2?" << endln;
    return 0;
  }

  int iData[1];
  double dData[4];
  
  int numData = 1;
  if (OPS_GetIntInput(&numData, iData) != 0) {
    opserr << "WARNING invalid tag for hystereticBackbone Bilinear" << endln;
    return 0;
  }

  numData = 4;
  if (OPS_GetDoubleInput(&numData, dData) != 0) {
    opserr << "WARNING invalid data for hystereticBackbone Bilinear" << endln;
    return 0;
  }

  theBackbone = new TrilinearBackbone(iData[0], dData[0], dData[1], dData[2], dData[3]);
  if (theBackbone == 0) {
    opserr << "WARNING could not create BilinearBackbone\n";
    return 0;
  }

  return theBackbone;
}

TrilinearBackbone::TrilinearBackbone (int tag, double x1, double y1, 
				      double x2, double y2, double x3, double y3):
  HystereticBackbone(tag,BACKBONE_TAG_Trilinear), e1(x1), e2(x2), e3(x3),
  s1(y1), s2(y2), s3(y3)
{
  bool error = false;
  
  if (e1 <= 0.0)
    error = true;
  
  if (e2 <= e1)
    error = true;
  
  if (e3 <= e2)
    error = true;
  
  if (error)
    opserr << "TrilinearBackbone::TrilinearBackbone -- input backbone is not unique (one-to-one)" << endln;
  
  E1 = s1/e1;
  E2 = (s2-s1)/(e2-e1);
  E3 = (s3-s2)/(e3-e2);
}

TrilinearBackbone::TrilinearBackbone (int tag, double x1, double y1, 
				      double x2, double y2):
  HystereticBackbone(tag,BACKBONE_TAG_Trilinear), e1(x1), e3(x2), s1(y1), s3(y2)
{
  bool error = false;
  
  if (e1 <= 0.0)
    error = true;
  
  if (e3 <= e1)
    error = true;
  
  if (error)
    opserr << "TrilinearBackbone::TrilinearBackbone -- input backbone is not unique (one-to-one)" << endln;
  
  e2 = 0.5*(e1+e3);
  s2 = 0.5*(s1+s3);
  
  E1 = s1/e1;
  E2 = (s3-s1)/(e3-e1);
  E3 = E2;
}

TrilinearBackbone::TrilinearBackbone():
  HystereticBackbone(0,BACKBONE_TAG_Trilinear), 
  E1(0.0), E2(0.0), E3(0.0),
  e1(0.0), e2(0.0), e3(0.0), 
  s1(0.0), s2(0.0), s3(0.0)
{

}

TrilinearBackbone::~TrilinearBackbone()
{
  
}

double
TrilinearBackbone::getTangent (double strain)
{
  if (strain < e1)
    return E1;
  else if (strain < e2)
    return E2;
  else if (strain < e3)
    return E3;
  else
    return E1*1.0e-9;
}

double
TrilinearBackbone::getStress (double strain)
{
  if (strain < e1)
    return E1*strain;
  else if (strain < e2)
    return s1 + E2*(strain-e1);
  else if (strain < e3)
    return s2 + E3*(strain-e2);
  else
    return s3;
}

double
TrilinearBackbone::getEnergy (double strain)
{
  if (strain < e1)
    return 0.5*strain*strain*E1;
  else if (strain < e2)
    return 0.5*e1*e1*E1 + 0.5*(getStress(strain)+s1)*(strain-e1);
  else if (strain < e3)
    return 0.5*e1*e1*E1 + 0.5*(s1+s2)*(e2-e1) + 0.5*(getStress(strain)+s2)*(strain-e2);
  else
    return 0.5*e1*e1*E1 + 0.5*(s1+s2)*(e2-e1) + 0.5*(s2+s3)*(e3-e2) + (getStress(strain)-s3)*(strain-e3);
}

double
TrilinearBackbone::getYieldStrain(void)
{
  return e1;
}

HystereticBackbone*
TrilinearBackbone::getCopy(void)
{
  TrilinearBackbone *theCopy = 
    new TrilinearBackbone (this->getTag(), e1, s1, e2, s2, e3, s3);
  
  theCopy->E1 = E1;
  theCopy->E2 = E2;
  theCopy->E3 = E3;
  
  return theCopy;
}

void
TrilinearBackbone::Print(OPS_Stream &s, int flag)
{
  s << "TrilinearBackbone, tag: " << this->getTag() << endln;
  s << "\te1: " << e1 << ", s1: " << s1 << endln;
  s << "\te2: " << e2 << ", s2: " << s2 << endln;
  s << "\te3: " << e3 << ", s3: " << s3 << endln;
  s << "\t\tE1: " << E1 << endln;
  s << "\t\tE2: " << E2 << endln;
  s << "\t\tE3: " << E3 << endln;
}

int
TrilinearBackbone::setVariable (char *argv)
{
  if (strcmp(argv,"yieldStrain") == 0)
    return 1;
  else
    return -1;
}

int
TrilinearBackbone::getVariable (int varID, double &theValue)
{
  switch (varID) {
  case 1:
    theValue = e1;
    return 1;
  default:
    return -1;
  }
}

int
TrilinearBackbone::sendSelf(int commitTag, Channel &theChannel)
{
  int res = 0;
  
  static Vector data(10);
  
  data(0) = this->getTag();
  data(1) = e1;
  data(2) = e2;
  data(3) = e3;
  data(4) = s1;
  data(5) = s2;
  data(6) = s3;
  data(7) = E1;
  data(8) = E2;
  data(9) = E3;
  
  res += theChannel.sendVector(this->getDbTag(), commitTag, data);
  if (res < 0) {
    opserr << "TrilinearBackbone::sendSelf -- could not send Vector" << endln;
    return res;
  }
  
  return res;
}

int
TrilinearBackbone::recvSelf(int commitTag, Channel &theChannel, 
			    FEM_ObjectBroker &theBroker)
{
  int res = 0;
  
  static Vector data(10);
  
  res = theChannel.recvVector(this->getDbTag(), commitTag, data);
  if (res < 0) {
    opserr << "TrilinearBackbone::recvSelf -- could not receive Vector" << endln;
    return res;
  }
  
  this->setTag(int(data(0)));
  e1 = data(1);
  e2 = data(2);
  e3 = data(3);
  s1 = data(4);
  s2 = data(5);
  s3 = data(6);
  E1 = data(7);
  E2 = data(8);
  E3 = data(9);
  
  return res;
}
