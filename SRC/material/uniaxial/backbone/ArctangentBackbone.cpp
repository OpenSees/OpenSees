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
// $Date: 2008-11-09 06:03:59 $
// $Source: /usr/local/cvs/OpenSees/SRC/material/uniaxial/backbone/ArctangentBackbone.cpp,v $

// Written: MHS
// Created: Aug 2000
//
// Description: This file contains the implementation of 
// ArctangentBackbone, which is a continuous function given
// by K1*atan(K2*strain); as developed by Ranzo and Petrangeli (1998)

#include <ArctangentBackbone.h>
#include <Vector.h>
#include <Channel.h>

#include <math.h>

#include <elementAPI.h>

void * OPS_ADD_RUNTIME_VPV(OPS_ArctangentBackbone)
{
  HystereticBackbone *theBackbone = 0;

  if (OPS_GetNumRemainingInputArgs() < 4) {
    opserr << "Invalid number of args, want: hystereticBackbone Arctangent tag? K1? gamma? alpha?" << endln;
    return 0;
  }

  int iData[1];
  double dData[3];
  
  int numData = 1;
  if (OPS_GetIntInput(&numData, iData) != 0) {
    opserr << "WARNING invalid tag for hystereticBackbone Arctangent" << endln;
    return 0;
  }

  numData = 3;
  if (OPS_GetDoubleInput(&numData, dData) != 0) {
    opserr << "WARNING invalid data for hystereticBackbone Arctangent" << endln;
    return 0;
  }

  theBackbone = new ArctangentBackbone(iData[0], dData[0], dData[1], dData[2]);
  if (theBackbone == 0) {
    opserr << "WARNING could not create ArctangentBackbone\n";
    return 0;
  }

  return theBackbone;
}

ArctangentBackbone::ArctangentBackbone(int tag, double k1, double gy, double a):
  HystereticBackbone(tag,BACKBONE_TAG_Arctangent),
  K1(k1), K2(0.0), gammaY(gy), alpha(a)
{
  if (gammaY == 0.0)
    opserr << "ArctangentBackbone::ArctangentBackbone -- gammaY is zero" << endln;
  
  gammaY = fabs(gammaY);
  alpha = fabs(alpha);
  
  K2 = tan(alpha)/gammaY;
}

ArctangentBackbone::ArctangentBackbone():
  HystereticBackbone(0,BACKBONE_TAG_Arctangent),
  K1(0.0), K2(0.0), gammaY(0.0), alpha(0.0)
{
  
}

ArctangentBackbone::~ArctangentBackbone()
{
  
}

double
ArctangentBackbone::getTangent (double strain)
{
  return K1/(1+pow(K2*strain,2))*K2;
}

double
ArctangentBackbone::getStress (double strain)
{
  return K1*atan(K2*strain);
}

double
ArctangentBackbone::getEnergy (double strain)
{
  double u = K2*strain;
  
  return K1/K2*(u*atan(u) - 0.5*log(1+u*u));
}

double
ArctangentBackbone::getYieldStrain(void)
{
  return gammaY;
}

HystereticBackbone*
ArctangentBackbone::getCopy(void)
{
  ArctangentBackbone *theCopy =
    new ArctangentBackbone (this->getTag(), K1, gammaY, alpha);
  
  return theCopy;
}

void
ArctangentBackbone::Print(OPS_Stream &s, int flag)
{
  s << "ArctangentBackbone, tag: " << this->getTag() << endln;
  s << "\tK1: " << K1 << endln;
  s << "\tgammaY: " << gammaY << endln;
  s << "\talpha: " << alpha << endln;
}

int
ArctangentBackbone::setVariable (char *argv)
{
  return -1;
}

int
ArctangentBackbone::getVariable (int varID, double &theValue)
{
  return -1;
}

int
ArctangentBackbone::sendSelf(int commitTag, Channel &theChannel)
{
  int res = 0;
  
  static Vector data(5);
  
  data(0) = this->getTag();
  data(1) = K1;
  data(2) = gammaY;
  data(3) = alpha;
  data(4) = K2;
  
  res += theChannel.sendVector(this->getDbTag(), commitTag, data);
  if (res < 0) {
    opserr << "ArctangentBackbone::sendSelf -- could not send Vector" << endln;

    return res;
  }
  
  return res;
}

int
ArctangentBackbone::recvSelf(int commitTag, Channel &theChannel, 
			     FEM_ObjectBroker &theBroker)
{
  int res = 0;
  
  static Vector data(5);
  
  res += theChannel.recvVector(this->getDbTag(), commitTag, data);
  if (res < 0) {
    opserr << "ArctangentBackbone::recvSelf -- could not receive Vector" << endln;

    return res;
  }
  
  this->setTag(int(data(0)));
  K1 = data(1);
  gammaY = data(2);
  alpha = data(3);
  K2 = data(4);
  
  return res;
}
