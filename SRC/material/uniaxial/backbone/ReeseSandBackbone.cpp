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
// $Date: 2008-11-10 22:40:58 $
// $Source: /usr/local/cvs/OpenSees/SRC/material/uniaxial/backbone/ReeseSandBackbone.cpp,v $

// Written: MHS
// Created: Aug 2000
//
// Description: This file contains the implementation of 
// ReeseSandBackbone

#include <ReeseSandBackbone.h>
#include <Vector.h>
#include <Channel.h>

#include <math.h>

ReeseSandBackbone::ReeseSandBackbone(int tag, double KX, 
				     double YM, double PM,
				     double YU, double PU):
  HystereticBackbone(tag,BACKBONE_TAG_ReeseSand),
  kx(KX), ym(YM), pm(PM), yu(YU), pu(PU)
{

}

ReeseSandBackbone::ReeseSandBackbone():
  HystereticBackbone(0,BACKBONE_TAG_ReeseSand),
  kx(0.0), ym(0.0), pm(0.0), yu(0.0), pu(0.0)
{
  
}

ReeseSandBackbone::~ReeseSandBackbone()
{
  
}

double
ReeseSandBackbone::getTangent (double strain)
{
  int signStrain = (strain > 0.0) ? 1 : -1;
  strain = signStrain*strain;

  double tangent = 0.0;

  double m = (pu-pm)/(yu-ym);
  double n = pm/(m*ym);
  double C = pm/pow(ym,1/n);

  double yk = pow(C/kx,n/(n-1));

  if (strain <= yk)
    tangent = kx;
  else if (strain <= ym)
    tangent = C/n*pow(strain,1/n-1);
  else if (strain <= yu)
    tangent = m;
  else
    tangent = 0.001*kx;

  return tangent;
}

double
ReeseSandBackbone::getStress (double strain)
{
  int signStrain = (strain > 0.0) ? 1 : -1;
  strain = signStrain*strain;

  double stress = 0.0;

  double m = (pu-pm)/(yu-ym);
  double n = pm/(m*ym);
  double C = pm/pow(ym,1/n);

  double yk = pow(C/kx,n/(n-1));

  if (strain <= yk)
    stress = kx*strain;
  else if (strain <= ym)
    stress = C*pow(strain,1/n);
  else if (strain <= yu)
    stress = pm + m*(strain-ym);
  else
    stress = pu;

  return signStrain*stress;
}

double
ReeseSandBackbone::getEnergy (double strain)
{
  return 0.0;
}

double
ReeseSandBackbone::getYieldStrain(void)
{
  return 0.0;
}

HystereticBackbone*
ReeseSandBackbone::getCopy(void)
{
  ReeseSandBackbone *theCopy =
    new ReeseSandBackbone (this->getTag(), kx, ym, pm, yu, pu);
  
  return theCopy;
}

void
ReeseSandBackbone::Print(OPS_Stream &s, int flag)
{
  s << "ReeseSandBackbone, tag: " << this->getTag() << endln;
  s << "\tkx: " << kx << endln;
  s << "\tym: " << ym << endln;
  s << "\tpm: " << pm << endln;
  s << "\tyu: " << yu << endln;
  s << "\tpu: " << pu << endln;
}

int
ReeseSandBackbone::setVariable (char *argv)
{
  return -1;
}

int
ReeseSandBackbone::getVariable (int varID, double &theValue)
{
  return -1;
}

int
ReeseSandBackbone::sendSelf(int commitTag, Channel &theChannel)
{
  int res = 0;
  
  static Vector data(6);
  
  data(0) = this->getTag();
  data(1) = kx;
  data(2) = ym;
  data(3) = pm;
  data(4) = yu;
  data(5) = pu;
  
  res += theChannel.sendVector(this->getDbTag(), commitTag, data);
  if (res < 0) {
    opserr << "ReeseSandBackbone::sendSelf -- could not send Vector" << endln;

    return res;
  }
  
  return res;
}

int
ReeseSandBackbone::recvSelf(int commitTag, Channel &theChannel, 
			     FEM_ObjectBroker &theBroker)
{
  int res = 0;
  
  static Vector data(6);
  
  res += theChannel.recvVector(this->getDbTag(), commitTag, data);
  if (res < 0) {
    opserr << "ReeseSandBackbone::recvSelf -- could not receive Vector" << endln;

    return res;
  }
  
  this->setTag(int(data(0)));
  kx = data(1);
  ym = data(2);
  pm = data(3);
  yu = data(4);
  pu = data(5);
  
  return res;
}
