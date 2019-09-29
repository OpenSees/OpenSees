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
// Created: Aug 2000
//
// Description: This file contains the implementation of 
// LinearCappedBackbone, which is a linear cap given by a
// cap deformation, slope, and residual strength.

#include <LinearCappedBackbone.h>
#include <Vector.h>
#include <ID.h>
#include <Channel.h>
#include <FEM_ObjectBroker.h>

LinearCappedBackbone::LinearCappedBackbone(int tag, HystereticBackbone &backbone,
					   double def, double slope, double res):
  HystereticBackbone(tag,BACKBONE_TAG_LinearCapped),
  eCap(def), E(slope), eRes(0.0), sRes(res)
{
  theBackbone = backbone.getCopy();

  if (theBackbone == 0)
    opserr << "LinearCappedBackbone::LinearCappedBackbone -- failed to get copy of backbone" << endln;
  
  sCap = theBackbone->getStress(eCap);
  
  if (E != 0.0)
    eRes = eCap + (sRes-sCap)/E;
  else
    eRes = eCap + (sRes-sCap)*1.0e9;
}

LinearCappedBackbone::LinearCappedBackbone():
  HystereticBackbone(0,BACKBONE_TAG_LinearCapped),
  theBackbone(0), eCap(0.0), E(0.0), eRes(0.0), sRes(0.0)
{
  
}

LinearCappedBackbone::~LinearCappedBackbone()
{
  if (theBackbone)
    delete theBackbone;
}

double
LinearCappedBackbone::getTangent (double strain)
{
  if (strain < eCap)
    return theBackbone->getTangent(strain);
  else if (strain < eRes)
    return E;
  else
    return 1.0e-9;
}

double
LinearCappedBackbone::getStress (double strain)
{
  if (strain < eCap)
    return theBackbone->getStress(strain);
  else if (strain < eRes)
    return sCap + E*(strain-eCap);
  else
    return sRes;
}

double
LinearCappedBackbone::getEnergy (double strain)
{
  if (strain < eCap)
    return theBackbone->getEnergy(strain);
  else
    return theBackbone->getEnergy(eCap) +
      0.5*(sCap + this->getStress(strain))*(strain-eCap);
}

double
LinearCappedBackbone::getYieldStrain(void)
{
  double e = theBackbone->getYieldStrain();
  
  return (e < eCap) ? e : eCap;
}

HystereticBackbone*
LinearCappedBackbone::getCopy(void)
{
  LinearCappedBackbone *theCopy = 
    new LinearCappedBackbone (this->getTag(), *theBackbone, eCap, E, sRes);
  
  return theCopy;
}

void
LinearCappedBackbone::Print (OPS_Stream &s, int flag)
{
  s << "LinearCappedBackbone, tag: " << this->getTag() << endln;
  s << "\tbackbone: " << theBackbone->getTag() << endln;
  s << "\teCap: " << eCap << endln;
  s << "\tE: " << E << endln;
  s << "\tsRes: " << sRes << endln;
  s << "\teRes: " << eRes << endln;
}

int
LinearCappedBackbone::setVariable (char *argv)
{
  return -1;
}

int
LinearCappedBackbone::getVariable (int varID, double &theValue)
{
  return -1;
}

int
LinearCappedBackbone::sendSelf(int cTag, Channel &theChannel)
{
  int res = 0;
  
  static Vector data(6);
  
  data(0) = this->getTag();
  data(1) = eCap;
  data(2) = sCap;
  data(3) = E;
  data(4) = eRes;
  data(5) = sRes;
  
  res += theChannel.sendVector(this->getDbTag(), cTag, data);
  if (res < 0) {
    opserr << "LinearCappedBackbone::sendSelf -- could not send Vector" << endln;
    return res;
  }
  
  static ID classTags(2);
  
  int clTag = theBackbone->getClassTag();
  int dbTag = theBackbone->getDbTag();
  
  classTags(0) = clTag;
  
  if (dbTag == 0) {
    dbTag = theChannel.getDbTag();
    if (dbTag != 0)
      theBackbone->setDbTag(dbTag);
  }
  
  classTags(1) = dbTag;
  
  res += theChannel.sendID(this->getDbTag(), cTag, classTags);
  if (res < 0) {
    opserr << "LinearCappedBackbone::sendSelf -- could not send ID" << endln;
    return res;
  }
  
  res += theBackbone->sendSelf(cTag, theChannel);
  if (res < 0) {
    opserr << "LinearCappedBackbone::sendSelf -- could not send HystereticBackbone" << endln;
    return res;
  }
  
  return res;
}

int
LinearCappedBackbone::recvSelf(int cTag, Channel &theChannel, 
			       FEM_ObjectBroker &theBroker)
{
  int res = 0;
  
  static Vector data(6);
  
  res += theChannel.recvVector(this->getDbTag(), cTag, data);
  if (res < 0) {
    opserr << "LinearCappedBackbone::recvSelf -- could not receive Vector" << endln;
    return res;
  }
  
  this->setTag((int)data(0));
  eCap = data(1);
  sCap = data(2);
  E = data(3);
  eRes = data(4);
  sRes = data(5);
  
  static ID classTags(2);
  
  res += theChannel.recvID(this->getDbTag(), cTag, classTags);
  if (res < 0) {
    opserr << "LinearCappedBackbone::recvSelf -- could not receive ID" << endln;
    return res;
  }
  
  // Check if the backbone is null; if so, get a new one
  if (theBackbone == 0) {
    //theBackbone = theBroker.getNewHystereticBackbone(classTags(0));
    if (theBackbone == 0) {
      opserr << "LinearCappedBackbone::recvSelf -- could not get a HystereticBackbone" << endln;
      return -1;
    }
  }
  // Check that the backbone is of the right type; if not, delete
  // the current one and get a new one of the right type
  if (theBackbone->getClassTag() != classTags(0)) {
    delete theBackbone;
    //theBackbone = theBroker.getNewHystereticBackbone(classTags(0));
    if (theBackbone == 0) {
      opserr << "LinearCappedBackbone::recvSelf -- could not get a HystereticBackbone" << endln;
      return -1;
    }
  }
  
  // Now, receive the backbone
  theBackbone->setDbTag(classTags(1));
  res += theBackbone->recvSelf(cTag, theChannel, theBroker);
  if (res < 0) {
    opserr << "LinearCappedBackbone::recvSelf -- could not receive HystereticBackbone" << endln;
    return res;
  }
  
  return res;
}
