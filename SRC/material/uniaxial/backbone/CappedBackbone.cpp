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
// CappedBackbone, which overlays two backbone curves.

#include <CappedBackbone.h>
#include <Vector.h>
#include <ID.h>
#include <Channel.h>
#include <FEM_ObjectBroker.h>

#include <math.h>

#include <G3Globals.h>

CappedBackbone::CappedBackbone (int tag, HystereticBackbone &backbone,
				HystereticBackbone &cap):
  HystereticBackbone(tag,BACKBONE_TAG_Capped), theBackbone(0), theCap(0), eCap(0.0)
{
  theBackbone = backbone.getCopy();
  
  if (theBackbone == 0)
    opserr << "CappedBackbone::CappedBackbone -- failed to get copy of backbone" << endln;
  
  theCap = cap.getCopy();
  
  if (theCap == 0)
    opserr << "CappedBackbone::CappedBackbone -- failed to get copy of cap" << endln;
  
  // Determine eCap, the root of h(x) := f(x) - g(x)
  // This is the strain at which the backbone and the cap intersect
  
  const double tol = 1.0e-14;
  eCap = 1.0e-14;	// Don't want to start at zero
  
  // Newton-Raphson
  double residual = theBackbone->getStress(eCap) - theCap->getStress(eCap);
  
  while (fabs(residual) > tol) {
    double tangent = theBackbone->getTangent(eCap) - theCap->getTangent(eCap);
    
    if (fabs(tangent) > tol)
      eCap -= residual/tangent;
    else {
      opserr << "CappedBackbone::CappedBackbone -- zero tangent encountered in finding cap strain" << endln;

      eCap = theBackbone->getYieldStrain();	// Or something else???
      break;
    }
    
    residual = theBackbone->getStress(eCap) - theCap->getStress(eCap);
  }
  
  /*
  // Regula-falsi ... slower, but converges
  double x1 = 0.0;
  double x2 = 0.1;
  
  double y1 = -1.0e-16;
  double y2 = theBackbone->getStress(x2) - theCap->getStress(x2) - 0.1;
  double y0;
  
  do {
  eCap = (x1*y2-x2*y1) / (y2-y1);
  y1 = theBackbone->getStress(x1) - theCap->getStress(x1);
  y2 = theBackbone->getStress(x2) - theCap->getStress(x2);
  
  y0 = theBackbone->getStress(eCap) - theCap->getStress(eCap);
  
  if (fabs(y0) < tol)	// Converged solution
  break;
  
  if (y1*y0 < 0.0)
  x2 = eCap;
  else
  x1 = eCap;
  } while (fabs(y1-y2) > tol);
  */
}

CappedBackbone::CappedBackbone():
  HystereticBackbone(0,BACKBONE_TAG_Capped),
  theBackbone(0), theCap(0), eCap(0.0)
{
  
}

CappedBackbone::~CappedBackbone()
{
  if (theBackbone)
    delete theBackbone;
  
  if (theCap)
    delete theCap;
}

double
CappedBackbone::getTangent (double strain)
{
  if (strain < eCap)
    return theBackbone->getTangent(strain);
  else
    return theCap->getTangent(strain);
}

double
CappedBackbone::getStress (double strain)
{
  if (strain < eCap)
    return theBackbone->getStress(strain);
  else
    return theCap->getStress(strain);
}

double
CappedBackbone::getEnergy (double strain)
{
  if (strain < eCap)
    return theBackbone->getEnergy(strain);
  else
    return theBackbone->getEnergy(eCap) + 
      theCap->getEnergy(strain) - theCap->getEnergy(eCap);
}

double
CappedBackbone::getYieldStrain(void)
{
  double eb = theBackbone->getYieldStrain();
  
  return (eCap < eb) ? eCap : eb;
}

HystereticBackbone*
CappedBackbone::getCopy(void)
{
  CappedBackbone *theCopy = 
    new CappedBackbone (this->getTag(), *theBackbone, *theCap);
  
  return theCopy;
}

void
CappedBackbone::Print(OPS_Stream &s, int flag)
{
  s << "CappedBackbone, tag: " << this->getTag() << endln;
  s << "\tBackbone: " << theBackbone->getTag() << endln;
  s << "\tCap: " << theCap->getTag() << endln;
}

int
CappedBackbone::setVariable (char *argv)
{
  return -1;
}

int
CappedBackbone::getVariable (int varID, double &theValue)
{
  return -1;
}

int
CappedBackbone::sendSelf(int cTag, Channel &theChannel)
{
  // Create and send Vector with CappedBackbone data
  static Vector data(2);
  data(0) = this->getTag();
  data(1) = eCap;
  
  int res = theChannel.sendVector(this->getDbTag(), cTag, data);
  if (res < 0) {
    opserr << "CappedBackbone::sendSelf -- could not send Vector" << endln;

    return res;
  }
  
  // Send backbone and cap class tags
  static ID classTags(4);
  classTags(0) = theBackbone->getClassTag();
  classTags(1) = theCap->getClassTag();
  
  int dbTag;
  
  dbTag = theBackbone->getDbTag();
  
  if (dbTag == 0) {
    dbTag = theChannel.getDbTag();
    if (dbTag != 0)
      theBackbone->setDbTag(dbTag);
  }
  
  classTags(2) = dbTag;
  
  dbTag = theCap->getDbTag();
  
  if (dbTag == 0) {
    dbTag = theChannel.getDbTag();
    if (dbTag != 0)
      theCap->setDbTag(dbTag);
  }
  
  classTags(3) = dbTag;
  
  res += theChannel.sendID(this->getDbTag(), cTag, classTags);
  if (res < 0) {
    opserr << "CappedBackbone::sendSelf -- could not send ID" << endln;

    return res;
  }
  
  // Ask backbone to send itself
  res += theBackbone->sendSelf(cTag, theChannel);
  if (res < 0) {
    opserr << "CappedBackbone::sendSelf -- could not send HystereticBackbone" << endln;

    return res;
  }
  
  // Ask cap to send itself
  res += theCap->sendSelf(cTag, theChannel);
  if (res < 0) {
    opserr << "CappedBackbone::sendSelf -- could not send HystereticBackbone" << endln;

    return res;
  }
  
  return res;
}

int
CappedBackbone::recvSelf(int cTag, Channel &theChannel, 
			 FEM_ObjectBroker &theBroker)
{
  // Create a Vector and receive CappedBackbone data
  static Vector data(2);
  int res = theChannel.recvVector(this->getDbTag(), cTag, data);
  if (res < 0) {
    opserr << "CappedBackbone::recvSelf -- could not receive Vector" << endln;

    return res;
  }
  
  this->setTag((int)data(0));
  eCap = data(1);
  
  // Receive the classTags for the backbone and cap
  static ID classTags(4);
  res += theChannel.recvID(this->getDbTag(), cTag, classTags);
  if (res < 0) {
    opserr << "CappedBackbone::recvSelf -- could not receive ID" << endln;
    
    return res;
  }
  
  // Check if the backbone is null; if so, get a new one
  if (theBackbone == 0) {
    //theBackbone = theBroker.getNewHystereticBackbone(classTags(0));
    if (theBackbone == 0) {
      opserr << "CappedBackbone::recvSelf -- could not receive HystereticBackbone" << endln;
      return -1;
    }
  }
  // Check that the backbone is of the right type; if not, delete
  // the current one and get a new one of the right type
  if (theBackbone->getClassTag() != classTags(0)) {
    delete theBackbone;
    //theBackbone = theBroker.getNewHystereticBackbone(classTags(0));
    if (theBackbone == 0) {
      opserr << "CappedBackbone::recvSelf -- could not get a HystereticBackbone" << endln;
      return -1;
    }
  }
  
  // Now, receive the backbone
  theBackbone->setDbTag(classTags(2));
  res += theBackbone->recvSelf(cTag, theChannel, theBroker);
  if (res < 0) {
    opserr << "CappedBackbone::recvSelf -- could not receive HystereticBackbone" << endln;
    return res;
  }
  
  // Check if the cap is null; if so, get a new one
  if (theCap == 0) {
    //theCap = theBroker.getNewHystereticBackbone(classTags(1));
    if (theCap == 0) {
      opserr << "CappedBackbone::recvSelf -- could not get a HystereticBackbone" << endln;
      return -1;
    }
  }
  // Check that the cap is of the right type; if not, delete
  // the current one and get a new one of the right type
  if (theCap->getClassTag() != classTags(1)) {
    delete theCap;
    //theCap = theBroker.getNewHystereticBackbone(classTags(1));
    if (theCap == 0) {
      opserr << "CappedBackbone::recvSelf -- could not get a HystereticBackbone" << endln;
      return -1;
    }
  }
  
  // Now, receive the cap
  theCap->setDbTag(classTags(3));
  res += theCap->recvSelf(cTag, theChannel, theBroker);
  if (res < 0) {
    opserr << "CappedBackbone::recvSelf -- could not receive HystereticBackbone" << endln;
    return res;
  }
  
  return res;
}
