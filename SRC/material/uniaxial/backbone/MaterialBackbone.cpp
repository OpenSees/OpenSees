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
// MaterialBackbone, which treats a uniaxial material as
// a hysteretic backbone by removing path dependency, i.e.
// commitState is never called on the uniaxial material.

#include <MaterialBackbone.h>
#include <UniaxialMaterial.h>
#include <ID.h>
#include <Channel.h>
#include <FEM_ObjectBroker.h>
#include <G3Globals.h>

#include <math.h>

MaterialBackbone::MaterialBackbone (int tag, UniaxialMaterial &material):
  HystereticBackbone(tag,BACKBONE_TAG_Material), theMaterial(0)
{
  theMaterial = material.getCopy();
  
  if (theMaterial == 0)
    opserr << "MaterialBackbone::MaterialBackbone -- failed to get copy of material" << endln;
}

MaterialBackbone::MaterialBackbone ():
  HystereticBackbone(0,BACKBONE_TAG_Material), theMaterial(0)
{

}

MaterialBackbone::~MaterialBackbone()
{
  if (theMaterial)
    delete theMaterial;
}

double
MaterialBackbone::getTangent (double strain)
{
  theMaterial->setTrialStrain(strain);
	
  return theMaterial->getTangent();
}

double
MaterialBackbone::getStress (double strain)
{
  theMaterial->setTrialStrain(strain);

  return theMaterial->getStress();
}

double
MaterialBackbone::getEnergy (double strain)
{
  const double incr = 1.0e-6;
  double energy = 0.0;
  
  // Mid-point integration
  for (double x = incr/2; x < strain; x += incr) {
    theMaterial->setTrialStrain(x);
    energy += theMaterial->getStress();
  }
  
  return energy*incr;
}

double
MaterialBackbone::getYieldStrain(void)
{
  return 0.002; // Fix later
}

HystereticBackbone*
MaterialBackbone::getCopy(void)
{
  MaterialBackbone *theCopy = 
    new MaterialBackbone (this->getTag(), *theMaterial);
  
  return theCopy;
}

void
MaterialBackbone::Print (OPS_Stream &s, int flag)
{
  s << "MaterialBackbone, tag: " << this->getTag() << endln;
  s << "\tmaterial: " << theMaterial->getTag() << endln;
}

int
MaterialBackbone::setVariable (char *argv)
{
  return -1;
}

int
MaterialBackbone::getVariable (int varID, double &theValue)
{
  return -1;
}

int
MaterialBackbone::sendSelf(int cTag, Channel &theChannel)
{
  int res = 0;
  
  static ID classTags(3);
  
  int clTag = theMaterial->getClassTag();
  int dbTag = theMaterial->getDbTag();
  
  classTags(0) = clTag;
  
  if (dbTag == 0) {
    dbTag = theChannel.getDbTag();
    if (dbTag != 0)
      theMaterial->setDbTag(dbTag);
  }
  
  classTags(1) = dbTag;
  classTags(2) = this->getTag();
  
  res += theChannel.sendID(this->getDbTag(), cTag, classTags);
  if (res < 0) {
    opserr << "MaterialBackbone::sendSelf -- could not send ID" << endln;
    return res;
  }
  
  res += theMaterial->sendSelf(cTag, theChannel);
  
  return res;
}

int
MaterialBackbone::recvSelf(int cTag, Channel &theChannel, 
			   FEM_ObjectBroker &theBroker)
{
  int res = 0;
  
  static ID classTags(3);
  
  res += theChannel.recvID(this->getDbTag(), cTag, classTags);
  if (res < 0) {
    opserr << "MaterialBackbone::recvSelf -- could not receive ID" << endln;
    return res;
  }
  
  this->setTag(classTags(2));
  
  // Check if the material is null; if so, get a new one
  if (theMaterial == 0) {
    theMaterial = theBroker.getNewUniaxialMaterial(classTags(0));
    if (theMaterial == 0) {
      opserr << "MaterialBackbone::recvSelf -- could not get a UniaxialMaterial" << endln;
      return -1;
    }
  }
  // Check that the material is of the right type; if not, delete
  // the current one and get a new one of the right type
  if (theMaterial->getClassTag() != classTags(0)) {
    delete theMaterial;
    theMaterial = theBroker.getNewUniaxialMaterial(classTags(0));
    if (theMaterial == 0) {
      opserr << "MaterialBackbone::recvSelf -- could not get a UniaxialMaterial" << endln;
      return -1;
    }
  }
  
  // Now, receive the material
  theMaterial->setDbTag(classTags(1));
  res += theMaterial->recvSelf(cTag, theChannel, theBroker);
  if (res < 0) {
    opserr << "MaterialBackbone::recvSelf -- could not receive UniaxialMaterial" << endln;
    return res;
  }
  
  return res;
}
