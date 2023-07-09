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

#include <OPS_Globals.h>
#include <elementAPI.h>

void* OPS_MaterialBackbone()
{
    int argc = OPS_GetNumRemainingInputArgs();
    if (argc < 2) {
	opserr << "WARNING insufficient arguments\n";
	opserr << "Want: hystereticBackbone tag? matTag? <-compression>\n";
	return 0;
    }
      
    // tag, bbTag;
    int idata[2];
    int numdata = 2;
    if (OPS_GetIntInput(&numdata, idata) < 0) {
	opserr << "WARNING invalid tag\n";
	opserr << "Backbone material: " << idata[0] << "\n";
	return 0;
    }
		
    UniaxialMaterial *material = OPS_getUniaxialMaterial(idata[1]);
		
    if (material == 0) {
	opserr << "WARNING material does not exist\n";
	opserr << "material: " << idata[1]; 
	opserr << "\nhystereticBackbone Material: " << idata[0] << "\n";
	return 0;
    }

    bool compression = false;
    if (OPS_GetNumRemainingInputArgs() > 0) {
      const char *opt = OPS_GetString();
      if (strcmp(opt,"-compression") == 0 || strcmp(opt,"compression") == 0)
	compression = true;
    }
    
    return new MaterialBackbone(idata[0], *material, compression);
}

MaterialBackbone::MaterialBackbone (int tag, UniaxialMaterial &material, bool compression):
  HystereticBackbone(tag,BACKBONE_TAG_Material), theMaterial(0), sign(1)
{
  theMaterial = material.getCopy();
  
  if (theMaterial == 0)
    opserr << "MaterialBackbone::MaterialBackbone -- failed to get copy of material" << endln;

  if (compression)
    sign = -1;
}

MaterialBackbone::MaterialBackbone ():
  HystereticBackbone(0,BACKBONE_TAG_Material), theMaterial(0), sign(1)
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
  theMaterial->setTrialStrain(sign*strain);
	
  return theMaterial->getTangent();
}

double
MaterialBackbone::getStress (double strain)
{
  theMaterial->setTrialStrain(sign*strain);

  return sign*theMaterial->getStress();
}

double
MaterialBackbone::getEnergy (double strain)
{
  const double incr = 1.0e-6;
  double energy = 0.0;
  
  // Mid-point integration
  for (double x = incr/2; x < strain; x += incr) {
    theMaterial->setTrialStrain(sign*x);
    energy += sign*theMaterial->getStress();
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

  theCopy->sign = sign;
  
  return theCopy;
}

void
MaterialBackbone::Print (OPS_Stream &s, int flag)
{
  s << "MaterialBackbone, tag: " << this->getTag() << endln;
  s << "\tmaterial: " << theMaterial->getTag() << endln;
  s << "\tsign: " << sign << endln;
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
  
  static ID classTags(4);
  
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
  classTags(3) = sign;
  
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
  
  static ID classTags(4);
  
  res += theChannel.recvID(this->getDbTag(), cTag, classTags);
  if (res < 0) {
    opserr << "MaterialBackbone::recvSelf -- could not receive ID" << endln;
    return res;
  }
  
  this->setTag(classTags(2));
  sign = classTags(3);
  
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
