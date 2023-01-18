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
// $Date: 2008-11-06 15:27:11 $
// $Source: /usr/local/cvs/OpenSees/SRC/material/uniaxial/BackboneMaterial.cpp,v $

// Written: MHS
// Created: Aug 2000
//
// Description: This file contains the class definition for 
// BackboneMaterial.  BackboneMaterial uses a HystereticBackbone
// object to represent a path-independent uniaxial material.  Since
// it is path-independent, no state information is stored by
// BackboneMaterial.

#include <BackboneMaterial.h>
#include <HystereticBackbone.h>

#include <stdlib.h>
#include <ID.h>
#include <Channel.h>
#include <FEM_ObjectBroker.h>

#include <OPS_Globals.h>
#include <elementAPI.h>

void* OPS_BackboneMaterial()
{
    int argc = OPS_GetNumRemainingInputArgs();
    if (argc < 2) {
	opserr << "WARNING insufficient arguments\n";
	opserr << "Want: uniaxialMaterial Backbone tag? bbTag? <multiplier?>\n";
	return 0;
    }
      
    // tag, bbTag;
    int idata[2];
    int numdata = 2;
    if (OPS_GetIntInput(&numdata, idata) < 0) {
	opserr << "WARNING invalid tags\n";
	opserr << "Backbone material: " << idata[0] << "\n";
	return 0;
    }

    double multiplier = 1.0;
    if (argc > 2) {
      numdata = 1;
      if (OPS_GetDoubleInput(&numdata, &multiplier) < 0) {
	opserr << "WARNING invalid multiplier\n";
	opserr << "Backbone material: " << idata[0] << "\n";
	return 0;      
      }
    }
    
    HystereticBackbone *backbone = OPS_getHystereticBackbone(idata[1]);
		
    if (backbone == 0) {
	opserr << "WARNING backbone does not exist\n";
	opserr << "backbone: " << idata[1]; 
	opserr << "\nuniaxialMaterial Backbone: " << idata[0] << "\n";
	return 0;
    }
      
    return new BackboneMaterial(idata[0], *backbone, multiplier);
}

BackboneMaterial::BackboneMaterial(int tag, HystereticBackbone &backbone, double mult)
  :UniaxialMaterial(tag,MAT_TAG_Backbone), theBackbone(0), strain(0.0), multiplier(mult)
{
  theBackbone = backbone.getCopy();

  if (theBackbone == 0) {
    opserr <<  "BackboneMaterial::BackboneMaterial -- failed to get copy of material\n";
    exit(-1);
  }
}

BackboneMaterial::BackboneMaterial()
  :UniaxialMaterial(0,MAT_TAG_Backbone), theBackbone(0), strain(0.0), multiplier(1.0)
{

}

BackboneMaterial::~BackboneMaterial()
{
  if (theBackbone != 0)
    delete theBackbone;
}

int 
BackboneMaterial::setTrialStrain(double eps, double epsdot)
{
  strain = eps;

  return 0;
}

double 
BackboneMaterial::getStress(void)
{
  return multiplier*theBackbone->getStress(strain);
}


double 
BackboneMaterial::getTangent(void)
{
  return multiplier*theBackbone->getTangent(strain);
}

double 
BackboneMaterial::getDampTangent(void)
{
  return 0.0;
}

double 
BackboneMaterial::getInitialTangent(void)
{
  return multiplier*theBackbone->getTangent(0.0);
}

double 
BackboneMaterial::getStrain(void)
{
  return strain;
}

double 
BackboneMaterial::getStrainRate(void)
{
  return 0.0;
}

int 
BackboneMaterial::commitState(void)
{
  return 0; // commit nothing, path independent
}

int 
BackboneMaterial::revertToLastCommit(void)
{
  return 0;
}

int 
BackboneMaterial::revertToStart(void)
{
  return 0;
}

UniaxialMaterial *
BackboneMaterial::getCopy(void)
{
  BackboneMaterial *theCopy = 
    new BackboneMaterial(this->getTag(), *theBackbone, multiplier);
  
  theCopy->strain = strain;

  return theCopy;
}


int 
BackboneMaterial::sendSelf(int cTag, Channel &theChannel)
{
  int dbTag = this->getDbTag();

  static ID data(3);
  data(0) = this->getTag();
  data(1) = theBackbone->getClassTag();
  int bbDbTag = theBackbone->getDbTag();
  if (bbDbTag == 0) {
    bbDbTag = theChannel.getDbTag();
    if (bbDbTag != 0)
      theBackbone->setDbTag(bbDbTag);
  }
  data(2) = bbDbTag;
  if (theChannel.sendID(dbTag, cTag, data) < 0) {
    opserr << "BackboneMaterial::sendSelf -- could not send ID" << endln;
    return -1;
  }

  static Vector ddata(2);
  ddata(0) = strain;
  ddata(1) = multiplier;
  if (theChannel.sendVector(dbTag, cTag, ddata) < 0) {
    opserr << "BackboneMaterial::sendSelf -- could not send Vector" << endln;
    return -2;
  }
  
  if (theBackbone->sendSelf(cTag, theChannel) < 0) {
    opserr << "BackboneMaterial::sendSelf -- could not send HystereticBackbone" << endln;
    return -3;
  }
  
  return 0;
}

int 
BackboneMaterial::recvSelf(int cTag, Channel &theChannel, 
			       FEM_ObjectBroker &theBroker)
{
  int dbTag = this->getDbTag();

  static ID data(3);
  if (theChannel.recvID(dbTag, cTag, data) < 0) {
    opserr << "BackboneMaterial::recvSelf -- could not receive ID" << endln;
    return -1;
  }
  this->setTag(data(0));


  int bbClassTag = data(1);
  // Check if the material is null; if so, get a new one
  if (theBackbone == 0) {
    theBackbone = theBroker.getNewHystereticBackbone(bbClassTag);
    if (theBackbone == 0) {
      opserr << " BackboneMaterial::recvSelf -- could not get a HystereticBackbone" << endln;
      return -1;
    }
  }
  // Check that the material is of the right type; if not, delete
  // the current one and get a new one of the right type
  if (theBackbone->getClassTag() != bbClassTag) {
    delete theBackbone;
    theBackbone = theBroker.getNewHystereticBackbone(bbClassTag);
    if (theBackbone == 0) {
      opserr << "BackboneMaterial::recvSelf -- could not get a HystereticBackbone" << endln;
      return -1;
    }
  }

  static Vector ddata(2);
  if (theChannel.recvVector(dbTag, cTag, ddata) < 0) {
    opserr << "BackboneMaterial::recvSelf -- could not receive Vector" << endln;
    return -2;
  }
  strain = ddata(0);
  multiplier = ddata(1);
  
  // Now, receive the material
  theBackbone->setDbTag(data(2));
  if (theBackbone->recvSelf(cTag, theChannel, theBroker) < 0) {
    opserr << "BackboneMaterial::recvSelf -- could not receive HystereticBackbone" << endln;;
    return -3;
  }
  
  return 0;
}

void 
BackboneMaterial::Print(OPS_Stream &s, int flag)
{
  s << "BackboneMaterial tag: " << this->getTag() << endln;
  s << " multiplier = " << multiplier << endln;
  theBackbone->Print(s, flag);
}
