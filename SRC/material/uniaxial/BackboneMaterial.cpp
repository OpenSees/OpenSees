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

void * OPS_ADD_RUNTIME_VPV(OPS_Backbone)
{
    int argc = OPS_GetNumRemainingInputArgs() + 2;
    if (argc < 4) {
	opserr << "WARNING insufficient arguments\n";
	opserr << "Want: uniaxialMaterial Backbone tag? bbTag?\n";
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
		
    HystereticBackbone *backbone = OPS_getHystereticBackbone(idata[1]);
		
    if (backbone == 0) {
	opserr << "WARNING backbone does not exist\n";
	opserr << "backbone: " << idata[1]; 
	opserr << "\nuniaxialMaterial Backbone: " << idata[0] << "\n";
	return 0;
    }
      
    return new BackboneMaterial(idata[0], *backbone);
}

BackboneMaterial::BackboneMaterial(int tag, HystereticBackbone &backbone)
  :UniaxialMaterial(tag,MAT_TAG_Backbone), theBackbone(0), strain(0.0)
{
  theBackbone = backbone.getCopy();

  if (theBackbone == 0) {
    opserr <<  "BackboneMaterial::BackboneMaterial -- failed to get copy of material\n";
    exit(-1);
  }
}

BackboneMaterial::BackboneMaterial()
  :UniaxialMaterial(0,MAT_TAG_Backbone), theBackbone(0), strain(0.0)
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
  return theBackbone->getStress(strain);
}


double 
BackboneMaterial::getTangent(void)
{
  return theBackbone->getTangent(strain);
}

double 
BackboneMaterial::getDampTangent(void)
{
  return 0.0;
}

double 
BackboneMaterial::getInitialTangent(void)
{
  return theBackbone->getTangent(0.0);
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
    new BackboneMaterial(this->getTag(), *theBackbone);
  
  theCopy->strain = strain;

  return theCopy;
}


int 
BackboneMaterial::sendSelf(int cTag, Channel &theChannel)
{
  return -1;

  /*
  int res = 0;

  static ID classTags(3);

  int clTag = theBackbone->getClassTag();
  int dbTag = theBackbone->getDbTag();
	
  classTags(0) = clTag;

  if (dbTag == 0) {
    dbTag = theChannel.getDbTag();
    if (dbTag != 0)
      theBackbone->setDbTag(dbTag);
  }
  
  classTags(1) = dbTag;
  classTags(2) = this->getTag();
  
  res = theChannel.sendID(dbTag, cTag, classTags);
  if (res < 0) {
    opserr << "BackboneMaterial::sendSelf -- could not send ID\n";
    return res;
  }
  
  res = theBackbone->sendSelf(cTag, theChannel);
  if (res < 0) {
    opserr << "BackboneMaterial::sendSelf -- could not send UniaxialMaterial\n";
    return res;
  }
  
  return res;
  */
}

int 
BackboneMaterial::recvSelf(int cTag, Channel &theChannel, 
			       FEM_ObjectBroker &theBroker)
{
  return -1;

  /*
  int res = 0;
  
  static ID classTags(3);

  int dbTag = this->getDbTag();

  res = theChannel.recvID(dbTag, cTag, classTags);
  if (res < 0) {
    opserr << "BackboneMaterial::recvSelf -- could not receive ID\n";
    return res;
  }

  this->setTag(int(classTags(2)));

  // Check if the material is null; if so, get a new one
  if (theBackbone == 0) {
    theBackbone = theBroker.getNewUniaxialMaterial(classTags(0));
    if (theBackbone == 0) {
      opserr << " BackboneMaterial::recvSelf -- could not get a UniaxialMaterial\n";
      return -1;
    }
  }
  // Check that the material is of the right type; if not, delete
  // the current one and get a new one of the right type
  if (theBackbone->getClassTag() != classTags(0)) {
    delete theBackbone;
    theBackbone = theBroker.getNewUniaxialMaterial(classTags(0));
    if (theBackbone == 0) {
      opserr << "BackboneMaterial::recvSelf -- could not get a UniaxialMaterial\n";
      exit(-1);
    }
  }
  
  // Now, receive the material
  theBackbone->setDbTag(classTags(1));
  res += theBackbone->recvSelf(cTag, theChannel, theBroker);
  if (res < 0) {
    opserr << "BackboneMaterial::recvSelf -- could not receive UniaxialMaterial\n";
    return res;
  }
  
  return res;
  */
}

void 
BackboneMaterial::Print(OPS_Stream &s, int flag)
{
  s << "BackboneMaterial tag: " << this->getTag() << endln;
  theBackbone->Print(s, flag);
}
