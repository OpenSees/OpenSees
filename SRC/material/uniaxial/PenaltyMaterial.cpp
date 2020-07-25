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
                                                                        
// $Revision: 1.6 $
// $Date: 2010-09-11 00:50:53 $
// $Source: /usr/local/cvs/OpenSees/SRC/material/uniaxial/PenaltyMaterial.cpp,v $

// Written: MHS
// Created: Aug 2001
//
// Description: This file contains the class definition for 
// PenaltyMaterial.  PenaltyMaterial wraps a UniaxialMaterial
// and imposes min and max strain limits.

#include <stdlib.h>
#include <string.h>

#include <PenaltyMaterial.h>
#include <ID.h>
#include <Channel.h>
#include <FEM_ObjectBroker.h>

#include <OPS_Globals.h>

#include <elementAPI.h>
#define OPS_Export 

OPS_Export void *
OPS_PenaltyMaterial(void)
{
  // Pointer to a uniaxial material that will be returned
  UniaxialMaterial *theMaterial = 0;
  UniaxialMaterial *theOtherMaterial = 0;
  int    iData[2];

  int argc = OPS_GetNumRemainingInputArgs();
  if (argc < 2) {
    opserr << "WARNING insufficient args, uniaxialMaterial Penalty $tag $otherTag $penalty" << endln;
    return 0;
  }

  int numData = 2;
  if (OPS_GetIntInput(&numData, iData) < 0) {
    opserr << "WARNING invalid uniaxialMaterial Penalty $tag $otherTag $penalty" << endln;
    return 0;
  }

  theOtherMaterial = OPS_GetUniaxialMaterial(iData[1]);
  if (theOtherMaterial == 0) {
    opserr << "WARNING invalid otherTag uniaxialMaterial Penalty tag: " << iData[0] << endln;
    return 0;	
  }

  double penalty = 0.0;
  numData = 1;
  if (OPS_GetDouble(&numData,&penalty) < 0) {
    opserr << "WARNING invalid input uniaxialMaterial Penalty tag: " << iData[0] << endln;
    return 0;
  }

  // Parsing was successful, allocate the material
  theMaterial = new PenaltyMaterial(iData[0], *theOtherMaterial, penalty);

  if (theMaterial == 0) {
    opserr << "WARNING could not create uniaxialMaterial of type PenaltyMaterial\n";
    return 0;
  }

  return theMaterial;
}

PenaltyMaterial::PenaltyMaterial(int tag, UniaxialMaterial &material, double mult)
  :UniaxialMaterial(tag,MAT_TAG_Penalty), theMaterial(0), penalty(mult)
{
  theMaterial = material.getCopy();

  if (theMaterial == 0) {
    opserr <<  "PenaltyMaterial::PenaltyMaterial -- failed to get copy of material\n";
    exit(-1);
  }
}

PenaltyMaterial::PenaltyMaterial()
  :UniaxialMaterial(0,MAT_TAG_Penalty), theMaterial(0), penalty(0.0)
{

}

PenaltyMaterial::~PenaltyMaterial()
{
  if (theMaterial)
    delete theMaterial;
}

int 
PenaltyMaterial::setTrialStrain(double strain, double strainRate)
{
  return theMaterial->setTrialStrain(strain, strainRate);
}


int 
PenaltyMaterial::setTrialStrain(double strain, double temp, double strainRate)
{
  return theMaterial->setTrialStrain(strain, temp, strainRate);
}


double 
PenaltyMaterial::getStress(void)
{
  return theMaterial->getStress() + penalty*theMaterial->getStrain();
}

double 
PenaltyMaterial::getTangent(void)
{
  return theMaterial->getTangent() + penalty;
}

double 
PenaltyMaterial::getDampTangent(void)
{
  return theMaterial->getDampTangent();
}



double 
PenaltyMaterial::getStrain(void)
{
  return theMaterial->getStrain();
}

double 
PenaltyMaterial::getStrainRate(void)
{
  return theMaterial->getStrainRate();
}

int 
PenaltyMaterial::commitState(void)
{	
  return theMaterial->commitState();
}

int 
PenaltyMaterial::revertToLastCommit(void)
{
  return theMaterial->revertToLastCommit();
}

int 
PenaltyMaterial::revertToStart(void)
{
  return theMaterial->revertToStart();
}

UniaxialMaterial *
PenaltyMaterial::getCopy(void)
{
  PenaltyMaterial *theCopy = 
    new PenaltyMaterial(this->getTag(), *theMaterial, penalty);
        
  return theCopy;
}

int 
PenaltyMaterial::sendSelf(int cTag, Channel &theChannel)
{
  int dbTag = this->getDbTag();

  static ID dataID(3);
  dataID(0) = this->getTag();
  dataID(1) = theMaterial->getClassTag();
  int matDbTag = theMaterial->getDbTag();
  if ( matDbTag == 0) {
    matDbTag = theChannel.getDbTag();
    theMaterial->setDbTag(matDbTag);
  }
  dataID(2) = matDbTag;
  if (theChannel.sendID(dbTag, cTag, dataID) < 0) {
    opserr << "PenaltyMaterial::sendSelf() - failed to send the ID\n";
    return -1;
  }

  static Vector dataVec(1);
  dataVec(0) = penalty;

  if (theChannel.sendVector(dbTag, cTag, dataVec) < 0) {
    opserr << "PenaltyMaterial::sendSelf() - failed to send the Vector\n";
    return -2;
  }

  if (theMaterial->sendSelf(cTag, theChannel) < 0) {
    opserr << "PenaltyMaterial::sendSelf() - failed to send the Material\n";
    return -3;
  }

  return 0;
}

int 
PenaltyMaterial::recvSelf(int cTag, Channel &theChannel, 
			 FEM_ObjectBroker &theBroker)
{
  int dbTag = this->getDbTag();

  static ID dataID(3);
  if (theChannel.recvID(dbTag, cTag, dataID) < 0) {
    opserr << "PenaltyMaterial::recvSelf() - failed to get the ID\n";
    return -1;
  }
  this->setTag(int(dataID(0)));

  // as no way to change material, don't have to check classTag of the material 
  if (theMaterial == 0) {
    int matClassTag = int(dataID(1));
    theMaterial = theBroker.getNewUniaxialMaterial(matClassTag);
    if (theMaterial == 0) {
      opserr << "PenaltyMaterial::recvSelf() - failed to create Material with classTag " 
	   << dataID(0) << endln;
      return -2;
    }
  }
  theMaterial->setDbTag(dataID(2));

  static Vector dataVec(1);
  if (theChannel.recvVector(dbTag, cTag, dataVec) < 0) {
    opserr << "PenaltyMaterial::recvSelf() - failed to get the Vector\n";
    return -3;
  }

  penalty = dataVec(0);

  if (theMaterial->recvSelf(cTag, theChannel, theBroker) < 0) {
    opserr << "PenaltyMaterial::recvSelf() - failed to get the Material\n";
    return -4;
  }
  return 0;
}

void 
PenaltyMaterial::Print(OPS_Stream &s, int flag)
{
  s << "PenaltyMaterial tag: " << this->getTag() << endln;
  s << "\tMaterial: " << theMaterial->getTag() << endln;
  s << "\tPenalty: " << penalty << endln;
}
