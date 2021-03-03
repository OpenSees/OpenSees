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
// $Source: /usr/local/cvs/OpenSees/SRC/material/uniaxial/MultiplierMaterial.cpp,v $

// Written: MHS
// Created: Aug 2001
//
// Description: This file contains the class definition for 
// MultiplierMaterial. MultiplierMaterial wraps a UniaxialMaterial
// and multiplies the stress and tangent of the wrapped UniaxialMaterial
// object by a factor. This wrapper can be used to apply overstrength
// factors to materials and p-y multipliers for shadowing effects in
// pile groups. An example is shown below with a multiplier of 0.8.
// and imposes min and max strain limits.
 
#include <stdlib.h>
#include <string.h>

#include <MultiplierMaterial.h>
#include <ID.h>
#include <Channel.h>
#include <FEM_ObjectBroker.h>
#include <Information.h>
#include <Parameter.h>

#include <OPS_Globals.h>

#include <elementAPI.h>
#define OPS_Export 

OPS_Export void *
OPS_MultiplierMaterial(void)
{
  // Pointer to a uniaxial material that will be returned
  UniaxialMaterial *theMaterial = 0;
  UniaxialMaterial *theOtherMaterial = 0;
  int    iData[2];

  int argc = OPS_GetNumRemainingInputArgs();
  if (argc < 2) {
    opserr << "WARNING insufficient args, uniaxialMaterial Multiplier $tag $otherTag $multiplier" << endln;
    return 0;
  }

  int numData = 2;
  if (OPS_GetIntInput(&numData, iData) < 0) {
    opserr << "WARNING invalid uniaxialMaterial Multiplier $tag $otherTag $multiplier" << endln;
    return 0;
  }

  theOtherMaterial = OPS_GetUniaxialMaterial(iData[1]);
  if (theOtherMaterial == 0) {
    opserr << "WARNING invalid otherTag uniaxialMaterial Multiplier tag: " << iData[0] << endln;
    return 0;	
  }

  double multiplier = 1.0;
  numData = 1;
  if (OPS_GetDouble(&numData,&multiplier) < 0) {
    opserr << "WARNING invalid input uniaxialMaterial Multiplier tag: " << iData[0] << endln;
    return 0;
  }

  // Parsing was successful, allocate the material
  theMaterial = new MultiplierMaterial(iData[0], *theOtherMaterial, multiplier);

  if (theMaterial == 0) {
    opserr << "WARNING could not create uniaxialMaterial of type MultiplierMaterial\n";
    return 0;
  }

  return theMaterial;
}

MultiplierMaterial::MultiplierMaterial(int tag, UniaxialMaterial &material, double mult)
  :UniaxialMaterial(tag,MAT_TAG_Multiplier), theMaterial(0), multiplier(mult), parameterID(0)
{
  theMaterial = material.getCopy();

  if (theMaterial == 0) {
    opserr <<  "MultiplierMaterial::MultiplierMaterial -- failed to get copy of material\n";
    exit(-1);
  }
}

MultiplierMaterial::MultiplierMaterial()
  :UniaxialMaterial(0,MAT_TAG_Multiplier), theMaterial(0), multiplier(1.0), parameterID(0)
{

}

MultiplierMaterial::~MultiplierMaterial()
{
  if (theMaterial)
    delete theMaterial;
}

int 
MultiplierMaterial::setTrialStrain(double strain, double strainRate)
{
  return theMaterial->setTrialStrain(strain, strainRate);
}


int 
MultiplierMaterial::setTrialStrain(double strain, double temp, double strainRate)
{
  return theMaterial->setTrialStrain(strain, temp, strainRate);
}


double 
MultiplierMaterial::getStress(void)
{
  return multiplier*theMaterial->getStress();
}

double 
MultiplierMaterial::getTangent(void)
{
  return multiplier*theMaterial->getTangent();
}

double 
MultiplierMaterial::getDampTangent(void)
{
  return multiplier*theMaterial->getDampTangent();
}



double 
MultiplierMaterial::getStrain(void)
{
  return theMaterial->getStrain();
}

double 
MultiplierMaterial::getStrainRate(void)
{
  return theMaterial->getStrainRate();
}

int 
MultiplierMaterial::commitState(void)
{	
  return theMaterial->commitState();
}

int 
MultiplierMaterial::revertToLastCommit(void)
{
  return theMaterial->revertToLastCommit();
}

int 
MultiplierMaterial::revertToStart(void)
{
  return theMaterial->revertToStart();
}

UniaxialMaterial *
MultiplierMaterial::getCopy(void)
{
  MultiplierMaterial *theCopy = 
    new MultiplierMaterial(this->getTag(), *theMaterial, multiplier);
        
  return theCopy;
}

int 
MultiplierMaterial::sendSelf(int cTag, Channel &theChannel)
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
    opserr << "MultiplierMaterial::sendSelf() - failed to send the ID\n";
    return -1;
  }

  static Vector dataVec(1);
  dataVec(0) = multiplier;

  if (theChannel.sendVector(dbTag, cTag, dataVec) < 0) {
    opserr << "MultiplierMaterial::sendSelf() - failed to send the Vector\n";
    return -2;
  }

  if (theMaterial->sendSelf(cTag, theChannel) < 0) {
    opserr << "MultiplierMaterial::sendSelf() - failed to send the Material\n";
    return -3;
  }

  return 0;
}

int 
MultiplierMaterial::recvSelf(int cTag, Channel &theChannel, 
			 FEM_ObjectBroker &theBroker)
{
  int dbTag = this->getDbTag();

  static ID dataID(3);
  if (theChannel.recvID(dbTag, cTag, dataID) < 0) {
    opserr << "MultiplierMaterial::recvSelf() - failed to get the ID\n";
    return -1;
  }
  this->setTag(int(dataID(0)));

  // as no way to change material, don't have to check classTag of the material 
  if (theMaterial == 0) {
    int matClassTag = int(dataID(1));
    theMaterial = theBroker.getNewUniaxialMaterial(matClassTag);
    if (theMaterial == 0) {
      opserr << "MultiplierMaterial::recvSelf() - failed to create Material with classTag " 
	   << dataID(0) << endln;
      return -2;
    }
  }
  theMaterial->setDbTag(dataID(2));

  static Vector dataVec(1);
  if (theChannel.recvVector(dbTag, cTag, dataVec) < 0) {
    opserr << "MultiplierMaterial::recvSelf() - failed to get the Vector\n";
    return -3;
  }

  multiplier = dataVec(0);

  if (theMaterial->recvSelf(cTag, theChannel, theBroker) < 0) {
    opserr << "MultiplierMaterial::recvSelf() - failed to get the Material\n";
    return -4;
  }
  return 0;
}

void 
MultiplierMaterial::Print(OPS_Stream &s, int flag)
{
  s << "MultiplierMaterial tag: " << this->getTag() << endln;
  s << "\tMaterial: " << theMaterial->getTag() << endln;
  s << "\tMultiplier: " << multiplier << endln;
}

int
MultiplierMaterial::setParameter(const char **argv, int argc, Parameter &param)
{
  if (strcmp(argv[0],"multiplier") == 0) {
    param.setValue(multiplier);
    return param.addObject(1,this);
  }
  return theMaterial->setParameter(argv, argc, param);
}

int
MultiplierMaterial::updateParameter(int parameterID, Information &info)
{
  switch (parameterID) {
  case -1:
    return -1;
  case 1:
    this->multiplier = info.theDouble;
    break;
  }

  return 0;
}

int
MultiplierMaterial::activateParameter(int paramID)
{
  parameterID = paramID;

  return 0;
}

double
MultiplierMaterial::getStressSensitivity(int gradIndex, bool conditional)
{
  // dsig = dF*sigma + F*dsigma
  if (parameterID == 1)
    return theMaterial->getStress(); // dF*sigma where dF=1
  else
    return multiplier*theMaterial->getStressSensitivity(gradIndex,conditional);
}

double
MultiplierMaterial::getStrainSensitivity(int gradIndex)
{
  return theMaterial->getStrainSensitivity(gradIndex);
}

double
MultiplierMaterial::getInitialTangentSensitivity(int gradIndex)
{
  if (parameterID == 1)
    return theMaterial->getInitialTangent();
  else
    return multiplier*theMaterial->getInitialTangentSensitivity(gradIndex);
}

double
MultiplierMaterial::getDampTangentSensitivity(int gradIndex)
{
  if (parameterID == 1)
    return theMaterial->getDampTangent();
  else
    return multiplier*theMaterial->getDampTangentSensitivity(gradIndex);
}

double
MultiplierMaterial::getRhoSensitivity(int gradIndex)
{
  return theMaterial->getRhoSensitivity(gradIndex);
}

int
MultiplierMaterial::commitSensitivity(double strainGradient, int gradIndex, int numGrads)
{
  return theMaterial->commitSensitivity(strainGradient, gradIndex, numGrads);
}
