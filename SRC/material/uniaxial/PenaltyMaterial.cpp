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
// PenaltyMaterial. PenaltyMaterial adds a small stiffness to
// its wrapped UniaxialMaterial object. This wrapper can help you
// avoid a singular stiffness due to perfect plasticity and is a
// downsized approach to putting the wrapped material in parallel
// with an ElasticMaterial.

#include <stdlib.h>
#include <string.h>

#include <PenaltyMaterial.h>
#include <ID.h>
#include <Channel.h>
#include <FEM_ObjectBroker.h>
#include <Information.h>
#include <Parameter.h>

#include <OPS_Globals.h>

#include <elementAPI.h>

void *
OPS_PenaltyMaterial(void)
{
  // Pointer to a uniaxial material that will be returned
  UniaxialMaterial *theMaterial = 0;
  UniaxialMaterial *theOtherMaterial = 0;
  int    iData[2];

  int argc = OPS_GetNumRemainingInputArgs();
  if (argc < 2) {
    opserr << "WARNING insufficient args, uniaxialMaterial Penalty $tag $otherTag $penalty <-noStress>" << endln;
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

  bool addStress = true;
  if (OPS_GetNumRemainingInputArgs() > 0) {
    const char *argvLoc = OPS_GetString();
    if (strcmp(argvLoc,"-noStress") == 0) {
      addStress = false;
    }
  }
  
  // Parsing was successful, allocate the material
  theMaterial = new PenaltyMaterial(iData[0], *theOtherMaterial, penalty, addStress);

  if (theMaterial == 0) {
    opserr << "WARNING could not create uniaxialMaterial of type PenaltyMaterial\n";
    return 0;
  }

  return theMaterial;
}

PenaltyMaterial::PenaltyMaterial(int tag, UniaxialMaterial &material, double mult, bool addSig)
  :UniaxialMaterial(tag,MAT_TAG_Penalty), theMaterial(0), penalty(mult), addStress(addSig), parameterID(0)
{
  theMaterial = material.getCopy();

  if (theMaterial == 0) {
    opserr <<  "PenaltyMaterial::PenaltyMaterial -- failed to get copy of material\n";
    //exit(-1);
  }
}

PenaltyMaterial::PenaltyMaterial()
  :UniaxialMaterial(0,MAT_TAG_Penalty), theMaterial(0), penalty(0.0), addStress(true), parameterID(0)
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
  if (theMaterial)
    return theMaterial->setTrialStrain(strain, strainRate);
  else
    return -1;
}


int 
PenaltyMaterial::setTrialStrain(double strain, double temp, double strainRate)
{
  if (theMaterial)
    return theMaterial->setTrialStrain(strain, temp, strainRate);
  else
    return -1;
}


double 
PenaltyMaterial::getStress(void)
{
  if (theMaterial) {
    double sig = theMaterial->getStress();
    if (addStress)
      sig += penalty*theMaterial->getStrain();
    return sig;
  }
  else
    return 0.0;
}

double 
PenaltyMaterial::getTangent(void)
{
  if (theMaterial)
    return theMaterial->getTangent() + penalty;
  else
    return 0.0;
}

double 
PenaltyMaterial::getDampTangent(void)
{
  if (theMaterial)
    return theMaterial->getDampTangent();
  else
    return 0.0;
}

double 
PenaltyMaterial::getStrain(void)
{
  if (theMaterial)
    return theMaterial->getStrain();
  else
    return 0.0;
}

double 
PenaltyMaterial::getStrainRate(void)
{
  if (theMaterial)  
    return theMaterial->getStrainRate();
  else
    return 0.0;
}

int 
PenaltyMaterial::commitState(void)
{
  if (theMaterial)
    return theMaterial->commitState();
  else
    return -1;
}

int 
PenaltyMaterial::revertToLastCommit(void)
{
  if (theMaterial)
    return theMaterial->revertToLastCommit();
  else
    return -1;
}

int 
PenaltyMaterial::revertToStart(void)
{
  if (theMaterial)  
    return theMaterial->revertToStart();
  else
    return -1;
}

UniaxialMaterial *
PenaltyMaterial::getCopy(void)
{
  PenaltyMaterial *theCopy = 0;
  if (theMaterial)
    theCopy = new PenaltyMaterial(this->getTag(), *theMaterial, penalty, addStress);
        
  return theCopy;
}

int 
PenaltyMaterial::sendSelf(int cTag, Channel &theChannel)
{
  if (theMaterial == 0) {
    opserr << "PenaltyMaterial::sendSelf() - theMaterial is null, nothing to send\n";
    return -1;
  }
  
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

  static Vector dataVec(2);
  dataVec(0) = penalty;
  dataVec(1) = addStress ? 1.0 : 0.0;
  
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
	   << dataID(1) << endln;
      return -2;
    }
  }
  theMaterial->setDbTag(dataID(2));

  static Vector dataVec(2);
  if (theChannel.recvVector(dbTag, cTag, dataVec) < 0) {
    opserr << "PenaltyMaterial::recvSelf() - failed to get the Vector\n";
    return -3;
  }

  penalty = dataVec(0);
  addStress = (dataVec(1) != 0.0) ? true : false;
  
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
  if (theMaterial)
    s << "\tMaterial: " << theMaterial->getTag() << endln;
  else
    s << "\tMaterial is NULL" << endln;  
  s << "\tPenalty: " << penalty << endln;
}

int
PenaltyMaterial::setParameter(const char **argv, int argc, Parameter &param)
{
  if (strcmp(argv[0],"penalty") == 0) {
    param.setValue(penalty);
    return param.addObject(1,this);
  }
  if (theMaterial)
    return theMaterial->setParameter(argv, argc, param);
  else
    return -1;
}

int
PenaltyMaterial::updateParameter(int parameterID, Information &info)
{
  switch (parameterID) {
  case -1:
    return -1;
  case 1:
    this->penalty = info.theDouble;
    break;
  }

  return 0;
}

int
PenaltyMaterial::activateParameter(int paramID)
{
  parameterID = paramID;

  return 0;
}

double
PenaltyMaterial::getStressSensitivity(int gradIndex, bool conditional)
{
  if (theMaterial == 0)
    return 0.0;
  
  // dsig = dsigma + dalpha*strain < + alpha*dstrain>
  if (parameterID == 1)
    return theMaterial->getStrain(); // dalpha*strain where dalpha=1
  else
    return theMaterial->getStressSensitivity(gradIndex,conditional);
}

double
PenaltyMaterial::getStrainSensitivity(int gradIndex)
{
  if (theMaterial)
    return theMaterial->getStrainSensitivity(gradIndex);
  else
    return 0.0;
}

double
PenaltyMaterial::getInitialTangentSensitivity(int gradIndex)
{
  if (parameterID == 1)
    return 1.0;
  else {
    if (theMaterial)
      return theMaterial->getInitialTangentSensitivity(gradIndex);
    else
      return 0.0;
  }
}

double
PenaltyMaterial::getDampTangentSensitivity(int gradIndex)
{
  if (theMaterial)
    return theMaterial->getDampTangentSensitivity(gradIndex);
  else
    return 0.0;
}

double
PenaltyMaterial::getRhoSensitivity(int gradIndex)
{
  if (theMaterial)
    return theMaterial->getRhoSensitivity(gradIndex);
  else
    return 0.0;
}

int
PenaltyMaterial::commitSensitivity(double strainGradient, int gradIndex, int numGrads)
{
  if (theMaterial)
    return theMaterial->commitSensitivity(strainGradient, gradIndex, numGrads);
  else
    return -1;
}
