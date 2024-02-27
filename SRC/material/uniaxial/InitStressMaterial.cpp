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
                                                                        
// $Revision: 1.1 $
// $Date: 2010-09-16 00:03:55 $
// $Source: /usr/local/cvs/OpenSees/SRC/material/uniaxial/InitStressMaterial.cpp,v $

// Written: fmk
// Created: Sep 2010

#include <stdlib.h>
#include <math.h>

#include <InitStressMaterial.h>
#include <ID.h>
#include <Channel.h>
#include <FEM_ObjectBroker.h>
#include <Information.h>
#include <Parameter.h>

#include <OPS_Globals.h>

#include <elementAPI.h>
#define OPS_Export 

OPS_Export void *
OPS_InitStressMaterial(void)
{
  // Pointer to a uniaxial material that will be returned
  UniaxialMaterial *theMaterial = 0;
  UniaxialMaterial *theOtherMaterial = 0;

  int    iData[2];
  double dData[1];
  int numData = 2;
  if (OPS_GetIntInput(&numData, iData) != 0) {
    opserr << "WARNING invalid uniaxialMaterial InitStressMaterial $tag $otherTag $sig0" << endln;
    return 0;
  }

  theOtherMaterial = OPS_GetUniaxialMaterial(iData[1]);
  if (theOtherMaterial == 0) {
    opserr << "Could not find material with tag: " << iData[1] << "uniaxialMaterial InitStress $tag $otherTag $sig0" << endln;
    return 0;	
  }

  numData = 1;
  if (OPS_GetDoubleInput(&numData, dData) != 0) {
    opserr << "Invalid Args want: uniaxialMaterial InitStress $tag $otherTag $sig0" << endln;
    return 0;	
  }

  // Parsing was successful, allocate the material
  theMaterial = new InitStressMaterial(iData[0], *theOtherMaterial, dData[0]);

  if (theMaterial == 0) {
    opserr << "WARNING could not create uniaxialMaterial of type InitStressMaterial\n";
    return 0;
  }

  return theMaterial;
}

int
InitStressMaterial::findInitialStrain(void)
{
  // determine the initial strain
  double tol = fabs(sigInit)*1e-12;
  double dSig = sigInit;
  double tStrain = 0.0, tStress = 0.0;
  int count = 0;

  do {
    count++;
    double K = theMaterial->getTangent();
    double dStrain = dSig/K;
    tStrain += dStrain;
    theMaterial->setTrialStrain(tStrain);
    tStress = theMaterial->getStress();
    dSig = sigInit-tStress;
  } while ((fabs(tStress-sigInit) > tol) && (count <= 100));

  epsInit = tStrain;

  if ((fabs(tStress-sigInit) < tol)) 
    theMaterial->setTrialStrain(epsInit);
  else {
    opserr << "WARNING: InitStressMaterial - could not find initStrain to within tol for material: " << theMaterial->getTag();
    opserr << " wanted sigInit: " << sigInit << " using tStress: " << theMaterial->getStress() << endln;
    return -1;
  }

  return 0;
}

InitStressMaterial::InitStressMaterial(int tag, 
				       UniaxialMaterial &material,
				       double sigini)
  :UniaxialMaterial(tag,MAT_TAG_InitStress), theMaterial(0),
   epsInit(0.0), sigInit(sigini)
{
  theMaterial = material.getCopy();

  if (theMaterial == 0) {
    opserr <<  "InitStressMaterial::InitStressMaterial -- failed to get copy of material\n";
    exit(-1);
  }

  if (this->findInitialStrain() == 0)
    theMaterial->commitState();
}

InitStressMaterial::InitStressMaterial()
  :UniaxialMaterial(0,MAT_TAG_InitStress), theMaterial(0),
   epsInit(0.0)
{

}

InitStressMaterial::~InitStressMaterial()
{
  if (theMaterial)
    delete theMaterial;
}

int 
InitStressMaterial::setTrialStrain(double strain, double strainRate)
{
  return theMaterial->setTrialStrain(strain+epsInit, strainRate);
}

double 
InitStressMaterial::getStress(void)
{
  return theMaterial->getStress();
}

double 
InitStressMaterial::getTangent(void)
{
  return theMaterial->getTangent();  
}

double 
InitStressMaterial::getDampTangent(void)
{
  return theMaterial->getDampTangent();
}

double 
InitStressMaterial::getStrain(void)
{
  return theMaterial->getStrain();
}

double 
InitStressMaterial::getStrainRate(void)
{
  return theMaterial->getStrainRate();
}

int 
InitStressMaterial::commitState(void)
{	
  return theMaterial->commitState();
}

int 
InitStressMaterial::revertToLastCommit(void)
{
  return theMaterial->revertToLastCommit();
}

int 
InitStressMaterial::revertToStart(void)
{
  int res = 0;
  res = theMaterial->revertToStart();
  res += theMaterial->setTrialStrain(epsInit);
  res += theMaterial->commitState();
  return res;
}

UniaxialMaterial *
InitStressMaterial::getCopy(void)
{
  InitStressMaterial *theCopy = 
    new InitStressMaterial(this->getTag(), *theMaterial, sigInit);
        
  return theCopy;
}

int 
InitStressMaterial::sendSelf(int cTag, Channel &theChannel)
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
    opserr << "InitStressMaterial::sendSelf() - failed to send the ID\n";
    return -1;
  }

  static Vector dataVec(1);
  dataVec(0) = epsInit;

  if (theChannel.sendVector(dbTag, cTag, dataVec) < 0) {
    opserr << "InitStressMaterial::sendSelf() - failed to send the Vector\n";
    return -2;
  }

  if (theMaterial->sendSelf(cTag, theChannel) < 0) {
    opserr << "InitStressMaterial::sendSelf() - failed to send the Material\n";
    return -3;
  }

  return 0;
}

int 
InitStressMaterial::recvSelf(int cTag, Channel &theChannel, 
			 FEM_ObjectBroker &theBroker)
{
  int dbTag = this->getDbTag();

  static ID dataID(3);
  if (theChannel.recvID(dbTag, cTag, dataID) < 0) {
    opserr << "InitStressMaterial::recvSelf() - failed to get the ID\n";
    return -1;
  }
  this->setTag(dataID(0));

  // as no way to change material, don't have to check classTag of the material 
  if (theMaterial == 0) {
    int matClassTag = dataID(1);
    theMaterial = theBroker.getNewUniaxialMaterial(matClassTag);
    if (theMaterial == 0) {
      opserr << "InitStressMaterial::recvSelf() - failed to create Material with classTag " 
	   << dataID(0) << endln;
      return -2;
    }
  }
  theMaterial->setDbTag(dataID(2));

  static Vector dataVec(1);
  if (theChannel.recvVector(dbTag, cTag, dataVec) < 0) {
    opserr << "InitStressMaterial::recvSelf() - failed to get the Vector\n";
    return -3;
  }

  epsInit = dataVec(0);
  
  if (theMaterial->recvSelf(cTag, theChannel, theBroker) < 0) {
    opserr << "InitStressMaterial::recvSelf() - failed to get the Material\n";
    return -4;
  }
  return 0;
}

void 
InitStressMaterial::Print(OPS_Stream &s, int flag)
{
    if (flag == OPS_PRINT_PRINTMODEL_MATERIAL) {
        s << "InitStressMaterial tag: " << this->getTag() << endln;
        s << "\tMaterial: " << theMaterial->getTag() << endln;
        s << "\tInitial stress: " << sigInit << endln;
        s << "\tInitial strain: " << epsInit << endln;
    }

    if (flag == OPS_PRINT_PRINTMODEL_JSON) {
        s << "\t\t\t{";
        s << "\"name\": \"" << this->getTag() << "\", ";
        s << "\"type\": \"InitStressMaterial\", ";
        s << "\"Material\": " << theMaterial->getTag() << ", ";
        s << "\"initialStress\": " << sigInit << ", "; 
        s << "\"initialStrain\": " << epsInit << "}";
    }
}

int 
InitStressMaterial::setParameter(const char **argv, int argc, Parameter &param)
{
  if (strcmp(argv[0],"sig0") == 0 || strcmp(argv[0],"f0") == 0 || strcmp(argv[0],"F0") == 0) {
    param.setValue(sigInit);
    return param.addObject(1, this);
  }
  return theMaterial->setParameter(argv, argc, param);
}

int
InitStressMaterial::updateParameter(int parameterID, Information &info)
{
  switch (parameterID) {
  case -1:
    return -1;
  case 1:
    this->sigInit = info.theDouble;
    this->findInitialStrain();
    break;
  }

  return 0;
}

double
InitStressMaterial::getStressSensitivity(int gradIndex, bool conditional)
{
  return theMaterial->getStressSensitivity(gradIndex, conditional);
}

double
InitStressMaterial::getInitialTangentSensitivity(int gradIndex)
{
  return theMaterial->getInitialTangentSensitivity(gradIndex);
}

int
InitStressMaterial::commitSensitivity(double strainGradient, 
				      int gradIndex, int numGrads)
{
  return theMaterial->commitSensitivity(strainGradient, gradIndex, numGrads);
}
