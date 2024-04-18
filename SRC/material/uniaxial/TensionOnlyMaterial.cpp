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
// $Source: /usr/local/cvs/OpenSees/SRC/material/uniaxial/TensionOnlyMaterial.cpp,v $

// Written: MHS
// Created: Aug 2001
//
// Description: This file contains the class definition for 
// TensionOnlyMaterial.  TensionOnlyMaterial wraps a UniaxialMaterial
// and imposes min and max strain limits.

#include <stdlib.h>
#include <string.h>

#include <TensionOnlyMaterial.h>
#include <ID.h>
#include <Channel.h>
#include <FEM_ObjectBroker.h>

#include <OPS_Globals.h>

#include <elementAPI.h>
#define OPS_Export 

OPS_Export void *
OPS_TensionOnlyMaterial(void)
{
  // Pointer to a uniaxial material that will be returned
  UniaxialMaterial *theMaterial = 0;
  UniaxialMaterial *theOtherMaterial = 0;
  double minStrain = -1.0e16;
  double maxStrain = 1.0e16;
  int    iData[2];

  int argc = OPS_GetNumRemainingInputArgs();
  if (argc < 2) {
    opserr << "WARNING invalid uniaxialMaterial TensionOnly $tag $otherTag <-min $minStrain> <-max $maxStrain>" << endln;
    return 0;
  }

  int numData = 2;
  if (OPS_GetIntInput(&numData, iData) != 0) {
    opserr << "WARNING invalid uniaxialMaterial TensionOnly $tag $otherTag" << endln;
    return 0;
  }

  theOtherMaterial = OPS_GetUniaxialMaterial(iData[1]);
  if (theOtherMaterial == 0) {
    opserr << "WARNING invalid otherTag uniaxialMaterial TensionOnly tag: " << iData[0] << endln;
    return 0;	
  }

  argc = OPS_GetNumRemainingInputArgs();  
  while (argc > 1) {
    //char argvLoc[10];
    const char *argvLoc = OPS_GetString();
    /*    if (OPS_GetString(argvLoc, 10) != 0) {
      opserr << "WARNING invalid string option uniaxialMaterial TensionOnly tag: " << iData[0] << endln;
      return 0;
    }
    */
    numData = 1;

    if ((strcmp(argvLoc, "-min") == 0) || (strcmp(argvLoc, "-Min") == 0) || (strcmp(argvLoc, "-MIN") == 0)) {
      if (OPS_GetDouble(&numData, &minStrain) != 0) {      
	opserr << "WARNING invalid min value  uniaxialMaterial TensionOnly tag: " << iData[0] << endln;	
	return 0;
      }
    } else if ((strcmp(argvLoc, "-max") == 0) || (strcmp(argvLoc, "-Max") == 0) || (strcmp(argvLoc, "-MAX") == 0)) {
      if (OPS_GetDouble(&numData, &maxStrain) != 0) {      
	opserr << "WARNING invalid min value  uniaxialMaterial TensionOnly tag: " << iData[0] << endln;  
	return 0;
      }
    } else {
      opserr << "WARNING invalid option:" << argvLoc << " uniaxialMaterial TensionOnly tag: " << iData[0] << endln;  
      return 0;
    }
    
    argc = OPS_GetNumRemainingInputArgs();
  }

  // Parsing was successful, allocate the material
  theMaterial = new TensionOnlyMaterial(iData[0], *theOtherMaterial);

  if (theMaterial == 0) {
    opserr << "WARNING could not create uniaxialMaterial of type TensionOnlyMaterial\n";
    return 0;
  }

  return theMaterial;
}

TensionOnlyMaterial::TensionOnlyMaterial(int tag, UniaxialMaterial &material)
  :UniaxialMaterial(tag,MAT_TAG_TensionOnly), theMaterial(0)
{
  theMaterial = material.getCopy();

  if (theMaterial == 0) {
    opserr <<  "TensionOnlyMaterial::TensionOnlyMaterial -- failed to get copy of material\n";
    exit(-1);
  }
}

TensionOnlyMaterial::TensionOnlyMaterial()
  :UniaxialMaterial(0,MAT_TAG_TensionOnly), theMaterial(0)
{

}

TensionOnlyMaterial::~TensionOnlyMaterial()
{
  if (theMaterial)
    delete theMaterial;
}

int 
TensionOnlyMaterial::setTrialStrain(double strain, double strainRate)
{
  return theMaterial->setTrialStrain(strain, strainRate);
}


int 
TensionOnlyMaterial::setTrialStrain(double strain, double temp, double strainRate)
{
  return theMaterial->setTrialStrain(strain, temp, strainRate);
}


double 
TensionOnlyMaterial::getStress(void)
{
  double f = theMaterial->getStress();
  if (f < 0.0)
    return 0.0001*f;
  else
    return f;
}

double 
TensionOnlyMaterial::getTangent(void)
{
  double E = theMaterial->getTangent();
  double f = theMaterial->getStress();
  if (f < 0.0)
    return 0.0001*E;
  else
    return E;
}

double 
TensionOnlyMaterial::getDampTangent(void)
{
  double D = theMaterial->getDampTangent();
  double f = theMaterial->getStress();
  if (f < 0.0)
    return 0.0001*D;
  else
    return D;
}

double 
TensionOnlyMaterial::getStrain(void)
{
  return theMaterial->getStrain();
}

double 
TensionOnlyMaterial::getStrainRate(void)
{
  return theMaterial->getStrainRate();
}

int 
TensionOnlyMaterial::commitState(void)
{
  double f = theMaterial->getStress();
  if (f >= 0.0)
    return theMaterial->commitState();
  else
    return 0;
}

int 
TensionOnlyMaterial::revertToLastCommit(void)
{
  return theMaterial->revertToLastCommit();
}

int 
TensionOnlyMaterial::revertToStart(void)
{
  return theMaterial->revertToStart();
}

UniaxialMaterial *
TensionOnlyMaterial::getCopy(void)
{
  TensionOnlyMaterial *theCopy = 
    new TensionOnlyMaterial(this->getTag(), *theMaterial);
        
  return theCopy;
}

int 
TensionOnlyMaterial::sendSelf(int cTag, Channel &theChannel)
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
    opserr << "TensionOnlyMaterial::sendSelf() - failed to send the ID\n";
    return -1;
  }

  if (theMaterial->sendSelf(cTag, theChannel) < 0) {
    opserr << "TensionOnlyMaterial::sendSelf() - failed to send the Material\n";
    return -3;
  }

  return 0;
}

int 
TensionOnlyMaterial::recvSelf(int cTag, Channel &theChannel, 
			 FEM_ObjectBroker &theBroker)
{
  int dbTag = this->getDbTag();

  static ID dataID(3);
  if (theChannel.recvID(dbTag, cTag, dataID) < 0) {
    opserr << "TensionOnlyMaterial::recvSelf() - failed to get the ID\n";
    return -1;
  }
  this->setTag(int(dataID(0)));

  // as no way to change material, don't have to check classTag of the material 
  if (theMaterial == 0) {
    int matClassTag = int(dataID(1));
    theMaterial = theBroker.getNewUniaxialMaterial(matClassTag);
    if (theMaterial == 0) {
      opserr << "TensionOnlyMaterial::recvSelf() - failed to create Material with classTag " 
	   << dataID(0) << endln;
      return -2;
    }
  }
  theMaterial->setDbTag(dataID(2));

  if (theMaterial->recvSelf(cTag, theChannel, theBroker) < 0) {
    opserr << "TensionOnlyMaterial::recvSelf() - failed to get the Material\n";
    return -4;
  }
  return 0;
}

void 
TensionOnlyMaterial::Print(OPS_Stream &s, int flag)
{
    if (flag == OPS_PRINT_PRINTMODEL_MATERIAL) {
        s << "TensionOnlyMaterial, tag: " << this->getTag() << endln;
        s << "  material: " << theMaterial->getTag() << endln;
    }
    
    if (flag == OPS_PRINT_PRINTMODEL_JSON) {
        s << "\t\t\t{";
        s << "\"name\": \"" << this->getTag() << "\", ";
        s << "\"type\": \"TensionOnlMaterial\", ";
        s << "\"material\": \"" << theMaterial->getTag() << "\", ";
    }
}

int
TensionOnlyMaterial::setParameter(const char **argv, int argc, Parameter &param)
{
  return theMaterial->setParameter(argv, argc, param);
}

int
TensionOnlyMaterial::updateParameter(int parameterID, Information &info)
{
  return 0;
}

double
TensionOnlyMaterial::getStressSensitivity(int gradIndex, bool conditional)
{
  double f = theMaterial->getStress();
  if (f < 0.0)
    return 0.0;
  else 
    return theMaterial->getStressSensitivity(gradIndex, conditional);
}

double
TensionOnlyMaterial::getStrainSensitivity(int gradIndex)
{
  return theMaterial->getStrainSensitivity(gradIndex);
}

double
TensionOnlyMaterial::getInitialTangentSensitivity(int gradIndex)
{
  return theMaterial->getInitialTangentSensitivity(gradIndex);
}

double
TensionOnlyMaterial::getDampTangentSensitivity(int gradIndex)
{
  double f = theMaterial->getStress();
  if (f < 0.0)
    return 0.0;
  else
    return theMaterial->getDampTangentSensitivity(gradIndex);
}

double
TensionOnlyMaterial::getRhoSensitivity(int gradIndex)
{
  return theMaterial->getRhoSensitivity(gradIndex);
}

int   
TensionOnlyMaterial::commitSensitivity(double strainGradient,
				       int gradIndex, int numGrads)
{
  double f = theMaterial->getStress();
  if (f >= 0.0)
    return theMaterial->commitSensitivity(strainGradient, gradIndex, numGrads);
  else
    return 0;
}
