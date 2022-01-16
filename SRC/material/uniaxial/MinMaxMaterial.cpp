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
// $Source: /usr/local/cvs/OpenSees/SRC/material/uniaxial/MinMaxMaterial.cpp,v $

// Written: MHS
// Created: Aug 2001
//
// Description: This file contains the class definition for 
// MinMaxMaterial.  MinMaxMaterial wraps a UniaxialMaterial
// and imposes min and max strain limits.

#include <stdlib.h>
#include <string.h>

#include <MinMaxMaterial.h>
#include <ID.h>
#include <Channel.h>
#include <FEM_ObjectBroker.h>

#include <OPS_Globals.h>

#include <elementAPI.h>
#define OPS_Export 

OPS_Export void * OPS_ADD_RUNTIME_VPV(OPS_MinMaxMaterial)
{
  // Pointer to a uniaxial material that will be returned
  UniaxialMaterial *theMaterial = 0;
  UniaxialMaterial *theOtherMaterial = 0;
  double minStrain = -1.0e16;
  double maxStrain = 1.0e16;
  int    iData[2];

  int argc = OPS_GetNumRemainingInputArgs();
  if (argc < 2) {
    opserr << "WARNING invalid uniaxialMaterial MinMaxMaterial $tag $otherTag <-min $minStrain> <-max $maxStrain>" << endln;
    return 0;
  }

  int numData = 2;
  if (OPS_GetIntInput(&numData, iData) != 0) {
    opserr << "WARNING invalid uniaxialMaterial MinMaxMaterial $tag $otherTag" << endln;
    return 0;
  }

  theOtherMaterial = OPS_GetUniaxialMaterial(iData[1]);
  if (theOtherMaterial == 0) {
    opserr << "WARNING invalid otherTag uniaxialMaterial MinMax tag: " << iData[0] << endln;
    return 0;	
  }

  argc = OPS_GetNumRemainingInputArgs();  
  while (argc > 1) {
    //char argvLoc[10];
    const char *argvLoc = OPS_GetString();
    /*    if (OPS_GetString(argvLoc, 10) != 0) {
      opserr << "WARNING invalid string option uniaxialMaterial MinMax tag: " << iData[0] << endln;
      return 0;
    }
    */
    numData = 1;

    if ((strcmp(argvLoc, "-min") == 0) || (strcmp(argvLoc, "-Min") == 0) || (strcmp(argvLoc, "-MIN") == 0)) {
      if (OPS_GetDouble(&numData, &minStrain) != 0) {      
	opserr << "WARNING invalid min value  uniaxialMaterial MinMax tag: " << iData[0] << endln;	
	return 0;
      }
    } else if ((strcmp(argvLoc, "-max") == 0) || (strcmp(argvLoc, "-Max") == 0) || (strcmp(argvLoc, "-MAX") == 0)) {
      if (OPS_GetDouble(&numData, &maxStrain) != 0) {      
	opserr << "WARNING invalid min value  uniaxialMaterial MinMax tag: " << iData[0] << endln;  
	return 0;
      }
    } else {
      opserr << "WARNING invalid option:" << argvLoc << " uniaxialMaterial MinMax tag: " << iData[0] << endln;  
      return 0;
    }
    
    argc = OPS_GetNumRemainingInputArgs();
  }

  // Parsing was successful, allocate the material
  theMaterial = new MinMaxMaterial(iData[0], *theOtherMaterial, minStrain, maxStrain);

  if (theMaterial == 0) {
    opserr << "WARNING could not create uniaxialMaterial of type MinMaxMaterial\n";
    return 0;
  }

  return theMaterial;
}

MinMaxMaterial::MinMaxMaterial(int tag, UniaxialMaterial &material,
			       double min, double max)
  :UniaxialMaterial(tag,MAT_TAG_MinMax), theMaterial(0),
   minStrain(min), maxStrain(max), Tfailed(false), Cfailed(false)
{
  theMaterial = material.getCopy();

  if (theMaterial == 0) {
    opserr <<  "MinMaxMaterial::MinMaxMaterial -- failed to get copy of material\n";
    exit(-1);
  }
}

MinMaxMaterial::MinMaxMaterial()
  :UniaxialMaterial(0,MAT_TAG_MinMax), theMaterial(0),
   minStrain(0.0), maxStrain(0.0), Tfailed(false), Cfailed(false)
{

}

MinMaxMaterial::~MinMaxMaterial()
{
  if (theMaterial)
    delete theMaterial;
}

int 
MinMaxMaterial::setTrialStrain(double strain, double strainRate)
{
  if (Cfailed)
    return 0;
  
  if (strain >= maxStrain || strain <= minStrain) {
    Tfailed = true;
    return 0;
  }
  else {
    Tfailed = false;
    return theMaterial->setTrialStrain(strain, strainRate);
  }
}


int 
MinMaxMaterial::setTrialStrain(double strain, double temp, double strainRate)
{
  if (Cfailed)
    return 0;
  
  if (strain >= maxStrain || strain <= minStrain) {
    Tfailed = true;
    return 0;
  }
  else {
    Tfailed = false;
    return theMaterial->setTrialStrain(strain, temp, strainRate);
  }
}


double 
MinMaxMaterial::getStress(void)
{
  if (Tfailed)
    return 0.0;
  else
    return theMaterial->getStress();
}

double 
MinMaxMaterial::getTangent(void)
{
  if (Tfailed)
    //return 0.0;
    return 1.0e-8*theMaterial->getInitialTangent();
  else
    return theMaterial->getTangent();
}

double 
MinMaxMaterial::getDampTangent(void)
{
  if (Tfailed)
    return 0.0;
  else
    return theMaterial->getDampTangent();
}



double 
MinMaxMaterial::getStrain(void)
{
  return theMaterial->getStrain();
}

double 
MinMaxMaterial::getStrainRate(void)
{
  return theMaterial->getStrainRate();
}

int 
MinMaxMaterial::commitState(void)
{	
  Cfailed = Tfailed;

  // Check if failed at current step
  if (Tfailed)
    return 0;
  else
    return theMaterial->commitState();
}

int 
MinMaxMaterial::revertToLastCommit(void)
{
  // Check if failed at last step
  if (Cfailed)
    return 0;
  else
    return theMaterial->revertToLastCommit();
}

int 
MinMaxMaterial::revertToStart(void)
{
  Cfailed = false;
  Tfailed = false;
  
  return theMaterial->revertToStart();
}

UniaxialMaterial *
MinMaxMaterial::getCopy(void)
{
  MinMaxMaterial *theCopy = 
    new MinMaxMaterial(this->getTag(), *theMaterial, minStrain, maxStrain);
        
  theCopy->Cfailed = Cfailed;
  theCopy->Tfailed = Tfailed;
  
  return theCopy;
}

int 
MinMaxMaterial::sendSelf(int cTag, Channel &theChannel)
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
    opserr << "MinMaxMaterial::sendSelf() - failed to send the ID\n";
    return -1;
  }

  static Vector dataVec(3);
  dataVec(0) = minStrain;
  dataVec(1) = maxStrain;
  if (Cfailed == true)
    dataVec(2) = 1.0;
  else
    dataVec(2) = 0.0;

  if (theChannel.sendVector(dbTag, cTag, dataVec) < 0) {
    opserr << "MinMaxMaterial::sendSelf() - failed to send the Vector\n";
    return -2;
  }

  if (theMaterial->sendSelf(cTag, theChannel) < 0) {
    opserr << "MinMaxMaterial::sendSelf() - failed to send the Material\n";
    return -3;
  }

  return 0;
}

int 
MinMaxMaterial::recvSelf(int cTag, Channel &theChannel, 
			 FEM_ObjectBroker &theBroker)
{
  int dbTag = this->getDbTag();

  static ID dataID(3);
  if (theChannel.recvID(dbTag, cTag, dataID) < 0) {
    opserr << "MinMaxMaterial::recvSelf() - failed to get the ID\n";
    return -1;
  }
  this->setTag(int(dataID(0)));

  // as no way to change material, don't have to check classTag of the material 
  if (theMaterial == 0) {
    int matClassTag = int(dataID(1));
    theMaterial = theBroker.getNewUniaxialMaterial(matClassTag);
    if (theMaterial == 0) {
      opserr << "MinMaxMaterial::recvSelf() - failed to create Material with classTag " 
	   << dataID(0) << endln;
      return -2;
    }
  }
  theMaterial->setDbTag(dataID(2));

  static Vector dataVec(3);
  if (theChannel.recvVector(dbTag, cTag, dataVec) < 0) {
    opserr << "MinMaxMaterial::recvSelf() - failed to get the Vector\n";
    return -3;
  }

  minStrain = dataVec(0);
  maxStrain = dataVec(1);
  
  if (dataVec(2) == 1.0)
    Cfailed = true;
  else
    Cfailed = false;

  Tfailed = Cfailed;

  if (theMaterial->recvSelf(cTag, theChannel, theBroker) < 0) {
    opserr << "MinMaxMaterial::recvSelf() - failed to get the Material\n";
    return -4;
  }
  return 0;
}

void 
MinMaxMaterial::Print(OPS_Stream &s, int flag)
{
    if (flag == OPS_PRINT_PRINTMODEL_MATERIAL) {
        s << "MinMaxMaterial, tag: " << this->getTag() << endln;
        s << "  material: " << theMaterial->getTag() << endln;
        s << "  min strain: " << minStrain << endln;
        s << "  max strain: " << maxStrain << endln;
    }
    
    if (flag == OPS_PRINT_PRINTMODEL_JSON) {
        s << "\t\t\t{";
        s << "\"name\": \"" << this->getTag() << "\", ";
        s << "\"type\": \"MinMaxMaterial\", ";
        s << "\"material\": \"" << theMaterial->getTag() << "\", ";
        s << "\"epsMin\": " << minStrain << ", ";
        s << "\"epsMax\": " << maxStrain << "}";
    }
}

int
MinMaxMaterial::setParameter(const char **argv, int argc, Parameter &param)
{
  //
  // I suppose epsMin and epsMax could be parameters, but not for now -- MHS
  //
  return theMaterial->setParameter(argv, argc, param);
}
  
int
MinMaxMaterial::updateParameter(int parameterID, Information &info)
{
  return 0;
}

double
MinMaxMaterial::getStressSensitivity(int gradIndex, bool conditional)
{
  if (Cfailed)
    return 0.0;
  else 
    return theMaterial->getStressSensitivity(gradIndex, conditional);
}

double
MinMaxMaterial::getStrainSensitivity(int gradIndex)
{
  return theMaterial->getStrainSensitivity(gradIndex);
}

double
MinMaxMaterial::getInitialTangentSensitivity(int gradIndex)
{
  return theMaterial->getInitialTangentSensitivity(gradIndex);
}

double
MinMaxMaterial::getDampTangentSensitivity(int gradIndex)
{
  return theMaterial->getDampTangentSensitivity(gradIndex);
}

double
MinMaxMaterial::getRhoSensitivity(int gradIndex)
{
  return theMaterial->getRhoSensitivity(gradIndex);
}

int   
MinMaxMaterial::commitSensitivity(double strainGradient,
				  int gradIndex, int numGrads)
{
  if (Cfailed)
    return 0;
  else
    return theMaterial->commitSensitivity(strainGradient, gradIndex, numGrads);
}
