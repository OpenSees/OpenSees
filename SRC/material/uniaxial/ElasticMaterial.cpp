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
                                                                        
// $Revision: 1.10 $
// $Date: 2008-08-26 16:30:55 $
// $Source: /usr/local/cvs/OpenSees/SRC/material/uniaxial/ElasticMaterial.cpp,v $
                                                                        
                                                                        
// Written: fmk 
// Created: 07/98
// Revision: A
//
// Description: This file contains the class implementation for 
// ElasticMaterial. 
//
// What: "@(#) ElasticMaterial.C, revA"

#include <ElasticMaterial.h>
#include <Vector.h>
#include <Channel.h>
#include <Information.h>
#include <Parameter.h>
#include <string.h>

#include <OPS_Globals.h>

#include <elementAPI.h>

void *
OPS_ElasticMaterial(void)
{
  // Pointer to a uniaxial material that will be returned
  UniaxialMaterial *theMaterial = 0;

  if (OPS_GetNumRemainingInputArgs() < 2) {
    opserr << "Invalid #args,  want: uniaxialMaterial Elastic tag? E? <eta?> <Eneg?> ... " << endln;
    return 0;
  }
  
  int iData[1];
  double dData[3];
  int numData = 1;
  if (OPS_GetIntInput(&numData, iData) != 0) {
    opserr << "WARNING invalid tag for uniaxialMaterial Elastic" << endln;
    return 0;
  }

  numData = OPS_GetNumRemainingInputArgs();

  if (numData >= 3) {
    numData = 3;
    if (OPS_GetDoubleInput(&numData, dData) != 0) {
      opserr << "Invalid data for uniaxial Elastic " << iData[0] << endln;
      return 0;	
    }
  } else if (numData >= 2) {
    numData = 2;
    if (OPS_GetDoubleInput(&numData, dData) != 0) {
      opserr << "Invalid data for uniaxial Elastic " << iData[0] << endln;
      return 0;
    }
    dData[2] = dData[0];
  } else {
    numData = 1;
    if (OPS_GetDoubleInput(&numData, dData) != 0) {
      opserr << "Invalid data for uniaxialMaterial Elastic " << iData[0] << endln;
      return 0;	
    }
    dData[1] = 0.0;
    dData[2] = dData[0];
  }

  // Parsing was successful, allocate the material
  theMaterial = new ElasticMaterial(iData[0], dData[0], dData[1], dData[2]);
  if (theMaterial == 0) {
    opserr << "WARNING could not create uniaxialMaterial of type ElasticMaterial" << endln;
    return 0;
  }

  return theMaterial;
}


ElasticMaterial::ElasticMaterial(int tag, double e, double et)
:UniaxialMaterial(tag,MAT_TAG_ElasticMaterial),
 trialStrain(0.0),  trialStrainRate(0.0),
 Epos(e), Eneg(e), eta(et), parameterID(0)
{

}


ElasticMaterial::ElasticMaterial(int tag, double ep, double et, double en)
:UniaxialMaterial(tag,MAT_TAG_ElasticMaterial),
 trialStrain(0.0),  trialStrainRate(0.0),
 Epos(ep), Eneg(en), eta(et), parameterID(0)
{

}


ElasticMaterial::ElasticMaterial()
:UniaxialMaterial(0,MAT_TAG_ElasticMaterial),
 trialStrain(0.0),  trialStrainRate(0.0),
 Epos(0.0), Eneg(0.0), eta(0.0), parameterID(0)
{

}


ElasticMaterial::~ElasticMaterial()
{
  // does nothing
}


int 
ElasticMaterial::setTrialStrain(double strain, double strainRate)
{
    trialStrain     = strain;
    trialStrainRate = strainRate;
    return 0;
}


int 
ElasticMaterial::setTrial(double strain, double &stress, double &tangent, double strainRate)
{
    trialStrain     = strain;
    trialStrainRate = strainRate;

    if (trialStrain >= 0.0) {
        stress = Epos*trialStrain + eta*trialStrainRate;
        tangent = Epos;
    } else {
        stress = Eneg*trialStrain + eta*trialStrainRate;
        tangent = Eneg;
    }

    return 0;
}


double 
ElasticMaterial::getStress(void)
{
    if (trialStrain >= 0.0)
        return Epos*trialStrain + eta*trialStrainRate;
    else
        return Eneg*trialStrain + eta*trialStrainRate;
}


double 
ElasticMaterial::getTangent(void)
{
    if (trialStrain >= 0.0)
        return Epos;
    else 
        return Eneg;
}


double 
ElasticMaterial::getInitialTangent(void)
{
  return Epos;
}


int 
ElasticMaterial::commitState(void)
{
  return 0;
}


int 
ElasticMaterial::revertToLastCommit(void)
{
  return 0;
}


int 
ElasticMaterial::revertToStart(void)
{
    trialStrain      = 0.0;
    trialStrainRate  = 0.0;
    return 0;
}


UniaxialMaterial *
ElasticMaterial::getCopy(void)
{
    ElasticMaterial *theCopy = new ElasticMaterial(this->getTag(),Epos,eta,Eneg);
    theCopy->trialStrain     = trialStrain;
    theCopy->trialStrainRate = trialStrainRate;
    theCopy->parameterID = parameterID;
    return theCopy;
}


int 
ElasticMaterial::sendSelf(int cTag, Channel &theChannel)
{
  int res = 0;
  static Vector data(5);
  data(0) = this->getTag();
  data(1) = Epos;
  data(2) = Eneg;
  data(3) = eta;
  data(4) = parameterID;
  res = theChannel.sendVector(this->getDbTag(), cTag, data);
  if (res < 0) 
    opserr << "ElasticMaterial::sendSelf() - failed to send data" << endln;

  return res;
}


int 
ElasticMaterial::recvSelf(int cTag, Channel &theChannel, 
			  FEM_ObjectBroker &theBroker)
{
  int res = 0;
  static Vector data(5);
  res = theChannel.recvVector(this->getDbTag(), cTag, data);
  
  if (res < 0) {
    opserr << "ElasticMaterial::recvSelf() - failed to receive data" << endln;
    Epos = Eneg = 0; 
    this->setTag(0);      
  }
  else {
    this->setTag(int(data(0)));
    Epos = data(1);
    Eneg = data(2);
    eta  = data(3);
    parameterID = (int)data(4);
  }
    
  return res;
}


void 
ElasticMaterial::Print(OPS_Stream &s, int flag)
{
	if (flag == OPS_PRINT_PRINTMODEL_MATERIAL) {
		s << "ElasticMaterial tag: " << this->getTag() << endln;
		s << "  Epos: " << Epos << " Eneg: " << Eneg << " eta: " << eta << endln;
	}
    
	if (flag == OPS_PRINT_PRINTMODEL_JSON) {
		s << "\t\t\t{";
		s << "\"name\": \"" << this->getTag() << "\", ";
		s << "\"type\": \"ElasticMaterial\", ";
		s << "\"Epos\": " << Epos << ", ";
		s << "\"Eneg\": " << Eneg << ", ";
		s << "\"eta\": " << eta << "}";
	}
}


int
ElasticMaterial::setParameter(const char **argv, int argc, Parameter &param)
{

  if (strcmp(argv[0],"E") == 0) {
    param.setValue(Epos);
    return param.addObject(1, this);
  }
  if (strcmp(argv[0],"Epos") == 0) {
    param.setValue(Epos);
    return param.addObject(2, this);
  }
  if (strcmp(argv[0],"Eneg") == 0) {
    param.setValue(Eneg);
    return param.addObject(3, this);
  }
  if (strcmp(argv[0],"eta") == 0) {
    param.setValue(eta);
    return param.addObject(4, this);
  }
  return -1;
}


int 
ElasticMaterial::updateParameter(int parameterID, Information &info)
{
  switch(parameterID) {
  case 1:
    Epos = info.theDouble;
    Eneg = info.theDouble;
    return 0;
  case 2:
    Epos = info.theDouble;
    return 0;
  case 3:
    Eneg = info.theDouble;
    return 0;
  case 4:
    eta = info.theDouble;
    return 0;
  default:
    return -1;
  }
}


int
ElasticMaterial::activateParameter(int paramID)
{
  parameterID = paramID;

  return 0;
}


double
ElasticMaterial::getStressSensitivity(int gradIndex, bool conditional)
{
  if (parameterID == 1)
    return trialStrain;
  if (parameterID == 2 && trialStrain >= 0.0)
    return trialStrain;
  if (parameterID == 3 && trialStrain < 0.0)
    return trialStrain;
  if (parameterID == 4)
    return trialStrainRate;

  return 0.0;
}


double
ElasticMaterial::getTangentSensitivity(int gradIndex)
{
  if (parameterID == 1)
    return 1.0;
  if (parameterID == 2 && trialStrain >= 0.0)
    return 1.0;
  if (parameterID == 3 && trialStrain < 0.0)
    return 1.0;

  return 0.0;
}


double
ElasticMaterial::getInitialTangentSensitivity(int gradIndex)
{
  if (parameterID == 1)
    return 1.0;
  if (parameterID == 2)
    return 1.0;

  return 0.0;
}


int
ElasticMaterial::commitSensitivity(double strainGradient,
				   int gradIndex, int numGrads)
{
  // Nothing to commit ... path independent
  return 0;
}
