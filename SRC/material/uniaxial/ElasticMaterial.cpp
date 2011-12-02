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
                                                                        
// $Revision: 1.7 $
// $Date: 2007-05-22 18:18:14 $
// $Source: /usr/local/cvs/OpenSees/SRC/material/uniaxial/ElasticMaterial.cpp,v $
                                                                        
                                                                        
// File: ~/material/ElasticMaterial.C
//
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

#include <OPS_Globals.h>

ElasticMaterial::ElasticMaterial(int tag, double e, double et)
:UniaxialMaterial(tag,MAT_TAG_ElasticMaterial),
 trialStrain(0.0),  trialStrainRate(0.0),
 E(e), eta(et), parameterID(0)
{

}

ElasticMaterial::ElasticMaterial()
:UniaxialMaterial(0,MAT_TAG_ElasticMaterial),
 trialStrain(0.0),  trialStrainRate(0.0),
 E(0.0), eta(0.0), parameterID(0)
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

    stress = E*trialStrain + eta*trialStrainRate;
    tangent = E;

    return 0;
}

double 
ElasticMaterial::getStress(void)
{
    return E*trialStrain + eta*trialStrainRate;
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
    ElasticMaterial *theCopy = new ElasticMaterial(this->getTag(),E,eta);
    theCopy->trialStrain     = trialStrain;
    theCopy->trialStrainRate = trialStrainRate;
    return theCopy;
}

int 
ElasticMaterial::sendSelf(int cTag, Channel &theChannel)
{
  int res = 0;
  static Vector data(3);
  data(0) = this->getTag();
  data(1) = E;
  data(2) = eta;

  res = theChannel.sendVector(this->getDbTag(), cTag, data);
  if (res < 0) 
    opserr << "ElasticMaterial::sendSelf() - failed to send data\n";

  return res;
}

int 
ElasticMaterial::recvSelf(int cTag, Channel &theChannel, 
			  FEM_ObjectBroker &theBroker)
{
  int res = 0;
  static Vector data(3);
  res = theChannel.recvVector(this->getDbTag(), cTag, data);
  
  if (res < 0) {
      opserr << "ElasticMaterial::recvSelf() - failed to receive data\n";
      E = 0; 
      this->setTag(0);      
  }
  else {
    this->setTag(data(0));
    E   = data(1);
    eta = data(2);
  }
    
  return res;
}

void 
ElasticMaterial::Print(OPS_Stream &s, int flag)
{
    s << "Elastic tag: " << this->getTag() << endln;
    s << "  E: " << E << " eta: " << eta << endln;
}

int
ElasticMaterial::setParameter(const char **argv, int argc, Parameter &param)
{
  if (argc < 1)
    return 0;

  if (strcmp(argv[0],"E") == 0)
    return param.addObject(1, this);

  else if (strcmp(argv[0],"eta") == 0)
    return param.addObject(2, this);

  else
    return 0;
}

int 
ElasticMaterial::updateParameter(int parameterID, Information &info)
{
  switch(parameterID) {
  case 1:
    E = info.theDouble;
    return 0;
  case 2:
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
ElasticMaterial::getStressSensitivity(int gradNumber, bool conditional)
{
  if (parameterID == 1)
    return trialStrain;
  else if (parameterID == 2)
    return trialStrainRate;
  else
    return 0.0;
}

double
ElasticMaterial::getInitialTangentSensitivity(int gradNumber)
{
  return 0.0;
}

int
ElasticMaterial::commitSensitivity(double strainGradient,
				   int gradNumber, int numGrads)
{
  // Nothing to commit ... path independent
  return 0.0;
}
