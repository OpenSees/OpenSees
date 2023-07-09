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
                                                                        
// $Revision: 1.3 $
// $Date: 2008-08-26 16:30:23 $
// $Source: /usr/local/cvs/OpenSees/SRC/material/uniaxial/ElasticBDMaterial.cpp,v $
                                                                        
                                                                        
// File: ~/material/ElasticBDMaterial.C
//
// Written: fmk 
// Created: 07/98
// Revision: A
//
// Description: This file contains the class implementation for 
// ElasticBDMaterial. 
//
// What: "@(#) ElasticBDMaterial.C, revA"

#include <ElasticBDMaterial.h>
#include <Vector.h>
#include <Channel.h>
#include <Information.h>
#include <Parameter.h>

#include <string.h>
#include <OPS_Globals.h>

ElasticBDMaterial::ElasticBDMaterial(int tag, double e, double bb, double dd, double aa, double et)
:UniaxialMaterial(tag,MAT_TAG_ElasticBDMaterial),
 trialStrain(0.0),  trialStrainRate(0.0),
 E(e), b(bb), d(dd), a(aa), eta(et), parameterID(0)
{

}

ElasticBDMaterial::ElasticBDMaterial()
:UniaxialMaterial(0,MAT_TAG_ElasticBDMaterial),
 trialStrain(0.0),  trialStrainRate(0.0),
 E(0.0), b(0.0), d(0.0), a(0.0), eta(0.0), parameterID(0)
{

}

ElasticBDMaterial::~ElasticBDMaterial()
{
  // does nothing
}

int 
ElasticBDMaterial::setTrialStrain(double strain, double strainRate)
{
    trialStrain     = strain;
    trialStrainRate = strainRate;
    return 0;
}


int 
ElasticBDMaterial::setTrial(double strain, double &stress, double &tangent, double strainRate)
{
    trialStrain     = strain;
    trialStrainRate = strainRate;

    stress = E*b*d*a*trialStrain + eta*trialStrainRate;
    tangent = E*b*d*a;

    return 0;
}

double 
ElasticBDMaterial::getStress(void)
{
    return E*b*d*a*trialStrain + eta*trialStrainRate;
}


int 
ElasticBDMaterial::commitState(void)
{
    return 0;
}


int 
ElasticBDMaterial::revertToLastCommit(void)
{
    return 0;
}


int 
ElasticBDMaterial::revertToStart(void)
{
    trialStrain      = 0.0;
    trialStrainRate  = 0.0;
    return 0;
}

UniaxialMaterial *
ElasticBDMaterial::getCopy(void)
{
    ElasticBDMaterial *theCopy = new ElasticBDMaterial(this->getTag(),E,b,d,a,eta);
    theCopy->trialStrain     = trialStrain;
    theCopy->trialStrainRate = trialStrainRate;
    return theCopy;
}

int 
ElasticBDMaterial::sendSelf(int cTag, Channel &theChannel)
{
  int res = 0;
  static Vector data(3);
  data(0) = this->getTag();
  data(1) = E;
  data(2) = eta;

  res = theChannel.sendVector(this->getDbTag(), cTag, data);
  if (res < 0) 
    opserr << "ElasticBDMaterial::sendSelf() - failed to send data\n";

  return res;
}

int 
ElasticBDMaterial::recvSelf(int cTag, Channel &theChannel, 
			  FEM_ObjectBroker &theBroker)
{
  int res = 0;
  static Vector data(3);
  res = theChannel.recvVector(this->getDbTag(), cTag, data);
  
  if (res < 0) {
      opserr << "ElasticBDMaterial::recvSelf() - failed to receive data\n";
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
ElasticBDMaterial::Print(OPS_Stream &s, int flag)
{
    s << "ElasticBD tag: " << this->getTag() << endln;
    s << "  E: " << E << " eta: " << eta << endln;
    s << "  b: " << b << endln;
    s << "  d: " << d << endln;
    s << "  a: " << a << endln;
}

int
ElasticBDMaterial::setParameter(const char **argv, int argc, Parameter &param)
{

  if (strcmp(argv[0],"E") == 0)
    return param.addObject(1, this);

  else if (strcmp(argv[0],"eta") == 0)
    return param.addObject(2, this);

  else if (strcmp(argv[0],"b") == 0)
    return param.addObject(3, this);

  else if (strcmp(argv[0],"d") == 0)
    return param.addObject(4, this);

  else if (strcmp(argv[0],"a") == 0)
    return param.addObject(5, this);

  return -1;
}

int 
ElasticBDMaterial::updateParameter(int parameterID, Information &info)
{
  switch(parameterID) {
  case 1:
    E = info.theDouble;
    return 0;
  case 2:
    eta = info.theDouble;
    return 0;
  case 3:
    b = info.theDouble;
    return 0;
  case 4:
    d = info.theDouble;
    return 0;
  case 5:
    a = info.theDouble;
    return 0;
  default:
    return -1;
  }
}

int
ElasticBDMaterial::activateParameter(int paramID)
{
  parameterID = paramID;

  return 0;
}

double
ElasticBDMaterial::getStressSensitivity(int gradIndex, bool conditional)
{
  if (parameterID == 1)
    return b*d*a*trialStrain;
  else if (parameterID == 2)
    return trialStrainRate;
  else if (parameterID == 3)
    return E*d*a*trialStrain;
  else if (parameterID == 4)
    return E*b*a*trialStrain;
  else if (parameterID == 5)
    return E*b*d*trialStrain;
  else
    return 0.0;
}

double
ElasticBDMaterial::getInitialTangentSensitivity(int gradIndex)
{
  return 0.0;
}

int
ElasticBDMaterial::commitSensitivity(double strainGradient,
				     int gradIndex, int numGrads)
{
  // Nothing to commit ... path independent
  return 0;
}
