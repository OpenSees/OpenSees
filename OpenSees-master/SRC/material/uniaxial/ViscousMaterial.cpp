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
// $Date: 2009-03-05 00:52:38 $
// $Source: /usr/local/cvs/OpenSees/SRC/material/uniaxial/ViscousMaterial.cpp,v $
                                                                        
// Written: Mehrdad Sasani 
// Created: June 2000
// Revision: A
//
// Description: This file contains the class implementation for 
// ViscousMaterial. 

#include <math.h>

#include <ViscousMaterial.h>
#include <Vector.h>
#include <Channel.h>
#include <Information.h>
#include <Parameter.h>
#include <string.h>

#include <OPS_Globals.h>


#include <elementAPI.h>

void *
OPS_ViscousMaterial(void)
{
  // Pointer to a uniaxial material that will be returned
  UniaxialMaterial *theMaterial = 0;

  int numArgs = OPS_GetNumRemainingInputArgs();
  if (numArgs < 3 || numArgs > 4) {
    opserr << "Invalid #args,  want: uniaxialMaterial Viscous tag? C? alpha? <minVel?> ... " << endln;
    return 0;
  }
  
  int iData[1];
  double dData[3];
  dData[2] = 1.0e-11; // setting default minVel

  int numData = 1;
  if (OPS_GetIntInput(&numData, iData) != 0) {
    opserr << "WARNING invalid tag for uniaxialMaterial Viscous" << endln;
    return 0;
  }

  numData = numArgs-1;
  if (OPS_GetDoubleInput(&numData, dData) != 0) {
    opserr << "Invalid data for uniaxial Viscous " << iData[0] << endln;
    return 0;	
  }

  // Parsing was successful, allocate the material
  theMaterial = new ViscousMaterial(iData[0], dData[0], dData[1], dData[2]);
  if (theMaterial == 0) {
    opserr << "WARNING could not create uniaxialMaterial of type Viscous\n";
    return 0;
  }

  return theMaterial;
}


ViscousMaterial::ViscousMaterial(int tag, double c, double a, double minV)
:UniaxialMaterial(tag,MAT_TAG_Viscous),
 trialRate(0.0), trialStrain(0.0), C(c), Alpha(a), minVel(minV), commitStrain(0.0), commitRate(0.0)
{
    if (Alpha < 0.0) {
      opserr << "ViscousMaterial::ViscousMaterial -- Alpha < 0.0, setting to 1.0\n";
      Alpha = 1.0;
    }
    
    minVel = fabs(minVel);
    if (minVel == 0.0) {
      opserr << "ViscousMaterial::ViscousMaterial -- minVel == 0.0, setting to 1.0e-21\n";
      minVel = 1.0e-21;
    }
}

ViscousMaterial::ViscousMaterial()
:UniaxialMaterial(0,MAT_TAG_Viscous),
 trialRate(0.0), trialStrain(0.0), C(0.0), Alpha(0.0), minVel(1e-11), commitStrain(0.0), commitRate(0.0)
{

}

ViscousMaterial::~ViscousMaterial()
{
  // does nothing
}

int 
ViscousMaterial::setTrialStrain(double strain, double strainRate)
{
  trialStrain = strain;
  trialRate = strainRate;

    return 0;
}

double 
ViscousMaterial::getStress(void)
{
    double stress = 0.0;
    double absRate = fabs(trialRate);

    if (absRate > minVel)
      stress = C*pow(absRate, Alpha);
    else
      stress = C*pow(minVel, Alpha);

    stress = C*pow(absRate, Alpha);    

    if (trialRate < 0.0)
        return -stress;
    else
        return  stress;
}

double 
ViscousMaterial::getTangent(void)
{
    return 0.0;
}

double 
ViscousMaterial::getInitialTangent(void)
{
    return 0.0;
}

double
ViscousMaterial::getDampTangent(void)
{
     double absRate = fabs(trialRate);

    if (absRate < minVel)
      return Alpha*C*pow(minVel,Alpha-1.0);
    else
      return Alpha*C*pow(absRate,Alpha-1.0);	
}


double 
ViscousMaterial::getStrain(void)
{
    return trialStrain;
}

double 
ViscousMaterial::getStrainRate(void)
{
    return trialRate;
}

int 
ViscousMaterial::commitState(void)
{
  commitStrain = trialStrain;
  commitRate = trialRate;
  return 0;
}

int 
ViscousMaterial::revertToLastCommit(void)
{
  trialStrain = commitStrain;
  trialRate = commitRate;
  return 0;
}

int 
ViscousMaterial::revertToStart(void)
{
    trialRate = 0.0;
    trialStrain = 0.0;
    commitStrain = 0.0;
    commitRate = 0.0;
    return 0;
}

UniaxialMaterial *
ViscousMaterial::getCopy(void)
{
  ViscousMaterial *theCopy = new ViscousMaterial(this->getTag(),C,Alpha, minVel);

    theCopy->trialRate = trialRate;

    return theCopy;
}

int 
ViscousMaterial::sendSelf(int cTag, Channel &theChannel)
{
  int res = 0;
  static Vector data(6);
  data(0) = this->getTag();
  data(1) = C;
  data(2) = Alpha;
  data(3) = trialRate;
  data(4) = minVel;
  data(5) = trialStrain;
  res = theChannel.sendVector(this->getDbTag(), cTag, data);
  if (res < 0) 
    opserr << "ViscousMaterial::sendSelf() - failed to send data\n";

  return res;
}

int 
ViscousMaterial::recvSelf(int cTag, Channel &theChannel, 
			       FEM_ObjectBroker &theBroker)
{
  int res = 0;
  static Vector data(6);
  res = theChannel.recvVector(this->getDbTag(), cTag, data);
  
  if (res < 0) {
      opserr << "ViscousMaterial::recvSelf() - failed to receive data\n";
      C = 0; 
      this->setTag(0);      
  }
  else {
    this->setTag((int)data(0));
    C = data(1);
    Alpha = data(2);
    trialRate = data(3);
    minVel = data(4);
    trialStrain = data(5);
    commitStrain = trialStrain;
    commitRate = trialRate;
  }
    
  return res;
}

void 
ViscousMaterial::Print(OPS_Stream &s, int flag)
{
    s << "Viscous tag: " << this->getTag() << endln;
    s << "  C: " << C << endln;
    s << "  Alpha: " << Alpha << endln;
    s << "  minVel: " << minVel << endln;
}



int
ViscousMaterial::setParameter(const char **argv, int argc, Parameter &param)
{

  if (strcmp(argv[0],"C") == 0) {
    param.setValue(C);
    return param.addObject(1, this);
  }
  if (strcmp(argv[0],"Alpha") == 0) {
    param.setValue(Alpha);
    return param.addObject(2, this);
  }
  if (strcmp(argv[0],"minVel") == 0) {
    param.setValue(minVel);
    return param.addObject(3, this);
  }
  return -1;
}


int 
ViscousMaterial::updateParameter(int parameterID, Information &info)
{
  switch(parameterID) {
  case 1:
    C = info.theDouble;
    return 0;
  case 2:
    Alpha = info.theDouble;
    return 0;
  case 3:
    minVel = info.theDouble;
    return 0;
  default:
    return -1;
  }
}


