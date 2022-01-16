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
                                                                        
// $Revision: 1.5 $
// $Date: 2008-04-14 22:01:20 $
// $Source: /usr/local/cvs/OpenSees/SRC/material/uniaxial/Elastic2Material.cpp,v $
                                                                        
                                                                        
// Written: ZHY
//
// Description: This file contains the class implementation for 
// Elastic2Material. 
//
// What: "@(#) Elastic2Material.C, revA"

#include <Elastic2Material.h>
#include <Vector.h>
#include <Channel.h>
#include <Parameter.h>
#include <string.h>

#include <OPS_Globals.h>
#include <elementAPI.h>

int Elastic2Material::zeroE = 0;

void * OPS_ADD_RUNTIME_VPV(OPS_Elastic2)
{
    int argc = OPS_GetNumRemainingInputArgs() + 2;
    
    if (argc < 4 || argc > 5) {
	opserr << "WARNING invalid number of arguments\n";
	opserr << "Want: uniaxialMaterial Elastic tag? E? <eta?>\n";
	return 0;
    }    

    int tag;
    int numdata = 1;
    if (OPS_GetIntInput(&numdata, &tag) < 0) {
	opserr << "WARNING invalid uniaxialMaterial Elastic tag\n";
	return 0;
    }

    // E, eta
    double data[2] = {0.0, 0.0};
    numdata = OPS_GetNumRemainingInputArgs();
    if (numdata > 2) numdata = 2;
    if (OPS_GetDoubleInput(&numdata, data) < 0) {
	opserr << "WARNING invalid double inputs\n";
	return 0;
    }
	
    // Parsing was successful, allocate the material
    return new Elastic2Material(tag, data[0], data[1]);
}

Elastic2Material::Elastic2Material(int tag, double e, double et)
:UniaxialMaterial(tag,MAT_TAG_Elastic2Material),
 trialStrain(0.0),  trialStrainRate(0.0),
  E(e), eta(et), initialStrain(99999.99)
{

}

Elastic2Material::Elastic2Material()
:UniaxialMaterial(0,MAT_TAG_Elastic2Material),
 trialStrain(0.0),  trialStrainRate(0.0),
  E(0.0), eta(0.0), initialStrain(99999.99)
{

}

Elastic2Material::~Elastic2Material()
{
  // does nothing
}

int 
Elastic2Material::setTrialStrain(double strain, double strainRate)
{
  if (initialStrain == 99999.99) initialStrain = strain;
  trialStrain     = strain - initialStrain;
  trialStrainRate = strainRate;
  return 0;
}


int 
Elastic2Material::setTrial(double strain, double &stress, double &tangent, double strainRate)
{
  if (initialStrain == 99999.99) initialStrain = strain;
  trialStrain     = strain - initialStrain;
  trialStrainRate = strainRate;
  
  stress = E*trialStrain + eta*trialStrainRate;
  tangent = E;
  
  if (zeroE==1) {
    stress = eta*trialStrainRate;
    tangent = 0;
  }
  
  return 0;
}

double 
Elastic2Material::getStress(void)
{
  if (zeroE==1) return eta*trialStrainRate;
  return E*trialStrain + eta*trialStrainRate;
}


int 
Elastic2Material::commitState(void)
{
    return 0;
}


int 
Elastic2Material::revertToLastCommit(void)
{
    return 0;
}


int 
Elastic2Material::revertToStart(void)
{
    initialStrain = 99999.99;
    trialStrain      = 0.0;
    trialStrainRate  = 0.0;
    return 0;
}

UniaxialMaterial *
Elastic2Material::getCopy(void)
{
    Elastic2Material *theCopy = new Elastic2Material(this->getTag(),E,eta);
    theCopy->trialStrain     = trialStrain;
    theCopy->trialStrainRate = trialStrainRate;
	theCopy->initialStrain = 99999.99;
    return theCopy;
}

int 
Elastic2Material::sendSelf(int cTag, Channel &theChannel)
{
  int res = 0;
  static Vector data(4);
  data(0) = this->getTag();
  data(1) = E;
  data(2) = eta;
  data(3) = zeroE;

  res = theChannel.sendVector(this->getDbTag(), cTag, data);
  if (res < 0) 
    opserr << "Elastic2Material::sendSelf() - failed to send data\n";

  return res;
}

int 
Elastic2Material::recvSelf(int cTag, Channel &theChannel, 
			       FEM_ObjectBroker &theBroker)
{
  int res = 0;
  static Vector data(4);
  res = theChannel.recvVector(this->getDbTag(), cTag, data);
  
  if (res < 0) {
      opserr << "Elastic2Material::recvSelf() - failed to receive data\n";
      E = 0; 
      this->setTag(0);      
  }
  else {
    this->setTag(data(0));
    E   = data(1);
    eta = data(2);
    zeroE = data(3);
  }
    
  return res;
}

void 
Elastic2Material::Print(OPS_Stream &s, int flag)
{
    s << "Elastic tag: " << this->getTag() << endln;
    s << "  E: " << E << " eta: " << eta << endln;
}

int
Elastic2Material::setParameter(const char **argv, int argc, Parameter &param)
{
  if (strcmp(argv[0],"zeroE") == 0)
    return param.addObject(3, this);
  else if (strcmp(argv[0],"E") == 0)
    return param.addObject(1, this);
  
  else if (strcmp(argv[0],"eta") == 0)
    return param.addObject(2, this);
  
  else if (strcmp(argv[0],"zeroE") == 0)
    return param.addObject(3, this);
  
  return -1;
}
  
int 
Elastic2Material::updateParameter(int parameterID, Information &info)
{
  switch(parameterID) {
  case -1:
    return -1;
  case 1:
    E = info.theDouble;
    return 0;
  case 2:
    eta = info.theDouble;
    return 0;
  case 3:
    zeroE = info.theInt;
    return 0;
  default:
    return -1;
  }

  return 0;
}

