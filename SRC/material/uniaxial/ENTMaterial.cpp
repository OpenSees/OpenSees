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
// $Date: 2003-02-25 23:33:38 $
// $Source: /usr/local/cvs/OpenSees/SRC/material/uniaxial/ENTMaterial.cpp,v $
                                                                        
                                                                        
// File: ~/material/ENTMaterial.C
//
// Written: fmk 
// Created: 07/98
// Revision: A
//
// Description: This file contains the class implementation for 
// ENTMaterial. 
//
// What: "@(#) ENTMaterial.C, revA"

#include <ENTMaterial.h>
#include <Vector.h>
#include <Channel.h>
#include <Information.h>

ENTMaterial::ENTMaterial(int tag, double e)
:UniaxialMaterial(tag,MAT_TAG_ENTMaterial),
 E(e), trialStrain(0.0) 
{

}

ENTMaterial::ENTMaterial()
:UniaxialMaterial(0,MAT_TAG_ENTMaterial),
 E(0.0), trialStrain(0.0)
{

}

ENTMaterial::~ENTMaterial()
{
  // does nothing
}

int 
ENTMaterial::setTrialStrain(double strain, double strainRate)
{
    trialStrain = strain;
    return 0;
}

double 
ENTMaterial::getStrain(void)
{
	return trialStrain;
}

double 
ENTMaterial::getStress(void)
{
	if (trialStrain < 0.0)
		return E*trialStrain;
	else
		return 0.0;
}

double 
ENTMaterial::getTangent(void)
{
	if (trialStrain <= 0.0)
		return E;
	else
		return 0.0;
}

int 
ENTMaterial::commitState(void)
{
    return 0;
}

int 
ENTMaterial::revertToLastCommit(void)
{
    return 0;
}

int 
ENTMaterial::revertToStart(void)
{
    return 0;
}

UniaxialMaterial *
ENTMaterial::getCopy(void)
{
    ENTMaterial *theCopy = new ENTMaterial(this->getTag(),E);
    theCopy->trialStrain = trialStrain;
    return theCopy;
}

int 
ENTMaterial::sendSelf(int cTag, Channel &theChannel)
{
  int res = 0;
  static Vector data(2);
  data(0) = this->getTag();
  data(1) = E;
  res = theChannel.sendVector(this->getDbTag(), cTag, data);
  if (res < 0) 
    opserr << "ENTMaterial::sendSelf() - failed to send data\n";

  return res;
}

int 
ENTMaterial::recvSelf(int cTag, Channel &theChannel, 
			       FEM_ObjectBroker &theBroker)
{
  int res = 0;
  static Vector data(2);
  res = theChannel.recvVector(this->getDbTag(), cTag, data);
  
  if (res < 0) {
      opserr << "ENTMaterial::recvSelf() - failed to receive data\n";
      E = 0; 
      this->setTag(0);      
  }
  else {
    this->setTag((int)data(0));
    E = data(1);
  }
    
  return res;
}

void 
ENTMaterial::Print(OPS_Stream &s, int flag)
{
    s << "ENTMaterial, tag: " << this->getTag() << endln;
    s << "  E: " << E << endln;
}

int
ENTMaterial::setParameter(const char **argv, int argc, Information &info)
{
	if (strcmp(argv[0],"E") == 0) {
		info.theType = DoubleType;
		return 1;
	}
	else
		return -1;
}

int 
ENTMaterial::updateParameter(int parameterID, Information &info)
{
	switch(parameterID) {
	case -1:
		return -1;
	case 1:
		E = info.theDouble;
		return 0;
	default:
		return -1;
	}
}
