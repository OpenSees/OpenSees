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
                                                                        
// $Revision: 1.1.1.1 $
// $Date: 2000-09-15 08:23:22 $
// $Source: /usr/local/cvs/OpenSees/SRC/material/uniaxial/ViscousMaterial.cpp,v $
                                                                        
// Written: Mehrdad Sasani 
// Created: June 2000
// Revision: A
//
// Description: This file contains the class implementation for 
// ViscousMaterial. 

#include <ViscousMaterial.h>
#include <Vector.h>
#include <Channel.h>

ViscousMaterial::ViscousMaterial(int tag, double c, double alpha)
:UniaxialMaterial(tag,MAT_TAG_Viscous),
 commitStrain(0.0), trialStrain(0.0), C(c), Alpha(alpha)
{

}

ViscousMaterial::ViscousMaterial()
:UniaxialMaterial(0,MAT_TAG_Viscous),
 commitStrain(0.0), C(0.0), Alpha(0.0)
{

}

ViscousMaterial::~ViscousMaterial()
{
  // does nothing

}

int 
ViscousMaterial::setTrialStrain(double velocity, double velocityRate)
{
    trialStrain = velocity;

    return 0;
}

double 
ViscousMaterial::getStress(void)
{
	int sgn = (trialStrain < 0.0) ? -1 : 1;

	double minvel = 1.e-11;

	if (fabs(trialStrain)<minvel) {
		// return Alpha*C*pow((10e-12),(Alpha-1));
	    return C*trialStrain/minvel*C*pow(minvel,Alpha);
	}
     
	return C*sgn*pow(fabs(trialStrain),Alpha);
}

double 
ViscousMaterial::getTangent(void)
{
	int sgn = (trialStrain < 0.0) ? -1 : 1;

	double minvel = 1.e-11;

    if (fabs(trialStrain)<minvel)
		return Alpha*C*pow((10e-12),(Alpha-1));
	else
		return Alpha*C*pow(fabs(trialStrain),(Alpha-1));
}

double
ViscousMaterial::getDampTangent(void)
{
	return this->getTangent();
}

double 
ViscousMaterial::getSecant(void)
{
	int sgn = (trialStrain < 0.0) ? -1 : 1;

	double minvel = 1.e-11;

    if (fabs(trialStrain)<minvel)
	    return C*pow(minvel,(Alpha-1));
	else
		return C*pow(fabs(trialStrain),(Alpha-1));
}

double 
ViscousMaterial::getStrain(void)
{
    return 0.0;
}

double 
ViscousMaterial::getStrainRate(void)
{
    return trialStrain;
}

int 
ViscousMaterial::commitState(void)
{
    commitStrain = trialStrain;
    return 0;
}

int 
ViscousMaterial::revertToLastCommit(void)
{
    trialStrain = commitStrain;
    return 0;
}

int 
ViscousMaterial::revertToStart(void)
{
    commitStrain = 0.0;
    trialStrain = 0.0;
    return 0;
}

UniaxialMaterial *
ViscousMaterial::getCopy(void)
{
    ViscousMaterial *theCopy = new ViscousMaterial(this->getTag(),C,Alpha);
    theCopy->commitStrain = commitStrain;
    return theCopy;
}


int 
ViscousMaterial::sendSelf(int cTag, Channel &theChannel)
{
  int res = 0;
  static Vector data(4);
  data(0) = this->getTag();
  data(1) = C;
  data(2) = Alpha;
  data(3) = commitStrain;
  res = theChannel.sendVector(this->getDbTag(), cTag, data);
  if (res < 0) 
    cerr << "ViscousMaterial::sendSelf() - failed to send data\n";

  return res;
}

int 
ViscousMaterial::recvSelf(int cTag, Channel &theChannel, 
			       FEM_ObjectBroker &theBroker)
{
  int res = 0;
  static Vector data(4);
  res = theChannel.recvVector(this->getDbTag(), cTag, data);
  
  if (res < 0) {
      cerr << "ViscousMaterial::recvSelf() - failed to receive data\n";
      C = 0; 
      this->setTag(0);      
  }
  else {
    this->setTag((int)data(0));
    C = data(1);
	Alpha = data(2);
    commitStrain = data(3);
    trialStrain = commitStrain;
  }
    
  return res;
}

void 
ViscousMaterial::Print(ostream &s, int flag)
{
    s << "Viscous tag: " << this->getTag() << endl;
    s << "  C: " << C << endl;
    s << "  Alpha: " << Alpha << endl;
}


