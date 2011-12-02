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
// $Source: /usr/local/cvs/OpenSees/SRC/material/uniaxial/HardeningMaterial.cpp,v $
                                                                        
                                                                        
// File: ~/material/HardeningMaterial.C
//
// Written: MHS
// Created: May 2000
// Revision: A
//
// Description: This file contains the class implementation for 
// HardeningMaterial. 
//
// What: "@(#) HardeningMaterial.C, revA"

#include <HardeningMaterial.h>
#include <Vector.h>
#include <Channel.h>
#include <math.h>

HardeningMaterial::HardeningMaterial(int tag, double e, double s,
				 double k, double h)
:UniaxialMaterial(tag,MAT_TAG_Hardening),
 E(e), sigmaY(s), K(k), H(h)
{
    this->revertToStart();
	this->revertToLastCommit();
}

HardeningMaterial::HardeningMaterial()
:UniaxialMaterial(0,MAT_TAG_Hardening),
 E(0.0), sigmaY(0.0), K(0.0), H(0.0)
{
	this->revertToStart();
	this->revertToLastCommit();
}

HardeningMaterial::~HardeningMaterial()
{
    // does nothing
}

int 
HardeningMaterial::setTrialStrain (double strain, double strainRate)
{
    // Set total strain
    Tstrain = strain;
    
    // Elastic trial stress
    Tstress = E * (Tstrain-CplasticStrain);
    
    // Compute trial stress relative to committed back stress
    double xsi = Tstress - CbackStress;

    // Compute yield criterion
    double f = fabs(xsi) - (sigmaY + K*Chardening);
    
    // Elastic step ... no updates required
    if (f <= 0.0)
    {
	// Set trial tangent
	Ttangent = E;
    }
    // Plastic step ... perform return mapping algorithm
    else
    {
	// Compute consistency parameter
	double dGamma = f / (E+K+H);

	// Find sign of xsi
	int sign = (xsi < 0) ? -1 : 1;

	// Bring trial stress back to yield surface
	Tstress -= dGamma*E*sign;
	
	// Update plastic strain
	TplasticStrain = CplasticStrain + dGamma*sign;
	
	// Update back stress
	TbackStress = CbackStress + dGamma*H*sign;
	
	// Update internal hardening variable
	Thardening = Chardening + dGamma;
	
	// Set trial tangent
	Ttangent = E*(H+K) / (E+H+K);
    }

    return 0;
}

double 
HardeningMaterial::getStress(void)
{
    return Tstress;
}

double 
HardeningMaterial::getTangent(void)
{
    return Ttangent;
}

double
HardeningMaterial::getSecant (void)
{
    if (Tstrain != 0.0)
	return Tstress/Tstrain;
    else
	return E;
}

double 
HardeningMaterial::getStrain(void)
{
    return Tstrain;
}

int 
HardeningMaterial::commitState(void)
{
    // Commit trial state variables
    Cstrain = Tstrain;
    CplasticStrain = TplasticStrain;
    CbackStress = TbackStress;
    Chardening = Thardening;

    // Do not need to commit the trial stress or trial tangent
    
    return 0;
}

int 
HardeningMaterial::revertToLastCommit(void)
{
    // Set trial state to last committed state
    Tstrain = Cstrain;
    TplasticStrain = CplasticStrain;
    TbackStress = CbackStress;
    Thardening = Chardening;

    // Recompute trial stress and trial tangent
	this->setTrialStrain(Cstrain);

    return 0;
}

int 
HardeningMaterial::revertToStart(void)
{
    // Reset committed state variables
    Cstrain = 0.0;
    CplasticStrain = 0.0;
    CbackStress = 0.0;
    Chardening = 0.0;

    // Reset the trial state
    this->revertToLastCommit();    
    
    return 0;
}

UniaxialMaterial *
HardeningMaterial::getCopy(void)
{
    HardeningMaterial *theCopy =
	new HardeningMaterial(this->getTag(), E, sigmaY, K, H);

    // Copy trial state variables
    theCopy->Tstrain = Tstrain;
    theCopy->Tstress = Tstress;
    theCopy->Ttangent = Ttangent;
    theCopy->TplasticStrain = TplasticStrain;
    theCopy->TbackStress = TbackStress;
    theCopy->Thardening = Thardening;
    
    // Copy committed state variables
    theCopy->Cstrain = Cstrain;
    theCopy->CplasticStrain = CplasticStrain;
    theCopy->CbackStress = CbackStress;
    theCopy->Chardening = Chardening;
    
    return theCopy;
}

int 
HardeningMaterial::sendSelf(int cTag, Channel &theChannel)
{
  int res = 0;
  
  static Vector data(9);
  
  data(0) = this->getTag();
  data(1) = E;
  data(2) = sigmaY;
  data(3) = K;
  data(4) = H;
  data(5) = Cstrain;
  data(6) = CplasticStrain;
  data(7) = CbackStress;
  data(8) = Chardening;
  
  res = theChannel.sendVector(this->getDbTag(), cTag, data);
  if (res < 0) 
    cerr << "HardeningMaterial::sendSelf() - failed to send data\n";

  return res;
}

int 
HardeningMaterial::recvSelf(int cTag, Channel &theChannel, 
			       FEM_ObjectBroker &theBroker)
{
  int res = 0;
  
  static Vector data(9);
  res = theChannel.recvVector(this->getDbTag(), cTag, data);
  
  if (res < 0) {
      cerr << "HardeningMaterial::recvSelf() - failed to receive data\n";
      E = 0; 
      this->setTag(0);      
  }
  else {
    this->setTag((int)data(0));
    E = data(1);
    sigmaY = data(2);
    K = data(3);
    H = data(4);
    Cstrain = data(5);
    CplasticStrain = data(6);
    CbackStress = data(7);
    Chardening = data(8);

    // Set the trial state variables
    revertToLastCommit();
  }
    
  return res;
}

void 
HardeningMaterial::Print(ostream &s, int flag)
{
    s << "HardeningMaterial, tag: " << this->getTag() << endl;
    s << "  E: " << E << endl;
    s << "  sigmaY: " << sigmaY << endl;
    s << "  K: " << K << endl;
    s << "  H: " << H << endl;
}


