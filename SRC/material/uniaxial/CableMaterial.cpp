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
                                                                        
// $Revision: 1.8 $
// $Date: 2008-04-14 21:17:42 $
// $Source: /usr/local/cvs/OpenSees/SRC/material/uniaxial/CableMaterial.cpp,v $
                                                                        
// Written: Charles Chadwell 
// Created: 07/01
//
// Description: This file contains the class definition for 
// CableMaterial. CableMaterial provides the abstraction
// of an elastic uniaxial material,
//
// The input parameters are the Prestress, E, Effective Self Weight (gravity component of 
// Weight per volume transverse to the cable), and Length of the cable.
//
// The cable Force Displacement is converted to Stress Strain material for use 
// with the truss element.  The stress strain ranges from slack (large strain at zero 
// stress) to taught (linear with modulus E).  The material has no history and is not
// path dependent.
//
//
// What: "@(#) CableMaterial.cpp, revA"

#include <CableMaterial.h>
#include <Vector.h>
#include <Channel.h>
#include <Information.h>

#include <elementAPI.h>

void * OPS_ADD_RUNTIME_VPV(OPS_CableMaterial)
{
  // Pointer to a uniaxial material that will be returned
  UniaxialMaterial *theMaterial = 0;

  int numArgs = OPS_GetNumRemainingInputArgs();
  if (numArgs != 5) {
    opserr << "Invalid # args, want: uniaxialMaterial Cable tag? $presetress $E $effUnitWeight $Lelement \n";
    return 0;
  }
  
  int iData[1];
  double dData[4];

  int numData = 1;
  if (OPS_GetIntInput(&numData, iData) != 0) {
    opserr << "WARNING invalid tag for uniaxialMaterial Cable" << endln;
    return 0;
  }

  numData = 4;
  if (OPS_GetDoubleInput(&numData, dData) != 0) {
    opserr << "Invalid data for uniaxial Cable " << iData[0] << endln;
    return 0;	
  }

  // Parsing was successful, allocate the material
  theMaterial = new CableMaterial(iData[0], dData[0], dData[1], dData[2], dData[3]);
  if (theMaterial == 0) {
    opserr << "WARNING could not create uniaxialMaterial of type Cable\n";
    return 0;
  }

  return theMaterial;
}


CableMaterial::CableMaterial(int tag, double PRESTRESS, double e, double UNIT_WEIGHT_EFF, double L_Element)
:UniaxialMaterial(tag,MAT_TAG_CableMaterial),
 Ps(PRESTRESS), E(e), Mue(UNIT_WEIGHT_EFF), 
 L(L_Element), trialStrain(0.0)
{

}

CableMaterial::CableMaterial()
:UniaxialMaterial(0,MAT_TAG_CableMaterial),
 Ps(0.0), E(0.0), Mue(0.0), 
 L(0.0), trialStrain(0.0)
{

}

CableMaterial::~CableMaterial()
{
  // does nothing
}

int 
CableMaterial::setTrialStrain(double strain, double strainRate)
{
  trialStrain     = strain;

  double derivE, derivG;

  // Check if out side the range of inportance 
  double testStress, dP, curstrain, e0;
  int i = 0;
    
  // Perameters for bisection
  double L_bound = 0, U_bound, middle = 0;
	
	
  if (trialStrain < 0) 
    U_bound = Ps;
  else {
    U_bound = E * trialStrain + Ps;
    testStress = U_bound;
  }
    
  // Check if slack in cable has been taken out and it is a bar
  e0 = Mue*Mue*L*L/(24*Ps*Ps) - Ps/E;
  if (trialStrain > 0 && abs(trialStrain - evalStress((trialStrain - e0)*E)) < 10e-9) 
    trialStress =  (trialStrain - e0)*E;

  // Check if all slack
  if (trialStrain < - Ps/E*10.0) 
    trialStress =  0.0; 

  // if stress is in between then do iterations -- Bisection
  dP = U_bound - L_bound;
    
  while (abs(dP)/U_bound > 0.00000001 && i < 100) {

    middle = .5 * (U_bound + L_bound);
    curstrain = evalStress(middle);
        
    if (curstrain <= trialStrain) {
      L_bound = middle;
    } else {
      U_bound = middle;
    }
    dP = U_bound - L_bound;
    i++;
  }
	
  // if it did not converge - return near zero stress
  if (i == 100) 
    trialStress =  0.0;
  else 
    trialStress = middle;

  if (trialStress <= 0.0) 
    trialTangent = 0.0;
	
  // Elastic Part
  derivE = 1 / E * (1. - Mue * Mue * L * L / (24. * trialStress * trialStress) * (1. - 2. * Ps / trialStress));
  // Geometric Part
  derivG = 1 / 12. * Mue * Mue * L * L / (trialStress * trialStress * trialStress);
  
  if (derivE + derivG != 0.0)
    trialTangent =  1.0 / (derivE + derivG);
  else 
    trialTangent =  1e-8;
  
  return 0;
}


int 
CableMaterial::setTrial(double strain, double &stress, double &tangent, double strainRate)
{
  this->setTrialStrain(strain, strainRate);

  // set the return values
  stress = trialStress;
  tangent = trialTangent;

  return 0;
}

double 
CableMaterial::getStress(void)
{
  return trialStress;

}

double 
CableMaterial::evalStress(double stress)
{
    double strainG, strainE;
    
    // Should never be zero or less than zero
    if (stress <= 0) {return -10;}
		
    // Elastic Part
    strainE = 1 / E * (stress - Ps) * (1 + Mue * Mue * L * L / (24 * stress));

    // Geometric Part
    strainG = 1 / 24 * Mue * Mue * L * L * (1 / (Ps * Ps) - 1 / (stress * stress));

    return strainE + strainG;
}


double 
CableMaterial::getTangent(void) 
{

  return trialTangent;

};
 
int 
CableMaterial::commitState(void)
{
    return 0;
}

int 
CableMaterial::revertToLastCommit(void)
{
    return 0;
}

int 
CableMaterial::revertToStart(void)
{
  trialStrain = 0.0;

  return 0;
}

UniaxialMaterial *
CableMaterial::getCopy(void)
{
    CableMaterial *theCopy = new CableMaterial(this->getTag(), Ps, E, Mue, L);
    theCopy->trialStrain = trialStrain;
    return theCopy;
}

int 
CableMaterial::sendSelf(int cTag, Channel &theChannel)
{
  int res = 0;
  static Vector data(5);
  data(0) = this->getTag();
  data(1) = Ps;
  data(2) = E;
  data(3) = Mue;
  data(4) = L;

  res = theChannel.sendVector(this->getDbTag(), cTag, data);
  if (res < 0) 
    opserr << "CableMaterial::sendSelf() - failed to send data\n";

  return res;
}

int 
CableMaterial::recvSelf(int cTag, Channel &theChannel, 
			       FEM_ObjectBroker &theBroker)
{
  int res = 0;
  static Vector data(5);
  res = theChannel.recvVector(this->getDbTag(), cTag, data);
  
  if (res < 0) {
      opserr << "CableMaterial::recvSelf() - failed to receive data\n";
      E = 0; 
      this->setTag(0);      
      return res;
  }

  this->setTag(data(0));
  Ps   = data(1);
  E = data(2);
  Mue     = data(3);
  L = data(4);
    
  return res;
}

void 
CableMaterial::Print(OPS_Stream &s, int flag)
{
   s << "CableMaterial tag: " << this->getTag() << endln;
    s << "  E: " << E << " Prestress: " << Ps << endln;
}

//int
//CableMaterial::setParameter(const char **argv, int argc, Information &info)
//{
//	if (strcmp(argv[0],"E") == 0) {
//		info.theType = DoubleType;
///		return 1;
//	}
//	else if (strcmp(argv[0],"eta") == 0) {
//		info.theType = DoubleType;
//		return 2;
//	}
//	else
//		return -1;
//}
//
//int 
//CableMaterial::updateParameter(int parameterID, Information &info)
//{
//	switch(parameterID) {
//	case -1:
//		return -1;
//	case 1:
//		E = info.theDouble;
//		return 0;
//	case 2:
//		eta = info.theDouble;
//		return 0;
//	default:
//		return -1;
//	}
//}

 
double 
CableMaterial::abs(double value)
{
	if (value < 0) return -value;
	else return value;
}
