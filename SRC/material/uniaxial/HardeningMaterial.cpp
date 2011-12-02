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

// $Revision: 1.11 $
// $Date: 2003-03-11 03:49:27 $
// $Source: /usr/local/cvs/OpenSees/SRC/material/uniaxial/HardeningMaterial.cpp,v $

// Written: MHS
// Created: May 2000
// Revision: A
//
// Description: This file contains the class implementation for 
// HardeningMaterial. 

#include <HardeningMaterial.h>
#include <Vector.h>
#include <Channel.h>
#include <math.h>
#include <Matrix.h>
#include <Information.h>
#include <float.h>

HardeningMaterial::HardeningMaterial(int tag, double e, double s,
				     double hi, double hk, double n)
:UniaxialMaterial(tag,MAT_TAG_Hardening),
 E(e), sigmaY(s), Hiso(hi), Hkin(hk), eta(n)
{
// AddingSensitivity:BEGIN /////////////////////////////////////
	parameterID = 0;
	SHVs = 0;
// AddingSensitivity:END //////////////////////////////////////

	// Initialize variables
    this->revertToStart();
}

HardeningMaterial::HardeningMaterial()
:UniaxialMaterial(0,MAT_TAG_Hardening),
 E(0.0), sigmaY(0.0), Hiso(0.0), Hkin(0.0), eta(0.0)
{
// AddingSensitivity:BEGIN /////////////////////////////////////
	parameterID = 0;
	SHVs = 0;
// AddingSensitivity:END //////////////////////////////////////

	// Initialize variables
	this->revertToStart();
}

HardeningMaterial::~HardeningMaterial()
{
// AddingSensitivity:BEGIN /////////////////////////////////////
	if (SHVs != 0) 
		delete SHVs;
// AddingSensitivity:END //////////////////////////////////////
}

int 
HardeningMaterial::setTrialStrain (double strain, double strainRate)
{

    if (fabs(Tstrain-strain) < DBL_EPSILON)
      return 0;

    // Set total strain
    Tstrain = strain;
    
    // Elastic trial stress
    Tstress = E * (Tstrain-CplasticStrain);

    // Compute trial stress relative to committed back stress
    double xsi = Tstress - CbackStress;

    // Compute yield criterion
    double f = fabs(xsi) - (sigmaY + Hiso*Chardening);

    // Elastic step ... no updates required
//    if (f <= -DBL_EPSILON * E)
    if (f <= 1.0e-8)
    {
	// Set trial tangent
	Ttangent = E;
    }

    // Plastic step ... perform return mapping algorithm
    else {
      double etadt = 0.0;

      if (eta != 0.0 || ops_Dt != 0)
	  etadt = eta/ops_Dt;

      // Compute consistency parameter
      double dGamma = f / (E+Hiso+Hkin+etadt);

      // Find sign of xsi
      int sign = (xsi < 0) ? -1 : 1;

      // Bring trial stress back to yield surface
      Tstress -= dGamma*E*sign;
	
      // Update plastic strain
      TplasticStrain = CplasticStrain + dGamma*sign;
	
      // Update back stress
      TbackStress = CbackStress + dGamma*Hkin*sign;
	
      // Update internal hardening variable
      Thardening = Chardening + dGamma;
	
      // Set trial tangent
      Ttangent = E*(Hkin+Hiso+etadt) / (E+Hkin+Hiso+etadt);
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
HardeningMaterial::getStrain(void)
{
    return Tstrain;
}

int 
HardeningMaterial::commitState(void)
{
    // Commit trial history variables
    CplasticStrain = TplasticStrain;
    CbackStress = TbackStress;
    Chardening = Thardening;
    
    return 0;
}

int 
HardeningMaterial::revertToLastCommit(void)
{
  return 0;
}

int 
HardeningMaterial::revertToStart(void)
{
    // Reset committed history variables
    CplasticStrain = 0.0;
    CbackStress = 0.0;
    Chardening = 0.0;

    // Reset trial history variables
    TplasticStrain = 0.0;
    TbackStress = 0.0;
    Thardening = 0.0;

	// Initialize state variables
	Tstrain = 0.0;
	Tstress = 0.0;
	Ttangent = E;

// AddingSensitivity:BEGIN /////////////////////////////////
	if (SHVs != 0) 
		SHVs->Zero();
// AddingSensitivity:END //////////////////////////////////

    return 0;
}

UniaxialMaterial *
HardeningMaterial::getCopy(void)
{
    HardeningMaterial *theCopy =
	new HardeningMaterial(this->getTag(), E, sigmaY, Hiso, Hkin, eta);

    // Copy committed history variables
    theCopy->CplasticStrain = CplasticStrain;
    theCopy->CbackStress = CbackStress;
    theCopy->Chardening = Chardening;

    // Copy trial history variables
    theCopy->TplasticStrain = TplasticStrain;
    theCopy->TbackStress = TbackStress;
    theCopy->Thardening = Thardening;

    // Copy trial state variables
    theCopy->Tstrain = Tstrain;
    theCopy->Tstress = Tstress;
    theCopy->Ttangent = Ttangent;
    
    return theCopy;
}

int 
HardeningMaterial::sendSelf(int cTag, Channel &theChannel)
{
  int res = 0;
  
  static Vector data(12);
  
  data(0) = this->getTag();
  data(1) = E;
  data(2) = sigmaY;
  data(3) = Hiso;
  data(4) = Hkin;
  data(5) = eta;
  data(6) = CplasticStrain;
  data(7) = CbackStress;
  data(8) = Chardening;
  data(9) = Tstrain;
  data(10) = Tstress;
  data(11) = Ttangent;
  
  res = theChannel.sendVector(this->getDbTag(), cTag, data);
  if (res < 0) 
    opserr << "HardeningMaterial::sendSelf() - failed to send data\n";

  return res;
}

int 
HardeningMaterial::recvSelf(int cTag, Channel &theChannel, 
			       FEM_ObjectBroker &theBroker)
{
  int res = 0;
  
  static Vector data(12);
  res = theChannel.recvVector(this->getDbTag(), cTag, data);
  
  if (res < 0) {
      opserr << "HardeningMaterial::recvSelf() - failed to receive data\n";
      E = 0; 
      this->setTag(0);      
  }
  else {
    this->setTag((int)data(0));
    E = data(1);
    sigmaY = data(2);
    Hiso = data(3);
    Hkin = data(4);
    eta = data(5);
    CplasticStrain = data(6);
    CbackStress = data(7);
    Chardening = data(8);
    Tstrain = data(9);
    Tstress = data(10);
    Ttangent = data(11);
  }
    
  return res;
}

void 
HardeningMaterial::Print(OPS_Stream &s, int flag)
{
    s << "HardeningMaterial, tag: " << this->getTag() << endln;
    s << "  E: " << E << endln;
    s << "  sigmaY: " << sigmaY << endln;
    s << "  Hiso: " << Hiso << endln;
    s << "  Hkin: " << Hkin << endln;
    s << "  eta: " << eta << endln;
}


// AddingSensitivity:BEGIN ///////////////////////////////////
int
HardeningMaterial::setParameter(const char **argv, int argc, Information &info)
{
	if (argc < 1)
		return -1;

	if (strcmp(argv[0],"sigmaY") == 0 || strcmp(argv[0],"fy") == 0) {
		info.theType = DoubleType;
		return 1;
	}
	if (strcmp(argv[0],"E") == 0) {
		info.theType = DoubleType;
		return 2;
	}
	if (strcmp(argv[0],"H_kin") == 0) {
		info.theType = DoubleType;
		return 3;
	}
	if (strcmp(argv[0],"H_iso") == 0) {
		info.theType = DoubleType;
		return 4;
	}
	else
		opserr << "WARNING: Could not set parameter in HardeningMaterial. " << endln;
                
	return -1;
}

int
HardeningMaterial::updateParameter(int parameterID, Information &info)
{
	switch (parameterID) {
	case -1:
		return -1;
	case 1:
		this->sigmaY = info.theDouble;
		break;
	case 2:
		this->E = info.theDouble;
		break;
	case 3:
		this->Hkin = info.theDouble;
		break;
	case 4:
		this->Hiso = info.theDouble;
		break;
	default:
		return -1;
	}

	return 0;
}



int
HardeningMaterial::activateParameter(int passedParameterID)
{
	parameterID = passedParameterID;

	return 0;
}




double
HardeningMaterial::getStressSensitivity(int gradNumber, bool conditional)
{

	// First set values depending on what is random
	double SigmaYSensitivity = 0.0;
	double ESensitivity = 0.0;
	double HkinSensitivity = 0.0;
	double HisoSensitivity = 0.0;

	if (parameterID == 1) {  // sigmaY
		SigmaYSensitivity = 1.0;
	}
	else if (parameterID == 2) {  // E
		ESensitivity = 1.0;
	}
	else if (parameterID == 3) {  // Hkin
		HkinSensitivity = 1.0;
	}
	else if (parameterID == 4) {  // Hiso
		HisoSensitivity = 1.0;
	}
	else {
		// Nothing random here, but may have to return something in any case
	}

	// Then pick up history variables for this gradient number
	double CplasticStrainSensitivity = 0.0;
	double CbackStressSensitivity	 = 0.0;
	double ChardeningSensitivity	 = 0.0;
	if (SHVs != 0) {
		CplasticStrainSensitivity = (*SHVs)(0,(gradNumber-1));
		CbackStressSensitivity	 = (*SHVs)(1,(gradNumber-1));
		ChardeningSensitivity	 = (*SHVs)(2,(gradNumber-1));
	}

	double sensitivity;
	// Elastic step
	if ( fabs(TplasticStrain-CplasticStrain) < 1.0e-8) { 

		sensitivity = ESensitivity*(Tstrain-CplasticStrain)-E*CplasticStrainSensitivity;

	}

	// Plastic step
	else { 

		double myTstress = E * (Tstrain-CplasticStrain);

		double xsi = myTstress - CbackStress;
		
		double f = fabs(xsi) - (sigmaY + Hiso*Chardening);
		
		double TstressSensitivity = ESensitivity*(Tstrain-CplasticStrain)-E*CplasticStrainSensitivity;
		
		int sign = (xsi < 0) ? -1 : 1;
		
		double dGamma = f / (E+Hiso+Hkin);
		
		double fSensitivity = (TstressSensitivity-CbackStressSensitivity)*sign
			- SigmaYSensitivity - HisoSensitivity*Chardening - Hiso*ChardeningSensitivity;
		
		double dGammaSensitivity = 
			(fSensitivity*(E+Hkin+Hiso)-f*(ESensitivity+HkinSensitivity+HisoSensitivity))
			/((E+Hkin+Hiso)*(E+Hkin+Hiso));
		
		sensitivity = (TstressSensitivity-dGammaSensitivity*E*sign-dGamma*ESensitivity*sign);
	}

	return sensitivity;
}



double
HardeningMaterial::getInitialTangentSensitivity(int gradNumber)
{
	// For now, assume that this is only called for initial stiffness 
	if (parameterID == 2) {
		return 1.0; 
	}
	else {
		return 0.0;
	}
}


int
HardeningMaterial::commitSensitivity(double TstrainSensitivity, int gradNumber, int numGrads)
{
	if (SHVs == 0) {
		SHVs = new Matrix(3,numGrads);
	}

	// First set values depending on what is random
	double SigmaYSensitivity = 0.0;
	double ESensitivity = 0.0;
	double HkinSensitivity = 0.0;
	double HisoSensitivity = 0.0;

	if (parameterID == 1) {  // sigmaY
		SigmaYSensitivity = 1.0;
	}
	else if (parameterID == 2) {  // E
		ESensitivity = 1.0;
	}
	else if (parameterID == 3) {  // Hkin
		HkinSensitivity = 1.0;
	}
	else if (parameterID == 4) {  // Hiso
		HisoSensitivity = 1.0;
	}
	else {
		// Nothing random here, but may have to save SHV's in any case
	}

	// Then pick up history variables for this gradient number
	double CplasticStrainSensitivity= (*SHVs)(0,(gradNumber-1));
	double CbackStressSensitivity	= (*SHVs)(1,(gradNumber-1));
	double ChardeningSensitivity	= (*SHVs)(2,(gradNumber-1));

	// Elastic step
	if ( fabs(TplasticStrain-CplasticStrain) < 1.0e-8) { 
		// No changes in the sensitivity history variables
	}

	// Plastic step
	else { 

		double myTstress = E * (Tstrain-CplasticStrain);

		double xsi = myTstress - CbackStress;
		
		double f = fabs(xsi) - (sigmaY + Hiso*Chardening);

		double TstressSensitivity = ESensitivity*(Tstrain-CplasticStrain)
			+ E*(TstrainSensitivity-CplasticStrainSensitivity);

		int sign = (xsi < 0) ? -1 : 1;

		double dGamma = f / (E+Hiso+Hkin);

		double fSensitivity = (TstressSensitivity-CbackStressSensitivity)*sign
			- SigmaYSensitivity - HisoSensitivity*Chardening - Hiso*ChardeningSensitivity;

		double dGammaSensitivity = 
			(fSensitivity*(E+Hkin+Hiso)-f*(ESensitivity+HkinSensitivity+HisoSensitivity))
			/((E+Hkin+Hiso)*(E+Hkin+Hiso));

		(*SHVs)(0,(gradNumber-1)) += dGammaSensitivity*sign;
		(*SHVs)(1,(gradNumber-1)) += dGammaSensitivity*Hkin*sign + dGamma*HkinSensitivity*sign;
		(*SHVs)(2,(gradNumber-1)) += dGammaSensitivity;

	}

	return 0;
}
// AddingSensitivity:END /////////////////////////////////////////////

