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

// $Revision: 1.15 $
// $Date: 2008-08-26 16:32:23 $
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
#include <Matrix.h>
#include <Information.h>
#include <Parameter.h>
#include <string.h>

#include <math.h>
#include <float.h>
#include <elementAPI.h>

void* OPS_HardeningMaterial()
{
    int numdata = OPS_GetNumRemainingInputArgs();
    if (numdata < 5) {
	opserr << "WARNING insufficient arguments\n";
	opserr << "Want: uniaxialMaterial Hardening tag? E? sigmaY? H_iso? H_kin? <eta?>" << endln;
	return 0;
    }

    int tag;
    numdata = 1;
    if (OPS_GetIntInput(&numdata,&tag) < 0) {
	opserr << "WARNING: failed to read tag\n";
	return 0;
    }

    double data[4];
    numdata = 4;
    if (OPS_GetDoubleInput(&numdata,data)) {
	opserr << "WARING: failed to read data\n";
	return 0;
    }

    double eta = 0.0;
    numdata = OPS_GetNumRemainingInputArgs();
    if (numdata > 0) {
	numdata = 1;
	if (OPS_GetDouble(&numdata,&eta)<0) {
	    opserr << "WARNING: failed to read eta\n";
	    return 0;
	}
    }

    UniaxialMaterial* mat = new HardeningMaterial(tag,data[0],data[1],data[2],data[3],eta);
    if (mat == 0) {
	opserr << "WARNING: failed to create Hardeningmaterial material\n";
	return 0;
    }

    return mat;
}

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
    double xsi = Tstress - Hkin*CplasticStrain;

    // Compute yield criterion
    double f = fabs(xsi) - (sigmaY + Hiso*Chardening);

    // Elastic step ... no updates required
    if (f <= -DBL_EPSILON * E) {
    //if (f <= 1.0e-8) {
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
    Chardening = 0.0;

    // Reset trial history variables
    TplasticStrain = 0.0;
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
    theCopy->Chardening = Chardening;

    // Copy trial history variables
    theCopy->TplasticStrain = TplasticStrain;
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
  
  static Vector data(11);
  
  data(0) = this->getTag();
  data(1) = E;
  data(2) = sigmaY;
  data(3) = Hiso;
  data(4) = Hkin;
  data(5) = eta;
  data(6) = CplasticStrain;
  data(7) = Chardening;
  data(8) = Tstrain;
  data(9) = Tstress;
  data(10) = Ttangent;
  
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
  
  static Vector data(11);
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
    Chardening = data(7);
    Tstrain = data(8);
    Tstress = data(9);
    Ttangent = data(10);
  }
    
  return res;
}

void 
HardeningMaterial::Print(OPS_Stream &s, int flag)
{
	if (flag == OPS_PRINT_PRINTMODEL_MATERIAL) {
		s << "HardeningMaterial, tag: " << this->getTag() << endln;
		s << "  E: " << E << endln;
		s << "  sigmaY: " << sigmaY << endln;
		s << "  Hiso: " << Hiso << endln;
		s << "  Hkin: " << Hkin << endln;
		s << "  eta: " << eta << endln;
	}
    
	if (flag == OPS_PRINT_PRINTMODEL_JSON) {
		s << "\t\t\t{";
		s << "\"name\": \"" << this->getTag() << "\", ";
		s << "\"type\": \"HardeningMaterial\", ";
		s << "\"E\": " << E << ", ";
		s << "\"fy\": " << sigmaY << ", ";
		s << "\"Hiso\": " << Hiso << ", ";
		s << "\"Hkin\": " << Hkin << ", ";
		s << "\"eta\": " << eta << "}";
	}
}


// AddingSensitivity:BEGIN ///////////////////////////////////
int
HardeningMaterial::setParameter(const char **argv, int argc, Parameter &param)
{
  if (strcmp(argv[0],"sigmaY") == 0 || strcmp(argv[0],"fy") == 0 || strcmp(argv[0],"Fy") == 0) {
    param.setValue(sigmaY);
    return param.addObject(1, this);
  }
  if (strcmp(argv[0],"E") == 0) {
    param.setValue(E);
    return param.addObject(2, this);
  }
  if (strcmp(argv[0],"H_kin") == 0 || strcmp(argv[0],"Hkin") == 0) {
    param.setValue(Hkin);
    return param.addObject(3, this);
  }
  if (strcmp(argv[0],"H_iso") == 0 || strcmp(argv[0],"Hiso") == 0) {
    param.setValue(Hiso);
    return param.addObject(4, this);
  }
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
HardeningMaterial::getStressSensitivity(int gradIndex, bool conditional)
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
	double ChardeningSensitivity	 = 0.0;
	if (SHVs != 0 && gradIndex < SHVs->noCols()) {
		CplasticStrainSensitivity = (*SHVs)(0,gradIndex);
		ChardeningSensitivity	 = (*SHVs)(1,gradIndex);
	}

	// Elastic trial stress
	double Tstress = E * (Tstrain-CplasticStrain);

	// Compute trial stress relative to committed back stress
	double xsi = Tstress - Hkin*CplasticStrain;

	// Compute yield criterion
	double f = fabs(xsi) - (sigmaY + Hiso*Chardening);

	double sensitivity;

	// Elastic step ... no updates required
	if (f <= -DBL_EPSILON * E) {
	//if (f <= 1.0e-8) {

		sensitivity = ESensitivity*(Tstrain-CplasticStrain)-E*CplasticStrainSensitivity;

	}

	// Plastic step
	else { 

		double TstressSensitivity = ESensitivity*(Tstrain-CplasticStrain)-E*CplasticStrainSensitivity;

		int sign = (xsi < 0) ? -1 : 1;

		double dGamma = f / (E+Hiso+Hkin);
		
		double CbackStressSensitivity = (HkinSensitivity*CplasticStrain + Hkin*CplasticStrainSensitivity);

		double fSensitivity = (TstressSensitivity-CbackStressSensitivity)*sign
			- SigmaYSensitivity - HisoSensitivity*Chardening - Hiso*ChardeningSensitivity;
		
		double dGammaSensitivity = 
			(fSensitivity*(E+Hkin+Hiso)-f*(ESensitivity+HkinSensitivity+HisoSensitivity))
			/((E+Hkin+Hiso)*(E+Hkin+Hiso));
		//double dGammaSensitivity = fSensitivity/(E+Hkin+Hiso);
		
		sensitivity = (TstressSensitivity-dGammaSensitivity*E*sign-dGamma*ESensitivity*sign);
		//sensitivity = TstressSensitivity-dGammaSensitivity*E*sign;
	}

	return sensitivity;
}


double
HardeningMaterial::getTangentSensitivity(int gradIndex)
{
  if (parameterID < 2 || parameterID > 4)
    return 0.0;


  // Elastic trial stress
  double Tstress = E * (Tstrain-CplasticStrain);
  
  // Compute trial stress relative to committed back stress
  double xsi = Tstress - Hkin*CplasticStrain;
  
  // Compute yield criterion
  double f = fabs(xsi) - (sigmaY + Hiso*Chardening);
  
  // Elastic step ... no updates required
  if (f <= -DBL_EPSILON * E) {

    if (parameterID == 2)
      return 1.0; 
  }
  
  // Plastic step
  else { 

    double EHK = E + Hiso + Hkin;
    double EHK2 = EHK*EHK;

    if (parameterID == 2)  // E
      return (EHK*(Hkin+Hiso)-E*(Hkin+Hiso)) / EHK2;
    else if (parameterID == 3)  // Hkin
      return (EHK*E          -E*(Hkin+Hiso)) / EHK2;
    else if (parameterID == 4)  // Hiso
      return (EHK*E          -E*(Hkin+Hiso)) / EHK2;
  }

  return 0.0;
}


double
HardeningMaterial::getInitialTangentSensitivity(int gradIndex)
{
  // For now, assume that this is only called for initial stiffness 
  if (parameterID == 2)
    return 1.0; 
  else
    return 0.0;
}


int
HardeningMaterial::commitSensitivity(double TstrainSensitivity, int gradIndex, int numGrads)
{
	if (SHVs == 0) {
		SHVs = new Matrix(2,numGrads);
	}

	if (gradIndex >= SHVs->noCols()) {
	  //opserr << gradIndex << ' ' << SHVs->noCols() << endln;
	  return 0;
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
	double CplasticStrainSensitivity= (*SHVs)(0,gradIndex);
	double ChardeningSensitivity	= (*SHVs)(1,gradIndex);

	// Elastic trial stress
	double Tstress = E * (Tstrain-CplasticStrain);

	// Compute trial stress relative to committed back stress
	double xsi = Tstress - Hkin*CplasticStrain;

	// Compute yield criterion
	double f = fabs(xsi) - (sigmaY + Hiso*Chardening);

	// Elastic step ... no updates required
	if (f <= -DBL_EPSILON * E) {
	//if (f <= 1.0e-8) {
		// No changes in the sensitivity history variables
	}

	// Plastic step
	else { 

		double TstressSensitivity = ESensitivity*(Tstrain-CplasticStrain)
			+ E*(TstrainSensitivity-CplasticStrainSensitivity);

		int sign = (xsi < 0) ? -1 : 1;
		//f = 0.0;
		double dGamma = f / (E+Hiso+Hkin);

		double CbackStressSensitivity = (HkinSensitivity*CplasticStrain + Hkin*CplasticStrainSensitivity);

		double fSensitivity = (TstressSensitivity-CbackStressSensitivity)*sign
			- SigmaYSensitivity - HisoSensitivity*Chardening - Hiso*ChardeningSensitivity;

	    double dGammaSensitivity = (fSensitivity*(E+Hkin+Hiso)-f*(ESensitivity+HkinSensitivity+HisoSensitivity))/((E+Hkin+Hiso)*(E+Hkin+Hiso));
		//double dGammaSensitivity = fSensitivity/(E+Hkin+Hiso);

		(*SHVs)(0,gradIndex) += dGammaSensitivity*sign;
		(*SHVs)(1,gradIndex) += dGammaSensitivity;
	}

	return 0;
}
// AddingSensitivity:END /////////////////////////////////////////////

