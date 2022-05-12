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
                                                                        
// $Revision: 1.00 $
// $Date: 2022-Apr-21 12:15:00 $
// $Source: /usr/local/cvs/OpenSees/SRC/material/uniaxial/HystereticAsym.cpp,v $
                                                                        
// Written: Salvatore Sessa Mail: salvatore.sessa2@unina.it
// Created: 04-2022
// Revision: A
//
// Description: This file contains the class implementation for 
// HystereticAsym. 
//
// What: "@(#) HystereticAsym.C, revA"


#include <HystereticAsym.h>
#include <Vector.h>
#include <Matrix.h>
#include <Channel.h>
#include <Information.h>
#include <Parameter.h>

#include <math.h>
#include <float.h>

#include <elementAPI.h>
#include <MaterialResponse.h>



void *
OPS_HystereticAsym()
{
	UniaxialMaterial *theMaterial = 0;
	int numArgs = OPS_GetNumRemainingInputArgs();
	if (numArgs < 7) {
		opserr << "WARNING insufficient arguments\n";
		opserr << "Want: uniaxialMaterial HystereticAsym tag? ka? kb? fo? b1? b2? gamma? <-alpha>" << endln;
		return 0;
	}


	int iData[1];
	double dData[6];
	


	int numData = 1;
	if (OPS_GetIntInput(&numData, iData) != 0) {
		opserr << "WARNING invalid tag for uniaxialMaterial HystereticAsym" << endln;
		return 0;
	}

	numData = 6;
	if (OPS_GetDoubleInput(&numData, dData) != 0) {
		opserr << "Invalid data for uniaxial HystereticAsym " << iData[0] << endln;
		return 0;
	}
	if (dData[0] <= 0.0) {
		opserr << "uniaxialMaterial HystereticAsym ka must be positive" << endln;
		return 0;
	}
	if (dData[1] >= dData[0]) {
		opserr << "uniaxialMaterial HystereticAsym kb must be < ka" << endln;
		return 0;
	}
	if (dData[2] <= 0.0) {
		opserr << "uniaxialMaterial HystereticAsym fo and/or alpha must be positive" << endln;
		return 0;
	}
	

	if (numArgs == 8) {
		const char* flag = OPS_GetString();

		if (strcmp(flag, "-alpha") == 0) {
			double newFbar = (dData[0] - dData[1]) * 0.5 / dData[2];
			dData[2] = newFbar;
		}
		else {
			opserr << "WARNING unknown parameter flag, ignoring it for uniaxialMaterial HystereticAsym " << iData[0] << endln;
		}
	}


	// Parsing was successful, allocate the material
	theMaterial = new HystereticAsym(iData[0], dData[0], dData[1], dData[2], dData[3], dData[4], dData[5]);
	if (theMaterial == 0) {
		opserr << "WARNING could not create uniaxialMaterial of type HystereticAsym\n";
		return 0;
	}

	return theMaterial;

}

HystereticAsym::HystereticAsym(int tag, double K1, double K2, double FBAR, double B1, double B2, double G)
   :UniaxialMaterial(tag, MAT_TAG_HystereticAsym),
   ka(K1), kb(K2), fo(FBAR), b1(B1), b2(B2), g(G)
{
   // Sets all history and state variables to initial values

	// Parameters:
	a = (ka - kb) * 0.5 / fo;
	
   // State variables
   Cstrain = 0.0;
   Cstress = 0.0;
   Tstrain = 0.0;
   Tstress = 0.0;

   // History variables
   st = 1.0;

   
   Ttangent = 0.5*(ka-kb);
   InitTangent = Ttangent;

   Ctangent = Ttangent;

   
// AddingSensitivity:BEGIN /////////////////////////////////////
	parameterID = 0;
	SHVs = 0;
// AddingSensitivity:END //////////////////////////////////////
}

HystereticAsym::HystereticAsym():UniaxialMaterial(0, MAT_TAG_HystereticAsym),
 ka(0.0), kb(0.0), fo(0.0), b1(0.0), b2(0.0), g(0.0)
{

// AddingSensitivity:BEGIN /////////////////////////////////////
	parameterID = 0;
	SHVs = 0;
// AddingSensitivity:END //////////////////////////////////////

}

HystereticAsym::~HystereticAsym ()
{
// AddingSensitivity:BEGIN /////////////////////////////////////
	if (SHVs != 0) 
		delete SHVs;
// AddingSensitivity:END //////////////////////////////////////
}

int HystereticAsym::setTrialStrain (double strain, double strainRate)
{
   // Determine change in strain from last converged state

   Tstrain = strain;
   dStrain = strain - Cstrain;
   st = signum(dStrain);	

   fe = -Cstrain * (b1 + b2) + exp(b1 * Cstrain) - exp(-b2 * Cstrain);
   delta2 = (fe - Cstress) * exp(-g * Cstrain) + kb * Cstrain + st * fo;
   uj = Cstrain * st + 1.0 / a * log(st * a * delta2 / (ka - kb));
   delta3 = kb * Tstrain - st * (ka - kb) / a * exp(-a*(Tstrain * st - uj)) + fo * st;
   Tstress = -Tstrain*(b1+b2) + exp(b1*Tstrain) - exp(-b2*Tstrain) + delta3*exp(g*Tstrain);
   
   double Ddelta3Dut = kb + (ka - kb) * exp(-a * (Tstrain*st-uj));
   Ttangent = -1.0*(b1 + b2) + b1*exp(b1 * Tstrain) +b2* exp(-b2 * Tstrain) + exp(g * Tstrain) * (g*delta3 + Ddelta3Dut);
   

   return 0;
}



double HystereticAsym::getStrain ()
{
   return Tstrain;
}

double HystereticAsym::getStress ()
{
   return Tstress;
}

double HystereticAsym::getTangent ()
{
   return Ttangent;
}

double HystereticAsym::getInitialTangent()
{
	InitTangent = 0.5 * (ka - kb);
	return InitTangent;
}

int HystereticAsym::commitState ()
{
   // State variables
   Cstrain = Tstrain;
   Cstress = Tstress;
   Ctangent = Ttangent;

   return 0;
}

int HystereticAsym::revertToLastCommit ()
{

   // State variables
   Tstrain = Cstrain;
   Tstress = Cstress;
   Ttangent = Ctangent;

   return 0;
}

int HystereticAsym::revertToStart ()
{
	// State variables
	Cstrain = 0.0;
	Cstress = 0.0;
	Tstrain = 0.0;
	Tstress = 0.0;

	// History variables
	st = 1.0;
	Ttangent = InitTangent;
	Ctangent = Ttangent;

   

// AddingSensitivity:BEGIN /////////////////////////////////
	if (SHVs != 0) 
		SHVs->Zero();
// AddingSensitivity:END //////////////////////////////////

   return 0;
}

UniaxialMaterial* HystereticAsym::getCopy ()
{
   HystereticAsym* theCopy = new HystereticAsym(this->getTag(), ka, kb, fo, b1, b2, g);

   // State variables
   theCopy->Cstrain = Cstrain;
   theCopy->Cstress = Cstress;
   theCopy->Ctangent = Ctangent;

   theCopy->Tstrain = Tstrain;
   theCopy->Tstress = Tstress;
   theCopy->Ttangent = Ttangent;

   return theCopy;
}

int HystereticAsym::sendSelf (int commitTag, Channel& theChannel)
{
   int res = 0;
   static Vector data(11);
   data(0) = this->getTag();

   // Material properties
   data(1) = ka;
   data(2) = kb;
   data(3) = fo;
   data(4) = b1;
   data(5) = b2;
   data(6) = g;
   
   // State variables from last converged state
   data(7) = Cstrain;
   data(8) = Cstress;
   data(9) = Ctangent;
   data(10) = InitTangent;

   // New variables


   // Data is only sent after convergence, so no trial variables
   // need to be sent through data vector

   res = theChannel.sendVector(this->getDbTag(), commitTag, data);
   if (res < 0) 
      opserr << "HystereticAsym::sendSelf() - failed to send data\n";

   return res;
}

int HystereticAsym::recvSelf (int commitTag, Channel& theChannel,
                                FEM_ObjectBroker& theBroker)
{
   int res = 0;
   static Vector data(11);
   res = theChannel.recvVector(this->getDbTag(), commitTag, data);
  
   if (res < 0) {
      opserr << "HystereticAsym::recvSelf() - failed to receive data\n";
      this->setTag(0);      
   }
   else {
      this->setTag(int(data(0)));

   ka = data(1);
   kb = data(2);
   fo = data(3);
   b1 = data(4);
   b2 = data(5);
   g = data(6);
   
   // State variables from last converged state
   Cstrain = data(7);
   Cstress = data(8);
   Ctangent = data(9);
   InitTangent = data(10);

   // New variables


      // Copy converged state values into trial values
      Tstrain = Cstrain;
      Tstress = Cstress;
      Ttangent = Ctangent;
   }
    
   return res;
}


/*
Response*
HystereticAsym::setResponse(const char **argv, int argc,
	OPS_Stream &theOutput)
{
	Response *theResponse = 0;
	double data = 0.0;
	if (strcmp(argv[0], "stress") == 0) {
		
		theResponse = new MaterialResponse(this, 100, data);
	}
	else if (strcmp(argv[0], "strain") == 0) {
		theResponse = new MaterialResponse(this, 101, data);
	}
	else if (strcmp(argv[0], "tangent") == 0) {
		theResponse = new MaterialResponse(this, 102, data);
	}
	else
		return this->UniaxialMaterial::setResponse(argv, argc, theOutput);

	return theResponse;
}

int
HystereticAsym::getResponse(int responseID, Information &matInfo)
{
	if (responseID == 100) {
		matInfo.theDouble = this->getStress();
	}
	else if (responseID == 101) {
		matInfo.theDouble = this->getStrain();
	}
	else if (responseID == 102) {
		matInfo.theDouble = this->getTangent();
	}
	else
		return this->UniaxialMaterial::getResponse(responseID, matInfo);

	return 0;
}
*/

void HystereticAsym::Print (OPS_Stream& s, int flag)
{
   s << "HystereticAsym tag: " << this->getTag() << endln;
   s << " strain: " << this->getStrain() << endln;
   s << " stress: " << this->getStress() << endln;
   s << " tangent: " << this->getTangent() << endln;
}




// AddingSensitivity:BEGIN ///////////////////////////////////
int
HystereticAsym::setParameter(const char **argv, int argc, Parameter &param)
{

  if (strcmp(argv[0],"ka") == 0 )
    return param.addObject(1, this);
  
  if (strcmp(argv[0],"kb") == 0)
    return param.addObject(2, this);
  
  if (strcmp(argv[0],"fo") == 0)
    return param.addObject(3, this);
  
  if (strcmp(argv[0],"b1") == 0)
    return param.addObject(4, this);
  
  if (strcmp(argv[0],"b2") == 0)
    return param.addObject(5, this);

  if (strcmp(argv[0], "gamma") == 0)
	  return param.addObject(6, this);
  
  return -1;
}



int
HystereticAsym::updateParameter(int parameterID, Information &info)
{
	switch (parameterID) {
	case -1:
		return -1;
	case 1:
		this->ka = info.theDouble;
		break;
	case 2:
		this->kb = info.theDouble;
		break;
	case 3:
		this->fo = info.theDouble;
		break;
	case 4:
		this->b1 = info.theDouble;
		break;
	case 5:
		this->b2 = info.theDouble;
		break;
	case 6:
		this->g = info.theDouble;
		break;
	default:
		return -1;
	}


	
	Ttangent = InitTangent;	// Initial stiffness

	return 0;
}




int
HystereticAsym::activateParameter(int passedParameterID)
{
	parameterID = passedParameterID;

	return 0;
}

double 
HystereticAsym::getStrainSensitivity(int gradNumber) {
	if (SHVs == 0) {
		opserr << "warning:HystereticAsym::getStrainsSensitivity, SHVs =0 " << endln;
		return 0.0;
	}
	double  Tsensitivity = (*SHVs)(0, (gradNumber));
	return Tsensitivity;

}

double
HystereticAsym::getStressSensitivity(int gradIndex, bool conditional)
{
	double gradient = 0.0;
	
	// Initialize return value
	


	// Pick up sensitivity history variables
	Duc = 0.0;
	Dfc = 0.0;
	Dut = 0.0;
	if (SHVs != 0) {
		Duc = (*SHVs)(0,gradIndex);
		Dfc = (*SHVs)(1,gradIndex);
	}


	// Assign values to parameter derivatives (depending on what's random)
	Dka = 0.0;
	Dkb = 0.0;
	Dfo = 0.0;
	Db1 = 0.0;
	Db2 = 0.0;
	Dg = 0.0;

	if (parameterID == 1) {
		Dka = 1.0;
	}
	else if (parameterID == 2) {
		Dkb = 1.0;
	}
	else if (parameterID == 3) {
		Dfo = 1.0;
	}
	else if (parameterID == 4) {
		Db1 = 1.0;
	}
	else if (parameterID == 5) {
		Db2 = 1.0;
	}
	else if (parameterID == 6) {
		Dg = 1.0;
	}

	st = signum(Tstrain-Cstrain);
	fe = -Cstrain * (b1 + b2) + exp(b1 * Cstrain) - exp(-b2 * Cstrain);
	delta2 = (fe - Cstress) * exp(-g * Cstrain) + kb * Cstrain + st * fo;
	uj = Cstrain * st + 1.0 / a * log(st * a * delta2 / (ka - kb));
	delta3 = kb * Tstrain - st * (ka - kb) / a * exp(-a * (Tstrain * st - uj)) + fo * st;

	double Da = 0.25/fo/fo*(2*fo*(Dka-Dkb)-2*Dfo*(ka-kb));
	double Dfe = -Cstrain * (Db1 + Db2) - (b1 + b2) * Duc + (Cstrain * Db1 + b1 * Duc) * exp(b1 * Cstrain) + (Cstrain * b2 + b2 * Duc) * exp(-b2 * Cstrain);
	double Dd2 = (Dfe - Dfc - Dg * Cstrain - g * Duc) * exp(-g * Cstrain) + Dkb * Cstrain + kb * Duc + Dfo * st;
	double Duj = Duc * st + (ka - kb) / delta2 / a / a * ((Da * delta2 + a * Dd2) / (ka - kb) - a * delta2 * (Dka - Dkb) / (ka - kb) / (ka - kb)) - Da / a / a * log(st * a * delta2 / (ka - kb));
	double Dd3 = Dkb * Tstrain + kb * Dut - st / a * (Dka - Dkb - (ka - kb) * (Da * (1 / a + Tstrain * st - uj) + a * (Dut * st - Duj))) * exp(-a * (Tstrain * st - uj)) + Dfo * st;

	gradient = -Tstrain*(Db1+Db2) - (b1+b2)*Dut + (Db1*Tstrain + b1*Dut)*exp(b1 * Tstrain) + (Db2 * Tstrain + b2 * Dut)*exp(-b2 * Tstrain) + Dd3 * exp(g * Tstrain) + delta3 * (Dg * Tstrain + g * Dut)*exp(g * Tstrain);


	
	
	return gradient;
}




double
HystereticAsym::getInitialTangentSensitivity(int gradIndex)
{
	
	// Initialize return value
	double gradient = 0.0;


	// Pick up sensitivity history variables
	Duc = 0.0;
	Dfc = 0.0;
	Dut = 0.0;


	// Assign values to parameter derivatives (depending on what's random)
	Dka = 0.0;
	Dkb = 0.0;
	Dfo = 0.0;
	Db1 = 0.0;
	Db2 = 0.0;
	Dg = 0.0;

	if (parameterID == 1) {
		Dka = 1.0;
	}
	else if (parameterID == 2) {
		Dkb = 1.0;
	}
	else if (parameterID == 3) {
		Dfo = 1.0;
	}
	else if (parameterID == 4) {
		Db1 = 1.0;
	}
	else if (parameterID == 5) {
		Db2 = 1.0;
	}
	else if (parameterID == 6) {
		Dg = 1.0;
	}


	st = 1;

	gradient = 0.5*(Dka+Dkb);

	return gradient;
}


int
HystereticAsym::commitSensitivity(double TstrainSensitivity, int gradIndex, int numGrads)
{
	if (SHVs == 0) {
		SHVs = new Matrix(2,numGrads);
	}


	// Initialize return value
	double gradient = 0.0;


	// Pick up sensitivity history variables
	Duc = 0.0;
	Dfc = 0.0;
	Dut = TstrainSensitivity;
	if (SHVs != 0) {
		Duc = (*SHVs)(0, gradIndex);
		Dfc = (*SHVs)(1, gradIndex);
	}


	// Assign values to parameter derivatives (depending on what's random)
	Dka = 0.0;
	Dkb = 0.0;
	Dfo = 0.0;
	Db1 = 0.0;
	Db2 = 0.0;
	Dg = 0.0;

	if (parameterID == 1) {
		Dka = 1.0;
	}
	else if (parameterID == 2) {
		Dkb = 1.0;
	}
	else if (parameterID == 3) {
		Dfo = 1.0;
	}
	else if (parameterID == 4) {
		Db1 = 1.0;
	}
	else if (parameterID == 5) {
		Db2 = 1.0;
	}
	else if (parameterID == 6) {
		Dg = 1.0;
	}


	st = signum(Tstrain - Cstrain);
	fe = -Cstrain * (b1 + b2) + exp(b1 * Cstrain) - exp(-b2 * Cstrain);
	delta2 = (fe - Cstress) * exp(-g * Cstrain) + kb * Cstrain + st * fo;
	uj = Cstrain * st + 1.0 / a * log(st * a * delta2 / (ka - kb));
	delta3 = kb * Tstrain - st * (ka - kb) / a * exp(-a * (Tstrain * st - uj)) + fo * st;

	double Da = 0.25 / fo / fo * (2 * fo * (Dka - Dkb) - 2 * Dfo * (ka - kb));
	double Dfe = -Cstrain * (Db1 + Db2) - (b1 + b2) * Duc + (Cstrain * Db1 + b1 * Duc) * exp(b1 * Cstrain) + (Cstrain * b2 + b2 * Duc) * exp(-b2 * Cstrain);
	double Dd2 = (Dfe - Dfc - Dg * Cstrain - g * Duc) * exp(-g * Cstrain) + Dkb * Cstrain + kb * Duc + Dfo * st;
	double Duj = Duc * st + (ka - kb) / delta2 / a / a * ((Da * delta2 + a * Dd2) / (ka - kb) - a * delta2 * (Dka - Dkb) / (ka - kb) / (ka - kb)) - Da / a / a * log(st * a * delta2 / (ka - kb));
	double Dd3 = Dkb * Tstrain + kb * Dut - st / a * (Dka - Dkb - (ka - kb) * (Da * (1 / a + Tstrain * st - uj) + a * (Dut * st - Duj))) * exp(-a * (Tstrain * st - uj)) + Dfo * st;

	gradient = -Tstrain * (Db1 + Db2) - (b1 + b2) * Dut + (Db1 * Tstrain + b1 * Dut) * exp(b1 * Tstrain) + (Db2 * Tstrain + b2 * Dut) * exp(-b2 * Tstrain) + Dd3 * exp(g * Tstrain) + delta3 * (Dg * Tstrain + g * Dut) * exp(g * Tstrain);


	// Commit 
	(*SHVs)(0,gradIndex) = TstrainSensitivity;
	(*SHVs)(1,gradIndex) = gradient;
	
	return 0;
}

// AddingSensitivity:END /////////////////////////////////////////////


double 
HystereticAsym::signum(double value)
{
	if (value > 0.0) {
		return 1.0;
	}
	else if (value == 0){
		return 1.0;
	}
	else {
		return -1.0;
	}
}