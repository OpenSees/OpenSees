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
// $Source: /usr/local/cvs/OpenSees/SRC/material/uniaxial/HystereticSmooth.cpp,v $

// Written: Salvatore Sessa Mail: salvatore.sessa2@unina.it
// Created: 04-2022
// Revision: A
//
// Description: This file contains the class implementation for 
// HystereticSmooth. 
//
// What: "@(#) HystereticSmooth.C, revA"


#include <HystereticSmooth.h>
#include <Vector.h>
#include <Matrix.h>
#include <Channel.h>
#include <Information.h>
#include <Parameter.h>

#include <math.h>
#include <float.h>

#include <elementAPI.h>
#include <MaterialResponse.h>



void*
OPS_HystereticSmooth()
{
	UniaxialMaterial* theMaterial = 0;
	int numArgs = OPS_GetNumRemainingInputArgs();
	if (numArgs < 5) {
		opserr << "WARNING insufficient arguments\n";
		opserr << "Want: uniaxialMaterial HystereticSmooth tag? ka? kb? fbar? beta? <-alpha>" << endln;
		return 0;
	}


	int iData[1];
	double dData[4];



	int numData = 1;
	if (OPS_GetIntInput(&numData, iData) != 0) {
		opserr << "WARNING invalid tag for uniaxialMaterial HystereticSmooth" << endln;
		return 0;
	}

	numData = 4;
	if (OPS_GetDoubleInput(&numData, dData) != 0) {
		opserr << "Invalid data for uniaxial HystereticSmooth " << iData[0] << endln;
		return 0;
	}
	if (dData[0] <= 0.0) {
		opserr << "uniaxialMaterial HystereticSmooth ka must be positive" << endln;
		return 0;
	}
	if (dData[1] >= dData[0]) {
		opserr << "uniaxialMaterial HystereticSmooth kb must be < ka" << endln;
		return 0;
	}
	if (dData[2] <= 0.0) {
		opserr << "uniaxialMaterial HystereticSmooth fbar and/or alpha must be positive" << endln;
		return 0;
	}

	if (numArgs == 6) {
		const char* flag = OPS_GetString();

		if (strcmp(flag, "-alpha") == 0) {
			double newFbar = (dData[0] - dData[1]) * 0.5 / dData[2];
			dData[2] = newFbar;
		}
		else {
			opserr << "WARNING unknown parameter flag, ignoring it for uniaxialMaterial HystereticSmooth " << iData[0] << endln;
	}
	}

	// Parsing was successful, allocate the material
	theMaterial = new HystereticSmooth(iData[0], dData[0], dData[1], dData[2], dData[3]);
	if (theMaterial == 0) {
		opserr << "WARNING could not create uniaxialMaterial of type HystereticSmooth\n";
		return 0;
	}

	return theMaterial;

}

HystereticSmooth::HystereticSmooth(int tag, double K1, double K2, double A, double B)
	:UniaxialMaterial(tag, MAT_TAG_HystereticSmooth),
	ka(K1), kb(K2), Fbar(A), beta(B)
{
	// Sets all history and state variables to initial values

	 // Parameters:
	kab = ka - kb;
	

	a = (ka - kb) * 0.5 / Fbar;
	

	InitTangent = kb - 2 * beta + a  * Fbar;

	// History variables

	st = 1.0;




	// State variables
	Cstrain = 0.0;
	Cstress = 0.0;
	Ctangent = InitTangent;


	Tstrain = 0.0;
	Tstress = 0.0;
	Ttangent = Ctangent;


	// AddingSensitivity:BEGIN /////////////////////////////////////
	parameterID = 0;
	SHVs = 0;
	// AddingSensitivity:END //////////////////////////////////////
}

HystereticSmooth::HystereticSmooth() :UniaxialMaterial(0, MAT_TAG_HystereticSmooth),
ka(0.0), kb(0.0), Fbar(0.0), beta(0.0)
{

	// AddingSensitivity:BEGIN /////////////////////////////////////
	parameterID = 0;
	SHVs = 0;
	// AddingSensitivity:END //////////////////////////////////////

}

HystereticSmooth::~HystereticSmooth()
{
	// AddingSensitivity:BEGIN /////////////////////////////////////
	if (SHVs != 0)
		delete SHVs;
	// AddingSensitivity:END //////////////////////////////////////
}

int HystereticSmooth::setTrialStrain(double strain, double strainRate)
{
	// Determine change in strain from last converged state

	Tstrain = strain;
	dStrain = strain - Cstrain;
	st = signum(dStrain);

	
	uj = Cstrain * st + 1.0 / a * log(a * st / kab * (-2 * beta * Cstrain + exp(beta * Cstrain) - exp(-beta * Cstrain) +kb*Cstrain + Fbar * st - Cstress));
	arg1 = -a * (Tstrain * st - uj);
	Tstress = -2 * beta * Tstrain + exp(beta * Tstrain) - exp(-1.0 * beta * Tstrain) + kb * Tstrain - st * kab / a * exp(arg1) + Fbar*st;
	Ttangent = kb - 2.0 * beta + beta * exp(beta * Tstrain) + beta * exp(-beta * Tstrain) + kab * exp(arg1);



	return 0;
}



double HystereticSmooth::getStrain()
{
	return Tstrain;
}

double HystereticSmooth::getStress()
{
	return Tstress;
}

double HystereticSmooth::getTangent()
{
	return Ttangent;
}
double HystereticSmooth::getInitialTangent()
{
	return InitTangent;
}

int HystereticSmooth::commitState()
{
	// State variables
	Cstrain = Tstrain;
	Cstress = Tstress;
	Ctangent = Ttangent;

	return 0;
}

int HystereticSmooth::revertToLastCommit()
{

	// State variables
	Tstrain = Cstrain;
	Tstress = Cstress;
	Ttangent = Ctangent;

	return 0;
}

int HystereticSmooth::revertToStart()
{
	uj = - 1.0 / a * log(a / kab * Fbar );



	// State variables
	Cstrain = 0.0;
	Cstress = 0.0;
	Ctangent = InitTangent;


	Tstrain = 0.0;
	Tstress = 0.0;
	Ttangent = Ctangent;



	// AddingSensitivity:BEGIN /////////////////////////////////
	if (SHVs != 0)
		SHVs->Zero();
	// AddingSensitivity:END //////////////////////////////////

	return 0;
}

UniaxialMaterial* HystereticSmooth::getCopy()
{
	HystereticSmooth* theCopy = new HystereticSmooth(this->getTag(), ka, kb, Fbar, beta);

	// State variables
	theCopy->Cstrain = Cstrain;
	theCopy->Cstress = Cstress;
	theCopy->Ctangent = Ctangent;

	theCopy->Tstrain = Tstrain;
	theCopy->Tstress = Tstress;
	theCopy->Ttangent = Ttangent;

	return theCopy;
}

int HystereticSmooth::sendSelf(int commitTag, Channel& theChannel)
{
	int res = 0;
	static Vector data(9);
	data(0) = this->getTag();

	// Material properties
	data(1) = ka;
	data(2) = kb;
	data(3) = a;
	data(4) = Fbar;
	data(5) = beta;

	// State variables from last converged state
	data(6) = Cstrain;
	data(7) = Cstress;
	data(8) = Ctangent;



	// Data is only sent after convergence, so no trial variables
	// need to be sent through data vector

	res = theChannel.sendVector(this->getDbTag(), commitTag, data);
	if (res < 0)
		opserr << "HystereticSmooth::sendSelf() - failed to send data\n";

	return res;
}

int HystereticSmooth::recvSelf(int commitTag, Channel& theChannel,
	FEM_ObjectBroker& theBroker)
{
	int res = 0;
	static Vector data(9);
	res = theChannel.recvVector(this->getDbTag(), commitTag, data);

	if (res < 0) {
		opserr << "HystereticSmooth::recvSelf() - failed to receive data\n";
		this->setTag(0);
	}
	else {
		this->setTag(int(data(0)));

		ka = data(1);
		kb = data(2);
		a = data(3);
		Fbar = data(4);
		beta = data(5);

		// State variables from last converged state
		Cstrain = data(6);
		Cstress = data(7);
		Ctangent = data(8);


		// Copy converged state values into trial values
		Tstrain = Cstrain;
		Tstress = Cstress;
		Ttangent = Ctangent;

		//auxiliary variables
		kab = ka - kb;
		a = (ka - kb) * 0.5 / Fbar;
		uj = -Cstrain * st - 1.0 / a * log(a * st / kab * (-2 * beta * Cstrain + exp(beta * Cstrain) - exp(-beta * Cstrain) + Fbar * st - Cstress));
		InitTangent = kb - 2 * beta + a * Fbar;
		st = 1.0;
	}

	return res;
}



void HystereticSmooth::Print(OPS_Stream& s, int flag)
{
	s << "HystereticSmooth tag: " << this->getTag() << endln;
	s << " strain: " << this->getStrain() << endln;
	s << " stress: " << this->getStress() << endln;
	s << " tangent: " << this->getTangent() << endln;
}




// AddingSensitivity:BEGIN ///////////////////////////////////
int
HystereticSmooth::setParameter(const char** argv, int argc, Parameter& param)
{

	if (strcmp(argv[0], "ka") == 0)
		return param.addObject(1, this);

	if (strcmp(argv[0], "kb") == 0)
		return param.addObject(2, this);

	if (strcmp(argv[0], "fbar") == 0)
		return param.addObject(3, this);

	if (strcmp(argv[0], "beta") == 0)
		return param.addObject(4, this);

	return -1;
}



int
HystereticSmooth::updateParameter(int parameterID, Information& info)
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
		this->Fbar= info.theDouble;
		break;
	case 4:
		this->beta = info.theDouble;
		break;
	default:
		return -1;
	}

	Ttangent = InitTangent;          // Initial stiffness

	return 0;
}




int
HystereticSmooth::activateParameter(int passedParameterID)
{
	parameterID = passedParameterID;

	return 0;
}

double
HystereticSmooth::getStrainSensitivity(int gradNumber) {
	if (SHVs == 0) {
		opserr << "warning:HystereticSmooth::getStrainsSensitivity, SHVs =0 " << endln;
		return 0.0;
	}
	double  Tsensitivity = (*SHVs)(0, (gradNumber));
	return Tsensitivity;

}

double
HystereticSmooth::getStressSensitivity(int gradIndex, bool conditional)
{
	double gradient = 0.0;

	// Initialize return value



	// Pick up sensitivity history variables
	Duc = 0.0;
	Dfc = 0.0;
	Dut = 0.0;
	if (SHVs != 0) {
		Duc = (*SHVs)(0, gradIndex);
		Dfc = (*SHVs)(1, gradIndex);
	}


	// Assign values to parameter derivatives (depending on what's random)
	Dka = 0.0;
	Dkb = 0.0;
	Dfbar = 0.0;
	Da = 0.0;
	Dbeta = 0.0;

	if (parameterID == 1) {
		Dka = 1.0;
	}
	else if (parameterID == 2) {
		Dkb = 1.0;
	}
	else if (parameterID == 3) {
		Dfbar = 1.0;
	}
	else if (parameterID == 4) {
		Dbeta = 1.0;
	}



	// This is the Tstress derivative with respect to Tstrain
	uj = Cstrain * st + 1.0 / a * log(a * st / kab * (-2 * beta * Cstrain + exp(beta * Cstrain) - exp(-beta * Cstrain) + kb * Cstrain + Fbar * st - Cstress));
	Da = (2 * Fbar * (Dka - Dkb) - 2 * Dfbar * (ka - kb)) / (4 * Fbar*Fbar);
	double arg = a * st / (ka - kb) * (2 * beta * Cstrain + exp(beta * Cstrain) - exp(-beta * Cstrain) + kb * Cstrain + Fbar * st - Cstress);
	double Darg = Da / a * arg - (Dka - Dkb) / (ka - kb) * arg + a * st / (ka - kb) * ((Dbeta * Cstrain + beta * Duc) * (2 + exp(beta * Cstrain) + exp(-beta * Cstrain)) + Dkb * Cstrain + kb * Duc + Dfbar * st - Dfc);
	double Duj = Duc * st - Darg/arg/a + Da*log(arg)/a/a;
	gradient = (Dbeta * Tstrain + beta * Dut) * (exp(beta * Tstrain) + exp(-1.0 * beta * Tstrain) - 2) + Dkb * Tstrain + kb * Dut + (-st * ((Dka - Dkb) * a - (ka - kb) * Da) / a / a + st * (ka - kb) / a * (Da * (Tstrain * st - uj) + a * (Dut * st - Duj))) * exp(-a * (Tstrain * st - uj)) + Dfbar * st;

	
	return gradient;

	
}




double
HystereticSmooth::getInitialTangentSensitivity(int gradIndex)
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
	Dfbar = 0.0;
	Dbeta = 0.0;

	if (parameterID == 1) {
		Dka = 1.0;
	}
	else if (parameterID == 2) {
		Dkb = 1.0;
	}
	else if (parameterID == 3) {
		Dfbar = 1.0;
	}
	else if (parameterID == 4) {
		Dbeta = 1.0;
	}
	


	Da = (2 * Fbar * (Dka - Dkb) - 2 * Dfbar * (ka - kb)) / (4 * Fbar * Fbar);
	gradient = Dkb - 2 * Dbeta + Da * Fbar + a*Dfbar ;

		return gradient;
}


int
HystereticSmooth::commitSensitivity(double TstrainSensitivity, int gradIndex, int numGrads)
{
	if (SHVs == 0) {
		SHVs = new Matrix(2, numGrads);
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
	Dfbar = 0.0;
	Dbeta = 0.0;

	if (parameterID == 1) {
		Dka = 1.0;
	}
	else if (parameterID == 2) {
		Dkb = 1.0;
	}
	else if (parameterID == 3) {
		Dfbar = 1.0;
	}
	else if (parameterID == 4) {
		Dbeta = 1.0;
	}


	// This is the Tstress derivative with respect to Tstrain

	uj = Cstrain * st + 1.0 / a * log(a * st / kab * (-2 * beta * Cstrain + exp(beta * Cstrain) - exp(-beta * Cstrain) + kb * Cstrain + Fbar * st - Cstress));
	Da = (2 * Fbar * (Dka - Dkb) - 2 * Dfbar * (ka - kb)) / (4 * Fbar * Fbar);
	double arg = a * st / (ka - kb) * (2 * beta * Cstrain + exp(beta * Cstrain) - exp(-beta * Cstrain) + kb * Cstrain + Fbar * st - Cstress);
	double Darg = Da / a * arg - (Dka - Dkb) / (ka - kb) * arg + a * st / (ka - kb) * ((Dbeta * Cstrain + beta * Duc) * (2 + exp(beta * Cstrain) + exp(-beta * Cstrain)) + Dkb * Cstrain + kb * Duc + Dfbar * st - Dfc);
	double Duj = Duc * st - Darg / arg / a + Da * log(arg) / a / a;
	gradient = (Dbeta * Tstrain + beta * Dut) * (exp(beta * Tstrain) + exp(-1.0 * beta * Tstrain) - 2) + Dkb * Tstrain + kb * Dut + (-st * ((Dka - Dkb) * a - (ka - kb) * Da) / a / a + st * (ka - kb) / a * (Da * (Tstrain * st - uj) + a * (Dut * st - Duj))) * exp(-a * (Tstrain * st - uj)) + Dfbar * st;


	// Commit 
	(*SHVs)(0, gradIndex) = TstrainSensitivity;
	(*SHVs)(1, gradIndex) = gradient;
	
	return 0;
}

// AddingSensitivity:END /////////////////////////////////////////////


double
HystereticSmooth::signum(double value)
{
	if (value > 0.0) {
		return 1.0;
	}
	else if (value == 0) {
		return 1.0;
	}
	else {
		return -1.0;
	}
}