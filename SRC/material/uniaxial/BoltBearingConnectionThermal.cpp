
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

// $Revision: 1.10 $
// $Date: 2008-08-26 16:30:55 $
// $Source: /usr/local/cvs/OpenSees/SRC/material/uniaxial/BoltBearingConnectionThermal.cpp,v $

//
// Description: This file contains the class implementation for 
// PlateBearingConnectionThermal. 



#include <BoltBearingConnectionThermal.h>
#include <Vector.h>
#include <Channel.h>
#include <Information.h>
#include <Parameter.h>
#include <string.h>
#include <math.h> 
#include <OPS_Globals.h> 

#include <elementAPI.h>

void*
OPS_BoltBearingConnectionThermal(void)
{
	// Pointer to a uniaxial material that will be returned
	UniaxialMaterial* theMaterial = 0;
	if (OPS_GetNumRemainingInputArgs() < 5) {
		opserr << "Invalid #args,  want: uniaxialMaterial BoltBearingConnectionThermal tag? Fu? E? db? As? Temp?" << endln;
		return 0;
	}

	int iData[1];
	double dData[5];
	int numData = 1;
	if (OPS_GetIntInput(&numData, iData) != 0) {
		opserr << "WARNING invalid tag for uniaxialMaterial BoltBearingConnectionThermal" << endln;
		return 0;
	}

	numData = 5;
	if (OPS_GetDoubleInput(&numData, dData) != 0) {
		opserr << "WARNING invalid input for uniaxialMaterial BoltBearingConnectionThermal" << endln;
		return 0;
	}

	// Parsing was successful, allocate the material
	theMaterial = new BoltBearingConnectionThermal(iData[0], dData[0], dData[1], dData[2], dData[3], dData[4]);
	if (theMaterial == 0) {
		opserr << "WARNING could not create uniaxialMaterial of type BoltBearingConnectionThermal\n";
		return 0;
	}

	return theMaterial;
}


BoltBearingConnectionThermal::BoltBearingConnectionThermal(int tag, double F, double Enot, double dbolt, double Abolt, double temp)
:UniaxialMaterial(tag, MAT_TAG_BoltBearingConnectionThermal),
	trialStrain(0.0), Fu0(F), E0(Enot), db(dbolt), As(Abolt), Temp(temp), committedStress(0.0), committedStrain(0.0)
{


}

// Default Constructor: Behavoir for default object
BoltBearingConnectionThermal::BoltBearingConnectionThermal()
	:UniaxialMaterial(0, MAT_TAG_BoltBearingConnectionThermal),

	trialStrain(0.0), Fu0(0.0), E0(0.0), db(1.0), As(1.0), Temp(20), committedStress(0.0), committedStrain(0.0)	
{

}

// Destructor
BoltBearingConnectionThermal::~BoltBearingConnectionThermal()
{
	// does nothing
}

// Set trial srain =  disp
int
BoltBearingConnectionThermal::setTrialStrain(double strain, double strainrate)
{
	trialStrain = strain;
	return 0;
}


double
BoltBearingConnectionThermal::getStrain(void)
{
	return trialStrain;
}



/*
interpolate a 1 d array of pairs [x1,y1,x2,y2,...]
given xTarget find the interpolated y value.
1d array of pairs b/c its easier to test this way, in one malloc compared to many.
Author: github@lozzoc, may 17,2020;
*/
double Interpolate(double* array, int numpairs, double searchTerm) {

	int mid, low = 0, high = numpairs - 1;
	double x1, /*searchTerm,*/ x3, y1, y3;

	if (array[0] >= searchTerm) {
		return array[1];
	}
	else if (array[(numpairs - 1) * 2] <= searchTerm) {
		return array[numpairs * 2 - 1];
	}
	else {


		while (low <= high) {

			mid = (low + high) / 2;
			if (array[mid * 2] == searchTerm) {
				return array[mid * 2 + 1];

				// based on the current pointer [mid*2] we either have a
			}
			else if (array[mid * 2] < searchTerm) {
				// if we have sandwhich facing right
				if (array[2 * (mid + 1)] > searchTerm) {
					x1 = array[2 * mid];
					x3 = array[2 * mid + 2];
					y1 = array[2 * mid + 1];
					y3 = array[2 * mid + 3];
					return (searchTerm - x1) * (y3 - y1) / (x3 - x1) + y1;
				}
				else {
					low = mid + 1;
				}

			}
			else {
				// or we might have a sandwhich faceing left
				if (array[2 * (mid - 1)] < searchTerm) {
					x1 = array[2 * mid - 2];
					x3 = array[2 * mid];
					y1 = array[2 * mid - 1];
					y3 = array[2 * mid + 1];
					return (searchTerm - x1) * (y3 - y1) / (x3 - x1) + y1;
				}
				else {
					high = mid - 1;

				}
			}
		}
		//end while loop
		return  0.0;
	}

}

double
BoltBearingConnectionThermal::getStress(void)
{
	// Thermal Reduction factors based on Eurocode

// Ultimate strength of steel, Fu - (x, y) pairs - (temperature, reduction factor)
	double FuRedfactor[18] = { 20,1,100,0.968,200, 0.935,300, 0.903 ,400, 0.8,500, 0.624 ,600, 0.376,700, 0.184, 800, 0.08 };
	double Fu = 0.0; // initialize Fu variable
	Fu = Fu0 * Interpolate(FuRedfactor, 9, Temp);

	// Young's modulus of elasticity, E - (x, y) pairs - (temperature, reduction factor)
	double ERedfactor[18] = { 20,1,100, 1,200, 0.9,300, 0.8 ,400, 0.7,500, 0.6,600, 0.31,700, 0.13,800, 0.09 };
	double E = 0.0; // initialize E variable
	E = E0 * Interpolate(ERedfactor, 9, Temp);

	// Bolt reuction factor, Zeta  - (x, y) pairs - (temperature, reduction factor)
	double zeta = 0.0; // initialize Zeta variable
	double zetafactor[20] = { 20,2.5,100, 2.8,200, 2.0 ,300, 2.2,400, 2.0 ,500, 2.0,600, 1.3,700, 0.6,800, 0.7,900, 0.02 };
	zeta = Interpolate(zetafactor, 10, Temp);

	// Bolt hole 
	double d = 0.0;
	d = db + (25.4 / 16);	// d = diameter of bolt hole  = diameter of bolt + 1/16"

	// initialize variables
	double Fvrd = 0;
	double Kv = 0;
	double dStress = 0;
	double dStrain = 0;
	double G = 0;
	double stressMax = 0.0;
	double strain1 = 0.0;
	double Fv = 0.0;
	double delta = 0.0;

	const int NUM_STRAIN_TRIALS = 100;
	double stressStrain[NUM_STRAIN_TRIALS * 2];


	//double stressStrain[100][2];

	// initialize model parameters
	G = E / (2.0 * (1.3)); // shear modulus
	Fvrd = 0.6 * Fu * As; // max stress
	Kv = (0.15 * G * As) / db; // shear stiffness
	delta = (Fvrd / Kv) + zeta;

	double Nsteps = 100;
	if (abs(trialStrain) < delta) {

		for (int i = 0; i < Nsteps; i++) {
			Fv = Fvrd * ((double)i / 100);
			stressStrain[i * 2] = (Fv / Kv) + zeta * pow((Fv / Fvrd), 6); // Strain
			stressStrain[i * 2 + 1] = Fv; // Stress
		}

		// Interpolate Stress from Strain via look-up table "stressStrain"
		// double Interpolate(double* array, int numpairs, double searchTerm)
		trialStress = Interpolate(stressStrain, NUM_STRAIN_TRIALS, fabs(trialStrain));



		dStress =  trialStress - committedStress;
		dStrain = trialStrain - committedStrain;


		if (dStrain != 0.0)
			trialTangent = fabs(dStress / dStrain);
		else
			trialTangent = this->getInitialTangent();

		if (trialStrain < 0) {
			return (trialStress * -1);
		}
		else
			return (trialStress);

	}

	else if (abs(trialStrain) >= (delta) && abs(trialStrain) < db / 2) {

		if (trialStrain < 0) {
			return (Fvrd * -1);
		}
		else
			return (Fvrd);
	}
	else if (fabs(trialStrain) >= db / 2 && fabs(trialStrain) < db) {

		trialStress = Fvrd + (trialStrain - (db / 2)) * ((0.0 - Fvrd) / (db / 2));

		if (dStrain != 0.0)
			trialTangent = fabs((0.0 - Fvrd) / (db / 2));
		else
			trialTangent = this->getInitialTangent();

		if (trialStrain < 0) {
			return (trialStress * -1);
		}
		else
			return (trialStress);
	}
	else if (fabs(trialStrain) >= db) {
		return 0.0;
	}
}


double
BoltBearingConnectionThermal::getTangent(void)
{
	//Tangent at current state
	return trialTangent;
	//return this->getInitialTangent();
}


double
BoltBearingConnectionThermal::getInitialTangent(void)
{
	return 173499.7406;
	

}

// Convergence reached, commit values 

int
BoltBearingConnectionThermal::commitState(void)
{
	//committedStrain = trialStrain;
	//committedStress = trialStress;

	return 0;
}


int
BoltBearingConnectionThermal::revertToLastCommit(void)
{
	trialStrain = committedStrain;
	trialStress = committedStress;
	return 0;
}


int
BoltBearingConnectionThermal::revertToStart(void)
{
	trialStrain = 0.0;
	trialStress = 0.0;
	committedStrain = 0.0;
	committedStress = 0.0;
	return 0;
}

// copy for every instance of material in model
UniaxialMaterial*
BoltBearingConnectionThermal::getCopy(void)
{
	BoltBearingConnectionThermal* theCopy = new BoltBearingConnectionThermal(this->getTag(),Fu0, E0, db, As, Temp);
	theCopy->trialStrain = trialStrain;
	theCopy->trialStress = trialStress;
	theCopy->committedStrain = committedStrain;
	theCopy->committedStress = committedStress;
	return theCopy;
}


int
BoltBearingConnectionThermal::sendSelf(int cTag, Channel& theChannel)
{
	int res = 0;
	/*
	static Vector data(6);
	data(0) = this->getTag();
	data(1) = Epos;
	data(2) = Eneg;
	data(3) = eta;
	data(4) = committedStrain;
	data(5) = committedStress;
	res = theChannel.sendVector(this->getDbTag(), cTag, data);
	if (res < 0)
		opserr << "PlateBearingConnectionThermal::sendSelf() - failed to send data\n";
		*/
	return res;
}


int
BoltBearingConnectionThermal::recvSelf(int cTag, Channel& theChannel,
	FEM_ObjectBroker& theBroker)
{
	int res = 0;
	/*
	static Vector data(6);
	res = theChannel.recvVector(this->getDbTag(), cTag, data);

	if (res < 0) {
		opserr << "PlateBearingConnectionThermal::recvSelf() - failed to receive data\n";
		Epos = Eneg = 0;
		this->setTag(0);
	}
	else {
		this->setTag(int(data(0)));
		Epos = data(1);
		Eneg = data(2);
		eta = data(3);
		committedStrain = data(4);
		committedStrainRate = data(5);
		this->revertToLastCommit();
	}
	*/
	return res;
}


void
BoltBearingConnectionThermal::Print(OPS_Stream& s, int flag)
{
	s << "BoltBearingConnectionThermal tag: " << this->getTag() << endln;
	s <<  " Fu: " << Fu0 << " E: " << E0 << " db: " << db << " As: " << As << " temp: " << Temp << endln;
}


