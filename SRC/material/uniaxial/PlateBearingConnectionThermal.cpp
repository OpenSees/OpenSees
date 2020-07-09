
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
// $Source: /usr/local/cvs/OpenSees/SRC/material/uniaxial/PlateBearingConnectionThermal.cpp,v $
                                                                        
//
// Description: This file contains the class implementation for 
// PlateBearingConnectionThermal. 
 


#include <PlateBearingConnectionThermal.h>
#include <Vector.h>
#include <Channel.h>
#include <Information.h>
#include <Parameter.h>
#include <string.h>
#include <math.h> 
#include <OPS_Globals.h> 

#include <elementAPI.h>

void *
OPS_PlateBearingConnectionThermal(void)
{
  // Pointer to a uniaxial material that will be returned
  UniaxialMaterial *theMaterial = 0;
  if (OPS_GetNumRemainingInputArgs() < 6) {
    opserr << "Invalid #args,  want: uniaxialMaterial PlateBearingConnectionThermal tag? Fy? Fu? E? db? e2? t?" << endln;
    return 0;
  }
  
  int iData[1];
  double dData[6];
  int numData = 1;
  if (OPS_GetIntInput(&numData, iData) != 0) {
    opserr << "WARNING invalid tag for uniaxialMaterial PlateBearingConnectionThermal" << endln;
    return 0;
  }

  numData = 6;
  if (OPS_GetDoubleInput(&numData, dData) != 0) {
	  opserr << "WARNING invalid input for uniaxialMaterial PlateBearingConnectionThermal" << endln;
	  return 0;
  }

  // Parsing was successful, allocate the material
  theMaterial = new PlateBearingConnectionThermal(iData[0], dData[0], dData[1], dData[2], dData[3],  dData[4], dData[5]);
  if (theMaterial == 0) {
    opserr << "WARNING could not create uniaxialMaterial of type PlateBearingConnectionThermal\n";
    return 0;
  }

  return theMaterial;
}


PlateBearingConnectionThermal::PlateBearingConnectionThermal(int tag, double Fy, double Fu, double E, double d, double x, double thickness)
:UniaxialMaterial(tag,MAT_TAG_PlateBearingConnectionThermal),
 trialStrain(0.0), trialStress(0.0), trialTangent(0.0), committedStrain(0.0), committedStress(0.0), committedTangent(0.0),
 Fu0(Fu), Fy0(Fy), E0(E), e2(x), t(thickness), db(d), Temp(20.0)
{
	

	
}

// Default Constructor: Behavoir for default object
PlateBearingConnectionThermal::PlateBearingConnectionThermal()
:UniaxialMaterial(0,MAT_TAG_PlateBearingConnectionThermal),
 trialStrain(0.0), trialStress(0.0), trialTangent(0.0), committedStrain(0.0), committedStress(0.0), committedTangent(0.0),
 Fu0(0.0), Fy0(0.0), E0(0.0), e2(1.0), t(1.0), db(1.0), Temp(20.0)
{
	
}

// Destructor
PlateBearingConnectionThermal::~PlateBearingConnectionThermal()
{
	// does nothing
}

// Set trial srain =   disp
int
PlateBearingConnectionThermal::setTrialStrain(double strain, double strainrate)
{

	trialStrain = strain;

	return 0;
}


double
PlateBearingConnectionThermal::getStrain(void)
{
	return trialStrain;
}

double
PlateBearingConnectionThermal::getStress(void)
{
	// Thermal Reduction factors based on Eurocode 


	// Ω & ψ & phi
	double omegafactor;
	double psyfactor;
	double phifactor;

	if (e2 <= 2 * db) {
		//double omegafactor[9] = { 145.0, 180.0, 180.0 , 180.0, 170.0 , 130.0, 80.0, 45.0  , 20.0 };
		//double psyfactor[9] = { 2.1 , 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 1.8 };
		//double phifactor[9] = { 0.012, 0.008, 0.008  , 0.008 , 0.008 , 0.008 , 0.008 , 0.008 , 0.008 };

		omegafactor = 145.0;

		psyfactor = 2.1;

		phifactor = 0.012;

	}
	else if (e2 >= 3 * db && db <= 20) {
		//double omegafactor[9] = { 250 , 220 , 220  , 220 , 200 , 170, 110.0, 40.0  , 20.0 };
		//double psyfactor[9] = { 1.7 , 1.7  , 1.7   , 1.7  , 1.7  , 1.7 , 1.7 , 1.7   , 1.7 };
		//double phifactor[9] = { 0.008  , 0.008 , 0.008 , 0.008 , 0.008 , 0.008 , 0.008 , 0.007, 0.007 };

		omegafactor = 250.0;
		psyfactor = 1.7;
		phifactor = 0.008;

	}
	else if (e2 >= 3 * db && db == 24) {
		//double omegafactor[9] = { 250 , 250 , 250  , 250 , 200 , 170, 110.0, 40.0  , 20.0 };
		//double psyfactor[9] = { 1.7 , 1.7  , 1.7   , 1.7  , 1.7  , 1.7 , 1.7 , 1.7   , 1.7 };
		//double phifactor[9] = { 0.011  , 0.011  , 0.011  , 0.011 , 0.009 , 0.007 , 0.0055 , 0.0055, 0.001 };

		omegafactor = 250.0;

		psyfactor = 1.7;

		phifactor = 0.011;

	}
	else {
		opserr << "BoltBearingConnectionThermal " << this->getTag() << "Connection geometery not applicable" << endln;
	}


	
	double Fbrd = 0;
	double Kbr = 0;
	double Kb = 0;
	double Kv = 0;
	double Ki = 0;
	double delta = 0;
	double dStress = 0;
	double dstrain = 0;
	double G = 0;
	double maxstress = 0.0;
	double deltaMax = 0.0;
	double stressMax = 0.0;

	dstrain = trialStrain - committedStrain;
	opserr << "trialStrain = " << trialStrain << ", committedStrain = " << committedStrain << endln;
	// min ( e2 & 2.76*
	if (e2 < 2.76 * db) { Fbrd = e2 * Fu0 * t; }

	else { Fbrd = 2.76 * db * Fu0 * t; }

	G = E0 / (2.0 * (1.3));
	Kbr = omegafactor * t * Fy0 * pow((db / 25.4), 0.8);
	Kb = 32 * E0 * t * pow(((e2 / db) - 0.5), 3);
	Kv = 6.67 * G * t;
	Ki = 1 / ((1 / Kbr) + (1 / Kb) + (1 / Kv));

	deltaMax = ((db / 2) * (Ki / Fbrd));
	stressMax = (psyfactor * deltaMax / pow((1 + sqrt(deltaMax)), 2) - phifactor * deltaMax) * Fbrd;

	opserr << "trial strain = " << trialStrain << endln;
	if (trialStrain < 0.5 * db) {

		delta = (fabs(dstrain) * (Ki / Fbrd));
		opserr << "dstrain = " << dstrain << endln;
		dStress = (psyfactor * delta / pow((1 + sqrt(delta)), 2) - phifactor * delta) * Fbrd;
		opserr << "dstress = " << dStress << endln;
		opserr << "committedStress = " << committedStress << endln;
		trialStress = dStress + committedStress;
		opserr << "trial stress = " << trialStress << endln;

		if (dstrain != 0.0)
		  trialTangent = (dStress / dstrain);
		else
		  trialTangent = this->getInitialTangent();
		opserr << "trial tangent = " << trialTangent << endln;

		return (trialStress);
	}
	else if (trialStrain >= (db / 2) && trialStrain < db) {

		trialStress = stressMax + (trialStrain - (db / 2)) * ((0.0 - stressMax) / (db / 2));

		trialTangent = ((0.0 - stressMax) / (db / 2));

		return (trialStress);

	}
	else 

		return 0.0;
	

	//return 100 * trialStrain;
	

	
}


double
PlateBearingConnectionThermal::getTangent(void)
{
//Tangent at current state
	return trialTangent;
	//return this->getInitialTangent();
}


double
PlateBearingConnectionThermal::getInitialTangent(void)
{
  opserr << "Initial tangent = " << E0*db*t/100 << endln;
	return E0 * (db * t) / 100;
}

// Convergence reached, commit values 

int
PlateBearingConnectionThermal::commitState(void)
{
  //committedStrain = trialStrain;
  //committedStress = trialStress;

	return 0;
}


int
PlateBearingConnectionThermal::revertToLastCommit(void)
{
  trialStrain = committedStrain;
  trialStress = committedStress;
	return 0;
}


int
PlateBearingConnectionThermal::revertToStart(void)
{
	trialStrain = 0.0;
	trialStress = 0.0;
	committedStrain = 0.0;
	committedStress = 0.0;
	return 0;
}

// copy for every instance of material in model
UniaxialMaterial *
PlateBearingConnectionThermal::getCopy(void)
{
	PlateBearingConnectionThermal *theCopy = new PlateBearingConnectionThermal(this->getTag(), Fy0, Fu0, E0, db, e2, t);
	theCopy->trialStrain = trialStrain;
	theCopy->trialStress = trialStress;
	theCopy->committedStrain = committedStrain;
	theCopy->committedStress = committedStress;
	return theCopy;
}


int
PlateBearingConnectionThermal::sendSelf(int cTag, Channel &theChannel)
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
PlateBearingConnectionThermal::recvSelf(int cTag, Channel &theChannel,
	FEM_ObjectBroker &theBroker)
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
PlateBearingConnectionThermal::Print(OPS_Stream &s, int flag)
{
	s << "PlateBearingConnectionThermal tag: " << this->getTag() << endln;
	s << "  Fy: " << Fy0 << " Fu: " << Fu0 << " E: " << E0 << " e2: " << e2 << " db: " << db << " t: " << t << endln;
}


