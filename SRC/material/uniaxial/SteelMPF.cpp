// Code written/implemented by:	Kristijan Kolozvari (kkolozvari@fullerton.edu)
//								California State University, Fullerton 
//								Kutay Orakcal
//								Bogazici University, Istanbul, Turkey
//								John Wallace
//								University of California, Los Angeles
//
// Created: 07/2015
//
// Description: This file contains the class implementation for uniaxialMaterial 
// SteelMPF, which represents the well-known uniaxial constitutive nonlinear 
// hysteretic material model for steel proposed by Menegotto and Pinto (1973), 
// and extended by Filippou et al. (1983) to include isotropic strain hardening 
// effects.
//
// References:
// 1) Menegotto, M., and Pinto, P.E. (1973). Method of analysis of cyclically 
// loaded RC plane frames including changes in geometry and non-elastic behavior 
// of elements under normal force and bending. Preliminary Report IABSE, vol 13. 
// 2) Filippou F.C., Popov, E.P., and Bertero, V.V. (1983). "Effects of Bond 
// Deterioration on Hysteretic Behavior of Reinforced Concrete Joints". Report 
// EERC 83-19, Earthquake Engineering Research Center, University of California, Berkeley.
//
// Source: /usr/local/cvs/OpenSees/SRC/material/uniaxial/SteelMPF.cpp
//
// Rev: 1


#include <SteelMPF.h>
#include <Vector.h>
#include <Matrix.h>
#include <Channel.h>
#include <Information.h>
#include <math.h>
#include <float.h>

#include <elementAPI.h>
#define OPS_Export 

// Read input parameters and build the material
OPS_Export void *OPS_SteelMPF(void)
{
	// Pointer to a uniaxial material that will be returned                       
	UniaxialMaterial *theMaterial = 0;

	int numArgs = OPS_GetNumRemainingInputArgs();

	// Parse the script for material parameters
	if (numArgs != 9 && numArgs != 11) {
		opserr << "Incorrect # args, Want: uniaxialMaterial SteelMPF tag? sigyieldp? sigyieldn? E0? bp? bn? R0? a1? a2? <a3? a4?>";
		return 0;
	}

	int iData[1];
	double dData[10];
	dData[8]=0.01; // set a3 to constructor default
	dData[9]=7.0;  // set a4 to constructor default

	int numData = 1;
	if (OPS_GetIntInput(&numData, iData) != 0) {
		opserr << "WARNING invalid uniaxialMaterial FSAM tag" << endln;
		return 0;
	}

	numData = numArgs-1;
	if (OPS_GetDoubleInput(&numData, dData) != 0) {
		opserr << "Invalid data for uniaxial Hysteretic " << dData[0] << endln;
		return 0;
	}

	theMaterial = new SteelMPF(iData[0], dData[0], dData[1], dData[2], dData[3], dData[4], 
		dData[5], dData[6], dData[7], dData[8], dData[9]);

	return theMaterial;
}

// Default Constructor
SteelMPF::SteelMPF
	(int tag, double FYp, double FYn, double E, double Bp, double Bn, double Rr, double A1, double A2):
UniaxialMaterial(tag,MAT_TAG_SteelMPF),
	sigyieldp(FYp), sigyieldn(FYn), E0(E), bp(Bp), bn(Bn), R0(Rr), a1(A1), a2(A2), 
	a3(0.01), a4(7.0) // default values for strain hardening

{
	// Sets all history and state variables to initial values
	// TRIAL History variables
	sigr = 0.0;
	sig0 = 0.0;
	er = 0.0;
	e0 = 0.0;
	emax = 0.0;
	inc = 0;
	nloop = 0;

	// CONVERGED History variables
	sigrold = 0.0;
	sig0old = 0.0;
	erold = 0.0;
	e0old = 0.0;
	emaxold = 0.0;
	incold = 0;
	nloopold = 0;

	// TRIAL State variables
	def = 0.0;
	F = 0.0;
	stif = E0;
	R = R0;

	// CONVERGED State variables
	defold = 0.0;
	Fold = 0.0;
	stifold = E0; 
	Rold = R0; 

	// Yield Strain
	eyieldp = sigyieldp/E0; 
	eyieldn = sigyieldn/E0; 

}

// Constructor with optional strain hardening parameters a3 and a4
SteelMPF::SteelMPF
	(int tag, double FYp, double FYn, double E, double Bp, double Bn, double Rr, double A1, double A2, double A3, double A4) :
UniaxialMaterial(tag, MAT_TAG_SteelMPF),
	sigyieldp(FYp), sigyieldn(FYn), E0(E), bp(Bp), bn(Bn), R0(Rr), a1(A1), a2(A2), a3(A3), a4(A4)

{
	// Sets all history and state variables to initial values
	// TRIAL History variables
	sigr = 0.0;
	sig0 = 0.0;
	er = 0.0;
	e0 = 0.0;
	emax = 0.0;
	inc = 0;
	nloop = 0;

	// CONVERGED History variables
	sigrold = 0.0;
	sig0old = 0.0;
	erold = 0.0;
	e0old = 0.0;
	emaxold = 0.0;
	incold = 0;
	nloopold = 0;

	// TRIAL State variables
	def = 0.0;
	F = 0.0;
	stif = E0;
	R = R0;

	// CONVERGED State variables
	defold = 0.0;
	Fold = 0.0;
	stifold = E0;
	Rold = R0;

	// Yield Strain
	eyieldp = sigyieldp / E0;
	eyieldn = sigyieldn / E0;

}

SteelMPF::SteelMPF() :UniaxialMaterial(0, MAT_TAG_SteelMPF),
	sigyieldp(0.0), sigyieldn(0.0), E0(0.0), bp(0.0), bn(0.0), R0(0.0), a1(0.0), a2(0.0),  a3(0.0), a4(0.0)
{

}

// Destructor
SteelMPF::~SteelMPF ()
{

}

int SteelMPF::setTrialStrain(double strain, double strainRate)
{
	// Reset history variables to last converged state
	sigr = sigrold;
	sig0 = sig0old;
	er = erold;
	e0 = e0old;
	emax = emaxold;
	inc = incold;
	nloop = nloopold;

	// "dummy" step 
	if (strain == 0.0 && strainRate == 0.0)
	{
		return 0;
	}

	// Set trial strain
	def = strain;

	// Calculate the trial state given the trial strain
	determineTrialState(def);

	return 0;
}


int SteelMPF::setTrial (double strain, double &stress, double &tangent, double strainRate)
{
	// Reset history variables to last converged state
	sigr = sigrold;
	sig0 = sig0old;
	er = erold;
	e0 = e0old;
	emax = emaxold;
	inc = incold;
	nloop = nloopold;

	// "dummy" step 
	if (strain == 0.0 && strainRate == 0.0)
	{
		return 0;
	}

	// Set trial strain
	def = strain;

	// Calculate the trial state given the trial strain
	this->determineTrialState(def);

	stress = F;
	tangent = stif;

	return 0;
}

// Calculate trial stress and stiffness
void SteelMPF::determineTrialState(double dStrain)
{
	double e = def;
	double eold = defold;
	double sigold = Fold;
	double b;

	// Determine initial loading condition
	if (incold == 0) {

		if (e < 0.0) {
			inc = -1;
		} else {
			inc = 1;
		}

		er = 0;
		sigr = 0;

		if (inc == 1) {

			e0 = eyieldp; 
			sig0 = sigyieldp;
			b = bp; 

		} else {

			e0 = -eyieldn;
			sig0 = -sigyieldn;
			b = bn;
		}

		nloop = 1;
		R = R0;
		emax = eyieldn;

	}

	else {

		// Cyclic Stress - Strain model for Steel derived from Menegotto - Pinto Equation

		if (def > defold) {
			inc = 1;
		} else if (def < defold) {
			inc = -1;
		} else {
			inc = incold;
		}

		if (incold == 1) {

			if (e < eold) {

				er = eold;
				sigr = sigold;
				nloop = nloop + 1;

				double sigy = -sigyieldn;
				b = bn;

				if (fabs(er) > emaxold) {
					emax = fabs(er);
				} else {
					emax = emaxold;
				}

				double sigstar = -sigy*a3*(emax / eyieldn - a4);

				if (sigstar < 0) {
					sigstar = 0;
				}

				sigy = sigy - sigstar;

				e0 = (sigy*(1.0 - b) + E0*er - sigr) / (E0*(1.0 - b));
				sig0 = sigr + E0*(e0 - er);

				double em;

				if (nloop == 1) {
					em = e0;
				} else if (nloop == 2) {
					em = -e0old;
				} else {
					em = erold;
				}

				double ksi = fabs((em - e0) / eyieldp);
				R = R0 - a1*ksi / (a2 + ksi);

			} else {
				b = bp;
				er = erold;
				e0 = e0old;
				sigr = sigrold;
				sig0 = sig0old;
				R = Rold;
				emax = emaxold;
			}
		}

		if (incold == -1) {
			if (e > eold) {

				er = eold;
				sigr = sigold;
				nloop = nloop + 1;
				double sigy = sigyieldp;
				b = bp;

				if (fabs(er) > emaxold) {
					emax = fabs(er);
				} else {
					emax = emaxold;
				}

				double sigstar = sigy*a3*(emax / eyieldn - a4);

				if (sigstar < 0) {
					sigstar = 0;
				}

				sigstar = 0;

				sigy = sigy + sigstar;

				e0 = (sigy*(1.0 - b) + E0*er - sigr) / (E0*(1.0 - b));
				sig0 = sigr + E0*(e0 - er);

				double em;

				if (nloop == 1) {
					em = e0;
				} else if (nloop == 2) {
					em = -e0old;
				} else {
					em = erold;
				}

				double ksi = fabs((em - e0) / eyieldn);
				R = R0 - a1*ksi / (a2 + ksi);
			}

			else {

				b = bn;
				er = erold;
				e0 = e0old;
				sigr = sigrold;
				sig0 = sig0old;
				R = Rold;
				emax = emaxold;
			}
		}

	}

	double est = (e - er) / (e0 - er);
	double sigst = b*est + ((1.0 - b) / pow((1.0 + pow(est, R)), (1.0 / R))) * est;
	double sig = sigr + sigst*(sig0 - sigr);
	double Et = ((sig0 - sigr) / (e0 - er)) * (b + ((1.0 - b) / pow((1.0 + pow(est, R)), (1.0 / R))) * (1.0 - pow(est, R) / (1.0 + pow(est, R))));

	F = sig;
	stif = Et;

}

double SteelMPF::getStrain()
{
	return def;
}

double SteelMPF::getStress()
{
	return F;
}

double SteelMPF::getTangent()
{
	return stif;
}

int SteelMPF::commitState()
{
	// History variables
	sigrold = sigr;
	sig0old = sig0;
	erold = er;
	e0old = e0;
	emaxold = emax;
	incold = inc;
	nloopold = nloop;

	// State variables
	defold = def;
	Fold = F;
	stifold = stif;
	Rold = R;

	return 0;
}

int SteelMPF::revertToLastCommit()
{
	// Reset trial history variables to last committed state
	sigr = sigrold;
	sig0 = sig0old;
	er = erold;
	e0 = e0old;
	emax = emaxold;
	inc = incold;
	nloop = nloopold;

	// Reset trial state variables to last committed state
	def = defold;
	F = Fold;
	stif = stifold;
	R = Rold;

	return 0;
}

int SteelMPF::revertToStart()
{
	// TRIAL History variables
	sigr = 0.0;
	sig0 = 0.0;
	er = 0.0;
	e0 = 0.0;
	emax = 0.0;
	inc = 0;
	nloop = 0;

	// CONVERGED History variables
	sigrold = 0.0;
	sig0old = 0.0;
	erold = 0.0;
	e0old = 0.0;
	emaxold = 0.0;
	incold = 0; 
	nloopold = 0;

	// TRIAL State variables
	def = 0.0;
	F = 0.0;
	stif = E0;
	R = R0;

	// CONVERGED State variables
	defold = 0.0;
	Fold = 0.0;
	stifold = E0; 
	Rold = R0;

	eyieldp = sigyieldp / E0; 
	eyieldn = sigyieldn / E0;

	return 0;
}

UniaxialMaterial* SteelMPF::getCopy()
{
	SteelMPF* theCopy = new SteelMPF(this->getTag(), sigyieldp, sigyieldn, E0, bp, bn, R0, a1, a2, a3, a4); 

	// CONVERGED History variables
	theCopy->sigrold = sigrold;
	theCopy->sig0old = sig0old;
	theCopy->erold = erold;
	theCopy->e0old = e0old;
	theCopy->emaxold = emaxold;
	theCopy->incold = incold;
	theCopy->nloopold = nloopold;

	// TRIAL History variables
	theCopy->sigr = sigr;
	theCopy->sig0 = sig0;
	theCopy->er = er;
	theCopy->e0 = e0;
	theCopy->emax = emax;
	theCopy->inc = inc;
	theCopy->nloop = nloop;

	// CONVERGED State variables
	theCopy->defold = defold;
	theCopy->Fold = Fold;
	theCopy->stifold = stifold;
	theCopy->Rold = Rold;

	// TRIAL State variables
	theCopy->def = def;
	theCopy->F = F;
	theCopy->stif = stif;
	theCopy->R = R;

	return theCopy;
}

int SteelMPF::sendSelf (int commitTag, Channel& theChannel)
{
	int res = 0;
	static Vector data(22);

	data(0) = this->getTag();

	// Material properties
	data(1) = sigyieldp;
	data(2) = sigyieldn;
	data(3) = E0;
	data(4) = bp;
	data(5) = bn;
	data(6) = R0;
	data(7) = a1;
	data(8) = a2;
	data(9) = a3;
	data(10) = a4;

	// CONVERGED History variables
	data(11) = sigrold;
	data(12) = sig0old;
	data(13) = erold;
	data(14) = e0old;
	data(15) = emaxold;
	data(16) = incold;
	data(17) = nloopold;

	// CONVERGED State variables
	data(18) = defold;
	data(19) = Fold;
	data(20) = stifold;
	data(21) = Rold;

	// Data is only sent after convergence, so no trial variables
	// need to be sent through data vector

	res = theChannel.sendVector(this->getDbTag(), commitTag, data);
	if (res < 0) 
		opserr << "SteelMPF::sendSelf() - failed to send data\n";

	return res;
}

int SteelMPF::recvSelf (int commitTag, Channel& theChannel, FEM_ObjectBroker& theBroker)
{
	int res = 0;
	static Vector data(22);
	res = theChannel.recvVector(this->getDbTag(), commitTag, data);

	if (res < 0) {
		opserr << "SteelMPF::recvSelf() - failed to receive data\n";
		this->setTag(0);      
	}
	else {
		this->setTag(int(data(0)));

		// Material properties
		sigyieldp = data(1);
		sigyieldn = data(2);
		E0 = data(3);
		bp = data(4);
		bn = data(5);
		R0 = data(6);
		a1 = data(7);
		a2 = data(8);
		a3 = data(9);
		a4 = data(10);

		// CONVERGED History variables
		sigrold = data(11);
		sig0old = data(12);
		erold = data(13);
		e0old = data(14);
		emaxold = data(15);
		incold = data(16);
		nloopold = data(17);

		// Copy converged history values into trial values since data is only
		// sent (received) after convergence
		sigr = sigrold;
		sig0 = sig0old;
		er = erold;
		e0 = e0old;
		emax = emaxold;
		inc = incold;
		nloop = nloopold;

		// CONVERGED State variables
		defold = data(18);
		Fold = data(19);
		stifold = data(20);
		Rold = data(21);

		// Copy converged state values into trial values
		def = defold;
		F = Fold;
		stif = stifold;
		R = Rold;

	}

	return res; 
}

void SteelMPF::Print (OPS_Stream& s, int flag)
{
	s << "SteelMPF tag: " << this->getTag() << endln;
	s << " fyp = " << sigyieldp << " ";
	s << " fyn = " << sigyieldn << " ";
	s << "  E0 = " << E0 << " ";
	s << "  bp = " << bp << " ";
	s << "  bn = " << bn << " ";
	s << "   R = " << R << ", a1 = " << a1 << ", a2 = " << a2 << ", a3 = " << a3 << ", a4 = " << a4 << endln;

}
