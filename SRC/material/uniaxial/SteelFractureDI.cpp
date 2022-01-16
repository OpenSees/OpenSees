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

// $Revision: 1.0 $
// $Date: 2021-05-18 00:17:20 $
// $Source: /usr/local/cvs/OpenSees/SRC/material/uniaxial/SteelFractureDI.cpp,v $

// Written: Francisco A. Galvis and Wen-Yi Yen
// Created: 05/2021
//
// Description: This file contains the class implementation of SteelFractureDI. 
// SteelFractureDI is based on the source code for Steel02
//

#include <math.h>

#include <stdlib.h>
#include "SteelFractureDI.h"
#include <OPS_Globals.h>
#include <MaterialResponse.h> // for recorder
#include <float.h>
#include <Channel.h>


#include <elementAPI.h>
#include <OPS_Globals.h>


void * OPS_ADD_RUNTIME_VPV(OPS_SteelFractureDI)
{
	// Pointer to a uniaxial material that will be returned
	UniaxialMaterial *theMaterial = 0;

	int    iData[1];
	double dData[15];
	int numData = 1;

	if (OPS_GetIntInput(&numData, iData) != 0) {
		opserr << "WARNING invalid uniaxialMaterial SteelFractureDI tag" << endln;
		return 0;
	}

	numData = OPS_GetNumRemainingInputArgs();

	if (numData != 15) {
		opserr << "Invalid #args, want: uniaxialMaterial SteelFractureDI " << iData[0] <<
			" Fy? Fyc? E? b? R0? cR1? cR2? a1? a2? a3? a4? sigcr? m? sigmin? FI_lim?" << endln;
		return 0;
	}
	else {
		if (OPS_GetDoubleInput(&numData, dData) != 0) {
			opserr << "Invalid arggs: uniaxialMaterial SteelFractureDI " << iData[0] <<
				" Fy? FyC? E? b? R0? cR1? cR2? a1? a2? a3? a4? sigcr? m? sigmin? FI_lim?" << endln;
			return 0;
		}
		theMaterial = new SteelFractureDI(iData[0], dData[0], dData[1], dData[2],
			dData[3], dData[4], dData[5], dData[6],
			dData[7], dData[8], dData[9], dData[10], dData[11], dData[12], dData[13], dData[14]);
	}



	if (theMaterial == 0) {
		opserr << "WARNING could not create uniaxialMaterial of type SteelFractureDI Material\n";
		return 0;
	}

	return theMaterial;
}



SteelFractureDI::SteelFractureDI(int tag,
	double _Fy, double _FyC, double _E0, double _b,
	double _R0, double _cR1, double _cR2,
	double _a1, double _a2, double _a3, double _a4, double _sigcr, double _m, double _sigmin, double _FI_lim):
	UniaxialMaterial(tag, MAT_TAG_SteelFractureDI),
	//UniaxialMaterial(tag, 0),
	Fy(_Fy), FyC(_FyC), E0(_E0), b(_b), R0(_R0), cR1(_cR1), cR2(_cR2), a1(_a1), a2(_a2), a3(_a3), a4(_a4),
	sigcr(_sigcr), m(_m), sigmin(_sigmin), FI_lim(_FI_lim)
{
	konP = 0;
	kon = 0;
	eP = E0;
	epsP = 0.0;
	sigP = 0.0;
	sig = 0.0;
	eps = 0.0;
	e = E0;
	
	epsmaxP = Fy / E0;
	epsminP = -FyC / E0; ////////////////////////////////////////////
	epsplP = 0.0;
	epss0P = 0.0;
	sigs0P = 0.0;
	epssrP = 0.0;
	sigsrP = 0.0;

	// ************** added for fracture ***************
	epsCont = 0.0;
	konf = 0;
	konC = 0;

	epsContP = 0.0;
	eps_0 = 0.0;
	eps_1 = 0.0;
	eps_r = 0;
	konfP = 0;
	konCP = 0;
	// ************** added for fracture ***************
	// ************** added for DI ********************
	DIP = 0.0;
	isStartP = 1;
	sigPDIP = 0.0;
	slopePP = 0.0;
	sumTenPP = 0.0;
	sumCompPP = 0.0;

	DI = 0.0;
	isStart = 1;
	sigPDI = 0.0;
	slopeP = 0.0;
	sumTenP = 0.0;
	sumCompP = 0.0;

	// ************** added for DI ********************

}


SteelFractureDI::SteelFractureDI(void) :
	UniaxialMaterial(0, MAT_TAG_SteelFractureDI)
	// UniaxialMaterial(0, 0)
{
	konP = 0;
}

SteelFractureDI::~SteelFractureDI(void)
{
	// Does nothing
}

UniaxialMaterial*
SteelFractureDI::getCopy(void)
{
	//Steel02 *theCopy = new Steel02(this->getTag(), Fy, FyC, E0, b, R0, cR1, cR2, a1, a2, a3, a4, sigini);
	SteelFractureDI *theCopy = new SteelFractureDI(this->getTag(), Fy, FyC, E0, b, R0, cR1, cR2, a1, a2, a3, a4, sigcr, m, sigmin, FI_lim);

	return theCopy;
}

double
SteelFractureDI::getInitialTangent(void)
{
	return E0;
}

int
SteelFractureDI::setTrialStrain(double trialStrain, double strainRate)
{
	double Esh = b * E0;
	double epsy = Fy / E0;
	double epsyc = FyC / E0; ////////////////////////////////////////////
	
	eps = trialStrain;
	double deps = eps - epsP;

	epsmax = epsmaxP;
	epsmin = epsminP;
	epspl = epsplP;
	epss0 = epss0P;
	sigs0 = sigs0P;
	epsr = epssrP;
	sigr = sigsrP;
	kon = konP;
	// ************** added for fracture ***************
	epsCont = epsContP;
	eps_0 = eps_0P;
	eps_1 = eps_1P;
	eps_r = eps_rP;
	konf = konfP;
	konC = konCP;
	// ************** added for fracture ***************
	// ************** added for DI ********************
	DI = DIP;
	isStart = isStartP;
	sigPDI = sigPDIP;
	slopeP = slopePP;
	sumTenP = sumTenPP;
	sumCompP = sumCompPP;

	// ************** added for DI ********************

	if (kon == 0 || kon == 3) { // modified C-P. Lamarche 2006


		if (fabs(deps) < 10.0*DBL_EPSILON) {

			e = E0;
			//sig = sigini;                // modified C-P. Lamarche 2006
			sig = 0;
			kon = 3;                     // modified C-P. Lamarche 2006 flag to impose initial stess/strain
			return 0;

		}
		else {

			epsmax = epsy;
			epsmin = -epsyc; ////////////////////////////////////////////
			if (deps < 0.0) {
				kon = 2;
				epss0 = epsmin;
				sigs0 = -FyC; ////////////////////////////////////////////
				epspl = epsmin;
			}
			else {
				kon = 1;
				epss0 = epsmax;
				sigs0 = Fy;
				epspl = epsmax;
			}
		}
	}

	// in case of load reversal from negative to positive strain increment, 
	// update the minimum previous strain, store the last load reversal 
	// point and calculate the stress and strain (sigs0 and epss0) at the 
	// new intersection between elastic and strain hardening asymptote 
	// To include isotropic strain hardening shift the strain hardening 
	// asymptote by sigsft before calculating the intersection point 
	// Constants a3 and a4 control this stress shift on the tension side 

	if (kon == 2 && deps > 0.0) {


		kon = 1;
		epsr = epsP;
		sigr = sigP;
		//epsmin = min(epsP, epsmin);
		if (epsP < epsmin)
			epsmin = epsP;
		double d1 = (epsmax - epsmin) / (2.0*(a4 * epsy));
		double shft = 1.0 + a3 * pow(d1, 0.8);
		epss0 = (Fy * shft - Esh * epsy * shft - sigr + E0 * epsr) / (E0 - Esh); 
		sigs0 = Fy * shft + Esh * (epss0 - epsy * shft);
		epspl = epsmax;

	}
	else if (kon == 1 && deps < 0.0) {

		// update the maximum previous strain, store the last load reversal 
		// point and calculate the stress and strain (sigs0 and epss0) at the 
		// new intersection between elastic and strain hardening asymptote 
		// To include isotropic strain hardening shift the strain hardening 
		// asymptote by sigsft before calculating the intersection point 
		// Constants a1 and a2 control this stress shift on compression side 

		kon = 2;
		epsr = epsP;
		sigr = sigP;
		//      epsmax = max(epsP, epsmax);
		if (epsP > epsmax)
			epsmax = epsP;

		double d1 = (epsmax - epsmin) / (2.0*(a2 * epsyc));
		double shft = 1.0 + a1 * pow(d1, 0.8);
		epss0 = (-FyC * shft + Esh * epsyc * shft - sigr + E0 * epsr) / (E0 - Esh);
		sigs0 = -FyC * shft + Esh * (epss0 + epsyc * shft);
		epspl = epsmin;
	}


	// calculate current stress sig and tangent modulus E 
	if (kon != 4) { // non-fracture case
		// tension
		double xi = fabs((epspl - epss0) / epsy);
		if (deps < 0.0)
			// compresion
			xi = fabs((epspl - epss0) / epsyc);

		double R = R0 * (1.0 - (cR1*xi) / (cR2 + xi));
		double epsrat = (eps - epsr) / (epss0 - epsr);
		double dum1 = 1.0 + pow(fabs(epsrat), R);
		double dum2 = pow(dum1, (1 / R));

		sig = b * epsrat + (1.0 - b)*epsrat / dum2;
		sig = sig * (sigs0 - sigr) + sigr;

		e = b + (1.0 - b) / (dum1*dum2);
		e = e * (sigs0 - sigr) / (epss0 - epsr);

		// calculate DI
		calcDI(sigcr, m, sigmin, FI_lim, isStart, sig, sigPDI, DI, slopeP, sumTenP, sumCompP);

		// check if fractured
		if (DI >= FI_lim) {
			kon = 4;
			konf = 1;

			// Point at compression yield envelope
			eps_1 = (-FyC + Esh * epsyc - sigP + E0 * epsP) / (E0 - Esh);
			sig_1 = -FyC + Esh * (eps_1 + epsyc);

			// Points at horizontal axis
			eps_0 = epsP - sigP / E0;
			epsCont = 2 * eps_0 - eps_1;
			eps_r = epsCont;
			epsr = eps_r;
			sigr = 0;

			konC = 1;

			// After fracture goes to zero stress
			sig = 0.0;
			e = 0.0;
		}
		//else {
		//	DI = DI_cache;
		//}

	}
	else if (kon == 4) { // fractured case
		if (eps >= epsCont) { // not contacted
			sig = 0.0;
			e = 0.0;
			if (deps > 0) {
				konf = 2;
			}
			else {
				konf = 1;
			}
			konC = 0;
		}
		else if (eps < epsCont) { // contacted
			if (!konC) { // at first contact
				konC = 1;
				konf = 2;
			}
			if (konf == 2 && deps > 0) { // compression --> tension so update variables
				konf = 1; // konf = 1 means tensile now

				// Only update if it is out of the linear branch
				if (sig < 0.7*sig_1) {
					// New properties for unloading and reloading
					eps_0 = epsP - sigP / E0; // new interception is at reversal point minus elastic unloading
					sigs0 = 0; // new interception always at horizontal axis since flange fractured

					eps_1 = (-FyC + Esh * epsyc + E0 * eps_0) / (E0 - Esh); // Point at compression yield envelope
					sig_1 = -FyC + Esh * (eps_1 + epsyc);

					eps_r = 2 * eps_0 - eps_1; // Point at zero stress and e
				}
			}
			else if (konf == 1 && deps < 0) { // tension --> compression so update variables
				konf = 2; // konf == 2 means compression now				
			}

			// compute sig and e
			if (konf == 1) { // compute if is decompressing	
				// Functional form of first half of the compression load curve all the way until complete unload
				double R = 14;
				double epsrat = (eps - eps_1) / (eps_1 - eps_0);
				double dum1 = 1.0 + pow(fabs(epsrat), R);
				double dum2 = pow(dum1, (1 / R));

				// Functional form for b = 0 since tend to zero stress
				sig = sig_1 * (epsrat / dum2 + 1);
				e = sig_1 / (eps_1 - eps_0) * (1 / (dum1*dum2));

				if (eps > eps_r) {
					sig = 0;
					e = 0;
				}

			}
			else { // compute if is compressing
				double R = 14;
				if (eps >= (eps_0 + eps_1) / 2) {
					// Update strain for zero stress and e
					// Functional form of first half of the compression load curve (b = 0 since tends to horizontal axis)
					double epsrat = (eps - eps_1) / (eps_1 - eps_0);
					double dum1 = 1.0 + pow(fabs(epsrat), R);
					double dum2 = pow(dum1, (1 / R));

					sig = sig_1 * (epsrat / dum2 + 1);
					e = sig_1 / (eps_1 - eps_0) * (1 / (dum1*dum2));

					if (eps > eps_r) {
						sig = 0;
						e = 0;
					}
				}
				else {
					double epsrat = (eps - eps_0) / (eps_1 - eps_0);
					double dum1 = 1.0 + pow(fabs(epsrat), R);
					double dum2 = pow(dum1, (1 / R));

					sig = b * epsrat + (1.0 - b)*epsrat / dum2;
					sig = sig * sig_1;

					e = b + (1.0 - b) / (dum1*dum2);
					e = e * sig_1 / (eps_1 - eps_0);
				}
			}
		}
	}

	return 0;
}

void
SteelFractureDI::calcDI(double sigcr, double m, double sigmin, double FI_lim, int& isStart, double sig, double& sigPDI, double& DI, double& slopeP, double& sumTenP, double& sumCompP)
{
	// initialize variables ****
	double slope;
	double currSign;
	double sumComp;
	double sumTen;

	// Already fractured
	if (DI > FI_lim) {
		return;
	}

	// if is starting point
	if (isStart) {
		isStart = 0;
		sigPDI = sig;
		return;
	}

	// determine slope sign
	slope = sig - sigPDI;
	if (slope == 0) {
		currSign = returnSign(slopeP);
	}
	else {
		currSign = returnSign(slope);
	}

	// Accumulate compressive and tensile stresses
	if (fabs(sig) > sigmin) {
		// Accumulate only if stress exceeds a threshold

		if (currSign == 1 && sig > sigmin) {
			// Tensile excursion (only accumulates in actual tension)
			sumComp = sumCompP;
			sumTen = sumTenP + fabs(slope);
		}
		else {
			// Compressive excursion
			sumTen = sumTenP;
			if (sumCompP + fabs(slope) < sumTen) {
				// Only considers compression when there is damage to heal
				sumComp = sumCompP + fabs(slope);
			}
			else {
				sumComp = sumCompP;
			}
		}
	}
	else {
		// If below threshold keep same cumulative stresses
		sumComp = sumCompP;
		sumTen = sumTenP;
	}

	DI = (sumTen - sumComp * m) / sigcr;

	if (DI < 0) {
		DI = 0;
	}

	// update variables
	sigPDI = sig;
	slopeP = slope;
	sumCompP = sumComp;
	sumTenP = sumTen;
}


int
SteelFractureDI::returnSign(double v) {
	if (v < 0) return -1;
	if (v > 0) return 1;
	return 0;
}

double
SteelFractureDI::getStrain(void)
{
	return eps;
}

double
SteelFractureDI::getStress(void)
{
	return sig;
}

double
SteelFractureDI::getTangent(void)
{
	return e;
}

double
SteelFractureDI::getDI(void)
{
	return DI;
}

int
SteelFractureDI::commitState(void)
{
	epsminP = epsmin;
	epsmaxP = epsmax;
	epsplP = epspl;
	epss0P = epss0;
	sigs0P = sigs0;
	epssrP = epsr;
	sigsrP = sigr;
	konP = kon;

	eP = e;
	sigP = sig;
	epsP = eps;

	// ************** added for fracture ***************
	epsContP = epsCont;
	eps_0P = eps_0;
	eps_1P = eps_1;
	eps_rP = eps_r;
	konfP = konf;
	konCP = konC;
	// ************** added for fracture ***************
	// ************** added for DI ********************
	DIP = DI;
	isStartP = isStart;
	sigPDIP = sigPDI;
	slopePP = slopeP;
	sumTenPP = sumTenP;
	sumCompPP = sumCompP;
	// ************** added for DI ********************
	return 0;
}

int
SteelFractureDI::revertToLastCommit(void)
{
	epsmin = epsminP;
	epsmax = epsmaxP;
	epspl = epsplP;
	epss0 = epss0P;
	sigs0 = sigs0P;
	epsr = epssrP;
	sigr = sigsrP;
	kon = konP;

	e = eP;
	sig = sigP;
	eps = epsP;

	// ************** added for fracture ***************
	epsCont = epsContP;
	eps_0 = eps_0P;
	eps_1 = eps_0P;
	eps_r = eps_rP;
	konf = konfP;
	konC = konCP;
	// ************** added for fracture ***************
	// ************** added for DI ********************
	DI = DIP;
	isStart = isStartP;
	sigPDI = sigPDIP;
	slopeP = slopePP;
	sumTenP = sumTenPP;
	sumCompP = sumCompPP;
	// ************** added for DI ********************

	return 0;
}

int
SteelFractureDI::revertToStart(void)
{
	konP = 0;
	kon = 0;
	eP = E0;
	epsP = 0.0;
	sigP = 0.0;
	sig = 0.0;
	eps = 0.0;
	e = E0;

	epsmaxP = Fy / E0;
	epsminP = -epsmaxP;
	epsplP = 0.0;
	epss0P = 0.0;
	sigs0P = 0.0;
	epssrP = 0.0;
	sigsrP = 0.0;

	// ************** added for fracture ***************
	epsCont = 0.0;
	eps_0 = 0.0;
	eps_1 = 0.0;
	eps_r = 0.0;
	konf = 0.0;
	konC = 0.0;

	epsContP = 0.0;
	eps_0P = 0.0;
	eps_1P = 0.0;
	eps_rP = 0.0;
	konfP = 0.0;
	konCP = 0.0;
	// ************** added for fracture ***************
	// ************** added for DI ********************
	DIP = 0.0;
	isStartP = 1;
	sigPDIP = 0.0;
	slopePP = 0.0;
	sumTenPP = 0.0;
	sumCompPP = 0.0;

	DI = 0.0;
	isStart = 1;
	sigPDI = 0.0;
	slopeP = 0.0;
	sumTenP = 0.0;
	sumCompP = 0.0;
	// ************** added for DI ********************

	return 0;
}

int
SteelFractureDI::sendSelf(int commitTag, Channel &theChannel)
{
	static Vector data(39);

	data(0) = this->getTag();

	// Material properties
	data(1) = Fy;
	data(2) = FyC;
	data(3) = E0;
	data(4) = b;
	data(5) = R0;
	data(6) = cR1;
	data(7) = cR2;
	data(8) = a1;
	data(9) = a2;
	data(10) = a3;
	data(11) = a4;
	data(12) = sigcr;
	data(13) = m;
	data(14) = sigmin;
	data(15) = FI_lim;

	// History variables
	data(16) = konP;
	data(17) = eP;
	data(18) = epsP;
	data(19) = sigP;
	data(20) = epsmaxP;
	data(21) = epsminP;
	data(22) = epsplP;
	data(23) = epss0P;
	data(24) = sigs0P;
	data(25) = epssrP;
	data(26) = sigsrP;
	// ************** added for fracture ***************
	data(27) = epsContP;
	data(28) = eps_0P;
	data(29) = eps_1P;
	data(30) = eps_rP;
	data(31) = konfP;
	data(32) = konCP;
	// ************** added for fracture ***************
	// ************** added for DI ********************
	data(33) = DIP;
	data(34) = isStartP;
	data(35) = sigPDIP;
	data(36) = slopePP;
	data(37) = sumTenPP;
	data(38) = sumCompPP;
	// ************** added for DI ********************	

	if (theChannel.sendVector(this->getDbTag(), commitTag, data) < 0) {
		opserr << "SteelFractureDI::sendSelf() - failed to sendSelf\n";
		return -1;
	}
	// opserr << "sendSelf SteelFractureDI completed tag data(0)\n"; 
	return 0;
}

int
SteelFractureDI::recvSelf(int commitTag, Channel &theChannel,
	FEM_ObjectBroker &theBroker)
{
	static Vector data(39);

	// opserr << "recvSelf SteelFractureDI started\n";

	if (theChannel.recvVector(this->getDbTag(), commitTag, data) < 0) {
		opserr << "SteelFractureDI::recvSelf() - failed to recvSelf\n";
		return -1;
	}

	this->setTag(int(data(0)));

	// Material properties
	Fy = data(1);
	FyC = data(2);
	E0 = data(3);
	b = data(4);
	R0 = data(5);
	cR1 = data(6);
	cR2 = data(7);
	a1 = data(8);
	a2 = data(9);
	a3 = data(10);
	a4 = data(11);
	sigcr = data(12);
	m = data(13);
	sigmin = data(14);
	FI_lim = data(15);

	// History variables
	konP = data(16);
	eP = data(17);
	epsP = data(18);
	sigP = data(19);
	epsmaxP = data(20);
	epsminP = data(21);
	epsplP = data(22);
	epss0P = data(23);
	sigs0P = data(24);
	epssrP = data(25);
	sigsrP = data(26);
	// ************** added for fracture ***************
	epsContP = data(27);
	eps_0P = data(28);
	eps_1P = data(29);
	eps_rP = data(30);
	konfP = data(31);
	konCP = data(32);
	// ************** added for fracture ***************
	// ************** added for DI ********************
	DIP = data(33);
	isStartP = data(34);
	sigPDIP = data(35);
	slopePP = data(36);
	sumTenPP = data(37);
	sumCompPP = data(38);
	// ************** added for DI ********************


	// Transient variables
	kon = konP;
	sig = sigP;
	eps = epsP;
	e = E0;

	epsCont = epsContP;
	eps_0 = eps_0P;
	eps_1 = eps_1P;
	eps_r = eps_rP;
	konf = konfP;
	konC = konCP;

	DI = DIP;
	isStart = isStartP;
	sigPDI = sigPDIP;
	slopeP = slopePP;
	sumTenP = sumTenPP;
	sumCompP = sumCompPP;

	// opserr << "--recvSelf SteelFractureDI completed\n";

	return 0;
	
}

void
SteelFractureDI::Print(OPS_Stream &s, int flag)
{
	if (flag == OPS_PRINT_PRINTMODEL_MATERIAL) {
		//    s << "Steel02:(strain, stress, tangent) " << eps << " " << sig << " " << e << endln;
		s << "SteelFracture tag: " << this->getTag() << endln;
		s << "  Fy: " << Fy << ", ";
		s << "  FyC: " << FyC << ", ";
		s << "  E0: " << E0 << ", ";
		s << "   b: " << b << ", ";
		s << "  R0: " << R0 << ", ";
		s << " cR1: " << cR1 << ", ";
		s << " cR2: " << cR2 << ", ";
		s << "  a1: " << a1 << ", ";
		s << "  a2: " << a2 << ", ";
		s << "  a3: " << a3 << ", ";
		s << "  a4: " << a4 << ", ";
		s << "  sigcr: " << sigcr << ", ";
		s << "  m: " << m << ", ";
		s << "  sigmin: " << sigmin << ", ";
		s << "  FI_lim: " << FI_lim << ", ";
	}

	if (flag == OPS_PRINT_PRINTMODEL_JSON) {
		s << "\t\t\t{";
		s << "\"name\": \"" << this->getTag() << "\", ";
		s << "\"type\": \"SteelFracture\", ";
		s << "\"E\": " << E0 << ", ";
		s << "\"Fy\": " << Fy << ", ";
		s << "\"FyC\": " << FyC << ", ";
		s << "\"b\": " << b << ", ";
		s << "\"R0\": " << R0 << ", ";
		s << "\"cR1\": " << cR1 << ", ";
		s << "\"cR2\": " << cR2 << ", ";
		s << "\"a1\": " << a1 << ", ";
		s << "\"a2\": " << a2 << ", ";
		s << "\"a3\": " << a3 << ", ";
		s << "\"a4\": " << a4 << ", ";
		s << "  sigcr: " << sigcr << ", ";
		s << "  m: " << m << ", ";
		s << "  sigmin: " << sigmin << ", ";
		s << "  FI_lim: " << FI_lim << ", ";
	}
}

Response*
SteelFractureDI::setResponse(const char **argv, int argc, OPS_Stream &theOutput)
{
	if (argc == 0)
		return 0;

	Response *theResponse = 0;

	theOutput.tag("UniaxialMaterialOutput");
	theOutput.attr("matType", this->getClassType());
	theOutput.attr("matTag", this->getTag());


	// stress
	if (strcmp(argv[0], "stress") == 0) {
		theOutput.tag("ResponseType", "sigma11");
		theResponse = new MaterialResponse(this, 1, this->getStress());
	}
	// tangent
	else if (strcmp(argv[0], "tangent") == 0) {
		theOutput.tag("ResponseType", "C11");
		theResponse = new MaterialResponse(this, 2, this->getTangent());
	}

	// strain
	else if (strcmp(argv[0], "strain") == 0) {
		theOutput.tag("ResponseType", "eps11");
		theResponse = new MaterialResponse(this, 3, this->getStrain());
	}

	// strain
	else if ((strcmp(argv[0], "stressStrain") == 0) ||
		(strcmp(argv[0], "stressANDstrain") == 0)) {
		theOutput.tag("ResponseType", "sig11");
		theOutput.tag("ResponseType", "eps11");
		theResponse = new MaterialResponse(this, 4, Vector(2));
	}

	// damage 
	else if (strcmp(argv[0], "damage") == 0) {
		theResponse = new MaterialResponse(this, 5, this->getDI());
		theOutput.tag("ResponseType", "DI");
		// added 6/9/2006
	}

	else if (strcmp(argv[0], "failure") == 0) {
		int res = 0;
		theResponse = new MaterialResponse(this, 6, res);
		theOutput.tag("ResponseType", "Failure");
	}
	// end add


	theOutput.endTag();
	return theResponse;
}

int
SteelFractureDI::getResponse(int responseID, Information &matInfo)
{
	static Vector stressStrain(2);
	static Vector cyclesAndRange(6);

	// each subclass must implement its own stuff    
	switch (responseID) {
	case 1:
		matInfo.setDouble(this->getStress());
		return 0;

	case 2:
		matInfo.setDouble(this->getTangent());
		return 0;

	case 3:
		matInfo.setDouble(this->getStrain());
		return 0;

	case 4:
		stressStrain(0) = this->getStress();
		stressStrain(1) = this->getStrain();
		matInfo.setVector(stressStrain);
		return 0;

	case 5:
		matInfo.setDouble(this->getDI());
		return 0;

	case 6:
		if (DI > FI_lim)
			matInfo.setInt(1);
		else
			matInfo.setInt(0);
		return 0;

	default:
		return -1;

	}
}
