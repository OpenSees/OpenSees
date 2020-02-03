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
                                                                        
// $Source: /usr/local/cvs/OpenSees/SRC/material/uniaxial/SPSW02.cpp,v $
//
//
// 
//
// Written: S. A. Jalalli 03/2015
//Reference: S.A. Jalali and M. Banazadeh, "Development of a new deteriorating hysteresis model for seismic collapse assessment of thin steel plate shear walls"
//link to reference: http://www.sciencedirect.com/science/article/pii/S026382311630249X
// Purpose: This file contains the implementation for the SPSW02 class.
//
//////////////////////////////////////////////////////////////////////

#include <elementAPI.h>
#include "SPSW02.h"
#include <MaterialResponse.h>

#include <Vector.h>
#include <Channel.h>
#include <math.h>
#include <Information.h>


#ifdef _USRDLL
#define OPS_Export extern "C" _declspec(dllexport)
#elif _MACOSX
#define OPS_Export extern "C" __attribute__((visibility("default")))
#else
#define OPS_Export extern "C"
#endif

static int numSPSWcall = 0;

void*
OPS_SPSW02()
{
	if (numSPSWcall == 0) {
		opserr << "------ SPSW02 unaxialMaterial, Written by SAJalali @ Amirkabir University of Technology, Tehran, 2015-------\n";
		opserr << "------------------------------ Please Send Comments to: seyyed-jalali@aut.ac.ir-----------------------------\n";
		opserr<<  "-------Syntax:\n";
		opserr<<  "-------UniaxialMaterial SPSW02 tag ";
		opserr<<  "-------E0 b <-geom Fpy t h l> <-params Fts Fcs cmpUnldngEFac sigTEFac sigTFfac epsTFfac> -R $R -Damage epsPCFac pstCapEFac gama c resFac\n\n";
		opserr << "------------------------------------------------------------------------------------------------------------\n\n\n";
			numSPSWcall = 1;
	}
	int tag;
	double fpy, E0,  b, t, hs, R, l, gama, c, epsPCFac, pstCapEFac, resFac;
	double Fts, Fcs, cmpUnldngEFac, sigTEFac, sigTFfac, epsTFfac;
	fpy = 0;
	bool paramsSet = false;
	UniaxialMaterial *theMaterial = 0;
	int argc = OPS_GetNumRemainingInputArgs();
	int numData = 1;
	int curArg = 2;
	if (OPS_GetIntInput(&numData, &tag) != 0) {
		opserr << "WARNING invalid -tag";
		opserr << "uniaxialMaterial SPSW02: " << tag << endln;
		return 0;
	}
	curArg ++;

	if (OPS_GetDoubleInput (&numData, &E0) != 0) {
		opserr << "WARNING invalid -E0";
		opserr << "uniaxialMaterial SPSW02: " << tag << endln;
		return 0;
	}
	curArg ++;

	if (OPS_GetDoubleInput (&numData, &b) != 0) {
		opserr << "WARNING invalid -b";
		opserr << "uniaxialMaterial SPSW02: " << tag << endln;
		return 0;
	}
	curArg ++;

	const char* str = OPS_GetString();
	curArg ++;

	if (strcmp(str , "-geom") == 0)
	{
		if (OPS_GetDoubleInput (&numData, &fpy) != 0) {
			opserr << "WARNING invalid -Fts";
			opserr << "uniaxialMaterial SPSW02: " << tag << endln;
			return 0;
		}
		curArg ++;

		if (OPS_GetDoubleInput (&numData, &t) != 0) {
			opserr << "WARNING invalid -t";
			opserr << "uniaxialMaterial SPSW02: " << tag << endln;
			return 0;
		}
		curArg ++;

		if (OPS_GetDoubleInput (&numData, &hs) != 0) {
			opserr << "WARNING invalid -h";
			opserr << "uniaxialMaterial SPSW02: " << tag << endln;
			return 0;
		}
		curArg ++;

		if (OPS_GetDoubleInput (&numData, &l) != 0) {
			opserr << "WARNING invalid -l";
			opserr << "uniaxialMaterial SPSW02: " << tag << endln;
			return 0;
		}
		curArg ++;
	} else if (strcmp(str , "-params") == 0)
	{
		paramsSet = true;
		//Fts, Fcs, cmpUnldngEFac, sigTEFac, sigTFfac, epsTFfac
		if (OPS_GetDoubleInput (&numData, &Fts) != 0) {
			opserr << "WARNING invalid Fts";
			opserr << "uniaxialMaterial SPSW02: " << tag << endln;
			return 0;
		}
		curArg ++;

		if (OPS_GetDoubleInput (&numData, &Fcs) != 0) {
			opserr << "WARNING invalid Fcs";
			opserr << "uniaxialMaterial SPSW02: " << tag << endln;
			return 0;
		}
		curArg ++;

		if (OPS_GetDoubleInput (&numData, &cmpUnldngEFac) != 0) {
			opserr << "WARNING invalid cmpUnldngEFac";
			opserr << "uniaxialMaterial SPSW02: " << tag << endln;
			return 0;
		}
		curArg ++;

		if (OPS_GetDoubleInput (&numData, &sigTEFac) != 0) {
			opserr << "WARNING invalid sigTEFac";
			opserr << "uniaxialMaterial SPSW02: " << tag << endln;
			return 0;
		}
		curArg ++;

		if (OPS_GetDoubleInput (&numData, &sigTFfac) != 0) {
			opserr << "WARNING invalid sigTFfac";
			opserr << "uniaxialMaterial SPSW02: " << tag << endln;
			return 0;
		}
		curArg ++;

		if (OPS_GetDoubleInput (&numData, &epsTFfac) != 0) {
			opserr << "WARNING invalid epsTFfac";
			opserr << "uniaxialMaterial SPSW02: " << tag << endln;
			return 0;
		}
		curArg ++;
	}
	if (fpy == 0 && !paramsSet) {
		opserr << "WARNING at least one of -params or -geom options must be provided";
		opserr << "uniaxialMaterial SPSW02: " << tag << endln;
		return 0;
	}
	if (fpy != 0 && paramsSet) {
		opserr << "WARNING both -params and -geom options cannot be used at the same time";
		opserr << "uniaxialMaterial SPSW02: " << tag << endln;
		return 0;
	}
	R = 50;
	
	if (argc >= curArg + 1)
	{
		/*char[] str[10];
		OPS_GetString(str, 5);*/		//for OpenSees 2.4.5 and older
		str = OPS_GetString();
		curArg ++;
		if (strcmp(str , "-R") == 0)
		{
			if (OPS_GetDoubleInput (&numData, &R) != 0) {
				opserr << "WARNING invalid -R";
				opserr << "uniaxialMaterial SPSW02: " << tag << endln;
				return 0;
			}
			curArg ++;
		}
	}
	epsPCFac = 1.e20;
	pstCapEFac = b;
	gama = 10000;
	c = 1;
	resFac = 0.001;
	if (argc >= curArg + 1)
	{
		/*char[] str[10];
		OPS_GetString(str, 124);*/
		str = OPS_GetString();
		curArg++;
		if (strcmp(str , "-Damage") == 0 || strcmp(str , "-damage") == 0)
		{
			if (OPS_GetDoubleInput (&numData, &epsPCFac) != 0) {
				opserr << "WARNING invalid -epsPCFac";
				opserr << "uniaxialMaterial SPSW02: " << tag << endln;
				return 0;
			}
			curArg++;
	
			if (OPS_GetDoubleInput (&numData, &pstCapEFac) != 0) {
				opserr << "WARNING invalid -pstCapEFac";
				opserr << "uniaxialMaterial SPSW02: " << tag << endln;
				return 0;
			}
			curArg++;
	
			if (OPS_GetDoubleInput (&numData, &gama) != 0) {
				opserr << "WARNING invalid -gama";
				opserr << "uniaxialMaterial SPSW02: " << tag << endln;
				return 0;
			}
			curArg++;
	
			if (OPS_GetDoubleInput (&numData, &c) != 0) {
				opserr << "WARNING invalid -c";
				opserr << "uniaxialMaterial SPSW02: " << tag << endln;
				return 0;
			}
			curArg++;
	
			if (OPS_GetDoubleInput (&numData, &resFac) != 0) {
				opserr << "WARNING invalid -resFac";
				opserr << "uniaxialMaterial SPSW02: " << tag << endln;
				return 0;
			}
			curArg++;
		}
	}
	if (paramsSet)
		theMaterial = new SPSW02(tag, E0,  b, Fts, Fcs, cmpUnldngEFac, sigTEFac, sigTFfac, epsTFfac, R, epsPCFac, pstCapEFac, gama, c, resFac);
	else
		theMaterial = new SPSW02(tag, fpy, E0,  b, t, hs, l, R, epsPCFac, pstCapEFac, gama, c, resFac);
	//opserr<<"ok\n";

	return theMaterial;
}


SPSW02::SPSW02(int tag, double _fpy, double _E0, double _b, double _t, double _hs,
		double _l, double _R, double _epsPCFac, double _pstCapEFac, double _gama, double _c, double _resFac):
	UniaxialMaterial(tag, MAT_TAG_SPSW02), fpy(_fpy), E0(_E0), b(_b), t(_t), hs(_hs), l(_l),
	R(_R), epsPCFac(_epsPCFac), pstcpEFac(_pstCapEFac), gama(_gama), c(_c), resFac(_resFac)
{
	givenParams = false;
	excurEnerg = totalEnerg = beta = 0;
	excurEnergP = totalEnergP = betaP = 0;
	cmpUnldngEFac = 0.2;
	sigTEFac = 0.5;						//to be multiplied by Fcs
	sigTFfac = 0.2;						//to be multiplied by Fts
	epsTFfac = 0.5;
	this->Calc_sigcr(/*Fts, Fcs*/);
	FTS = Fts;
	FCS = Fcs;
	FailEnerg = gama * Fts * Fts/E0;
	epsmaxP	=	Fts / E0	;
	sigmaxP	=	Fts	;
	epss0P	=	0.0	;
	sigs0P	=	0.0	;
	epssrP	=	0.0	;
	sigsrP	=	0.0	;
	epsTFP	=	0.0	;
	plstrP	=	0.0	;
	konP	=	0.0	;
	eP		=	0.0	;
	sigP	=	0.0	;
	epsP	=	0.0	;
	sig		=	0.0	;
	eps		=	0.0	;
	e		=	0.0	;
}

SPSW02::SPSW02(int tag, double _E0, double _b, double _FTS, double _FCS, double _cmpUnldngEFac,
		double _sigTEFac, double _sigTFfac, double _epsTFfac, double _R, double _epsPCFac,
		double _pstCapEFac, double _gama, double _c, double _resFac):
		UniaxialMaterial(tag, MAT_TAG_SPSW02), E0(_E0), b(_b), FTS(_FTS), FCS(_FCS), cmpUnldngEFac(_cmpUnldngEFac),
		sigTEFac(_sigTEFac), sigTFfac(_sigTFfac), epsTFfac(_epsTFfac),
		R(_R), epsPCFac(_epsPCFac), pstcpEFac(_pstCapEFac), gama(_gama), c(_c), resFac(_resFac)
{
	givenParams = true;
	excurEnerg = totalEnerg = beta = 0;
	excurEnergP = totalEnergP = betaP = 0;
	Fts = FTS;
	Fcs = FCS;
	FailEnerg = gama * Fts * Fts/E0;
	epsmaxP	=	Fts / E0	;
	sigmaxP	=	Fts	;
	epss0P	=	0.0	;
	sigs0P	=	0.0	;
	epssrP	=	0.0	;
	sigsrP	=	0.0	;
	epsTFP	=	0.0	;
	plstrP	=	0.0	;
	konP	=	0.0	;
	eP		=	0.0	;
	sigP	=	0.0	;
	epsP	=	0.0	;
	sig		=	0.0	;
	eps		=	0.0	;
	e		=	0.0	;
}

SPSW02::SPSW02():
	UniaxialMaterial(0, MAT_TAG_SPSW02)
{
	konP = 0;
}

SPSW02::~SPSW02()
{

}

UniaxialMaterial* SPSW02::getCopy()
{
	SPSW02* theCopy;
	if (givenParams)
		theCopy = new SPSW02(this->getTag(), E0, b, FTS, FCS, cmpUnldngEFac, sigTEFac, sigTFfac, epsTFfac, R, epsPCFac, pstcpEFac, gama, c, resFac);
	else
		theCopy = new SPSW02(this->getTag(), fpy, E0, b, t, hs, l, R, epsPCFac, pstcpEFac, gama, c, resFac);
	return theCopy;
}

double SPSW02::getInitialTangent()
{
	return E0;
}

int SPSW02::setTrialStrain(double trialStrain, double strainRate)
{
	double resF = resFac * FTS;
	double Esh = b * E0;
	double epsyn = FCS / E0;
	double epsyp = FTS / E0;
	eps = trialStrain;
	double deps = eps - epsP;
	epsmax = epsmaxP;
	sigmax = sigmaxP;
	epss0 = epss0P;
	sigs0 = sigs0P;
	epsr = epssrP;
	sigr = sigsrP;
	epsTF = epsTFP;
	plstr = plstrP;
	kon = konP;
	double sigTF = sigTFfac*sigmax;			//sig at tension field redevelopment point
	double epsPC = epsPCFac * epsyp;

	//kon == 0 (value at start): elastic with (fabs(eps) < epsyn || eps < -epsyn)
	//kon == 11: comp. branch
		//kon == 12: reversal from comp. when sarts at eps >  plstr - 2.*Fcs/(cmpUnldngEFac*E0) (pinching branch excluded)
		//kon == 13: reversal from comp. when sarts at eps <  plstr - 2.*Fcs/(cmpUnldngEFac*E0)  (pinching branch included)
	//kon == 21: Menegotto-Pinto tension loading

	if (kon == 0)
	{
		if (fabs(eps) <= epsyn)
		{
			sig = eps * E0;
			e = E0;
		}
		else if (eps < -epsyn)
		{
			if (deps > 0)
				kon = 11;
			else
			{
				sig = -Fcs;
				e = 1.0e-15*E0;
			}
		}
		else if (eps > epsyn)
		{
			kon = 21;
			epsr = 0.0;
			sigr = 0.0;
			epss0 = epsyp;
			sigs0 = Fts;
		}
	}
	
	
	if (kon == 11)
	{
		if (deps > 0)
		{
			if (epsP > plstr - 2.*Fcs/(cmpUnldngEFac*E0))
				kon = 12;
			else
			{
				kon = 13;
				epsTF = plstr - epsTFfac*(plstr - (epsP + Fcs / (cmpUnldngEFac*E0))) + sigTF / E0;
				epsr = epsP;
				sigr = sigP;
				epss0 = epsr - (sigr - sigTEFac*Fcs) / (cmpUnldngEFac*E0);
				sigs0 = sigTEFac*Fcs;
			}
		}
		else
		{
			sig = -Fcs;
			e = 1.0e-15*E0;
		}
	}

	// compression reloading branch
	if (kon == 12)
	{
		sig = sigP + cmpUnldngEFac*E0*deps;
		e = cmpUnldngEFac*E0;
		double sigUnloading = E0 * (eps - plstr);
		if (sig <= -Fcs)
		{
			kon = 11;
			sig = -Fcs;
			e = 1.0e-15*E0;
		}
		else if (sig <= sigUnloading)		// check to see if reloading branch coincides with the previous unloading branch
		{
			kon = 21;
			sig = sigUnloading;
			e= E0;
			epsr = eps;
			sigr = sig;
			double ET = (Fts - sigr) / (epsmax - epsr);
			epss0 = (Fts - Esh * epsyp - sigr + ET * epsr) / (ET - Esh);
			sigs0 = Fts + Esh * (epss0 - epsyp);
			//this->MenegottoPinto(eps, Esh, R, sig, e);
			return 0;
		}
	}
	if (kon == 13)
	{
		if (eps > epsTF)
		{
			kon = 21;
			epsr = epsTF;
			sigr = sigTF;
			double ET = (Fts - sigr) / (epsmax - epsr);
			epss0 = (Fts - Esh * epsyp - sigr + ET * epsr) / (ET - Esh);
			sigs0 = Fts + Esh * (epss0 - epsyp);
			//this->MenegottoPinto(eps, Esh, R, sig, e);
		}
		else
		{
			double bT = (sigTF - sigs0) / (epsTF - epss0)/*/(0.2*E0)*/;
			double sigEnvlp, eEnvlp;
			this->MenegottoPinto(eps, bT, R, sigEnvlp, eEnvlp);
			sig = sigP + cmpUnldngEFac*E0*deps;
			e = cmpUnldngEFac*E0;
			if (sig <= -Fcs)
			{
				kon = 11; sig = -Fcs;
				e = 1.0e-15*E0;
			}
			//commented to handle non-convergence ny using the same branch in loading and unloading directions:
			else /*if (sig >= sigEnvlp)*/
			{
				sig = sigEnvlp;
				e = eEnvlp;
			}
		}
	}
	// tension branch
	if (kon == 21)
	{
		if (epsP > epsmax)
		{
			sigmax = sigP;
			epsmax = epsP;
		}
		if (epsmax > epsPC)
		{
			if (deps > 0)
			{
				
				if (eps < epsmax)
				{
					e = (sigmax - sigr) / (epsmax - epsr);
					sig = sigP + e * deps;
				}
				else
				{
					e = pstcpEFac * E0;
					sig = sigP + e * deps;
					if (sig < resF)
					{
						sig = resF;
						e = 1.0e-15*E0;
					}
				}
			}
			else
			{
				sig = sigP + E0 * deps;
				epsr = eps;
				sigr = sig;
				e = E0;
				if (sig <= -Fcs)
				{
					kon = 11;
					sig = -Fcs;
					e = 1.0e-15*E0;
					plstr = epsP - sigP/E0;
				}
			}
		}
		else
		{
			double sigEnvlp, eEnvlp;
			this->MenegottoPinto(eps, Esh, R, sigEnvlp, eEnvlp);
			sig = sigP + E0 * deps;
			e = E0;
			if (sig <= -Fcs)
			{
				kon = 11;
				sig = -Fcs;
				e = 1.0e-15*E0;
				plstr = epsP - sigP/E0;
			}
			else if (sig >= sigEnvlp)
			{
				sig = sigEnvlp;
				e = eEnvlp;
			}
		}
	}
	return 0;
}

void SPSW02::updateDamage()
{
	if ( (sig < 0.0 && sigP >= 0) || (sig == 0 && sigP > 0))
	{
		//submit damage and reset for new excursion
		double zeroSigEps = epsP - sigP/E0;
		double dE = 0.5 * sigP * (zeroSigEps - epsP);
		totalEnerg += dE;
		if (totalEnerg < 0) totalEnerg = 0.;
		if (gama > 9999)
		{
			return;
		}
		excurEnerg += dE;
		if (excurEnerg < 0) excurEnerg = 0.;
		beta = pow( excurEnerg / ( FailEnerg - totalEnerg) , c );
		if (beta > 0.999 || beta < 0)
		{
			opserr<< "\nSPSW02:"<< this->getTag()<< " WARNING! Maximum Energy Absorbance Capacity Reached\n"<< endln;
			beta = 0.999;

		}
		sigmaxP = (1. - beta)*sigmaxP + beta * resFac*FTS;
		Fts = (1. - beta)*Fts + beta * resFac*FTS;
		if (Fcs > Fts)
			Fcs = Fts;
		excurEnerg = 0.0;
	} 
	else if ( sig > 0.0 )
	{
		double dE = 0.5 * (sig + sigP) * (eps - epsP);
		excurEnerg += dE;
		totalEnerg += dE;
    }
}

double SPSW02::getStrain()
{
	return eps;
}

double SPSW02::getStress()
{
	return sig;
}

double SPSW02::getTangent()
{
	return e;
}

int SPSW02::commitState()
{
	epsmaxP = epsmax;
	sigmaxP = sigmax;
	epss0P = epss0;
	sigs0P = sigs0;
	epssrP = epsr;
	sigsrP = sigr;
	epsTFP = epsTF;
	plstrP = plstr;
	konP = kon;

	this->updateDamage();
	eP = e;
	sigP = sig;
	epsP = eps;
	totalEnergP = totalEnerg;
	excurEnergP = excurEnerg;
	
	return 0;
}

int SPSW02::revertToLastCommit()
{
	epsmax	=	epsmaxP	;
	sigmax	=	sigmaxP	;
	epss0	=	epss0P	;
	sigs0	=	sigs0P	;
	epsr	=	epssrP	;
	sigr	=	sigsrP	;
	epsTF	=	epsTFP	;
	plstr	=	plstrP	;
	kon		=	konP	;
	e		=	eP		;
	sig		=	sigP	;
	eps		=	epsP	;
	excurEnerg = excurEnergP;
	totalEnerg = totalEnergP;
	beta = betaP;
	
	return 0;
}

int SPSW02::revertToStart()
{
	opserr <<"revert called\n";
	excurEnerg = totalEnerg = beta = 0;
	excurEnergP = totalEnergP = betaP = 0;
	Fts = FTS;
	Fcs = FCS;
	if (!givenParams)
	{
		this->Calc_sigcr(/*Fts, Fcs*/);
		FTS = Fts;
		FCS = Fcs;
	}
	FailEnerg = gama * Fts * Fts/E0;
	epsmaxP	=	Fts / E0	;
	sigmaxP	=	Fts	;
	epss0P	=	0.0	;
	sigs0P	=	0.0	;
	epssrP	=	0.0	;
	sigsrP	=	0.0	;
	epsTFP	=	0.0	;
	plstrP	=	0.0	;
	konP	=	0.0	;
	eP		=	0.0	;
	sigP	=	0.0	;
	epsP	=	0.0	;
	sig		=	0.0	;
	eps		=	0.0	;
	e		=	0.0	;
	return 0;
}

int SPSW02::sendSelf(int commitTag, Channel & theChannel)
{
				
	int res = 0;
	static Vector data(38);
	data(0)  = this->getTag();
	data(1)  = t			;
	data(2)  = hs			;
	data(3)  = l			;
	data(4)  = fpy			;
	data(5)  = E0			;
	data(6)  = b			;
	data(7)  = R			;
	data(8)  = Fts			;
	data(9)  = Fcs			;
	data(10) = FTS			;
	data(11) = FCS			;
	data(12) = epsPCFac		;
	data(13) = pstcpEFac	;
	data(14) = gama			;
	data(15) = FailEnerg	;
	data(16) = c			;
	data(17) = resFac		;
	data(18) = givenParams	;
	data(19) = cmpUnldngEFac;
	data(20) = sigTEFac		;
	data(21) = sigTFfac		;
	data(22) = epsTFfac		;
	data(23) = epsmaxP		;
	data(24) = sigmaxP		;
	data(25) = epss0P		;
	data(26) = sigs0P		;
	data(27) = epssrP		;
	data(28) = sigsrP		;
	data(29) = epsTFP		;
	data(30) = plstrP		;
	data(31) = konP			;
	data(32) = epsP			;
	data(33) = sigP			;
	data(34) = eP			;
	data(35) = excurEnergP	;
	data(36) = totalEnergP	;
	data(37) = betaP		;
	
	// Data is only sent after convergence, so no trial variables
	// need to be sent through data vector

	res = theChannel.sendVector(this->getDbTag(), commitTag, data);
	if (res < 0)
		opserr << "SPSW02::sendSelf() - failed to send data\n";

	return res;
}

int SPSW02::recvSelf(int commitTag, Channel & theChannel, FEM_ObjectBroker & theBroker)
{
	int res = 0;
	static Vector data(38);
	res = theChannel.recvVector(this->getDbTag(), commitTag, data);

	if (res < 0) {
		opserr << "SPSW02::recvSelf() - failed to receive data\n";
		this->setTag(0);
		return res;
	}
	this->setTag(int(data(0)));
	t = data(1);
	hs = data(2);
	l = data(3);
	fpy = data(4);
	E0 = data(5);
	b = data(6);
	R = data(7);
	Fts = data(8);
	Fcs = data(9);
	FTS = data(10);
	FCS = data(11);
	epsPCFac = data(12);
	pstcpEFac = data(13);
	gama = data(14);
	FailEnerg = data(15);
	c = data(16);
	resFac = data(17);
	givenParams = data(18);
	cmpUnldngEFac = data(19);
	sigTEFac = data(20);
	sigTFfac = data(21);
	epsTFfac = data(22);
	epsmaxP = data(23);
	sigmaxP = data(24);
	epss0P = data(25);
	sigs0P = data(26);
	epssrP = data(27);
	sigsrP = data(28);
	epsTFP = data(29);
	plstrP = data(30);
	konP = data(31);
	epsP = data(32);
	sigP = data(33);
	eP = data(34);
	excurEnergP = data(35);
	totalEnergP = data(36);
	betaP = data(37);
	epsmax = epsmaxP;
	sigmax = sigmaxP;
	epss0 = epss0P;
	sigs0 = sigs0P;
	epsr = epssrP;
	sigr = sigsrP;
	epsTF = epsTFP;
	plstr = plstrP;
	kon = konP;
	eps = epsP;
	sig = sigP;
	e = eP;
	excurEnerg = excurEnergP;
	totalEnerg = totalEnergP;
	beta = betaP;
}

void SPSW02::Print(OPS_Stream & s, int flag)
{
	s << "SPSW02:(strain, stress, tangent)" << eps << " " << sig << " " << e << endln;
}

void SPSW02::Calc_sigcr(/*double & _Fts, double &_Fcs*/)
{
	// calculate plate buckling strengthh, sigcr
	double hsToL = hs/l;
	double ks = 5.6 + 8.98 / hsToL / hsToL;
	if (hsToL > 1.0)
		ks = 8.98 + 5.6 / hsToL / hsToL;
	
	double pi = 2 * asin (1.0);
	double nu = 0.3;
	double lToT = l/t;
	//double sigcr = ks * pi * pi * E0 / (12 * (1.0 - nu * nu) * lToT * lToT); 
	Fcs = ks * pi * pi * E0 / (12 * (1.0 - nu * nu) * lToT * lToT); 
	// plate clamped on four edges
	//_Fcs = sigcr * sin(2*alpha/180.0 * pi);									//buckling strength of compression strip
	Fts = pow( fpy * fpy - 0.75 * Fcs * Fcs, 0.5) - 0.5 * Fcs;		//yield strength of tensional strip acording to Von-Mises rule

	/*Fcs = t;
	_Fts = Fts;*/
	return;
}

void SPSW02::MenegottoPinto(double epsc, double E2, double R, double & sigc, double& ec)
{
	double E1 = (sigs0 - sigr) / (epss0 - epsr);
	double bT = E2/E1;
	// double bT = E2;
	double epsrat = (epsc-epsr) / (epss0 - epsr);
	double dum1 = 1.0 + pow (fabs(epsrat),R);
	double dum2 = pow(dum1, 1.0/R);
	sigc = bT * epsrat + (1.0 - bT) * epsrat / dum2;
	sigc = sigc * (sigs0 - sigr) + sigr;
	ec = bT + (1.0 - bT) / (dum1 * dum2);
	ec *= (sigs0 - sigr) / (epss0 - epsr);
}

double SPSW02::getEnergy()
{
	return totalEnergP;
}

//Response* 
//SPSW02::setResponse(const char **argv, int argc,
//			      OPS_Stream &theOutput)
//{
//	Response *theResponse = UniaxialMaterial::setResponse(argv, argc, theOutput);
//	if (theResponse != 0)
//		return theResponse;
//
//	if ( (strcmp(argv[0],"energy") == 0) ||
//		(strcmp(argv[0],"Energy") == 0) ) {
//    
//		theOutput.tag("UniaxialMaterialOutput");
//		theOutput.attr("matType", this->getClassType());
//		theOutput.attr("matTag", this->getTag());
//    
//		theOutput.tag("ResponseType", "energy");
//		theResponse =  new MaterialResponse(this, 6, totalEnergP);
//		theOutput.endTag();
//  }
//  
//  return theResponse;
//
//}
// 
//int 
//SPSW02::getResponse(int responseID, Information &matInfo)
//{
//	if (UniaxialMaterial::getResponse(responseID, matInfo) == 0)
//		return 0;
//	//opserr << "getResponse-2\n";
//	switch (responseID) {
//	case 6:
//		matInfo.setDouble(totalEnergP);
//		return 0;
//	default:      
//		return -1;
//	}
//}
//
