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
                                                                        
// $Revision: 1.1.1.1 $
// $Date: 2000-09-15 08:23:22 $
// $Source: /usr/local/cvs/OpenSees/SRC/material/uniaxial/HystereticMaterial.cpp,v $

// Written: MHS
// Created: July 2000
//
// Description: This file contains the implementation of 
// HystereticMaterial.  HystereticMaterial is
// a one-dimensional hysteretic model with pinching of both
// force and deformation, damage due to deformation and energy, and
// degraded unloading stiffness based on maximum ductility.  This
// is a modified implementation of Hyster2.f90 by Filippou.

#include <HystereticMaterial.h>
#include <G3Globals.h>
#include <math.h>

HystereticMaterial::HystereticMaterial(int tag,
			double m1p, double r1p, double m2p, double r2p, double m3p, double r3p,
			double m1n, double r1n, double m2n, double r2n, double m3n, double r3n,
			double px, double py, double d1, double d2, double b):
UniaxialMaterial(tag, MAT_TAG_Hysteretic),
mom1p(m1p), rot1p(r1p), mom2p(m2p), rot2p(r2p), mom3p(m3p), rot3p(r3p),
mom1n(m1n), rot1n(r1n), mom2n(m2n), rot2n(r2n), mom3n(m3n), rot3n(r3n),
pinchX(px), pinchY(py), damfc1(d1), damfc2(d2), beta(b)
{
	bool error = false;
	
	// Positive backbone parameters
	if (rot1p <= 0.0)
		error = true;

	if (rot2p <= rot1p)
		error = true;

	if (rot3p <= rot2p)
		error = true;

	// Negative backbone parameters
	if (rot1n >= 0.0)
		error = true;

	if (rot2n >= rot1n)
		error = true;

	if (rot3n >= rot2n)
		error = true;
	
	if (error)
		g3ErrorHandler->fatal("%s -- input backbone is not unique (one-to-one)",
		"HystereticMaterial::HystereticMaterial");

	energyA = 0.5 * (rot1p*mom1p + (rot2p-rot1p)*(mom2p+mom1p) + (rot3p-rot2p)*(mom3p+mom2p) +
		rot1n*mom1n + (rot2n-rot1n)*(mom2n+mom1n) * (rot3n-rot2n)*(mom3n+mom2n));

	// Set envelope slopes
	this->setEnvelope();

	// Initialize history variables
	this->revertToStart();
	this->revertToLastCommit();
}

HystereticMaterial::HystereticMaterial(int tag,
			double m1p, double r1p, double m2p, double r2p,
			double m1n, double r1n, double m2n, double r2n,
			double px, double py, double d1, double d2, double b):
UniaxialMaterial(tag, MAT_TAG_Hysteretic),
mom1p(m1p), rot1p(r1p), mom3p(m2p), rot3p(r2p),
mom1n(m1n), rot1n(r1n), mom3n(m2n), rot3n(r2n),
pinchX(px), pinchY(py), damfc1(d1), damfc2(d2), beta(b)
{
	bool error = false;
	
	// Positive backbone parameters
	if (rot1p <= 0.0)
		error = true;

	if (rot3p <= rot1p)
		error = true;

	// Negative backbone parameters
	if (rot1n >= 0.0)
		error = true;

	if (rot3n >= rot1n)
		error = true;

	if (error)
		g3ErrorHandler->fatal("%s -- input backbone is not unique (one-to-one)",
		"HystereticMaterial::HystereticMaterial");

	energyA = 0.5 * (rot1p*mom1p + (rot3p-rot1p)*(mom3p+mom1p) +
		rot1n*mom1n + (rot3n-rot1n)*(mom3n+mom1n));

	mom2p = 0.5*(mom1p+mom3p);
	mom2n = 0.5*(mom1n+mom3n);

	rot2p = 0.5*(rot1p+rot3p);
	rot2n = 0.5*(rot1n+rot3n);

	// Set envelope slopes
	this->setEnvelope();

	// Initialize history variables
	this->revertToStart();
	this->revertToLastCommit();
}

HystereticMaterial::HystereticMaterial():
UniaxialMaterial(0, MAT_TAG_Hysteretic),
mom1p(0.0), rot1p(0.0), mom2p(0.0), rot2p(0.0), mom3p(0.0), rot3p(0.0),
mom1n(0.0), rot1n(0.0), mom2n(0.0), rot2n(0.0), mom3n(0.0), rot3n(0.0),
pinchX(0.0), pinchY(0.0), damfc1(0.0), damfc2(0.0), beta(0.0)
{

}

HystereticMaterial::~HystereticMaterial()
{

}

int
HystereticMaterial::setTrialStrain(double strain, double strainRate)
{
	TrotMax = CrotMax;
	TrotMin = CrotMin;
	TenergyD = CenergyD;

	Tstrain = strain;
	double dStrain = Tstrain - Cstrain;

	TloadIndicator = CloadIndicator;
	
	if (TloadIndicator == 0)
		TloadIndicator = (dStrain < 0.0) ? 2 : 1;

	if (Tstrain >= CrotMax) {
		TrotMax = Tstrain;
		Ttangent = posEnvlpTangent(Tstrain);
		Tstress = posEnvlpStress(Tstrain);
	}
	else if (Tstrain <= CrotMin) {
		TrotMin = Tstrain;
		Ttangent = negEnvlpTangent(Tstrain);
		Tstress = negEnvlpStress(Tstrain);
	}
	else
		(dStrain < 0.0) ? negativeIncrement(dStrain) : positiveIncrement(dStrain);

	TenergyD = CenergyD + 0.5*(Cstress+Tstress)*dStrain;

	return 0;
}


double
HystereticMaterial::getStrain(void)
{
	return Tstrain;
}

double
HystereticMaterial::getStress(void)
{
	return Tstress;
}

double
HystereticMaterial::getTangent(void)
{
	return Ttangent;
}

double
HystereticMaterial::getSecant(void)
{
	return Ttangent;
}

void
HystereticMaterial::positiveIncrement(double dStrain)
{
	TrotNu = CrotNu;

	if (TrotMax < rot1p)
		TrotMax = rot1p;

	double k = pow(TrotMax/rot1p,beta);
	k = (k < 1.0) ? 1.0 : 1.0/k;

	if (TloadIndicator == 2) {
		TloadIndicator = 1;
		if (Cstress <= 0.0) {
			TrotNu = Cstrain - Cstress/(E1n*k);
			double energy = CenergyD - 0.5*Cstress/(E1n*k)*Cstress;
			double damfc = 1.0;
			if (CrotMin < rot1n) {
				damfc += damfc2*energy/energyA;

				if (Cstrain == CrotMin) {
					damfc += damfc1*(CrotMax/rot1p-1.0);
				}
			}

			TrotMax *= damfc;
		}
	}

	double maxmom = posEnvlpStress(TrotMax);
	double rotlim = negEnvlpRotlim(CrotMin);
	double rotrel = (rotlim > TrotNu) ? rotlim : TrotNu;

	double rotmp1 = rotrel + pinchY*(TrotMax-rotrel);
	double rotmp2 = TrotMax - (1.0-pinchY)*maxmom/(E1p*k);
	double rotch = rotmp1 + (rotmp2-rotmp1)*pinchX;

	double tmpmo1;
	double tmpmo2;

	if (Tstrain < TrotNu) {
		Ttangent = E1n*k;
		Tstress = Cstress + Ttangent*dStrain;
		if (Tstress >= 0.0) {
			Tstress = 0.0;
			Ttangent = 1.0e-9;
		}
	}

	else if (Tstrain >= TrotNu && Tstrain < rotch) {
		if (Tstrain <= rotrel) {
			Tstress = 0.0;
			Ttangent = 1.0e-9;
		}
		else {
			Ttangent = maxmom*pinchY/(rotch-rotrel);
			tmpmo1 = Cstress + E1p*k*dStrain;
			tmpmo2 = (Tstrain-rotrel)*Ttangent;
			if (tmpmo1 < tmpmo2) {
				Tstress = tmpmo1;
				Ttangent = E1p*k;
			}
			else {
				Tstress = tmpmo2;
			}
		}
	}

	else {
		Ttangent = (1.0-pinchY)*maxmom/(TrotMax-rotch);
		tmpmo1 = Cstress + E1p*k*dStrain;
		tmpmo2 = pinchY*maxmom + (Tstrain-rotch)*Ttangent;
		if (tmpmo1 < tmpmo2) {
			Tstress = tmpmo1;
			Ttangent = E1p*k;
		}
		else {
			Tstress = tmpmo2;
		}
	}
}

void
HystereticMaterial::negativeIncrement(double dStrain)
{
	TrotPu = CrotPu;

	if (TrotMin > rot1n)
		TrotMin = rot1n;

	double k = pow(TrotMin/rot1n,beta);
	k = (k < 1.0) ? 1.0 : 1.0/k;

	if (TloadIndicator == 1) {
		TloadIndicator = 2;
		if (Cstress >= 0.0) {
			TrotPu = Cstrain - Cstress/(E1p*k);
			double energy = CenergyD - 0.5*Cstress/(E1p*k)*Cstress;
			double damfc = 1.0;
			if (CrotMax > rot1p) {
				damfc += damfc2*energy/energyA;

				if (Cstrain == CrotMax) {
					damfc += damfc1*(CrotMin/rot1n-1.0);
				}
			}

			TrotMin *= damfc;
		}
	}

	double minmom = negEnvlpStress(TrotMin);
	double rotlim = posEnvlpRotlim(CrotMax);
	double rotrel = (rotlim < TrotPu) ? rotlim : TrotPu;

	double rotmp1 = rotrel + pinchY*(TrotMin-rotrel);
	double rotmp2 = TrotMin - (1.0-pinchY)*minmom/(E1n*k);
	double rotch = rotmp1 + (rotmp2-rotmp1)*pinchX;

	double tmpmo1;
	double tmpmo2;

	if (Tstrain > TrotPu) {
		Ttangent = E1p*k;
		Tstress = Cstress + Ttangent*dStrain;
		if (Tstress <= 0.0) {
			Tstress = 0.0;
			Ttangent = 1.0e-9;
		}
	}

	else if (Tstrain <= TrotPu && Tstrain > rotch) {
		if (Tstrain >= rotrel) {
			Tstress = 0.0;
			Ttangent = 1.0e-9;
		}
		else {
			Ttangent = minmom*pinchY/(rotch-rotrel);
			tmpmo1 = Cstress + E1n*k*dStrain;
			tmpmo2 = (Tstrain-rotrel)*Ttangent;
			if (tmpmo1 > tmpmo2) {
				Tstress = tmpmo1;
				Ttangent = E1n*k;
			}
			else {
				Tstress = tmpmo2;
			}
		}
	}

	else {
		Ttangent = (1.0-pinchY)*minmom/(TrotMin-rotch);
		tmpmo1 = Cstress + E1n*k*dStrain;
		tmpmo2 = pinchY*minmom + (Tstrain-rotch)*Ttangent;
		if (tmpmo1 > tmpmo2) {
			Tstress = tmpmo1;
			Ttangent = E1n*k;
		}
		else {
			Tstress = tmpmo2;
		}
	}
}

int
HystereticMaterial::commitState(void)
{
	CrotMax = TrotMax;
	CrotMin = TrotMin;
	CrotPu = TrotPu;
	CrotNu = TrotNu;
	CenergyD = TenergyD;
	CloadIndicator = TloadIndicator;

	Cstress = Tstress;
	Cstrain = Tstrain;

	return 0;
}

int
HystereticMaterial::revertToLastCommit(void)
{
	TrotMax = CrotMax;
	TrotMin = CrotMin;
	TrotPu = CrotPu;
	TrotNu = CrotNu;
	TenergyD = CenergyD;
	TloadIndicator = CloadIndicator;

	Tstress = Cstress;
	Tstrain = Cstrain;

	return 0;
}

int
HystereticMaterial::revertToStart(void)
{
	CrotMax = 0.0;
	CrotMin = 0.0;
	CrotPu = 0.0;
	CrotNu = 0.0;
	CenergyD = 0.0;
	CloadIndicator = 0;

	Cstress = 0.0;
	Cstrain = 0.0;

	Ttangent = E1p;

	return 0;
}

UniaxialMaterial*
HystereticMaterial::getCopy(void)
{
	HystereticMaterial *theCopy = new HystereticMaterial (this->getTag(),
		mom1p, rot1p, mom2p, rot2p, mom3p, rot3p,
		mom1n, rot1n, mom2n, rot2n, mom3n, rot3n,
		pinchX, pinchY, damfc1, damfc2, beta);

	theCopy->CrotMax = CrotMax;
	theCopy->CrotMin = CrotMin;
	theCopy->CrotPu = CrotPu;
	theCopy->CrotNu = CrotNu;
	theCopy->CenergyD = CenergyD;
	theCopy->CloadIndicator = CloadIndicator;
	theCopy->Cstress = Cstress;
	theCopy->Cstrain = Cstrain;

	theCopy->Ttangent = Ttangent;

	return theCopy;
}

int
HystereticMaterial::sendSelf(int commitTag, Channel &theChannel)
{
	return -1;
}

int
HystereticMaterial::recvSelf(int commitTag, Channel &theChannel, 
			FEM_ObjectBroker &theBroker)
{
	return -1;
}
    
void
HystereticMaterial::Print(ostream &s, int flag)
{
	s << "Hysteretic Material, tag: " << this->getTag() << endl;
	s << "mom1p: " << mom1p << endl;
	s << "rot1p: " << rot1p << endl;
	s << "E1p: " << E1p << endl;
	s << "mom2p: " << mom2p << endl;
	s << "rot2p: " << rot2p << endl;
	s << "E2p: " << E2p << endl;
	s << "mom3p: " << mom3p << endl;
	s << "rot3p: " << rot3p << endl;
	s << "E3p: " << E3p << endl;

	s << "mom1n: " << mom1n << endl;
	s << "rot1n: " << rot1n << endl;
	s << "E1n: " << E1n << endl;
	s << "mom2n: " << mom2n << endl;
	s << "rot2n: " << rot2n << endl;
	s << "E2n: " << E2n << endl;
	s << "mom3n: " << mom3n << endl;
	s << "rot3n: " << rot3n << endl;
	s << "E3n: " << E3n << endl;

	s << "pinchX: " << pinchX << endl;
	s << "pinchY: " << pinchY << endl;
	s << "damfc1: " << damfc1 << endl;
	s << "damfc2: " << damfc2 << endl;
	s << "energyA: " << energyA << endl;
	s << "beta: " << beta << endl;
}

void
HystereticMaterial::setEnvelope(void)
{
	E1p = mom1p/rot1p;
	E2p = (mom2p-mom1p)/(rot2p-rot1p);
	E3p = (mom3p-mom2p)/(rot3p-rot2p);

	E1n = mom1n/rot1n;
	E2n = (mom2n-mom1n)/(rot2n-rot1n);
	E3n = (mom3n-mom2n)/(rot3n-rot2n);
}

double
HystereticMaterial::posEnvlpStress(double strain)
{
	if (strain <= 0.0 || strain > posEnvlpRotlim(strain))
		return 0.0;
	else if (strain <= rot1p)
		return E1p*strain;
	else if (strain <= rot2p)
		return mom1p + E2p*(strain-rot1p);
	else
		return mom2p + E3p*(strain-rot2p);
}

double
HystereticMaterial::negEnvlpStress(double strain)
{
	if (strain >= 0.0 || strain < negEnvlpRotlim(strain))
		return 0.0;
	if (strain >= rot1n)
		return E1n*strain;
	if (strain >= rot2n)
		return mom1n + E2n*(strain-rot1n);

	return mom2n + E3n*(strain-rot2n);
}

double
HystereticMaterial::posEnvlpTangent(double strain)
{
	if (strain < 0.0 || strain > posEnvlpRotlim(strain))
		return 1.0e-9;
	else if (strain <= rot1p)
		return E1p;
	else if (strain <= rot2p)
		return E2p;
	else
		return E3p;
}

double
HystereticMaterial::negEnvlpTangent(double strain)
{
	if (strain > 0.0 || strain < negEnvlpRotlim(strain))
		return 1.0e-9;
	else if (strain >= rot1n)
		return E1n;
	else if (strain >= rot2n)
		return E2n;
	else
		return E3n;
}

double
HystereticMaterial::posEnvlpRotlim(double strain)
{
	if (strain <= rot1p)
		return POS_INF_STRAIN;
	if (strain > rot1p && strain <= rot2p && E2p < 0.0)
		return rot1p - mom1p/E2p;
	if (strain > rot2p && E3p < 0.0)
		return rot2p - mom2p/E3p;
	
	return POS_INF_STRAIN;
}

double
HystereticMaterial::negEnvlpRotlim(double strain)
{
	if (strain >= rot1n)
		return NEG_INF_STRAIN;
	if (strain < rot1n && strain >= rot2n && E2n < 0.0)
		return rot1n - mom1n/E2n;
	if (strain < rot2n && E3n < 0.0)
		return rot2n - mom2n/E3n;
	
	return NEG_INF_STRAIN;
}