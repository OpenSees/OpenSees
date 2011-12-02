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
// $Date: 2003-04-02 22:02:42 $
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
#include <stdlib.h>
#include <HystereticMaterial.h>
#include <OPS_Globals.h>
#include <math.h>
#include <float.h>
#include <Channel.h>

HystereticMaterial::HystereticMaterial(int tag,
			double m1p, double r1p, double m2p, double r2p, double m3p, double r3p,
			double m1n, double r1n, double m2n, double r2n, double m3n, double r3n,
			double px, double py, double d1, double d2, double b):
UniaxialMaterial(tag, MAT_TAG_Hysteretic),
pinchX(px), pinchY(py), damfc1(d1), damfc2(d2), beta(b),
mom1p(m1p), rot1p(r1p), mom2p(m2p), rot2p(r2p), mom3p(m3p), rot3p(r3p),
mom1n(m1n), rot1n(r1n), mom2n(m2n), rot2n(r2n), mom3n(m3n), rot3n(r3n)
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
	
	if (error) {
	  opserr << "HystereticMaterial::HystereticMaterial -- input backbone is not unique (one-to-one)\n";
	  exit(-1);
	}		

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
pinchX(px), pinchY(py), damfc1(d1), damfc2(d2), beta(b),
mom1p(m1p), rot1p(r1p), mom3p(m2p), rot3p(r2p),
mom1n(m1n), rot1n(r1n), mom3n(m2n), rot3n(r2n)
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

	if (error) {
	  opserr << "HystereticMaterial::HystereticMaterial -- input backbone is not unique (one-to-one)\n";
	  exit(-1);
	}

				      

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
pinchX(0.0), pinchY(0.0), damfc1(0.0), damfc2(0.0), beta(0.0),
mom1p(0.0), rot1p(0.0), mom2p(0.0), rot2p(0.0), mom3p(0.0), rot3p(0.0),
mom1n(0.0), rot1n(0.0), mom2n(0.0), rot2n(0.0), mom3n(0.0), rot3n(0.0)
{

}

HystereticMaterial::~HystereticMaterial()
{
	// Nothing to do
}

int
HystereticMaterial::setTrialStrain(double strain, double strainRate)
{
	TrotMax = CrotMax;
	TrotMin = CrotMin;
	TenergyD = CenergyD;
	TrotPu = CrotPu;
	TrotNu = CrotNu;

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
	else {
	  if (dStrain < 0.0)
	    negativeIncrement(dStrain);
	  else if (dStrain > 0.0)
	    positiveIncrement(dStrain);
	}

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

void
HystereticMaterial::positiveIncrement(double dStrain)
{
	double kn = pow(CrotMin/rot1n,beta);
	kn = (kn < 1.0) ? 1.0 : 1.0/kn;
	double kp = pow(CrotMax/rot1p,beta);
	kp = (kp < 1.0) ? 1.0 : 1.0/kp;

	if (TloadIndicator == 2) {
		TloadIndicator = 1;
		if (Cstress <= 0.0) {
			TrotNu = Cstrain - Cstress/(E1n*kn);
			double energy = CenergyD - 0.5*Cstress/(E1n*kn)*Cstress;
			double damfc = 1.0;
			if (CrotMin < rot1n) {
				damfc += damfc2*energy/energyA;

				if (Cstrain == CrotMin) {
					damfc += damfc1*(CrotMax/rot1p-1.0);
				}
			}

			TrotMax = CrotMax * damfc;
		}
	}

  TloadIndicator = 1;

	TrotMax = (TrotMax > rot1p) ? TrotMax : rot1p;

	double maxmom = posEnvlpStress(TrotMax);
	double rotlim = negEnvlpRotlim(CrotMin);
	double rotrel = (rotlim > TrotNu) ? rotlim : TrotNu;
	rotrel = TrotNu;
	if (negEnvlpStress(CrotMin) >= 0.0)
	  rotrel = rotlim;

	double rotmp1 = rotrel + pinchY*(TrotMax-rotrel);
	double rotmp2 = TrotMax - (1.0-pinchY)*maxmom/(E1p*kp);
	double rotch = rotmp1 + (rotmp2-rotmp1)*pinchX;

	double tmpmo1;
	double tmpmo2;

	if (Tstrain < TrotNu) {
		Ttangent = E1n*kn;
		Tstress = Cstress + Ttangent*dStrain;
		if (Tstress >= 0.0) {
			Tstress = 0.0;
			Ttangent = E1n*1.0e-9;
		}
	}

	else if (Tstrain >= TrotNu && Tstrain < rotch) {
		if (Tstrain <= rotrel) {
			Tstress = 0.0;
			Ttangent = E1p*1.0e-9;
		}
		else {
			Ttangent = maxmom*pinchY/(rotch-rotrel);
			tmpmo1 = Cstress + E1p*kp*dStrain;
			tmpmo2 = (Tstrain-rotrel)*Ttangent;
			if (tmpmo1 < tmpmo2) {
				Tstress = tmpmo1;
				Ttangent = E1p*kp;
			}
			else
				Tstress = tmpmo2;
		}
	}

	else {
		Ttangent = (1.0-pinchY)*maxmom/(TrotMax-rotch);
		tmpmo1 = Cstress + E1p*kp*dStrain;
		tmpmo2 = pinchY*maxmom + (Tstrain-rotch)*Ttangent;
		if (tmpmo1 < tmpmo2) {
			Tstress = tmpmo1;
			Ttangent = E1p*kp;
		}
		else
			Tstress = tmpmo2;
	}
}

void
HystereticMaterial::negativeIncrement(double dStrain)
{
	double kn = pow(CrotMin/rot1n,beta);
	kn = (kn < 1.0) ? 1.0 : 1.0/kn;
	double kp = pow(CrotMax/rot1p,beta);
	kp = (kp < 1.0) ? 1.0 : 1.0/kp;

	if (TloadIndicator == 1) {
		TloadIndicator = 2;
		if (Cstress >= 0.0) {
			TrotPu = Cstrain - Cstress/(E1p*kp);
			double energy = CenergyD - 0.5*Cstress/(E1p*kp)*Cstress;
			double damfc = 1.0;
			if (CrotMax > rot1p) {
				damfc += damfc2*energy/energyA;

				if (Cstrain == CrotMax) {
					damfc += damfc1*(CrotMin/rot1n-1.0);
				}
			}

			TrotMin = CrotMin * damfc;
		}
	}

  TloadIndicator = 2;

	TrotMin = (TrotMin < rot1n) ? TrotMin : rot1n;

	double minmom = negEnvlpStress(TrotMin);
	double rotlim = posEnvlpRotlim(CrotMax);
	double rotrel = (rotlim < TrotPu) ? rotlim : TrotPu;
	rotrel = TrotPu;
	if (posEnvlpStress(CrotMax) <= 0.0)
	  rotrel = rotlim;

	double rotmp1 = rotrel + pinchY*(TrotMin-rotrel);
	double rotmp2 = TrotMin - (1.0-pinchY)*minmom/(E1n*kn);
	double rotch = rotmp1 + (rotmp2-rotmp1)*pinchX;

	double tmpmo1;
	double tmpmo2;

	if (Tstrain > TrotPu) {
		Ttangent = E1p*kp;
		Tstress = Cstress + Ttangent*dStrain;
		if (Tstress <= 0.0) {
			Tstress = 0.0;
			Ttangent = E1p*1.0e-9;
		}
	}

	else if (Tstrain <= TrotPu && Tstrain > rotch) {
		if (Tstrain >= rotrel) {
			Tstress = 0.0;
			Ttangent = E1n*1.0e-9;
		}
		else {
			Ttangent = minmom*pinchY/(rotch-rotrel);
			tmpmo1 = Cstress + E1n*kn*dStrain;
			tmpmo2 = (Tstrain-rotrel)*Ttangent;
			if (tmpmo1 > tmpmo2) {
				Tstress = tmpmo1;
				Ttangent = E1n*kn;
			}
			else
				Tstress = tmpmo2;
		}
	}

	else {
		Ttangent = (1.0-pinchY)*minmom/(TrotMin-rotch);
		tmpmo1 = Cstress + E1n*kn*dStrain;
		tmpmo2 = pinchY*minmom + (Tstrain-rotch)*Ttangent;
		if (tmpmo1 > tmpmo2) {
			Tstress = tmpmo1;
			Ttangent = E1n*kn;
		}
		else
			Tstress = tmpmo2;
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

	Tstrain = 0;
	Tstress = 0;
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
  int res = 0;
  
  static Vector data(27);
  
  data(0) = this->getTag();
  data(1) = mom1p;
  data(2) = rot1p;
  data(3) = mom2p;
  data(4) = rot2p;
  data(5) = mom3p;
  data(6) = rot3p;
  data(7) = mom1n;
  data(8) = rot1n;
  data(9) = mom2n;
  data(10) = rot2n;
  data(11) = mom3n;
  data(12) = rot3n;
  data(13) = pinchX;
  data(14) = pinchY;
  data(15) = damfc1;
  data(16) = damfc2;
  data(17) = beta;
  data(18) = CrotMax;
  data(19) = CrotMin;
  data(20) = CrotPu;
  data(21) = CrotNu;
  data(22) = CenergyD;
  data(23) = CloadIndicator;
  data(24) = Cstress;
  data(25) = Cstrain;
  data(26) = Ttangent;

  res = theChannel.sendVector(this->getDbTag(), commitTag, data);
  if (res < 0) 
    opserr << "HystereticMaterial::sendSelf() - failed to send data\n";


  return res;
}

int
HystereticMaterial::recvSelf(int commitTag, Channel &theChannel, 
			FEM_ObjectBroker &theBroker)
{
  int res = 0;
  
  static Vector data(27);
  res = theChannel.recvVector(this->getDbTag(), commitTag, data);
  
  if (res < 0) {
      opserr << "HystereticMaterial::recvSelf() - failed to receive data\n";
      return res;
  }
  else {
    this->setTag((int)data(0));
    mom1p = data(1);
    rot1p = data(2);
    mom2p = data(3);
    rot2p = data(4);
    mom3p = data(5);
    rot3p = data(6);
    mom1n = data(7);
    rot1n = data(8);
    mom2n = data(9);
    rot2n = data(10);
    mom3n = data(11);
    rot3n = data(12);
    pinchX = data(13);
    pinchY = data(14);
    damfc1 = data(15);
    damfc2 = data(16);
    beta = data(17);
    CrotMax = data(18);
    CrotMin = data(19);
    CrotPu = data(20);
    CrotNu = data(21);
    CenergyD = data(22);
    CloadIndicator = int(data(23));
    Cstress = data(24);
    Cstrain = data(25);
    Ttangent = data(26);

    // set the trial values
    TrotMax = CrotMax;
    TrotMin = CrotMin;
    TrotPu = CrotPu;
    TrotNu = CrotNu;
    TenergyD = CenergyD;
    TloadIndicator = CloadIndicator;
    Tstress = Cstress;
    Tstrain = Cstrain;

  }

  return 0;
}
    
void
HystereticMaterial::Print(OPS_Stream &s, int flag)
{
	s << "Hysteretic Material, tag: " << this->getTag() << endln;
	s << "mom1p: " << mom1p << endln;
	s << "rot1p: " << rot1p << endln;
	s << "E1p: " << E1p << endln;
	s << "mom2p: " << mom2p << endln;
	s << "rot2p: " << rot2p << endln;
	s << "E2p: " << E2p << endln;
	s << "mom3p: " << mom3p << endln;
	s << "rot3p: " << rot3p << endln;
	s << "E3p: " << E3p << endln;

	s << "mom1n: " << mom1n << endln;
	s << "rot1n: " << rot1n << endln;
	s << "E1n: " << E1n << endln;
	s << "mom2n: " << mom2n << endln;
	s << "rot2n: " << rot2n << endln;
	s << "E2n: " << E2n << endln;
	s << "mom3n: " << mom3n << endln;
	s << "rot3n: " << rot3n << endln;
	s << "E3n: " << E3n << endln;

	s << "pinchX: " << pinchX << endln;
	s << "pinchY: " << pinchY << endln;
	s << "damfc1: " << damfc1 << endln;
	s << "damfc2: " << damfc2 << endln;
	s << "energyA: " << energyA << endln;
	s << "beta: " << beta << endln;
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
	if (strain <= 0.0)
		return 0.0;
	else if (strain <= rot1p)
		return E1p*strain;
	else if (strain <= rot2p)
		return mom1p + E2p*(strain-rot1p);
	else if (strain <= rot3p || E3p > 0.0)
		return mom2p + E3p*(strain-rot2p);
	else
		return mom3p;
}

double
HystereticMaterial::negEnvlpStress(double strain)
{
	if (strain >= 0.0)
		return 0.0;
	else if (strain >= rot1n)
		return E1n*strain;
	else if (strain >= rot2n)
		return mom1n + E2n*(strain-rot1n);
	else if (strain >= rot3n || E3n > 0.0)
		return mom2n + E3n*(strain-rot2n);
	else
		return mom3n;
}

double
HystereticMaterial::posEnvlpTangent(double strain)
{
	if (strain < 0.0)
		return E1p*1.0e-9;
	else if (strain <= rot1p)
		return E1p;
	else if (strain <= rot2p)
		return E2p;
	else if (strain <= rot3p || E3p > 0.0)
		return E3p;
	else
		return E1p*1.0e-9;
}

double
HystereticMaterial::negEnvlpTangent(double strain)
{
	if (strain > 0.0)
		return E1n*1.0e-9;
	else if (strain >= rot1n)
		return E1n;
	else if (strain >= rot2n)
		return E2n;
	else if (strain >= rot3n || E3n > 0.0)
		return E3n;
	else
		return E1n*1.0e-9;
}

double
HystereticMaterial::posEnvlpRotlim(double strain)
{
  double strainLimit = POS_INF_STRAIN;

  if (strain <= rot1p)
    return POS_INF_STRAIN;
  if (strain > rot1p && strain <= rot2p && E2p < 0.0)
    strainLimit = rot1p - mom1p/E2p;
  if (strain > rot2p && E3p < 0.0)
    strainLimit = rot2p - mom2p/E3p;

  if (strainLimit == POS_INF_STRAIN)
    return POS_INF_STRAIN;
  else if (posEnvlpStress(strainLimit) > 0)
    return POS_INF_STRAIN;
  else
    return strainLimit;
}

double
HystereticMaterial::negEnvlpRotlim(double strain)
{
  double strainLimit = NEG_INF_STRAIN;

  if (strain >= rot1n)
    return NEG_INF_STRAIN;
  if (strain < rot1n && strain >= rot2n && E2n < 0.0)
    strainLimit = rot1n - mom1n/E2n;
  if (strain < rot2n && E3n < 0.0)
    strainLimit = rot2n - mom2n/E3n;

  if (strainLimit == NEG_INF_STRAIN)
    return NEG_INF_STRAIN;
  else if (negEnvlpStress(strainLimit) < 0)
    return NEG_INF_STRAIN;
  else
    return strainLimit;
}
