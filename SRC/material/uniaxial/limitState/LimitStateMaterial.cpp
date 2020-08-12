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

// $Revision: 1.2 $
// $Date: 2006-09-05 22:32:05 $
// $Source: /usr/local/cvs/OpenSees/SRC/material/uniaxial/limitState/LimitStateMaterial.cpp,v $
                                                                        
// Written: KJE
// Created: July 2001
// Modified July 2002
// Based on HystereticMaterial by MHS
//
// Description: This file contains the implementation of 
// LimitStateMaterial.  LimitStateMaterial is
// a one-dimensional hysteretic model with pinching of both
// force and deformation, damage due to deformation and energy, and
// degraded unloading stiffness based on maximum ductility.  Based on 
// HystereticMaterial which is a modified implementation of Hyster2.f90 
// by Filippou.  
//
// Modified July 2002 by KJE to include the use of limitCurve to 
// detect failure in hysteretic material (see commitState).
// Option only available for 3 point specification. Behaves 
// same as HystereticMaterial if no limit curve specified.
//
// All code specific to LimitStateMaterial separated by ////////////////

#include <stdlib.h>
#include "LimitCurve.h"

#include <LimitStateMaterial.h>
#include <G3Globals.h>
#include <math.h>
#include <float.h>
#include <Channel.h>
#include <Vector.h>
#include <MaterialResponse.h>
#include <elementAPI.h>
#include <string>

void* OPS_LimiStateMaterial()
{
    UniaxialMaterial* mat = 0;

    int argc = OPS_GetNumRemainingInputArgs()+2;
    if (argc != 20 && argc != 19 && argc != 16 && argc != 15 && argc != 22 && argc != 23) {
	opserr << "WARNING insufficient arguments\n";
	opserr << "Want: uniaxialMaterial LimitState tag? mom1p? rot1p? mom2p? rot2p? mom3p? rot3p? "
	       << "\nmom1n? rot1n? mom2n? rot2n? mom3n? rot3n? pinchX? pinchY? damfc1? damfc2? beta? "
	       << "\n<curveTag? curveType?>";
	return 0;
    }

    int tag;
    int numdata = 1;
    if (OPS_GetIntInput(&numdata,&tag) < 0) {
	opserr << "WARNING invalid tag\n";
	return 0;
    }

    double sp12[4];
    numdata = 4;
    if (OPS_GetDoubleInput(&numdata,sp12) < 0) {
	opserr << "WARNING invalid double inputs\n";
	return 0;
    }

    double sp3[2];
    if (argc > 16) {
	numdata = 2;
	if (OPS_GetDoubleInput(&numdata,sp3) < 0) {
	    opserr << "WARNING invalid double inputs\n";
	    return 0;
	}
    }

    double sn12[4];
    numdata = 4;
    if (OPS_GetDoubleInput(&numdata,sn12) < 0) {
	opserr << "WARNING invalid double inputs\n";
	return 0;
    }

    double sn3[2];
    if (argc > 16) {
	numdata = 2;
	if (OPS_GetDoubleInput(&numdata,sn3) < 0) {
	    opserr << "WARNING invalid double inputs\n";
	    return 0;
	}
    }

    double data[4];
    numdata = 4;
    if (OPS_GetDoubleInput(&numdata,data) < 0) {
	opserr << "WARNING invalid double inputs\n";
	return 0;
    }

    double beta = 0.0;
    numdata = 1;
    if (argc == 20 || argc == 16 || argc >= 22 ) {
	if (OPS_GetDoubleInput(&numdata,&beta) < 0) {
	    opserr << "WARNING invalid beta\n";
	    return 0;
	}
    }

    int degrade = 0;
    if (argc == 22 || argc == 23) {

	double curveData[2];
	numdata = 2;
	if (OPS_GetDoubleInput(&numdata,curveData) < 0) {
	    opserr << "WARNING invalid int inputs\n";
	    return 0;
	}

	//LimitCurve *theCurve = GetLimitCurve(curveTag); //MRL Commented
	LimitCurve *theCurve = 0;//MRL Added
	theCurve = OPS_getLimitCurve(curveData[0]); //MRL Added
    
	if (theCurve == 0) {
	    opserr << "WARNING limit curve does not exist\n";
	    opserr << "limit curve: " << curveData[0]; 
	    opserr << "\nLimitStateMaterial: " << tag << "\n";
	    return 0;
	}
    
	if (argc == 23) {
	    numdata = 1;
	    if (OPS_GetIntInput(&numdata,&degrade) < 0) {
		opserr << "WARNING invalid degrade\n";
		return 0;
	    }
	} 
	mat = new LimitStateMaterial(tag,
				     sp12[0], sp12[1], sp12[2], sp12[3], sp3[0], sp3[1],
				     sn12[0], sn12[1], sn12[2], sn12[3], sn3[0], sn3[1],
				     data[0], data[1], data[2], data[3], beta,
				     *theCurve, curveData[1], degrade);
    }

    // Parsing was successful, allocate the material
    if (argc == 20 || argc == 19) {		
	mat = new LimitStateMaterial(tag,
				     sp12[0], sp12[1], sp12[2], sp12[3], sp3[0], sp3[1],
				     sn12[0], sn12[1], sn12[2], sn12[3], sn3[0], sn3[1],
				     data[0], data[1], data[2], data[3], beta);
	
    } else if (argc == 16 || argc == 15) {
	mat = new LimitStateMaterial(tag,
				     sp12[0], sp12[1], sp12[2], sp12[3],
				     sn12[0], sn12[1], sn12[2], sn12[3],
				     data[0], data[1], data[2], data[3], beta);
    }
    return mat;
}

LimitStateMaterial::LimitStateMaterial(int tag,
			double m1p, double r1p, double m2p, double r2p, double m3p, double r3p,
			double m1n, double r1n, double m2n, double r2n, double m3n, double r3n,
			double px, double py, double d1, double d2, double b):
UniaxialMaterial(tag, MAT_TAG_LimitState),
pinchX(px), pinchY(py), damfc1(d1), damfc2(d2), beta(b),
mom1p(m1p), rot1p(r1p), mom2p(m2p), rot2p(r2p), mom3p(m3p), rot3p(r3p),
mom1n(m1n), rot1n(r1n), mom2n(m2n), rot2n(r2n), mom3n(m3n), rot3n(r3n)
{
	constructorType = 1;

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
//		this->Print(cout);  // Commented out by Terje
		
//		g3ErrorHandler->fatal("%s -- input backbone is not unique (one-to-one) 3-point backbone",
//		"LimitStateMaterial::LimitStateMaterial");
	}

	// Store original input parameters
	pinchX_orig = px;
	pinchY_orig = py;
	damfc1_orig = d1;
	damfc2_orig = d2;
	beta_orig = b;
	mom1p_orig = m1p;
	rot1p_orig = r1p;
	mom2p_orig = m2p;
	rot2p_orig = r2p;
	mom3p_orig = m3p;
	rot3p_orig = r3p;
	mom1n_orig = m1n;
	rot1n_orig = r1n;
	mom2n_orig = m2n;
	rot2n_orig = r2n;
	mom3n_orig = m3n;
	rot3n_orig = r3n;

	energyA = 0.5 * (rot1p*mom1p + (rot2p-rot1p)*(mom2p+mom1p) + (rot3p-rot2p)*(mom3p+mom2p) +
		rot1n*mom1n + (rot2n-rot1n)*(mom2n+mom1n) * (rot3n-rot2n)*(mom3n+mom2n));
	
	// Set envelope slopes
	this->setEnvelope();

	// Initialize history variables
	this->revertToStart();
	this->revertToLastCommit();
	
	/////////////////
	// not using limit curve option
	curveType = 0;
	degrade = 0;
	/////////////////
}

LimitStateMaterial::LimitStateMaterial(int tag,
			double m1p, double r1p, double m2p, double r2p,
			double m1n, double r1n, double m2n, double r2n,
			double px, double py, double d1, double d2, double b):
UniaxialMaterial(tag, MAT_TAG_LimitState),
pinchX(px), pinchY(py), damfc1(d1), damfc2(d2), beta(b),
mom1p(m1p), rot1p(r1p), mom3p(m2p), rot3p(r2p),
mom1n(m1n), rot1n(r1n), mom3n(m2n), rot3n(r2n)
{
	constructorType = 2;
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
//		g3ErrorHandler->fatal("%s -- input backbone is not unique (one-to-one) 2-point backbone",
//		"LimitStateMaterial::LimitStateMaterial");
	}
	
	// Store original input parameters
	pinchX_orig = px;
	pinchY_orig = py;
	damfc1_orig = d1;
	damfc2_orig = d2;
	beta_orig = b;
	mom1p_orig = m1p;
	rot1p_orig = r1p;
	mom2p_orig = m2p;
	rot2p_orig = r2p;
	mom3p_orig = m2p;
	rot3p_orig = r2p;
	mom1n_orig = m1n;
	rot1n_orig = r1n;
	mom2n_orig = m2n;
	rot2n_orig = r2n;
	mom3n_orig = m2n;
	rot3n_orig = r2n;	
	
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

	//this->Print(cout);

	/////////////////
	// not using limit curve option
	curveType = 0;
	degrade = 0;
	/////////////////
}

LimitStateMaterial::LimitStateMaterial(int tag,
			double m1p, double r1p, double m2p, double r2p, double m3p, double r3p,
			double m1n, double r1n, double m2n, double r2n, double m3n, double r3n,
			double px, double py, double d1, double d2, double b, LimitCurve &curve, 
			int cType, int deg):
UniaxialMaterial(tag, MAT_TAG_LimitState),
mom1p(m1p), rot1p(r1p), mom2p(m2p), rot2p(r2p), mom3p(m3p), rot3p(r3p),
mom1n(m1n), rot1n(r1n), mom2n(m2n), rot2n(r2n), mom3n(m3n), rot3n(r3n),
pinchX(px), pinchY(py), damfc1(d1), damfc2(d2), beta(b), theCurve(0), curveType(cType), 
degrade(deg)
{
	constructorType = 3;
	theCurve = curve.getCopy();

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
//		this->Print(cout); Commented out by Terje

//		g3ErrorHandler->fatal("%s -- input backbone is not unique (one-to-one) limit curve option",
//		"LimitStateMaterial::LimitStateMaterial");
	}

	// Store original input parameters
	pinchX_orig = px;
	pinchY_orig = py;
	damfc1_orig = d1;
	damfc2_orig = d2;
	beta_orig = b;
	mom1p_orig = m1p;
	rot1p_orig = r1p;
	mom2p_orig = m2p;
	rot2p_orig = r2p;
	mom3p_orig = m3p;
	rot3p_orig = r3p;
	mom1n_orig = m1n;
	rot1n_orig = r1n;
	mom2n_orig = m2n;
	rot2n_orig = r2n;
	mom3n_orig = m3n;
	rot3n_orig = r3n;


	energyA = 0.5 * (rot1p*mom1p + (rot2p-rot1p)*(mom2p+mom1p) + (rot3p-rot2p)*(mom3p+mom2p) +
		rot1n*mom1n + (rot2n-rot1n)*(mom2n+mom1n) * (rot3n-rot2n)*(mom3n+mom2n));

	
	////////////////////
	// get copy of the limit state curve 
	//theCurve = curve.getCopy(); Terje Xmas

	if (theCurve == 0) {
//		g3ErrorHandler->fatal("WARNING LimitStateMaterial - could not get copy of limitCurve");
	////////////////////
	}
	
	// Set envelope slopes
	this->setEnvelope();

	/////////////////
	// get elastic slope
	Eelasticp = E1p;
	Eelasticn = E1n;
	//////////////////

	// Initialize history variables
	this->revertToStart();
	this->revertToLastCommit();

//	this->Print(cout); // commented out by Terje
}

LimitStateMaterial::LimitStateMaterial():
UniaxialMaterial(0, MAT_TAG_LimitState),
pinchX(0.0), pinchY(0.0), damfc1(0.0), damfc2(0.0), beta(0.0),
mom1p(0.0), rot1p(0.0), mom2p(0.0), rot2p(0.0), mom3p(0.0), rot3p(0.0),
mom1n(0.0), rot1n(0.0), mom2n(0.0), rot2n(0.0), mom3n(0.0), rot3n(0.0)
{
	constructorType = 4;

	// Store original input parameters
	pinchX_orig = 0.0;
	pinchY_orig = 0.0;
	damfc1_orig = 0.0;
	damfc2_orig = 0.0;
	beta_orig = 0.0;
	mom1p_orig = 0.0;
	rot1p_orig = 0.0;
	mom2p_orig = 0.0;
	rot2p_orig = 0.0;
	mom3p_orig = 0.0;
	rot3p_orig = 0.0;
	mom1n_orig = 0.0;
	rot1n_orig = 0.0;
	mom2n_orig = 0.0;
	rot2n_orig = 0.0;
	mom3n_orig = 0.0;
	rot3n_orig = 0.0;

	curveType = 0;
}

LimitStateMaterial::~LimitStateMaterial()
{
  ////////////////////
  if (curveType != 0)
    delete theCurve;
  ////////////////////
}

int
LimitStateMaterial::setTrialStrain(double strain, double strainRate)
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
LimitStateMaterial::getStrain(void)
{
	/////////////////
	// Return trail strain plus strain due to loss of axial load.
	// Ploss will be zero if no axial failure or not using AxialCurve.
	// Ploss is always positive.
	// E3 set to any number if not using limit curve, 
	// otherwise should be negative for axial curve.
	double strain;
	double E3;
	
	if (curveType != 0)
		E3 = theCurve->getDegSlope();
	else 
		E3 = 1.0;

	if (Tstrain < 0.0)
		strain = Tstrain + Ploss/E3;
	else
		strain = Tstrain - Ploss/E3;

 	return strain;
	//////////////////
}

double
LimitStateMaterial::getStress(void)
{
	/////////////////
	// Return trail stress minus the loss of axial load.
	// Ploss will be zero if no axial failure or not using AxialCurve.
	// Ploss is always positive.
	// For axial failure Tstress is negative in compression
	double stress;

	stress = Tstress + Ploss;

	return stress;
	/////////////////
}

double
LimitStateMaterial::getTangent(void)
{
	//////////////////
	// If on the limit state surface use degrading slope, 
	// but if beyond third corner point use approx zero slope 
	// (axial curve only)
	if (curveType == 1) 
	{
		double E3 = theCurve->getDegSlope();

		if (CstateFlag == 1 || CstateFlag == 2) {
			if (Tstrain > 0.0) {
				if (Tstrain > rot3p) {
					Ttangent = E1p*1.0e-9;
				} else {
					Ttangent = E3p;
				} 
			} else {
				if (Tstrain < rot3n) {
					Ttangent = E1p*1.0e-9;
				} else {
					Ttangent = E3n;
				} 
			}			
		}
	}
	///////////////////

	return Ttangent;
}

void
LimitStateMaterial::positiveIncrement(double dStrain)
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

//////////////////////////
	
	if (degrade == 1)
	{
		// limit TrotMax to maximum strain from opposite direction
		if (TrotMax < fabs(CrotMin))
			TrotMax = fabs(CrotMin);
	}
//////////////////////////

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
LimitStateMaterial::negativeIncrement(double dStrain)
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

//////////////////////////
	if (degrade == 1)
	{
		// limit TrotMax to maximum strain from opposite direction
		if (TrotMin > -1.0*(CrotMax))
			TrotMin = -1.0*(CrotMax);
	}
//////////////////////////

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
LimitStateMaterial::commitState(void)
{
	CrotMax = TrotMax;
	CrotMin = TrotMin;
	CrotPu = TrotPu;
	CrotNu = TrotNu;
	CenergyD = TenergyD;
	CloadIndicator = TloadIndicator;

	Cstress = Tstress;
	Cstrain = Tstrain;
	
	int result = 0;
	
	////////////////////
	// check element state if using limit curve option
	// and not beyond residual capacity (CstateFlag == 4)
	if (curveType != 0 && CstateFlag != 4)
	{
		////////////////////
		// Check state of element relative to the limit state surface.
		// Note that steps should be kept small to minimize error
		// caused by committed state being far beyond limit state surface
		int stateFlag = theCurve->checkElementState(Cstress);

		// If beyond limit state surface for first time,
		// get the new final slope and residual capacity 
		// for this LimitState material
		if (stateFlag == 1)
		{
//			g3ErrorHandler->warning("LimitStateMaterial tag %d - failure detected", this->getTag());

			// display warning if Cstrain = max strain experienced.
			if (Cstrain != CrotMax && Cstrain != CrotMin) {
//				g3ErrorHandler->warning("WARNING LimitStateMaterial - failure occurred while not at peak in displacement",
//						"History variables may not be correct. Need to reset history variables in commitState");
			}

			result += getNewBackbone(stateFlag); // get backbone in current direction
			
			// if not an axial curve, cause failure in both directions 
			if (curveType != 1) 
				result += mirrorBackbone(); 

//			this->Print(cout); // Commented out by Terje
		}

		// special functions for axial curve
		if (curveType == 1) {

			// If on surface, get axial load lost  
			if (stateFlag == 1 || stateFlag == 2 || stateFlag == 4) {
				Ploss += theCurve->getUnbalanceForce();
				opserr << "Axial load loss: " << Ploss << endln;
			}

			// if moving off surface then get new backbone with elastic 3rd slope
			if (CstateFlag == 2 || CstateFlag == 1) {
				if (stateFlag == 3) {
					result += getNewBackbone(stateFlag);
//					this->Print(cout); // Commented out by Terje
				}

			}
			
			// if moving onto surface then get new backbone with degrading slope
			if (CstateFlag == 3) {
				if (stateFlag == 2) {
					result += getNewBackbone(stateFlag);
//					this->Print(cout); Commented out by Terje
				}
			}

			// if forceSurface governed by residual capacity set new flat backbone
			// do not allow backbone to be changed again.
			if (stateFlag == 4) {
				result += getNewBackbone(stateFlag);
//				this->Print(cout); commented out by Terje
			}
		}

		// commit the current state if needed outside commitState 
		CstateFlag = stateFlag;

	}	//////////////////

	return 0;
}	

int
LimitStateMaterial::revertToLastCommit(void)
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
LimitStateMaterial::revertToStart(void)
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

	/////////////////
	// set committed stateFlag to prior to failure
	CstateFlag = 0;
	// set axial load loss to zero
	Ploss = 0.0;
	/////////////////

	// Do the same thin inthe limit curve (done by Terje)
	theCurve->revertToStart();


	// Store original input parameters
	pinchX = pinchX_orig;
	pinchY = pinchY_orig;
	damfc1 = damfc1_orig;
	damfc2 = damfc2_orig;
	beta = beta_orig;
	mom1p = mom1p_orig;
	rot1p = rot1p_orig;
	mom2p = mom2p_orig;
	rot2p = rot2p_orig;
	mom3p = mom3p_orig;
	rot3p = rot3p_orig;
	mom1n = mom1n_orig;
	rot1n = rot1n_orig;
	mom2n = mom2n_orig;
	rot2n = rot2n_orig;
	mom3n = mom3n_orig;
	rot3n = rot3n_orig;

	energyA = 0.5 * (rot1p*mom1p + (rot2p-rot1p)*(mom2p+mom1p) + (rot3p-rot2p)*(mom3p+mom2p) +
		rot1n*mom1n + (rot2n-rot1n)*(mom2n+mom1n) * (rot3n-rot2n)*(mom3n+mom2n));
	

	if (constructorType == 2) {
		mom2p = 0.5*(mom1p+mom3p);
		mom2n = 0.5*(mom1n+mom3n);

		rot2p = 0.5*(rot1p+rot3p);
		rot2n = 0.5*(rot1n+rot3n);
	}


	// Set envelope slopes
	this->setEnvelope();


	Eelasticp = E1p;
	Eelasticn = E1n;


	// Initialize history variables
	this->revertToLastCommit();


	return 0;
}

UniaxialMaterial*
LimitStateMaterial::getCopy(void)
{
	//////////////////
	if (curveType == 0) {
		LimitStateMaterial *theCopy = new LimitStateMaterial (this->getTag(),
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
		theCopy->CstateFlag = CstateFlag;

		return theCopy;
	}
	else {
		LimitStateMaterial *theCopy = new LimitStateMaterial (this->getTag(),
			mom1p, rot1p, mom2p, rot2p, mom3p, rot3p,
			mom1n, rot1n, mom2n, rot2n, mom3n, rot3n,
			pinchX, pinchY, damfc1, damfc2, beta, *theCurve, curveType, degrade);

		theCopy->CrotMax = CrotMax;
		theCopy->CrotMin = CrotMin;
		theCopy->CrotPu = CrotPu;
		theCopy->CrotNu = CrotNu;
		theCopy->CenergyD = CenergyD;
		theCopy->CloadIndicator = CloadIndicator;
		theCopy->Cstress = Cstress;
		theCopy->Cstrain = Cstrain;
		theCopy->Ttangent = Ttangent;
		theCopy->CstateFlag = CstateFlag;

		return theCopy;
	}
	///////////////////
}

int
LimitStateMaterial::sendSelf(int commitTag, Channel &theChannel)
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
    opserr << "LimitStateMaterial::sendSelf() - failed to send data\n";


  return res;
}

int
LimitStateMaterial::recvSelf(int commitTag, Channel &theChannel, 
			FEM_ObjectBroker &theBroker)
{
  int res = 0;
  
  static Vector data(27);
  res = theChannel.recvVector(this->getDbTag(), commitTag, data);
  
  if (res < 0) {
      opserr << "LimitStateMaterial::recvSelf() - failed to receive data\n";
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
LimitStateMaterial::Print(OPS_Stream &s, int flag)
{
	s << "LimitState Material, tag: " << this->getTag() << endln;
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
///////////////////
	s << "CstateFlag: " << CstateFlag << endln;
	s << "Cstress: " << Cstress << endln;
	s << "Cstrain: " << Cstrain << endln;
///////////////////
}

void
LimitStateMaterial::setEnvelope(void)
{
	E1p = mom1p/rot1p;
	E2p = (mom2p-mom1p)/(rot2p-rot1p);
	E3p = (mom3p-mom2p)/(rot3p-rot2p);

	E1n = mom1n/rot1n;
	E2n = (mom2n-mom1n)/(rot2n-rot1n);
	E3n = (mom3n-mom2n)/(rot3n-rot2n);
}

double
LimitStateMaterial::posEnvlpStress(double strain)
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
LimitStateMaterial::negEnvlpStress(double strain)
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
LimitStateMaterial::posEnvlpTangent(double strain)
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
LimitStateMaterial::negEnvlpTangent(double strain)
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
LimitStateMaterial::posEnvlpRotlim(double strain)
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
LimitStateMaterial::negEnvlpRotlim(double strain)
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

///////////////////
int
LimitStateMaterial::getNewBackbone(int flag)
{
	double k3 = theCurve->getDegSlope();
	double Fr = theCurve->getResForce();

	// set backbone with flat slope if at residual capacity
	if (flag == 4) {
		if (Cstress > 0.0) {
			mom3p = Cstress;
			rot3p = Cstrain;
			rot2p = (rot3p+rot1p)/2;
			mom2p = mom3p-(rot3p-rot2p)*(E1p*1.0e-9);
		} 
		else {
			mom3n = Cstress;
			rot3n = Cstrain;
			rot2n = (rot3n+rot1n)/2;
			mom2n = mom3n-(rot3n-rot2n)*(E1n*1.0e-9);
		}
	}
	else {
		// determine new corner points for post-failure envelope.
		//
		// set point of failure as second corner point
		if (Cstress > 0.0) {
			mom2p = Cstress;
			rot2p = Cstrain;
		} 
		else {
			mom2n = Cstress;
			rot2n = Cstrain;
		}
		// if not yet beyond first corner point,
		// set new first corner point on elastic slope,
		// otherwise leave first corner point as is.
		if (Cstrain <= rot1p && Cstrain >= rot1n) {
				if (Cstress > 0.0) {
					mom1p = mom2p/2.0;
					rot1p = mom1p/Eelasticp;
				}
				else {
					mom1n = mom2n/2.0;
					rot1n = mom1n/Eelasticn;
				}
		}
		
		// set new third corner point
		if (flag == 3 && curveType == 1) { //if coming off surface
			if (Cstress > 0.0) {
				mom3p = 10*mom2p;
				//rot3p = rot2p + (mom3p-mom2p)/Eelasticp;
				rot3p = rot2p + (mom3p-mom2p)/(Eelasticp*0.01);
			}
			else {
				mom3n = 10*mom2n;
				//rot3n = rot2n + (mom3n-mom2n)/Eelasticn;
				rot3n = rot2n + (mom3n-mom2n)/(Eelasticn*0.01);
			}
		}
		else {									//if hitting surface
			if (Cstress > 0.0) {
				mom3p = Fr;
				rot3p = rot2p + (mom3p-mom2p)/k3;
			}
			else {
				mom3n = -1*Fr;
				rot3n = rot2n + (mom3n-mom2n)/k3;
			}
		}
	}

	// Check backbone
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
	{
//		this->Print(cout); // Commented out by Terje
//		g3ErrorHandler->fatal("%s -- post-failure backbone is not unique (one-to-one)",
//			"LimitStateMaterial::getNewBackbone");
	}

	// recalculate energy for damfc2 parameter
	energyA = 0.5 * (rot1p*mom1p + (rot2p-rot1p)*(mom2p+mom1p) + (rot3p-rot2p)*(mom3p+mom2p) +
		rot1n*mom1n + (rot2n-rot1n)*(mom2n+mom1n) * (rot3n-rot2n)*(mom3n+mom2n));
	
	if (Cstress > 0.0) {
		E1p = mom1p/rot1p;
		E2p = (mom2p-mom1p)/(rot2p-rot1p);
		E3p = (mom3p-mom2p)/(rot3p-rot2p);
	}
	else {
		E1n = mom1n/rot1n;
		E2n = (mom2n-mom1n)/(rot2n-rot1n);
		E3n = (mom3n-mom2n)/(rot3n-rot2n);
	}

	if (Cstress > 0.0)
		return 1;
	else
		return -1;

}

int
LimitStateMaterial::mirrorBackbone(void)
{
	if (Cstress > 0.0) {
		E1n = E1p;
		E2n = E2p;
		E3n = E3p;
		mom1n = -mom1p;
		mom2n = -mom2p;
		mom3n = -mom3p;
		rot1n = -rot1p;
		rot2n = -rot2p;
		rot3n = -rot3p;
	} 
	else {
		E1p = E1n;
		E2p = E2n;
		E3p = E3n;
		mom1p = -mom1n;
		mom2p = -mom2n;
		mom3p = -mom3n;
		rot1p = -rot1n;
		rot2p = -rot2n;
		rot3p = -rot3n;
	}
	return 0;
}
///////////////////




int
LimitStateMaterial::setParameter(const char **argv, int argc, Parameter &param)
{
  return theCurve->setParameter(argv, argc, param);
}


Response* 
LimitStateMaterial::setResponse(const char **argv, int argc,
				OPS_Stream &theOutput)
{
  Response *theResponse = this->UniaxialMaterial::setResponse(argv, argc, theOutput);

  if (theResponse != 0)
    return theResponse;

  // stress
  if (strcmp(argv[0],"stateFlag") == 0) {
    theOutput.tag("UniaxialMaterialOutput");
    theOutput.attr("matType", this->getClassType());
    theOutput.attr("matTag", this->getTag());
    
    theOutput.tag("ResponseType", "stateFlag");
    theResponse =  new MaterialResponse(this, 101, CstateFlag*1.0);
    theOutput.endTag();
  }  

  return theResponse;
}


int 
LimitStateMaterial::getResponse(int responseID, Information &matInfo)
{
  if (responseID == 101) {
      matInfo.setDouble(CstateFlag);
      return 0;
  } else
    return this->UniaxialMaterial::getResponse(responseID, matInfo);
}
