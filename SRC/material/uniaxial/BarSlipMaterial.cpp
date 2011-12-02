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
// $Date: 2004-10-06 19:21:12 $
// $Source: /usr/local/cvs/OpenSees/SRC/material/uniaxial/BarSlipMaterial.cpp,v $
                                                                        
                                                                        
// Written: NM (nmitra@u.washington.edu) 
// Created: January 2002
// Updated: September 2004
//
// Description: This file contains the class implementation for 
// bar-slip material which is defined by 4 points on the positive and 
// negative envelopes and a bunch of damage parameters. The material accounts for
// 3 types of damage rules : Strength degradation, Stiffness degradation, 
// unloading stiffness degradation. 
// The theory for the material has been provided by Prof. Lowes
// Updates: new damage calculations and user interfaces


#include <BarSlipMaterial.h>
#include <math.h>
#include <float.h>
#include <OPS_Globals.h>

BarSlipMaterial::BarSlipMaterial(int tag,
				 double f, double fs, double es, double fsu,
				 double eh, double dbar, double ljoint, int n, double w, double d, 
				 int bsf, int typ):
  UniaxialMaterial(tag, MAT_TAG_BarSlip),
  tagMat(tag), bsflag(bsf), unit(0), type(typ), damage(1), width(w), depth(d),
  envlpPosStress(6), envlpPosStrain(6), envlpNegStress(6), envlpNegStrain(6),
  fc(f), fy(fs), Es(es), fu(fsu), Eh(eh), db(dbar), nbars(n), ld(ljoint), 
  eP(4,2), eN(4,2), envlpPosDamgdStress(6), envlpNegDamgdStress(6), state3Stress(4), 
  state3Strain(4), state4Stress(4), state4Strain(4)
{
  rDispP = 0.25; rForceP = 0.25; uForceP = 0.0;
  rDispN = 0.25; rForceN = 0.25; uForceN = 0.0;
  gammaK1 = 0.3; gammaK2 = 0.0; gammaK3 = 0.1; gammaK4 = 0.0; gammaKLimit = 0.4;
  gammaD1 = 0.6; gammaD2 = 0.0; gammaD3 = 0.2; gammaD4 = 0.0; gammaDLimit = 0.25;
  gammaF1 = 0.7; gammaF2 = 0.3; gammaF3 = 0.5; gammaF4 = 0.1; gammaFLimit = 0.0;
  gammaE = 10.0; 

	getBondStrength();
	getBarSlipEnvelope();
	createMaterial();
}

BarSlipMaterial::BarSlipMaterial(int tag,
				 double f, double fs, double es, double fsu,
				 double eh, double dbar, double ljoint, int n, double w, double d, 
				 int bsf, int typ, int dam, int unt):
  UniaxialMaterial(tag, MAT_TAG_BarSlip),
  tagMat(tag), bsflag(bsf), unit(unt), type(typ), damage(dam), width(w), depth(d),
  envlpPosStress(6), envlpPosStrain(6), envlpNegStress(6), envlpNegStrain(6),
  fc(f), fy(fs), Es(es), fu(fsu), Eh(eh), db(dbar), nbars(n), ld(ljoint), 
  eP(4,2), eN(4,2), envlpPosDamgdStress(6), envlpNegDamgdStress(6), state3Stress(4), 
  state3Strain(4), state4Stress(4), state4Strain(4)
{
  rDispP = 0.25; rForceP = 0.25; uForceP = 0.0;
  rDispN = 0.25; rForceN = 0.25; uForceN = 0.0;
  gammaK1 = 0.3; gammaK2 = 0.0; gammaK3 = 0.1; gammaK4 = 0.0; gammaKLimit = 0.4;
  gammaD1 = 0.6; gammaD2 = 0.0; gammaD3 = 0.2; gammaD4 = 0.0; gammaDLimit = 0.25;
  gammaF1 = 0.7; gammaF2 = 0.3; gammaF3 = 0.5; gammaF4 = 0.1; gammaFLimit = 0.0;
  gammaE = 10.0; 


	if (damage == 0)
	{
		gammaK1 = 0.0; gammaK2 = 0.0; gammaK3 = 0.0; gammaK4 = 0.0; gammaKLimit = 0.0;
		gammaD1 = 0.0; gammaD2 = 0.0; gammaD3 = 0.0; gammaD4 = 0.0; gammaDLimit = 0.0;
		gammaF1 = 0.0; gammaF2 = 0.0; gammaF3 = 0.0; gammaF4 = 0.0; gammaFLimit = 0.0;
	}
	if (damage == 1)
	{
		// activation of linear degradation of strength
		gammaF1 = 0.0; gammaF2 = 0.0; gammaF3 = 0.0; gammaF4 = 0.0; gammaFLimit = 0.0;
	}
	if (damage == 2)
	{
		// activation of logarithmic degradation of strength
		gammaF1 = 11.8986; gammaF2 = 0.0; gammaF3 = 3.9694; gammaF4 = 0.0; gammaFLimit = 0.85;
	}

	getBondStrength();
	getBarSlipEnvelope();
	createMaterial();
}


BarSlipMaterial::BarSlipMaterial():
  UniaxialMaterial(0, MAT_TAG_BarSlip),
  tagMat(0), bsflag(0), unit(0), type(0), damage(0), width(0.0), depth(0.0),
  envlpPosStress(6), envlpPosStrain(6), envlpNegStress(6), envlpNegStrain(6),
  fc(0.0), fy(0.0), Es(0.0), fu(0.0), Eh(0.0), db(0.0), nbars(0), ld(0.0),
  eP(4,2), eN(4,2), envlpPosDamgdStress(6), envlpNegDamgdStress(6), state3Stress(4), 
  state3Strain(4), state4Stress(4), state4Strain(4)
{
	// default constructor --- does nothing
}

BarSlipMaterial::~BarSlipMaterial()
{
	// destructor
//	fn->close();
//	fg->close();

}

void BarSlipMaterial::getBondStrength(void)
{
	if (fc <= 0.0)
		opserr << "WARNING : BAR-SLIP -- fc should be positive entry" << endln;

	if (unit == 0)  // unit has not been specified : either psi or MPa
	{
		if (fc >= 1000.0)
		{
			unit = 2;   // psi unit being used
		}
		else
		{
			unit = 1;   // MPa unit being used
		}
	}

	if (unit == 1)   //  "mpa" (N/sq. mm) system being used
	{
		if (bsflag == 1)  // consider bond strength flag --- weak
		{
			tauYT = 0.05*pow(fc,0.5);
			tauET = 1.8*pow(fc,0.5);
			tauEC = 2.2*pow(fc,0.5);
			tauYC = 3.7*pow(fc,0.5);
			tauR = 0.15*pow(fc,0.5);
		}
		else if (bsflag == 0) // bond strength -- strong
		{
			tauYT = 0.4*pow(fc,0.5);
			tauET = 1.8*pow(fc,0.5);
			tauEC = 2.2*pow(fc,0.5);
			tauYC = 3.7*pow(fc,0.5);
			tauR = 0.15*pow(fc,0.5);
		}
	}
	else if (unit == 2)  // "psi" (lb/sq. inch) system used
	{
		if (bsflag == 1)   // bond strength -- weak
		{
			tauYT = 0.6*pow(fc,0.5);
			tauET = 21*pow(fc,0.5);
			tauEC = 26*pow(fc,0.5);
			tauYC = 43*pow(fc,0.5);
			tauR =  1.8*pow(fc,0.5);
		}
		else if (bsflag == 0)  // bond strength -- strong
		{
			tauYT = 4.8*pow(fc,0.5);
			tauET = 21*pow(fc,0.5);
			tauEC = 26*pow(fc,0.5);
			tauYC = 43*pow(fc,0.5);
			tauR = 1.8*pow(fc,0.5);
		}
	}
	else if (unit == 3)  // "Pa" (N/sq. m) system being used
	{
		if (bsflag == 1)   // bond strength -- weak
		{
			tauYT = 1.58*pow(fc,0.5);
			tauET = 56.92*pow(fc,0.5);
			tauEC = 69.57*pow(fc,0.5);
			tauYC = 117*pow(fc,0.5);
			tauR =  4.74*pow(fc,0.5);
		}
		else if (bsflag == 0)  // bond strength -- strong
		{
			tauYT = 12.65*pow(fc,0.5);
			tauET = 56.92*pow(fc,0.5);
			tauEC = 69.57*pow(fc,0.5);
			tauYC = 117*pow(fc,0.5);
			tauR = 4.74*pow(fc,0.5);
		}
	}
	else if (unit == 4)  // "psf" (lb/sq.ft) system being used
	{
		if (bsflag == 1)   // bond strength -- weak
		{
			tauYT = 7.2*pow(fc,0.5);
			tauET = 252*pow(fc,0.5);
			tauEC = 312*pow(fc,0.5);
			tauYC = 516*pow(fc,0.5);
			tauR =  21.6*pow(fc,0.5);
		}
		else if (bsflag == 0)  // bond strength -- strong
		{
			tauYT = 57.6*pow(fc,0.5);
			tauET = 252*pow(fc,0.5);
			tauEC = 312*pow(fc,0.5);
			tauYC = 516*pow(fc,0.5);
			tauR =  21.6*pow(fc,0.5);
		}
	}
	else if (unit == 5)  // "ksi" (kilo-lb/sq.inch) system being used
	{
		if (bsflag == 1)   // bond strength -- weak
		{
			tauYT = 0.02*pow(fc,0.5);
			tauET = 0.66*pow(fc,0.5);
			tauEC = 0.82*pow(fc,0.5);
			tauYC = 1.36*pow(fc,0.5);
			tauR =  0.06*pow(fc,0.5);
		}
		else if (bsflag == 0)  // bond strength -- strong
		{
			tauYT = 0.15*pow(fc,0.5);
			tauET = 0.66*pow(fc,0.5);
			tauEC = 0.82*pow(fc,0.5);
			tauYC = 1.36*pow(fc,0.5);
			tauR =  0.06*pow(fc,0.5);
		}
	}
	else if (unit == 6)  // "ksf" (kilo-lb/sq.ft) system being used
	{
		if (bsflag == 1)   // bond strength -- weak
		{
			tauYT = 0.24*pow(fc,0.5);
			tauET = 7.92*pow(fc,0.5);
			tauEC = 9.84*pow(fc,0.5);
			tauYC = 16.32*pow(fc,0.5);
			tauR =  0.72*pow(fc,0.5);
		}
		else if (bsflag == 0)  // bond strength -- strong
		{
			tauYT = 1.8*pow(fc,0.5);
			tauET = 7.92*pow(fc,0.5);
			tauEC = 9.84*pow(fc,0.5);
			tauYC = 16.32*pow(fc,0.5);
			tauR =  0.72*pow(fc,0.5);
		}
	}
}

void BarSlipMaterial::getBarSlipEnvelope(void)
{
//fn = new FileStream();
//fn->setFile("BarSlipEnvelope.out", APPEND);
//fg = new FileStream();
//fg->setFile("BarSlipDamage.out", APPEND);

	double delta = 0.0;
	double del_ult = 0.0;
	if (unit == 1)  // MPa
	{ 
		delta = 3.0;
		 del_ult = 10.0;
	}
	else if (unit == 2 || unit == 5) //psi or ksi
	{
		delta = 3.0/25.43;
		 del_ult = 10.0/25.43;
	}
	else if (unit == 3)  // Pa
	{
		delta = 0.003;
		del_ult = 0.01;
	}
	else if (unit == 4 || unit == 6)  // psf or ksf
	{
		delta = 3.0/(25.43*12.0);
		del_ult = 10.0/(25.43*12.0);
	}

	double Ab = PI*pow(db,2)/4;
	double As = Ab*nbars;
	
	double frP = 0.0;
	double frN = 0.0;
	double geo = 0.0;
	double let = 0.0;
	double lec = 0.0;
	double lyc = 0.0;
	double lyt = 0.0;
	double k1 = 0.0;
	double k2 = 0.0;
	double k3 = 0.0;
	double le1 = 0.0;
	double le2 = 0.0;
	double unloadNForce = 0.0;
	double unloadPForce = 0.0;
	eP.Zero();
	eN.Zero();

	frP = tauR*ld*PI*db*As/Ab;
	frN = -tauR*ld*PI*db*As/Ab;

	geo = PI*db/Ab;
	let = fy/(tauET*geo);       // length of elastic part in tension
	lyt = (fu-fy)/(tauYT*geo);  // length of yielded part in tension
	lec = fy/(tauEC*geo);       // length of elastic part in compression
	lyc = (fu-fy)/(tauYC*geo);  //  length of yielded part in compression

	k1 = 2*Es*(tauET/fy)*geo*As;  // initial slope

	eP(0,0) = 0.5*fy*As/k1;
	eP(0,1) = 0.5*fy*As;
	eP(1,0) = fy*As/k1;
	eP(1,1) = fy*As;

	if (let+lyt<ld && bsflag ==0 )
	{
		k2 = (fu-fy)*As/(fy*lyt/Es + 0.5*geo*tauYT*pow(lyt,2)/Eh); // second slope after yield
	}
	else
	{
		le1 = fy/(tauET*geo);
		le2 = fy/(tauYT*geo);
		k2 = (fu-fy)*As/(geo*0.5*tauYT*(pow(le2,2)/Es - pow(le1,2)/Es + pow(lyt,2)/Eh) + fy*lyt/Es);
	}


	eP(2,0) = (delta<(fy*As/k1 + (fu-fy)*As/k2))? delta:(fy*As/k1 + (fu-fy)*As/k2);
	
	if (eP(2,0) == delta)
	{
		eP(2,1) = fy*As + (delta - fy*As/k1)*k2;
	}
	else
	{
		eP(2,1) = fu*As;
	}

	eP(3,0) = del_ult;
	eP(3,1) = eP(2,1) + (eP(2,1) - eP(1,1))*(eP(3,0) - eP(2,0))/(eP(2,0) - eP(1,0));

	gammaFLimit = 1.0 - frP/eP(2,1);

	double dth = depth;  
	double dd = 0.1*depth;  
	double j = 0.0;
	double beta = 0.85;
	double fcc = 0.0;
	
	if (unit == 2)
	{
		fcc = fc;
	}
	else if (unit == 1)
	{
		fcc = 145*fc;
	}
	else if (unit == 3)
	{
		fcc = 0.000145*fc;	
	}
	else if (unit == 4)
	{
		fcc = 0.00694*fc;
	}
	else if (unit == 5)
	{
		fcc = 1000*fc;
	}
	else if (unit == 6)
	{
		fcc = 6.94*fc;
	}

	double b1 = (fcc - 4000)*0.05/1000;
	if (b1 <= 0.0)
	{
		b1 = 0.0;
	}
	if (b1 >= 0.2)
	{
		b1 = 0.2;
	}
	
	double beta1 = 0.85 - b1;

	double cForce = 0.0;

	if ((type == 0) || (type == 1))  // member type -- beam bottom and beam top (or just beam)
	{
		j = 0.85;
	}
	else if (type ==2)  // member type -- column
	{
		j = 0.75;
	}

	cForce = 1 + (beta*fc*dth*width*2*(1-j)/(Es*As*0.003*beta1*(1-dd*beta1/(2*dth*(1-j)))));
	As = As*cForce;
	k1 = 2*Es*(tauEC/fy)*geo*As;


	eN(0,0) = -.5*fy*As/k1;
	eN(0,1) = -.5*fy*As;
	eN(1,0) = -fy*As/k1;
	eN(1,1) = -fy*As;

	if (lec+lyc < ld && bsflag ==0 )
	{
		k2 = (fu-fy)*As/(fy*lyc/Es + 0.5*geo*tauYC*pow(lyc,2)/Eh);
	}
	else
	{
		le1 = fy/(tauEC*geo);
		le2 = fy/(tauYC*geo);
		k2 = (fu-fy)*As/(geo*0.5*tauYC*(pow(le2,2)/Es - pow(le1,2)/Es + pow(lyc,2)/Eh) + fy*lyc/Es);
	}

	eN(2,0) = (delta<(fy*As/k1 + (fu-fy)*As/k2))? -delta:-(fy*As/k1 + (fu-fy)*As/k2);
	
	if (eN(2,0) == -delta)
	{
		eN(2,1) = -fy*As + (-delta + fy*As/k1)*k2;
	}
	else
	{
		eN(2,1) = -fu*As;
	}

	k3 = 1e-3*k1;
	eN(3,0) = -del_ult;
	eN(3,1) = eN(2,1) + k3*(eN(3,0)-eN(2,0));
		
	As = As/cForce;
	
	double tP = (ld<(lec+lyc)) ? ld:(lec+lyc);
	double tN = (ld<(let+lyt)) ? ld:(let+lyt);
	
	unloadPForce = tauR*PI*db*As*tP/Ab;
	unloadNForce = -tauR*PI*db*As*tN/Ab;

	uForceP = unloadPForce/eP(2,1);
	uForceN = unloadNForce/eN(2,1);

	double ratio = 2.0;  //(research underway to induce more pinching in the material)
	rForceP = (uForceP>uForceN) ? ratio*uForceP:ratio*uForceN;
	rForceN = (uForceP>uForceN) ? ratio*uForceP:ratio*uForceN;
//	rForceP = 0.25; rForceN = 0.25;

}

void BarSlipMaterial::createMaterial(void)
{
/*	 material = new Pinching4Material (tagMat, eP(0,1), eP(0,0), eP(1,1), eP(1,0), eP(2,1), eP(2,0), eP(3,1), eP(3,0),
		eN(0,1), eN(0,0), eN(1,1), eN(1,0), eN(2,1), eN(2,0), eN(3,1), eN(3,0), rDispP, rForceP, uForceP, rDispN,
		rForceN, uForceN, gammaK1, gammaK2, gammaK3, gammaK4, gammaKLimit, gammaD1, gammaD2, gammaD3, 
		gammaD4, gammaDLimit, gammaF1, gammaF2, gammaF3, gammaF4, gammaFLimit, gammaE, 0);
*/
	bool error = false;

	// check backbone parameters formed
	if (eP(0,0) <= 0.0 || eP(1,0) <= 0.0 || eP(2,0) <= 0.0 || eP(3,0) <= 0.0)
		error = true;

	if (eN(0,0) >= 0.0 || eN(1,0) >= 0.0 || eN(2,0) >= 0.0 || eN(3,0) >= 0.0)
		error = true;

	if (error) {
		opserr << "Error: -- input backbone not unique, BarSlipMaterial::BarSlipMaterial" << "\a";
	}

	envlpPosStress.Zero(); envlpPosStrain.Zero(); envlpNegStress.Zero(); envlpNegStrain.Zero(); 
	energyCapacity = 0.0; kunload = 0.0; elasticStrainEnergy = 0.0;

	// set envelope slopes
	SetEnvelope();

	envlpPosDamgdStress = envlpPosStress; envlpNegDamgdStress = envlpNegStress;
	state3Stress.Zero(); state3Strain.Zero(); state4Stress.Zero(); state4Strain.Zero();
	// Initialize history variables
	revertToStart();
	revertToLastCommit();


}
//*******************************************************************************
/*
int BarSlipMaterial::setTrialStrain(double strain, double strainrate)
{

	return material->setTrialStrain(strain, strainrate);
}

double BarSlipMaterial::getStrain(void)
{
	return material->getStrain();
}

double BarSlipMaterial::getStress(void)
{
	return material->getStress();
}

double BarSlipMaterial::getTangent(void)
{
	return material->getTangent();
}

double BarSlipMaterial::getInitialTangent(void)
{
	return material->getInitialTangent();
}
int BarSlipMaterial::commitState(void)					   
{
	return material->commitState();
}

int BarSlipMaterial::revertToLastCommit(void)
{
	
	return material->revertToLastCommit();
}

int BarSlipMaterial::revertToStart(void)
{
	return material->revertToStart();
}
*/
UniaxialMaterial* BarSlipMaterial::getCopy(void)
{
	BarSlipMaterial *theCopy = new BarSlipMaterial (this->getTag(), fc, fy, Es,
		fu, Eh, db, ld, nbars, width, depth, bsflag, type, damage, unit);
	
	theCopy->rDispN = rDispN;
	theCopy->rDispP = rDispP;
	theCopy->rForceN = rForceN;
	theCopy->rForceP = rForceP;
	theCopy->uForceN = uForceN;
	theCopy->uForceP = uForceP;

	// Trial state variables
	theCopy->Tstress = Tstress;
	theCopy->Tstrain = Tstrain;
	theCopy->Ttangent = Ttangent;

	// Coverged material history parameters
	theCopy->Cstate = Cstate;
	theCopy->Cstrain = Cstrain;
	theCopy->Cstress = Cstress;
	theCopy->CstrainRate = CstrainRate;

	theCopy->lowCstateStrain = lowCstateStrain;
	theCopy->lowCstateStress = lowCstateStress;
	theCopy->hghCstateStrain = hghCstateStrain;
	theCopy->hghCstateStress = hghCstateStress;
	theCopy->CminStrainDmnd = CminStrainDmnd;
	theCopy->CmaxStrainDmnd = CmaxStrainDmnd;
	theCopy->Cenergy = Cenergy;
	theCopy->CgammaK = CgammaK;
	theCopy->CgammaD = CgammaD;
	theCopy->CgammaF = CgammaF;
	theCopy->gammaKUsed = gammaKUsed;
	theCopy->gammaFUsed = gammaFUsed;

	// trial material history parameters
	theCopy->Tstate = Tstate;
	theCopy->dstrain = dstrain;
	theCopy->lowTstateStrain = lowTstateStrain;
	theCopy->lowTstateStress = lowTstateStress;
	theCopy->hghTstateStrain = hghTstateStrain;
	theCopy->hghTstateStress = hghTstateStress;
	theCopy->TminStrainDmnd = TminStrainDmnd;
	theCopy->TmaxStrainDmnd = TmaxStrainDmnd;
	theCopy->Tenergy = Tenergy;
	theCopy->TgammaK = TgammaK;
	theCopy->TgammaD = TgammaD;
	theCopy->TgammaF = TgammaF;

	// Strength and stiffness parameters
	theCopy->kElasticPos = kElasticPos;
	theCopy->kElasticNeg = kElasticNeg;
	theCopy->kElasticPosDamgd = kElasticPosDamgd;
	theCopy->kElasticNegDamgd = kElasticNegDamgd;
	theCopy->uMaxDamgd = uMaxDamgd;
	theCopy->uMinDamgd = uMinDamgd;

	for (int i = 0; i<6; i++)
	{
		theCopy->envlpPosStrain(i) = envlpPosStrain(i);
		theCopy->envlpPosStress(i) = envlpPosStress(i);
		theCopy->envlpNegStrain(i) = envlpNegStrain(i);
		theCopy->envlpNegStress(i) = envlpNegStress(i);
		theCopy->envlpNegDamgdStress(i) = envlpNegDamgdStress(i);
		theCopy->envlpPosDamgdStress(i) = envlpPosDamgdStress(i);
	}

	for (int j = 0; j<4; j++)
	{
		theCopy->state3Strain(j) = state3Strain(j);
		theCopy->state3Stress(j) = state3Stress(j);
		theCopy->state4Strain(j) = state4Strain(j);
		theCopy->state4Stress(j) = state4Stress(j);
	}

	theCopy->energyCapacity = energyCapacity;
	theCopy->kunload = kunload;
	theCopy->elasticStrainEnergy = elasticStrainEnergy;
	
	return theCopy;
}

int BarSlipMaterial::sendSelf(int commitTag, Channel &theChannel)
{
	return -1;
}

int BarSlipMaterial::recvSelf(int commitTag, Channel &theChannel,
							   FEM_ObjectBroker & theBroker)
{
	return -1;
}

void BarSlipMaterial::Print(OPS_Stream &s, int flag)		  
{
	s << "BarSlipMaterial, tag: " << this-> getTag() << endln;
	s << "Positive Envelope : " << eP << endln;
	s << "Negative Envelope : " << eN << endln;
}

//**********************************************************************************
void BarSlipMaterial::SetEnvelope(void)
{
	double kPos = eP(0,1)/eP(0,0);
	double kNeg = eN(0,1)/eN(0,0);

	double k = (kPos>kNeg) ? kPos:kNeg;
	double u = (eP(0,0)>-eN(0,0)) ? 1e-4*eP(0,0):-1e-4*eN(0,0);

	envlpPosStrain(0) = u; envlpPosStress(0) = u*k;
	envlpNegStrain(0) = -u; envlpNegStress(0) = -u*k;

	for (int i1 = 1; i1 < 5; i1++) {
		envlpPosStrain(i1) = eP(i1-1,0); 
		envlpPosStress(i1) = eP(i1-1,1);
		envlpNegStrain(i1) = eN(i1-1,0);
		envlpNegStress(i1) = eN(i1-1,1);
	}

	double k1 = (eP(3,1) - eP(2,1))/(eP(3,0) - eP(2,0));
	double k2 = (eN(3,1) - eN(2,1))/(eN(3,0) - eN(2,0));

	envlpPosStrain(5) = 1e+6*eP(3,0); envlpNegStrain(5) = 1e+6*eN(3,0);
	envlpPosStress(5) = (k1>0) ? eP(3,1) + k1*(envlpPosStrain(5)-envlpPosStrain(4)):envlpPosStress(4)*1.1;
	envlpNegStress(5) = (k2>0) ? eN(3,1) + k2*(envlpNegStrain(5)-envlpNegStrain(4)):envlpNegStress(4)*1.1;

	kElasticPos = envlpPosStress(1)/envlpPosStrain(1);
	kElasticNeg = envlpNegStress(1)/envlpNegStrain(1);

	double energyPos = 0.5*envlpPosStrain(0)*envlpPosStress(0);
	double energyNeg = 0.5*envlpNegStrain(0)*envlpNegStress(0);

	for (int i2 = 0; i2<4; i2++) {
		energyPos += 0.5*(envlpPosStress(i2)+envlpPosStress(i2+1))*(envlpPosStrain(i2+1) - envlpPosStrain(i2));
	}

	for (int i3 = 0; i3<4; i3++) {
		energyNeg += 0.5*(envlpNegStress(i3)+envlpNegStress(i3+1))*(envlpNegStrain(i3+1) - envlpNegStrain(i3));
	}

	double max_energy = (energyPos>energyNeg) ? energyPos:energyNeg;
	energyCapacity = gammaE*max_energy;
}

int BarSlipMaterial::setTrialStrain(double strain, double CstrainRate)
{
	Tstate = Cstate;
	Tenergy = Cenergy;
	Tstrain = strain;
	lowTstateStrain = lowCstateStrain;
	hghTstateStrain = hghCstateStrain;
	lowTstateStress = lowCstateStress;
	hghTstateStress = hghCstateStress;
	TminStrainDmnd = CminStrainDmnd;
	TmaxStrainDmnd = CmaxStrainDmnd;
	TgammaF = CgammaF;
	TgammaK = CgammaK; 
	TgammaD = CgammaD;

	dstrain = Tstrain - Cstrain;
	if (dstrain<1e-12 && dstrain>-1e-12){
		dstrain = 0.0;
	}


	// determine new state if there is a change in state
	getstate(Tstrain,dstrain);  

	switch (Tstate)
	{ 

	case 0:
		Ttangent = envlpPosStress(0)/envlpPosStrain(0);
		Tstress = Ttangent*Tstrain;
		break;
	case 1:
		Tstress = posEnvlpStress(strain);
		Ttangent = posEnvlpTangent(strain);
		break;
	case 2:
		Ttangent = negEnvlpTangent(strain);
		Tstress = negEnvlpStress(strain);
		break;
	case 3:
		kunload = (hghTstateStrain<0.0) ? kElasticNegDamgd:kElasticPosDamgd; 	
			state3Strain(0) = lowTstateStrain;
			state3Strain(3) = hghTstateStrain;
			state3Stress(0) = lowTstateStress;
			state3Stress(3) = hghTstateStress;

		getState3(state3Strain,state3Stress,kunload);
		Ttangent = Envlp3Tangent(state3Strain,state3Stress,strain);
		Tstress = Envlp3Stress(state3Strain,state3Stress,strain);
		
		//Print(opserr,1);
		break;
	case 4:
		kunload = (lowTstateStrain<0.0) ? kElasticNegDamgd:kElasticPosDamgd;
			state4Strain(0) = lowTstateStrain;
			state4Strain(3) = hghTstateStrain;
			state4Stress(0) = lowTstateStress;
			state4Stress(3) = hghTstateStress;

		getState4(state4Strain,state4Stress,kunload);
		Ttangent = Envlp4Tangent(state4Strain,state4Stress,strain);
		Tstress = Envlp4Stress(state4Strain,state4Stress,strain);
		break;
	}

	double denergy = 0.5*(Tstress+Cstress)*dstrain;
	elasticStrainEnergy = (Tstrain>0.0) ? 0.5*Tstress/kElasticPosDamgd*Tstress:0.5*Tstress/kElasticNegDamgd*Tstress;

	Tenergy = Cenergy + denergy;

	updateDmg(Tstrain);
	return 0;
}

double BarSlipMaterial::getStrain(void)
{
	return Tstrain;
}

double BarSlipMaterial::getStress(void)
{
	return Tstress;
}

double BarSlipMaterial::getTangent(void)
{
	return Ttangent;
}

double BarSlipMaterial::getInitialTangent(void)
{
	return envlpPosStress(0)/envlpPosStrain(0);
}

int BarSlipMaterial::commitState(void)					   
{
	Cstate = Tstate;

	if (dstrain>1e-12||dstrain<-(1e-12)) {
		CstrainRate = dstrain;}
	else {
		CstrainRate = TstrainRate;}

	lowCstateStrain = lowTstateStrain;
	lowCstateStress = lowTstateStress;
	hghCstateStrain = hghTstateStrain;
	hghCstateStress = hghTstateStress;
	CminStrainDmnd = TminStrainDmnd;
	CmaxStrainDmnd = TmaxStrainDmnd;
	Cenergy = Tenergy;

	Cstress = Tstress;
	Cstrain = Tstrain;

	CgammaK = TgammaK;
	CgammaD = TgammaD;
	CgammaF = TgammaF;
	
	// define adjusted strength and stiffness parameters
	kElasticPosDamgd = kElasticPos*(1 - gammaKUsed);
	kElasticNegDamgd = kElasticNeg*(1 - gammaKUsed);

	uMaxDamgd = TmaxStrainDmnd*(1 + CgammaD);   
	uMinDamgd = TminStrainDmnd*(1 + CgammaD);

	envlpPosDamgdStress = envlpPosStress*(1-gammaFUsed);
	envlpNegDamgdStress = envlpNegStress*(1-gammaFUsed);

//(*fn) << tagMat << "  " << envlpNegStrain(4) << "  " << envlpNegStrain(3) << "  " << envlpNegStrain(2) << "  " << envlpNegStrain(1) << "  " << envlpPosStrain(1) << "  " << envlpPosStrain(2) << "  " << envlpPosStrain(3) << "  " << envlpPosStrain(4) << "  " << envlpNegDamgdStress(4) << "  " << envlpNegDamgdStress(3) << "  " << envlpNegDamgdStress(2) << "  " << envlpNegDamgdStress(1) << "  " << envlpPosDamgdStress(1) << "  " << envlpPosDamgdStress(2) << "  " << envlpPosDamgdStress(3) << "  " << envlpPosDamgdStress(4) <<  endln;
//(*fg) << tagMat << "  " << gammaF1 << "  " << gammaF2 << " " << gammaF3 << " " << gammaF4 << " " << gammaK1 << "  " << gammaK2 << " " << gammaK3 << " " << gammaK4 << " "<< gammaD1 << "  " << gammaD2 << " " << gammaD3 << " " << gammaD4 << endln;
//(*fg) << tagMat << "  " << TgammaF << "  " << TgammaD << " " << TgammaK << " " << CgammaF << " " << CgammaD << "  " << CgammaK << endln;
	return 0;
}

int BarSlipMaterial::revertToLastCommit(void)
{
	Tstate = Cstate;

	TstrainRate = CstrainRate;

	lowTstateStrain = lowCstateStrain;
	lowTstateStress = lowCstateStress;
	hghTstateStrain = hghCstateStrain;
	hghTstateStress = hghCstateStress;
	TminStrainDmnd = CminStrainDmnd;
	TmaxStrainDmnd = CmaxStrainDmnd;
	Tenergy = Cenergy;

	Tstrain = Cstrain; Tstress = Cstress;

	TgammaD = CgammaD;
	TgammaK = CgammaK;
	TgammaF = CgammaF;

	return 0;
}

int BarSlipMaterial::revertToStart(void)
{
    Cstate = 0;
	Cstrain = 0.0;
	Cstress = 0.0;
	CstrainRate = 0.0;
	lowCstateStrain = envlpNegStrain(0);
	lowCstateStress = envlpNegStress(0);
	hghCstateStrain = envlpPosStrain(0);
	hghCstateStress = envlpPosStress(0);
	CminStrainDmnd = envlpNegStrain(1);
	CmaxStrainDmnd = envlpPosStrain(1);
	Cenergy = 0.0;
	CgammaK = 0.0;
	CgammaD = 0.0;
	CgammaF = 0.0;

	Ttangent = envlpPosStress(0)/envlpPosStrain(0);
	dstrain = 0.0;       
	gammaKUsed = 0.0;
	gammaFUsed = 0.0;
	
	kElasticPosDamgd = kElasticPos;
	kElasticNegDamgd = kElasticNeg;
	uMaxDamgd = CmaxStrainDmnd;
	uMinDamgd = CminStrainDmnd;

	return 0;
}

void BarSlipMaterial::getstate(double u,double du)
{
	int cid = 0;
	int cis = 0;
	int newState = 0;
	if (du*CstrainRate<=0.0){   
		cid = 1;
	}
	if (u<lowTstateStrain || u>hghTstateStrain || cid) {                
		if (Tstate == 0) {                                              
			if (u>hghTstateStrain) {
				cis = 1;
				newState = 1;
				lowTstateStrain = envlpPosStrain(0);
				lowTstateStress = envlpPosStress(0);
				hghTstateStrain = envlpPosStrain(5);
				hghTstateStress = envlpPosStress(5);
			}
			else if (u<lowTstateStrain){
				cis = 1;
				newState = 2;
				lowTstateStrain = envlpNegStrain(5);
				lowTstateStress = envlpNegStress(5);
				hghTstateStrain = envlpNegStrain(0);
				hghTstateStress = envlpNegStress(0);
			}
		}
		else if (Tstate==1 && du<0.0) {
			cis = 1;
			if (Cstrain>TmaxStrainDmnd) {
				TmaxStrainDmnd = u - du;
			}
			if (TmaxStrainDmnd<uMaxDamgd) {
				TmaxStrainDmnd = uMaxDamgd;
			}
			if (u<uMinDamgd) {
				newState = 2;
				gammaFUsed = CgammaF;     
				for (int i=0; i<=5; i++) {
					envlpNegDamgdStress(i) = envlpNegStress(i)*(1.0-gammaFUsed);
				}
				lowTstateStrain = envlpNegStrain(5);
				lowTstateStress = envlpNegStress(5);
				hghTstateStrain = envlpNegStrain(0);
				hghTstateStress = envlpNegStress(0);
			}
			else {
				newState = 3;
				lowTstateStrain = uMinDamgd;
				gammaFUsed = CgammaF;        
				for (int i=0; i<=5; i++) {
					envlpNegDamgdStress(i) = envlpNegStress(i)*(1.0-gammaFUsed);
				}
				lowTstateStress = negEnvlpStress(uMinDamgd); 
				hghTstateStrain = Cstrain;
				hghTstateStress = Cstress;
			}
			gammaKUsed = CgammaK;
			kElasticPosDamgd = kElasticPos*(1.0-gammaKUsed);
		}
		else if (Tstate ==2 && du>0.0){
			cis = 1;
			if (Cstrain<TminStrainDmnd) {
				TminStrainDmnd = Cstrain;
			}
			if (TminStrainDmnd>uMinDamgd) {
				TminStrainDmnd = uMinDamgd;
			}
			if (u>uMaxDamgd) {
				newState = 1;
				gammaFUsed = CgammaF;      
				for (int i=0; i<=5; i++) {
					envlpPosDamgdStress(i) = envlpPosStress(i)*(1.0-gammaFUsed);
				}
				lowTstateStrain = envlpPosStrain(0);
				lowTstateStress = envlpPosStress(0);
				hghTstateStrain = envlpPosStrain(5);
				hghTstateStress = envlpPosStress(5);
			}
			else {
				newState = 4;
				lowTstateStrain = Cstrain;
				lowTstateStress = Cstress;
				hghTstateStrain = uMaxDamgd;
				gammaFUsed = CgammaF;         
				for (int i=0; i<=5; i++) {
					envlpPosDamgdStress(i) = envlpPosStress(i)*(1.0-gammaFUsed);
				}
				hghTstateStress = posEnvlpStress(uMaxDamgd);
			}
			gammaKUsed = CgammaK;
			kElasticNegDamgd = kElasticNeg*(1.0-gammaKUsed);
		}
			else if (Tstate ==3) {
				if (u<lowTstateStrain){
					cis = 1;
					newState = 2;
					lowTstateStrain = envlpNegStrain(5);
					hghTstateStrain = envlpNegStrain(0);
					lowTstateStress = envlpNegDamgdStress(5);
					hghTstateStress = envlpNegDamgdStress(0);
				}
				else if (u>uMaxDamgd && du>0.0) {
					cis = 1;
					newState = 1;
					lowTstateStrain = envlpPosStrain(0);
					lowTstateStress = envlpPosStress(0);
					hghTstateStrain = envlpPosStrain(5);
					hghTstateStress = envlpPosStress(5);
				}
				else if (du>0.0) {
					cis = 1;
					newState = 4;
					lowTstateStrain = Cstrain;
					lowTstateStress = Cstress;
					hghTstateStrain = uMaxDamgd;
					gammaFUsed = CgammaF;
					for (int i=0; i<=5; i++) {
					envlpPosDamgdStress(i) = envlpPosStress(i)*(1.0-gammaFUsed);
				    }
					hghTstateStress = posEnvlpStress(uMaxDamgd);
					gammaKUsed = CgammaK;
					kElasticNegDamgd = kElasticNeg*(1.0-gammaKUsed); 
				}
			}
			else if (Tstate == 4){
				if (u>hghTstateStrain){
					cis = 1;
					newState = 1;
					lowTstateStrain = envlpPosStrain(0);
					lowTstateStress = envlpPosDamgdStress(0);
					hghTstateStrain = envlpPosStrain(5);
					hghTstateStress = envlpPosDamgdStress(5);
				}
				else if (u<uMinDamgd && du <0.0) {
					cis = 1;
					newState = 2;
					lowTstateStrain = envlpNegStrain(5);
					lowTstateStress = envlpNegDamgdStress(5);
					hghTstateStrain = envlpNegStrain(0);
					hghTstateStress = envlpNegDamgdStress(0);
				}
				else if (du<0.0) { 
					cis = 1;
					newState = 3;
					lowTstateStrain = uMinDamgd;
					gammaFUsed = CgammaF;         
					for (int i=0; i<=5; i++) {
					envlpNegDamgdStress(i) = envlpNegStress(i)*(1.0-gammaFUsed);
				     }
					lowTstateStress = negEnvlpStress(uMinDamgd);
					hghTstateStrain = Cstrain;
					hghTstateStress = Cstress;
					gammaKUsed = CgammaK;
					kElasticPosDamgd = kElasticPos*(1.0-gammaKUsed);
				}
			}
			}
			if (cis) {
				Tstate = newState;
				
			}
}

double BarSlipMaterial::posEnvlpStress(double u)
{
			double k = 0.0;
			int i = 0;
			double f = 0.0;
			while (k==0.0 && i<=4){   
				 
				 if (u<=envlpPosStrain(i+1)){
					 k = (envlpPosDamgdStress(i+1)-envlpPosDamgdStress(i))/(envlpPosStrain(i+1)-envlpPosStrain(i));
					 f = envlpPosDamgdStress(i) + (u-envlpPosStrain(i))*k;
				 }
				 i++;
			}
            

			if (k==0.0){
				k = (envlpPosDamgdStress(5) - envlpPosDamgdStress(4))/(envlpPosStrain(5) - envlpPosStrain(4));
				f = envlpPosDamgdStress(5) + k*(u-envlpPosStrain(5));
			}

			return f;

}

double BarSlipMaterial::posEnvlpTangent(double u)
{
			double k = 0.0;
			int i = 0;
			while (k==0.0 && i<=4){        
				 
				 if (u<=envlpPosStrain(i+1)){
					 k = (envlpPosDamgdStress(i+1)-envlpPosDamgdStress(i))/(envlpPosStrain(i+1)-envlpPosStrain(i));
			}
				 i++;
			}

			if (k==0.0){
				k = (envlpPosDamgdStress(5) - envlpPosDamgdStress(4))/(envlpPosStrain(5) - envlpPosStrain(4));
		   }

			return k;

}

double BarSlipMaterial::negEnvlpStress(double u)
{
			double k = 0.0;
			int i = 0;
			double f = 0.0;
			while (k==0.0 && i<=4){      				 
				 if (u>=envlpNegStrain(i+1)){
					 k = (envlpNegDamgdStress(i)-envlpNegDamgdStress(i+1))/(envlpNegStrain(i)-envlpNegStrain(i+1));
					 f = envlpNegDamgdStress(i+1)+(u-envlpNegStrain(i+1))*k;
				 }
				 i++;
			}

			if (k==0.0){
				k = (envlpNegDamgdStress(4) - envlpNegDamgdStress(5))/(envlpNegStrain(4)-envlpNegStrain(5));
				f = envlpNegDamgdStress(5) + k*(u-envlpNegStrain(5));
			}
			return f;
}

double BarSlipMaterial::negEnvlpTangent(double u)
{
			double k = 0.0;
			int i = 0;
			while (k==0.0 && i<=4){              
				 
				 if (u>=envlpNegStrain(i+1)){
					 k = (envlpNegDamgdStress(i)-envlpNegDamgdStress(i+1))/(envlpNegStrain(i)-envlpNegStrain(i+1));
					}
				 i++;
			}

			if (k==0.0){
				k = (envlpNegDamgdStress(4) - envlpNegDamgdStress(5))/(envlpNegStrain(4)-envlpNegStrain(5));
				}
			return k;

}

void BarSlipMaterial::getState3(Vector& state3Strain, Vector& state3Stress, double kunload)
		{

			double kmax = (kunload>kElasticNegDamgd) ? kunload:kElasticNegDamgd;

			if (state3Strain(0)*state3Strain(3) <0.0){
				// trilinear unload reload path expected, first define point for reloading
				state3Strain(1) = lowTstateStrain*rDispN;
				if (rForceN-uForceN > 1e-8) {
					state3Stress(1) = lowTstateStress*rForceN;
				}
				else {
					if (TminStrainDmnd < envlpNegStrain(3)) {
						double st1 = lowTstateStress*uForceN*(1.0+1e-6);
						double st2 = envlpNegDamgdStress(4)*(1.0+1e-6);
						state3Stress(1) = (st1<st2) ? st1:st2;
					}
					else {
						double st1 = envlpNegDamgdStress(3)*uForceN*(1.0+1e-6);
						double st2 = envlpNegDamgdStress(4)*(1.0+1e-6);
						state3Stress(1) = (st1<st2) ? st1:st2;
					}
				}
				// if reload stiffness exceeds unload stiffness, reduce reload stiffness to make it equal to unload stiffness
				if ((state3Stress(1)-state3Stress(0))/(state3Strain(1)-state3Strain(0)) > kElasticNegDamgd) {
					state3Strain(1) = lowTstateStrain + (state3Stress(1)-state3Stress(0))/kElasticNegDamgd;
				}
				// check that reloading point is not behind point 4 
				if (state3Strain(1)>state3Strain(3)) {
					// path taken to be a straight line between points 1 and 4
					double du = state3Strain(3) - state3Strain(0);
					double df = state3Stress(3) - state3Stress(0);
					state3Strain(1) = state3Strain(0) + 0.33*du;
					state3Strain(2) = state3Strain(0) + 0.67*du;
					state3Stress(1) = state3Stress(0) + 0.33*df;
					state3Stress(2) = state3Stress(0) + 0.67*df;
				}
				else {
					if (TminStrainDmnd < envlpNegStrain(3)) {
						state3Stress(2) = uForceN*envlpNegDamgdStress(4);
					}
					else {
						state3Stress(2) = uForceN*envlpNegDamgdStress(3);
					}
					state3Strain(2) = hghTstateStrain - (hghTstateStress-state3Stress(2))/kunload;

					if (state3Strain(2) > state3Strain(3)) {
						// point 3 should be along a line between 2 and 4
						double du = state3Strain(3) - state3Strain(1);
						double df = state3Stress(3) - state3Stress(1);
						state3Strain(2) = state3Strain(1) + 0.5*du;
						state3Stress(2) = state3Stress(1) + 0.5*df;
					}
					else if ((state3Stress(2) - state3Stress(1))/(state3Strain(2) - state3Strain(1)) > kmax) {
						// linear unload-reload path expected
						double du = state3Strain(3) - state3Strain(0);
						double df = state3Stress(3) - state3Stress(0);
						state3Strain(1) = state3Strain(0) + 0.33*du;
						state3Strain(2) = state3Strain(0) + 0.67*du;
						state3Stress(1) = state3Stress(0) + 0.33*df;
						state3Stress(2) = state3Stress(0) + 0.67*df;
					}
					else if ((state3Strain(2) < state3Strain(1))||((state3Stress(2)-state3Stress(1))/(state3Strain(2)-state3Strain(1))<0)) {
						if (state3Strain(2)<0.0) {
							// pt 3 should be along a line between 2 and 4
							double du = state3Strain(3)-state3Strain(1);
							double df = state3Stress(3)-state3Stress(1);
							state3Strain(2) = state3Strain(1) + 0.5*du;
							state3Stress(2) = state3Stress(1) + 0.5*df;
						}
						else if (state3Strain(1) > 0.0) {
							// pt 2 should be along a line between 1 and 3
							double du = state3Strain(2)-state3Strain(0);
							double df = state3Stress(2)-state3Stress(0);
							state3Strain(1) = state3Strain(0) + 0.5*du;
							state3Stress(1) = state3Stress(0) + 0.5*df;
						}
						else {
							double avgforce = 0.5*(state3Stress(2) + state3Stress(1));
							double dfr = 0.0;
							if (avgforce < 0.0){
								dfr = -avgforce/100;
							}
							else {
								dfr = avgforce/100;
							}
							double slope12 = (state3Stress(1) - state3Stress(0))/(state3Strain(1) - state3Strain(0));
							double slope34 = (state3Stress(3) - state3Stress(2))/(state3Strain(3) - state3Strain(2));
							state3Stress(1) = avgforce - dfr;
							state3Stress(2) = avgforce + dfr;
							state3Strain(1) = state3Strain(0) + (state3Stress(1) - state3Stress(0))/slope12;
							state3Strain(2) = state3Strain(3) - (state3Stress(3) - state3Stress(2))/slope34;
						}
					}
				}
			}
				else {
					// linear unload reload path is expected		 
					double du = state3Strain(3)-state3Strain(0);
					double df = state3Stress(3)-state3Stress(0);
					state3Strain(1) = state3Strain(0) + 0.33*du;
					state3Strain(2) = state3Strain(0) + 0.67*du;
					state3Stress(1) = state3Stress(0) + 0.33*df;
					state3Stress(2) = state3Stress(0) + 0.67*df;
				}
			
				
				double checkSlope = state3Stress(0)/state3Strain(0);
				double slope = 0.0;

				// final check
				int i = 0;
				while (i<3) {
					double du = state3Strain(i+1)-state3Strain(i);
					double df = state3Stress(i+1)-state3Stress(i);
					if (du<0.0 || df<0.0) {
						double du = state3Strain(3)-state3Strain(0);
						double df = state3Stress(3)-state3Stress(0);
						state3Strain(1) = state3Strain(0) + 0.33*du;
						state3Strain(2) = state3Strain(0) + 0.67*du;
						state3Stress(1) = state3Stress(0) + 0.33*df;
						state3Stress(2) = state3Stress(0) + 0.67*df;
						slope = df/du;
						i = 3;
					}
					if (slope > 1e-8 && slope < checkSlope) {
						state3Strain(1) = 0.0; state3Stress(1) = 0.0;
						state3Strain(2) = state3Strain(3)/2; state3Stress(2) = state3Stress(3)/2;
					} 
					i++;
				}
		// added jun 9th 2004 to prevent slope from getting negative
		if (state3Stress(2) <= state3Stress(1)) {
			state3Stress(1) = state3Stress(2)*1.02;
		}
}

void BarSlipMaterial::getState4(Vector& state4Strain,Vector& state4Stress, double kunload)
		{

			double kmax = (kunload>kElasticPosDamgd) ? kunload:kElasticPosDamgd;

			if (state4Strain(0)*state4Strain(3) <0.0){
				// trilinear unload reload path expected
				state4Strain(2) = hghTstateStrain*rDispP;
				if (uForceP==0.0){
					state4Stress(2) = hghTstateStress*rForceP;
				}
				else if (rForceP-uForceP > 1e-8) {
					state4Stress(2) = hghTstateStress*rForceP;
				}
				else {
					if (TmaxStrainDmnd > envlpPosStrain(3)) {
						double st1 = hghTstateStress*uForceP*(1.0+1e-6);
						double st2 = envlpPosDamgdStress(4)*(1.0+1e-6);
						state4Stress(2) = (st1>st2) ? st1:st2;
					}
					else {
						double st1 = envlpPosDamgdStress(3)*uForceP*(1.0+1e-6);
						double st2 = envlpPosDamgdStress(4)*(1.0+1e-6);
						state4Stress(2) = (st1>st2) ? st1:st2;
					}
				}
				// if reload stiffness exceeds unload stiffness, reduce reload stiffness to make it equal to unload stiffness
				if ((state4Stress(3)-state4Stress(2))/(state4Strain(3)-state4Strain(2)) > kElasticPosDamgd) {
					state4Strain(2) = hghTstateStrain - (state4Stress(3)-state4Stress(2))/kElasticPosDamgd;
				}
				// check that reloading point is not behind point 1 
				if (state4Strain(2)<state4Strain(0)) {
					// path taken to be a straight line between points 1 and 4
					double du = state4Strain(3) - state4Strain(0);
					double df = state4Stress(3) - state4Stress(0);
					state4Strain(1) = state4Strain(0) + 0.33*du;
					state4Strain(2) = state4Strain(0) + 0.67*du;
					state4Stress(1) = state4Stress(0) + 0.33*df;
					state4Stress(2) = state4Stress(0) + 0.67*df;
				}
				else {
					if (TmaxStrainDmnd > envlpPosStrain(3)) {
						state4Stress(1) = uForceP*envlpPosDamgdStress(4);
					}
					else {
						state4Stress(1) = uForceP*envlpPosDamgdStress(3);
					}
					state4Strain(1) = lowTstateStrain + (-lowTstateStress+state4Stress(1))/kunload;

					if (state4Strain(1) < state4Strain(0)) {
						// point 2 should be along a line between 1 and 3
						double du = state4Strain(2) - state4Strain(0);
						double df = state4Stress(2) - state4Stress(0);
						state4Strain(1) = state4Strain(0) + 0.5*du;
						state4Stress(1) = state4Stress(0) + 0.5*df;
					}
					else if ((state4Stress(2) - state4Stress(1))/(state4Strain(2) - state4Strain(1)) > kmax) {
						// linear unload-reload path expected
						double du = state4Strain(3) - state4Strain(0);
						double df = state4Stress(3) - state4Stress(0);
						state4Strain(1) = state4Strain(0) + 0.33*du;
						state4Strain(2) = state4Strain(0) + 0.67*du;
						state4Stress(1) = state4Stress(0) + 0.33*df;
						state4Stress(2) = state4Stress(0) + 0.67*df;
					}
					else if ((state4Strain(2) < state4Strain(1))||((state4Stress(2)-state4Stress(1))/(state4Strain(2)-state4Strain(1))<0)) {
						if (state4Strain(1)>0.0) {
							// pt 2 should be along a line between 1 and 3
							double du = state4Strain(2)-state4Strain(0);
							double df = state4Stress(2)-state4Stress(0);
							state4Strain(1) = state4Strain(0) + 0.5*du;
							state4Stress(1) = state4Stress(0) + 0.5*df;
						}
						else if (state4Strain(2) < 0.0) {
							// pt 2 should be along a line between 2 and 4
							double du = state4Strain(3)-state4Strain(1);
							double df = state4Stress(3)-state4Stress(1);
							state4Strain(2) = state4Strain(1) + 0.5*du;
							state4Stress(2) = state4Stress(1) + 0.5*df;
						}
						else {
							double avgforce = 0.5*(state4Stress(2) + state4Stress(1));
							double dfr = 0.0;
							if (avgforce < 0.0){
								dfr = -avgforce/100;
							}
							else {
								dfr = avgforce/100;
							}
							double slope12 = (state4Stress(1) - state4Stress(0))/(state4Strain(1) - state4Strain(0));
							double slope34 = (state4Stress(3) - state4Stress(2))/(state4Strain(3) - state4Strain(2));
							state4Stress(1) = avgforce - dfr;
							state4Stress(2) = avgforce + dfr;
							state4Strain(1) = state4Strain(0) + (state4Stress(1) - state4Stress(0))/slope12;
							state4Strain(2) = state4Strain(3) - (state4Stress(3) - state4Stress(2))/slope34;
						}
					}
				}
			}
				else {
					// linear unload reload path is expected
					double du = state4Strain(3)-state4Strain(0);
					double df = state4Stress(3)-state4Stress(0);
					state4Strain(1) = state4Strain(0) + 0.33*du;
					state4Strain(2) = state4Strain(0) + 0.67*du;
					state4Stress(1) = state4Stress(0) + 0.33*df;
					state4Stress(2) = state4Stress(0) + 0.67*df;
				}
			

				
				double checkSlope = state4Stress(0)/state4Strain(0);
				double slope = 0.0;

				// final check
				int i = 0;
				while (i<3) {
					double du = state4Strain(i+1)-state4Strain(i);
					double df = state4Stress(i+1)-state4Stress(i);
					if (du<0.0 || df<0.0) {
						double du = state4Strain(3)-state4Strain(0);
						double df = state4Stress(3)-state4Stress(0);
						state4Strain(1) = state4Strain(0) + 0.33*du;
						state4Strain(2) = state4Strain(0) + 0.67*du;
						state4Stress(1) = state4Stress(0) + 0.33*df;
						state4Stress(2) = state4Stress(0) + 0.67*df;
						slope = df/du;
						i = 3;
					}
					if (slope > 1e-8 && slope < checkSlope) {
						state4Strain(1) = 0.0; state4Stress(1) = 0.0;
						state4Strain(2) = state4Strain(3)/2; state4Stress(2) = state4Stress(3)/2;
					} 

					i++;
				}
		// added jun 9th 2004 to prevent getting a negative slope
		if (state4Stress(2) <= state4Stress(1)) {
			state4Stress(2) = state4Stress(1)*1.02;
		}
}

double BarSlipMaterial::Envlp3Tangent(Vector s3Strain, Vector s3Stress, double u)
			{
				double k = 0.0;
				int i = 0;
				while ((k==0.0||i<=2) && (i<=2)) 
				{
					if (u>= s3Strain(i)) {
						k = (s3Stress(i+1)-s3Stress(i))/(s3Strain(i+1)-s3Strain(i));
					}
					i++;
				}
				if (k==0.0) {
					if (u<s3Strain(0)) {
						i = 0;
					}
					else {
						i = 2;
					}
					k = (s3Stress(i+1)-s3Stress(i))/(s3Strain(i+1)-s3Strain(i));
					
				}
				return k;
}

double BarSlipMaterial::Envlp4Tangent(Vector s4Strain, Vector s4Stress, double u)
			{
				double k = 0.0;
				int i = 0;
				while ((k==0.0||i<=2) && (i<=2)) 
				{
					if (u>= s4Strain(i)) {
						k = (s4Stress(i+1)-s4Stress(i))/(s4Strain(i+1)-s4Strain(i));
					}
					i++;
				}
				if (k==0.0) {
					if (u<s4Strain(0)) {
						i = 0;
					}
					else {
						i = 2;
					}
					k = (s4Stress(i+1)-s4Stress(i))/(s4Strain(i+1)-s4Strain(i));
					
				}
				return k;
}

double BarSlipMaterial::Envlp3Stress(Vector s3Strain, Vector s3Stress, double u)
			{
				double k = 0.0;
				int i = 0;
				double f = 0.0;
				while ((k==0.0||i<=2) && (i<=2)) 
				{
					if (u>= s3Strain(i)) {
						k = (s3Stress(i+1)-s3Stress(i))/(s3Strain(i+1)-s3Strain(i));
						f = s3Stress(i)+(u-s3Strain(i))*k;
					}
					i++;
				}
				if (k==0.0) {
					if (u<s3Strain(0)) {
						i = 0;
					}
					else {
						i = 2;
					}
					k = (s3Stress(i+1)-s3Stress(i))/(s3Strain(i+1)-s3Strain(i));
					f = s3Stress(i)+(u-s3Strain(i))*k;
				}
				return f;
}

double BarSlipMaterial::Envlp4Stress(Vector s4Strain, Vector s4Stress, double u)
			{
				double k = 0.0;
				int i = 0;
				double f = 0.0;
				while ((k==0.0||i<=2) && (i<=2)) 
				{
					if (u>= s4Strain(i)) {
						k = (s4Stress(i+1)-s4Stress(i))/(s4Strain(i+1)-s4Strain(i));
						f = s4Stress(i)+(u-s4Strain(i))*k;
					}
					i++;
				}
				if (k==0.0) {
					if (u<s4Strain(0)) {
						i = 0;
					}
					else {
						i = 2;
					}
					k = (s4Stress(i+1)-s4Stress(i))/(s4Strain(i+1)-s4Strain(i));
					f = s4Stress(i)+(u-s4Strain(i))*k;
				}
				return f;
}	

void BarSlipMaterial::updateDmg(double strain)
{
	double umaxAbs = (TmaxStrainDmnd>-TminStrainDmnd) ? TmaxStrainDmnd:-TminStrainDmnd;
	double uultAbs = (envlpPosStrain(4)>-envlpNegStrain(4)) ? envlpPosStrain(4):-envlpNegStrain(4);
	if ((strain<uultAbs && strain>-uultAbs)&& Tenergy< energyCapacity)
	{
		TgammaK = gammaK1*pow((umaxAbs/uultAbs),gammaK3);
		TgammaD = gammaD1*pow((umaxAbs/uultAbs),gammaD3);
		if (damage == 2 || damage == 0) {
			TgammaF = gammaF1*pow((umaxAbs/uultAbs),gammaF3);
		}

		if (damage == 1) {
			if (umaxAbs >= envlpPosStrain(3)) {
				double a = gammaFLimit/0.7;
				double b = -gammaFLimit*3.0/7.0;
				double x = (umaxAbs/uultAbs);
				TgammaF = a*x + b;
			}
		}

		if (Tenergy>elasticStrainEnergy) {
			TgammaK = TgammaK + gammaK2*pow(((Tenergy-elasticStrainEnergy)/energyCapacity),gammaK4);
			TgammaD = TgammaD + gammaD2*pow(((Tenergy-elasticStrainEnergy)/energyCapacity),gammaD4);
			TgammaF = TgammaF + gammaF2*pow(((Tenergy-elasticStrainEnergy)/energyCapacity),gammaF4);
		}
		double kminP = (posEnvlpStress(TmaxStrainDmnd)/TmaxStrainDmnd);
		double kminN = (negEnvlpStress(TminStrainDmnd)/TminStrainDmnd);
		double kmin = (kminP/kElasticPos>kminN/kElasticNeg) ? kminP/kElasticPos:kminN/kElasticNeg;
		double gammaKLimEnv = (0.0>(1-kmin)) ? 0.0:(1-kmin);
		double k1 = (TgammaK<gammaKLimit) ? TgammaK:gammaKLimit;
		TgammaK = (k1<gammaKLimEnv) ? k1:gammaKLimEnv;
		TgammaD = (TgammaD<gammaDLimit) ? TgammaD:gammaDLimit;
		TgammaF = (TgammaF<gammaFLimit) ? TgammaF:gammaFLimit;
	}
	else if (strain<uultAbs && strain>-uultAbs) {
		double kminP = (posEnvlpStress(TmaxStrainDmnd)/TmaxStrainDmnd);
		double kminN = (negEnvlpStress(TminStrainDmnd)/TminStrainDmnd);
		double kmin = (kminP/kElasticPos>kminN/kElasticNeg) ? kminP/kElasticPos:kminN/kElasticNeg;
		double gammaKLimEnv = (0.0>(1-kmin)) ? 0.0:(1-kmin);
			
		TgammaK = (gammaKLimit<gammaKLimEnv) ? gammaKLimit:gammaKLimEnv;    
		TgammaD = gammaDLimit;
		TgammaF = gammaFLimit;
	}		
}
