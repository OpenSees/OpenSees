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
                                                                        
// $Revision: 1.1 $
// $Date: 2003-06-11 18:19:44 $
// $Source: /usr/local/cvs/OpenSees/SRC/material/uniaxial/BarSlipMaterial.cpp,v $
                                                                        
                                                                        
// Written: NM (nmitra@u.washington.edu) 
// Created: January 2002
//
// Description: This file contains the class implementation for 
// bar-slip material which is defined by 4 points on the positive and 
// negative envelopes and a bunch of damage parameters. The material accounts for
// 3 types of damage rules : Strength degradation, Stiffness degradation, 
// unloading stiffness degradation. The theory for the material has been provided by Prof. Lowes


#include <BarSlipMaterial.h>
#include <math.h>
#include <float.h>

BarSlipMaterial::BarSlipMaterial(int tag,
				 double f, double fs, double es, double fsu,
				 double eh, double dbar, double ljoint, int n, double w, double d, 
				 int bsf, int typ):
  UniaxialMaterial(tag, MAT_TAG_BarSlip),
  bsflag(bsf), type(typ), width(w), depth(d),
  fc(f), fy(fs), Es(es), fu(fsu), Eh(eh), db(dbar), nbars(n), ld(ljoint), 
  rDispP(0.25), rForceP(0.25), uForceP(0.0),
  rDispN(0.25), rForceN(0.25), uForceN(0.0),
  gammaK1(0.3), gammaK2(0.0), gammaK3(0.1), gammaK4(0.0), gammaKLimit(0.4),
  gammaD1(0.6), gammaD2(0.0), gammaD3(0.2), gammaD4(0.0), gammaDLimit(0.25),
  gammaF1(0.7), gammaF2(0.3), gammaF3(0.5), gammaF4(0.1), gammaFLimit(0.0),
  gammaE(10.0), eP(4,2), eN(4,2)
{
	if (fc <= 0)
		opserr << "WARNING : BAR-SLIP -- fc should be positive entry" << endln;

	if (fc >= 1000)
	{
		unit = 1;
	}
	else
	{
		unit = 0;
	}


	if (unit == 0)   // consider "mpa" system being used
	{
		if (bsflag == 1)  // consider bond strength flag --- weak
		{
			tauYT = 0.05*pow(fc,0.5);
			tauET = 1.8*pow(fc,0.5);
			tauEC = 2.2*pow(fc,0.5);
			tauYC = 3.6*pow(fc,0.5);
			tauR = 0.15*pow(fc,0.5);
		}
		else if (bsflag == 0) // bond strength -- strong
		{
			tauYT = 0.4*pow(fc,0.5);
			tauET = 1.8*pow(fc,0.5);
			tauEC = 2.2*pow(fc,0.5);
			tauYC = 3.6*pow(fc,0.5);
			tauR = 0.15*pow(fc,0.5);
		}
	}
	else if (unit == 1)  // "psi" system used
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

	this->getBarSlipEnvelope();

//	opserr << eP << endln;
//	opserr << eN << endln;
	 material = new Pinching4Material (tag, eP(0,1), eP(0,0), eP(1,1), eP(1,0), eP(2,1), eP(2,0), eP(3,1), eP(3,0),
		eN(0,1), eN(0,0), eN(1,1), eN(1,0), eN(2,1), eN(2,0), eN(3,1), eN(3,0), rDispP, rForceP, uForceP, rDispN,
		rForceN, uForceN, gammaK1, gammaK2, gammaK3, gammaK4, gammaKLimit, gammaD1, gammaD2, gammaD3, 
		gammaD4, gammaDLimit, gammaF1, gammaF2, gammaF3, gammaF4, gammaFLimit, gammaE);

}

BarSlipMaterial::BarSlipMaterial(int tag,
				 double f, double fs, double es, double fsu,
				 double eh, double dbar, double ljoint, int n, double w, double d, 
				 double mdp, double mfp, double msp,
				 double mdn, double mfn, double msn,
				 double gk1, double gk2, double gk3, double gk4, double gklim,
				 double gd1, double gd2, double gd3, double gd4, double gdlim,
				 double gf1, double gf2, double gf3, double gf4, double gflim, double ge, int bsf, int typ):
  UniaxialMaterial(tag, MAT_TAG_BarSlip),
  bsflag(bsf), type(typ), width(w), depth(d),
  fc(f), fy(fs), Es(es), fu(fsu), Eh(eh), db(dbar), nbars(n), ld(ljoint), 
  rDispP(mdp), rForceP(mfp), uForceP(msp),
  rDispN(mdn), rForceN(mfn), uForceN(msn),
  gammaK1(gk1), gammaK2(gk2), gammaK3(gk3), gammaK4(gk4), gammaKLimit(gklim),
  gammaD1(gd1), gammaD2(gd2), gammaD3(gd3), gammaD4(gd4), gammaDLimit(gdlim),
  gammaF1(gf1), gammaF2(gf2), gammaF3(gf3), gammaF4(gf4), gammaFLimit(gflim),
  gammaE(ge),  eP(4,2), eN(4,2)
{
	if (fc <= 0)
		opserr << "WARNING : BAR-SLIP -- fc should be positive entry" << endln;

	if (fc >= 1000)
	{
		unit = 1;
	}
	else
	{
		unit = 0;
	}


	if (unit == 0)   // consider "mpa" system being used
	{
		if (bsflag == 1)  // consider bond strength flag --- weak
		{
			tauYT = 0.05*pow(fc,0.5);
			tauET = 1.8*pow(fc,0.5);
			tauEC = 2.2*pow(fc,0.5);
			tauYC = 3.6*pow(fc,0.5);
			tauR = 0.15*pow(fc,0.5);
		}
		else if (bsflag == 0) // bond strength -- strong
		{
			tauYT = 0.4*pow(fc,0.5);
			tauET = 1.8*pow(fc,0.5);
			tauEC = 2.2*pow(fc,0.5);
			tauYC = 3.6*pow(fc,0.5);
			tauR = 0.15*pow(fc,0.5);
		}
	}
	else if (unit == 1)  // "psi" system used
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

	this->getBarSlipEnvelope();

//	opserr << eP << endln;
//	opserr << eN << endln;
	 material = new Pinching4Material (tag, eP(0,1), eP(0,0), eP(1,1), eP(1,0), eP(2,1), eP(2,0), eP(3,1), eP(3,0),
		eN(0,1), eN(0,0), eN(1,1), eN(1,0), eN(2,1), eN(2,0), eN(3,1), eN(3,0), rDispP, rForceP, uForceP, rDispN,
		rForceN, uForceN, gammaK1, gammaK2, gammaK3, gammaK4, gammaKLimit, gammaD1, gammaD2, gammaD3, 
		gammaD4, gammaDLimit, gammaF1, gammaF2, gammaF3, gammaF4, gammaFLimit, gammaE);

}


BarSlipMaterial::BarSlipMaterial():
  UniaxialMaterial(0, MAT_TAG_BarSlip),
  bsflag(0), type(0), width(0.0), depth(0.0),
  fc(0.0), fy(0.0), Es(0.0), fu(0.0), Eh(0.0), db(0.0), nbars(0), ld(0.0),
  tauET(0.0), tauYT(0.0), tauEC(0.0), tauYC(0.0), tauR(0.0),
  rDispP(0.0), rForceP(0.0), uForceP(0.0),
  rDispN(0.0), rForceN(0.0), uForceN(0.0),
  gammaK1(0.0), gammaK2(0.0), gammaK3(0.0), gammaK4(0.0), gammaKLimit(0.0),
  gammaD1(0.0), gammaD2(0.0), gammaD3(0.0), gammaD4(0.0), gammaDLimit(0.0),
  gammaF1(0.0), gammaF2(0.0), gammaF3(0.0), gammaF4(0.0), gammaFLimit(0.0),
  gammaE(0.0),  eP(4,2), eN(4,2)
{
	// default constructor --- does nothing
}

BarSlipMaterial::~BarSlipMaterial()
{
	// destructor
}

void BarSlipMaterial::getBarSlipEnvelope(void)
{
	double delta;
	if (unit == 0)  // MPa
	{ 
		delta = 3.0;
	}
	else if (unit == 1) //psi
	{
		delta = 3.0/25.43;
	}
	double Ab = PI*pow(db,2)/4;
	double As = Ab*nbars;
	
	double frP;
	double frN;
	double geo;
	double let;
	double lec;
	double lyc;
	double lyt;
	double k1;
	double k2;
	double k3;
	double le1;
	double le2;
	double unloadNForce;
	double unloadPForce;

	frP = tauR*ld*PI*db*As/Ab;
	frN = -tauR*ld*PI*db*As/Ab;

	geo = PI*db/Ab;
	let = fy/(tauET*geo);
	lyt = (fu-fy)/(tauYT*geo);
	lec = fy/(tauEC*geo);
	lyc = (fu-fy)/(tauYC*geo);

	k1 = 2*Es*(tauET/fy)*geo*As;

	eP(0,0) = 0.5*fy*As/k1;
	eP(0,1) = 0.5*fy*As;
	eP(1,0) = fy*As/k1;
	eP(1,1) = fy*As;

	if (let+lyt<ld && bsflag ==0 )
	{
		k2 = (fu-fy)*As/(fy*lyt/Es + 0.5*geo*tauYT*pow(lyt,2)/Eh);
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

	if ((type == 0) || (type == 1))    // for beam
	{
		k3 = -0.1*k1;
	}
	else if (type == 2)  // for column
	{
		k3 = -0.005*k1;
	}
	eP(3,1) = frP; // eP(2,1);    // intentionally changed to remove softening region (originally ---- frP;) // dtd 19th may
	eP(3,0) = eP(2,0) + (frP - eP(2,1))/k3;

	double dth = 0.9*depth;
	double j;
	double beta = 0.85;
	
	double cForce;

	if ((type == 1) || (type == 2))  // member type -- beam bottom and column
	{
		j = 0.85;
	}
	else if (type ==0)  // member type -- beam top
	{
		j = 0.75;
	}


	cForce = 1 + (beta*fc*dth*width/Es/As*2*(1-j)/.003/(1-.1/.9*beta/2/(1-j)));
	As = As*cForce;
	k1 = 2*Es*(tauET/fy)*geo*As;


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
	eN(3,0) = 1.1*eN(2,0);
	eN(3,1) = eN(2,1) + k3*(eN(3,0)-eN(2,0));

	As = As/cForce;
	
	double tP = (ld<(lec+lyc)) ? ld:(lec+lyc);
	double tN = (ld<(let+lyt)) ? ld:(let+lyt);
	
	unloadPForce = tauR*PI*db*As*tP/Ab;
	unloadNForce = -tauR*PI*db*As*tN/Ab;


	uForceP = unloadPForce/eP(2,1);
	uForceN = unloadNForce/eN(2,1);

}


//*******************************************************************************
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

UniaxialMaterial* BarSlipMaterial::getCopy(void)
{
	BarSlipMaterial *theCopy = new BarSlipMaterial (this->getTag(), fc, fy, Es,
		fu, Eh, db, ld, nbars, width, depth, bsflag, type);
	
	theCopy->fc = fc;
	theCopy->fy = fy;
	theCopy->Es = Es;
	theCopy->fu = fu;
	theCopy->Eh = Eh;
	theCopy->db = db;
	theCopy->ld = ld;
	theCopy->nbars = nbars;
	theCopy->width = width;
	theCopy->depth = depth;
	theCopy->bsflag = bsflag;
	theCopy->type = type;
	
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
