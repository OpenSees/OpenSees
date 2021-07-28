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
                                                                        
// $Revision: 1.7 $
// $Date: 2009/03/23 23:17:04 $
// $Source: /usr/local/cvs/OpenSees/PACKAGES/NewMaterial/cpp/Trilinwp.cpp,v $
                                                                        
#include <elementAPI.h>
#include <Trilinwp2.h>

#include <Vector.h>
#include <Channel.h>
#include <math.h>
#include <float.h>


void* OPS_Trilinwp2()
{
  // print out some KUDO's
 // if (numTrilinwp == 0) {
  //  opserr << "Trilineal with pinching unaxial material - Written by GST UNcuyo Copyright 2017 - Use at your Own Peril\n";
  //  numTrilinwp =1;
  //}

  // Pointer to a uniaxial material that will be returned
 // UniaxialMaterial *theMaterial = 0;

  //
  // parse the input line for the material parameters
  //

  int    iData[2];
  double dData[15];
  int numData; 
  int numDatatot;

  numData = 1;
  if (OPS_GetIntInput(&numData, &iData[0]) != 0) {
    opserr << "WARNING invalid uniaxialMaterial Trilinwp2 tag" << endln;
    return 0;
  }
  numDatatot=numData;

  numData = 15;
  if (OPS_GetDoubleInput(&numData, dData) != 0) {
    opserr << "WARNING invalid parameters\n";
    return 0;	
  }
  numDatatot=numDatatot+numData;
   numData = 1;
  if (OPS_GetIntInput(&numData, &iData[1]) != 0) {
    opserr << "WARNING invalid uniaxialMaterial Trilinwp2 type" << endln;
    return 0;
  }
  numDatatot=numDatatot+numData;
   if (numDatatot != 17 ) { 
    opserr << "Invalid Args want: uniaxialMaterial Trilinwp2 tag? Fcrp? dcrp? Fyp? dyp? Fup? dup? px? py? d1? d2? beta? Pt? Pb? Pc? Mb? itype?  ";            
    return 0;        
  }


   int itype=iData[1];

  // 
  // create a new material
  //

  UniaxialMaterial* mat = new Trilinwp2(iData[0], dData[0], dData[1], dData[2], dData[3], dData[4], dData[5], dData[6], dData[7], dData[8], dData[9], dData[10], dData[11], dData[12], dData[13], dData[14], iData[1]);


  if (mat == 0) {
    opserr << "WARNING could not create uniaxialMaterial of type Trilinwp2\n";
    return 0;
  }

  // return the material
  return mat;
}




Trilinwp2::Trilinwp2(int tag,
			double m1p, double r1p, double m2p, double r2p, double m3p, double r3p,
			double px, double py, double d1, double d2, double b, double ptn, double pbn,double pcn, double mbn, int it):
UniaxialMaterial(tag, MAT_TAG_Trilinwp2),
pinchX(px), pinchY(py), damfc1(d1), damfc2(d2), beta(b),
mom1p(m1p), rot1p(r1p), mom2p(m2p), rot2p(r2p), mom3p(m3p), rot3p(r3p),pt(ptn), pb(pbn),pc(pcn),mb(mbn),itype(it)
{
	bool error = false;
	// Positive backbone parameters
	if (rot1p <= 0.0)
		error = true;

	if (rot2p <= rot1p)
		error = true;

	if (rot3p <= rot2p)
		error = true;

	rot1n=-1.0*rot1p;
	rot2n=-1.0*rot2p;
	rot3n=-1.0*rot3p;
	mom1n=-1.0*mom1p;
	mom2n=-1.0*mom2p;
	mom3n=-1.0*mom3p;

	// Negative backbone parameters
	if (rot1n >= 0.0)
		error = true;

	if (rot2n >= rot1n)
		error = true;

	if (rot3n >= rot2n)
		error = true;
	
	if (error) {
	  opserr << "Trilinwp2::Trilinwp2 -- input backbone is not unique (one-to-one)\n";
	  exit(-1);
	}		
	
	mom1pi=mom1p;
	mom2pi=mom2p;
	mom3pi=mom3p;
	mom1ni=mom1n;
	mom2ni=mom2n;
	mom3ni=mom3n;
	rot1pi=rot1p;
	rot2pi=rot2p;
	rot3pi=rot3p;
	
	energyA = 0.5 * (rot1p*mom1p + (rot2p-rot1p)*(mom2p+mom1p) + (rot3p-rot2p)*(mom3p+mom2p) +
		rot1n*mom1n + (rot2n-rot1n)*(mom2n+mom1n) + (rot3n-rot2n)*(mom3n+mom2n));

	// Set envelope slopes
	this->setEnvelope();

	// Initialize history variables
	this->revertToStart();
	this->revertToLastCommit();

}



Trilinwp2::Trilinwp2():
UniaxialMaterial(0, MAT_TAG_Trilinwp2),
pinchX(0.0), pinchY(0.0), damfc1(0.0), damfc2(0.0), beta(0.0),
mom1p(0.0), rot1p(0.0), mom2p(0.0), rot2p(0.0), mom3p(0.0), rot3p(0.0),
pt(0.0), pb(0.0),pc(0.0),mb(0.0),itype(0)
{

}

Trilinwp2::~Trilinwp2()
{
	// Nothing to do
}

int
Trilinwp2::setTrialStrain(double strain, double strainRate)
{
 
//	opserr<<"m1p= "<<mom1p<<"  m1n= "<<mom1n<<endln;
	
if (TloadIndicator == 0 && strain == 0.0)
    return 0;
//double N=strainRate;   //Fuerza axial >0= traccion
//opserr << "N= "<< N<<endln;
//envolvente de flexion
double N=1.0*(strainRate);

if (itype == 1) {
	//flexure: ref.: ph.D. thesis Gonzalo Torrisi (2012)
		//mom1: crack
		//mom2: yield
		//mom3: nominal
	// Balanced moment Mb and Axial compression capacity Pc not used yet...
	double duct0 = rot3pi / rot2pi;
	double stiff1 = mom1pi / rot1pi;
	double stiff2 = (mom3pi - mom2pi) / (rot3pi - rot2pi);

	if (N > 0) {
		double duct = duct0;
	}
	else if (N < 0 && N > pb) {
		double duct = (1.0 - duct0) * N / pb + duct0;
	}
	else if (N < pb) {
		double duct = 1.1;
	}
	else {
		double duct = duct0;
	}
	duct = (duct > duct0) ? duct0 : duct;
	if (N > 0 && N < pt) {
		mom3p = mom3pi * (1 - (N * N / (pb * pb)) - N / (1.2 * pt));
		mom3p = (mom3p < 0) ? 0 : mom3p;
		mom3p = (mom3p > mom3pi) ? mom3pi : mom3p;
		mom2p = mom2pi * (1 - (N * N / (pb * pb)) - N / (1.2 * pt));
		mom2p = (mom2p < 0) ? 0 : mom2p;
		mom2p = (mom2p > mom2pi) ? mom2pi : mom2p;
		mom1p = mom2p / 10;
		rot1p = mom1p / stiff1;
		rot2p = rot2pi;
		rot1p = (rot1p > rot2p) ? 0.9 * rot2p : rot1p;
		rot3p = duct * rot2p;
	}
	else if (N > pt) {
		mom1p = mom1pi / 100.0;
		mom1n = mom1ni / 100.0;
		mom2p = mom2pi / 100.0;
		mom2n = mom2ni / 100.0;
		mom3p = mom3pi / 100.0;
		mom3n = mom3ni / 100.0;
	}
	else if (N<0 && N>pb/3) {
		mom3p = mom3pi * (1 - (N * N / (pb * pb)) - N / (1.2 * pt));
		mom3p = (mom3p < 0) ? 0 : mom3p;
		mom2p = mom2pi * (1 - (N * N / (pb * pb)) - N / (1.2 * pt));
		mom2p = (mom2p < 0) ? 0 : mom2p;
		mom1p = mom1pi * (1 - (N * N / (pb * pb)) - N / (1.2 * pt));
		mom1p = (mom1p < 0) ? 0 : mom1p;
		rot1p = mom1p / stiff1;
		if (rot1p > rot2pi) {
			rot1p = 0.9 * rot2pi;
		}
		rot2p = rot2pi;
		rot3p = rot2p * duct;
	}
	else if (N <=pb/3 && N >pb/2 ) {
		mom3p = mom3pi * (1 - (N * N / (pb * pb)) - N / (1.2 * pt));
		mom3p = (mom3p < 0) ? 0 : mom3p;
		mom2p = mom2pi * mom3p / mom3pi;
		mom1p = 0.9*mom2p;
		rot1p = mom1p/stiff1;
		rot2p = rot2pi;
		if (rot1p > rot2pi) {
			rot1p = 0.9 * rot2pi;
		}
		rot3p = rot2p * duct;
	}
	else if (N <= pb / 2 && N>pb) {
		mom3p = mom3pi * (1 - (N * N / (pb * pb)) - N / (1.2 * pt));
		mom3p = (mom3p < 0) ? 0 : mom3p;
		mom3p = (mom3p < 0) ? -1.0 * mom3p : mom3p;
		mom2p = 0.95 * mom3p;
		mom1p = 0.95 * mom2p;
		rot1p = mom1p/stiff1;
		if (rot1p > rot2pi) {
			rot1p = 0.9 * rot2pi;
		}
		rot2p = rot2pi;
		rot3p = duct * rot2p;
	}
	else if (N <= pb) {
		mom3p = mom3pi * (1 - (N * N / (pb * pb)) - N / (1.2 * pt));
		mom3p = (mom3p < 0) ? 0 : mom3p;
		mom2p = 0.95 * mom3p;
		mom1p = 0.95 * mom2p;
		rot1p = 0.95 * rot2pi;
		rot2p = rot2pi;
		rot3p = duct * rot2p;
	}
	else {
		mom1p = mom1pi;
		mom1n = mom1ni;
		mom2p = mom2pi;
		mom2n = mom2ni;
		mom3p = mom3pi;
		mom3n = mom3ni;
	}
	mom1n = -1.0 * mom1p;
	mom2n = -1.0 * mom2p;
	mom3n = -1.0 * mom3p;
	rot1n = -1.0 * rot1p;
	rot2n = -1.0 * rot2p;
	rot3n = -1.0 * rot3p;

}
else if (itype == 2) {
	double ksp = mom2pi / (rot2pi - rot1pi);
	if (N > 0 && N < pt) {
		double expo = 2.5;
		mom1p = mom1pi * (1.0 - pow(N / pt, expo));
		mom1n = mom1ni * (1.0 - pow(N / pt, expo));
		mom2p = mom2pi * (1.0 - pow(N / pt, expo));
		mom2n = mom2ni * (1.0 - pow(N / pt, expo));
		mom3p=mom3pi*(1.0-pow(N/pt,expo));
		mom3n=mom3ni*(1.0-pow(N/pt,expo));

	}
	else if (N > pt) {
		mom1p = mom1pi / 100.0;
		mom1n = mom1ni / 100.0;
		mom2p = mom2pi / 100.0;
		mom2n = mom2ni / 100.0;
		mom3p = mom3pi / 100.0;
		mom3n = mom3ni / 100.0;

	}
	else if (N<0 && N>pb) {
		mom1p = mom1pi * (1.0 + pow(N / pb, 2));
		mom1n = mom1ni * (1.0 + pow(N / pb, 2));
		mom2p = mom2pi * (1.0 + pow(N / pb, 2));
		mom2n = mom2ni * (1.0 + pow(N / pb, 2));
		mom3p=mom3pi*(1.0+pow(N/pb,2));
		mom3n=mom3ni*(1.0+pow(N/pb,2));

	}
	else if (N < 0 && N < pb) {
		mom1p = mom1pi * (1.0 + pow(N / pb, 2));
		mom1n = mom1ni * (1.0 + pow(N / pb, 2));
		mom2p = mom2pi * (1.0 + pow(N / pb, 2));
		mom2n = mom2ni * (1.0 + pow(N / pb, 2));
		mom3p=mom3pi*(1.0+pow(N/pb,2));
		mom3n=mom3ni*(1.0+pow(N/pb,2));
	}
	else {
		mom1p = mom1pi;
		mom1n = mom1ni;
		mom2p = mom2pi;
		mom2n = mom2ni;
		mom3p = mom3pi;
		mom3n = mom3ni;
	}
	rot2p = (mom2p / ksp) * (1 - 1.4 * N / pb) + rot1p;
    rot2n = (mom2n / ksp) * (1 - 1.4 * N / pb) + rot1n;
    rot3p = 2 * rot2p;
    rot3n = 2 * rot2n;
}
else {
mom1p = mom1pi;
mom1n = mom1ni;
mom2p = mom2pi;
mom2n = mom2ni;
mom3p = mom3pi;
mom3n = mom3ni;
}

this->setEnvelope();

  TrotMax = CrotMax;
  TrotMin = CrotMin;
  TenergyD = CenergyD;
  TrotPu = CrotPu;
  TrotNu = CrotNu;

  Tstrain = strain;
  double dStrain = Tstrain - Cstrain;

  if (fabs(dStrain) < DBL_EPSILON)
    return 0;
  
  TloadIndicator = CloadIndicator;
  
  if (TloadIndicator == 0)
    TloadIndicator = (dStrain < 0.0) ? 2 : 1;
  
  if (Tstrain >= CrotMax) {
    TrotMax = Tstrain;
    Ttangent = posEnvlpTangent(Tstrain);
    Tstress = posEnvlpStress(Tstrain);
	TloadIndicator=1;
  }
  else if (Tstrain <= CrotMin) {
    TrotMin = Tstrain;
    Ttangent = negEnvlpTangent(Tstrain);
    Tstress = negEnvlpStress(Tstrain);
	TloadIndicator=2;
  }
  else {
    if (dStrain < 0.0)
      negativeIncrement(dStrain);
    else if (dStrain > 0.0)
      positiveIncrement(dStrain);
  }
  
  TenergyD = CenergyD + 0.5*(Cstress+Tstress)*dStrain;


  //  if (this->getTag() == 40)
  //    opserr << "setTrial: " << Tstrain << " " << Ttangent << " " << Tstress << endln;
  
  return 0;
}


double
Trilinwp2::getStrain(void)
{
	return Tstrain;
}

double
Trilinwp2::getStress(void)
{
	return Tstress;
}


double
Trilinwp2::getTangent(void)
{
  return Ttangent;
}

void
Trilinwp2::positiveIncrement(double dStrain)
{
	double kn = pow(CrotMin/rot1n,beta);
	kn = (kn < 1.0) ? 1.0 : 1.0/kn;
	double kp = pow(CrotMax/rot1p,beta);
	kp = (kp < 1.0) ? 1.0 : 1.0/kp;

	if (TloadIndicator == 2) {
		TloadIndicator = 1;
		if (Cstress <= 0.0) {
			TrotNu = Cstrain - Cstress/(Eun*kn);
			double energy = CenergyD - 0.5*Cstress/(Eun*kn)*Cstress;
			double damfc = 0.0;
			if (CrotMin < rot1n) {
				damfc = damfc2*energy/energyA;
				damfc += damfc1*(CrotMin-rot1n)/rot1n;
			}

			TrotMax = CrotMax*(1.0+damfc);
		}
	}

  TloadIndicator = 1;

	TrotMax = (TrotMax > rot1p) ? TrotMax : rot1p;

	double maxmom = posEnvlpStress(TrotMax);
	double rotlim = negEnvlpRotlim(CrotMin);
	double rotrel = (rotlim > TrotNu) ? rotlim : TrotNu;

	// rotrel = TrotNu;
	// if (negEnvlpStress(CrotMin) >= 0.0)
	//    rotrel = rotlim;
	
	//	double rotmp1 = rotrel + pinchY*(TrotMax-rotrel);

	double rotmp2 = TrotMax - (1.0-pinchY)*maxmom/(Eup*kp);
	//double rotmp2 = TrotMax-(1-pinchY)*maxmom/Eup;
	//	double rotch = rotmp1 + (rotmp2-rotmp1)*pinchX;
	double rotch = rotrel + (rotmp2-rotrel)*pinchX;                   // changed on 7/11/2006

	double tmpmo1;
	double tmpmo2;

	if (Tstrain < TrotNu) {
		Ttangent = Eun*kn;
		Tstress = Cstress + Ttangent*dStrain;
		if (Tstress >= 0.0) {
			Tstress = 0.0;
			Ttangent = Eun*1.0e-9;
		}
	}

	else if (Tstrain >= TrotNu && Tstrain < rotch) {
		if (Tstrain <= rotrel) {
			Tstress = 0.0;
			Ttangent = Eup*1.0e-9;
		}
		else {
			Ttangent = maxmom*pinchY/(rotch-rotrel);
			tmpmo1 = Cstress + Eup*kp*dStrain;
			tmpmo2 = (Tstrain-rotrel)*Ttangent;
			if (tmpmo1 < tmpmo2) {
				Tstress = tmpmo1;
				Ttangent = Eup*kp;
			}
			else
				Tstress = tmpmo2;
		}
	}

	else {
		Ttangent = (1.0-pinchY)*maxmom/(TrotMax-rotch);
		tmpmo1 = Cstress + Eup*kp*dStrain;
		tmpmo2 = pinchY*maxmom + (Tstrain-rotch)*Ttangent;
		if (tmpmo1 < tmpmo2) {
			Tstress = tmpmo1;
			Ttangent = Eup*kp;
		}
		else
			Tstress = tmpmo2;
	}
}

void
Trilinwp2::negativeIncrement(double dStrain)
{
	double kn = pow(CrotMin/rot1n,beta);
	kn = (kn < 1.0) ? 1.0 : 1.0/kn;
	double kp = pow(CrotMax/rot1p,beta);
	kp = (kp < 1.0) ? 1.0 : 1.0/kp;

	if (TloadIndicator == 1) {
		TloadIndicator = 2;
		if (Cstress >= 0.0) {
			TrotPu = Cstrain - Cstress/(Eup*kp);
			double energy = CenergyD - 0.5*Cstress/(Eup*kp)*Cstress;
			double damfc = 0.0;
			if (CrotMax > rot1p) {
				damfc = damfc2*energy/energyA;
				damfc += damfc1*(CrotMax-rot1p)/rot1p;
			}

			TrotMin = CrotMin*(1.0+damfc);
		}
	}

  TloadIndicator = 2;

	TrotMin = (TrotMin < rot1n) ? TrotMin : rot1n;

	double minmom = negEnvlpStress(TrotMin);
	double rotlim = posEnvlpRotlim(CrotMax);
	double rotrel = (rotlim < TrotPu) ? rotlim : TrotPu;

	//rotrel = TrotPu;
	//if (posEnvlpStress(CrotMax) <= 0.0)
	//  rotrel = rotlim;

	//double rotmp1 = rotrel + pinchY*(TrotMin-rotrel);
	double rotmp2 = TrotMin - (1.0-pinchY)*minmom/(Eun*kn);
	//double rotmp2 = TrotMin-(1-pinchY)*minmom/Eun;	
	//double rotch = rotmp1 + (rotmp2-rotmp1)*pinchX;
	double rotch = rotrel + (rotmp2-rotrel)*pinchX;                   // changed on 7/11/2006

	double tmpmo1;
	double tmpmo2;

	if (Tstrain > TrotPu) {
		Ttangent = Eup*kp;
		Tstress = Cstress + Ttangent*dStrain;
		if (Tstress <= 0.0) {
			Tstress = 0.0;
			Ttangent = Eup*1.0e-9;
		}
	}

	else if (Tstrain <= TrotPu && Tstrain > rotch) {
		if (Tstrain >= rotrel) {
			Tstress = 0.0;
			Ttangent = Eun*1.0e-9;
		}
		else {
			Ttangent = minmom*pinchY/(rotch-rotrel);
			tmpmo1 = Cstress + Eun*kn*dStrain;
			tmpmo2 = (Tstrain-rotrel)*Ttangent;
			if (tmpmo1 > tmpmo2) {
				Tstress = tmpmo1;
				Ttangent = Eun*kn;
			}
			else
				Tstress = tmpmo2;
		}
	}

	else {
		Ttangent = (1.0-pinchY)*minmom/(TrotMin-rotch);
		tmpmo1 = Cstress + Eun*kn*dStrain;
		tmpmo2 = pinchY*minmom + (Tstrain-rotch)*Ttangent;
		if (tmpmo1 > tmpmo2) {
			Tstress = tmpmo1;
			Ttangent = Eun*kn;
		}
		else
			Tstress = tmpmo2;
	}
}

int
Trilinwp2::commitState(void)
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
Trilinwp2::revertToLastCommit(void)
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
Trilinwp2::revertToStart(void)
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

	mom1pi=mom1p;
	mom2pi=mom2p;
	mom3pi=mom3p;
	mom1ni=mom1n;
	mom2ni=mom2n;
	mom3ni=mom3n;
	return 0;
}


UniaxialMaterial*
Trilinwp2::getCopy(void)
{
	Trilinwp2 *theCopy = new Trilinwp2 (this->getTag(),
		mom1pi, rot1pi, mom2pi, rot2pi, mom3pi, rot3pi,
		pinchX, pinchY, damfc1, damfc2, beta,pt,pb,pc,mb, itype);

	theCopy->CrotMax = CrotMax;
	theCopy->CrotMin = CrotMin;
	theCopy->CrotPu = CrotPu;
	theCopy->CrotNu = CrotNu;
	theCopy->CenergyD = CenergyD;
	theCopy->CloadIndicator = CloadIndicator;
	theCopy->Cstress = Cstress;
	theCopy->Cstrain = Cstrain;
	theCopy->Ttangent = Ttangent;
	//theCopy->N = N;
	//theCopy->mom1p = mom1p;
	//theCopy->mom2p = mom2p;
	//theCopy->mom3p = mom3p;
	//theCopy->mom1n = mom1n;
	//theCopy->mom2n = mom2n;
	//theCopy->mom3n = mom3n;
	//theCopy->rot1p = rot1p;
	//theCopy->rot2p = rot2p;
	//theCopy->rot3p = rot3p;
	//theCopy->rot1n = rot1n;
	//theCopy->rot2n = rot2n;
	//theCopy->rot3n = rot3n;



	return theCopy;
}
int
Trilinwp2::sendSelf(int commitTag, Channel &theChannel)
{
  int res = 0;
  
  static Vector data(32);
  
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
  data(27) = pt;
  data(28) = pb;
  data(29) = pc;
  data(30) = mb;
  data(31)=itype;

  res = theChannel.sendVector(this->getDbTag(), commitTag, data);
  if (res < 0) 
    opserr << "Trilinwp2::sendSelf() - failed to send data\n";


  return res;
}

int
Trilinwp2::recvSelf(int commitTag, Channel &theChannel, 
			FEM_ObjectBroker &theBroker)
{
  int res = 0;
  
  static Vector data(32);
  res = theChannel.recvVector(this->getDbTag(), commitTag, data);
  
  if (res < 0) {
      opserr << "Trilinwp2::recvSelf() - failed to receive data\n";
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
	pt = data(27);
	pb = data(28);
	pc = data(29);
	mb = data(30);
	itype=data(31);
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

  // Set envelope slopes
  this->setEnvelope();
  
  return 0;
}

void
Trilinwp2::Print(OPS_Stream &s, int flag)
{
	s << "Trilineal with pinching material - GST(2017), tag: " << this->getTag() << endln;
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
	s << "Pt: " << pt << endln;	
	s << "Pb: " << pb << endln;
	s << "Pc: " << pc << endln;
	s << "Mb: " << mb << endln;
	s << "itype "<< itype<< endln;
	if(itype==1){
	s << "Type: Flexure" << endln;
	}
	else if(itype==2){
	s << "Type: Shear" << endln;
	}
	else {
	s << "Shear and flexure with no Axial Interaction" << endln;
	}
}

void
Trilinwp2::setEnvelope(void)
{
	E1p = mom1p/rot1p;
	E2p = (mom2p-mom1p)/(rot2p-rot1p);
	E3p = (mom3p-mom2p)/(rot3p-rot2p);

	E1n = mom1n/rot1n;
	E2n = (mom2n-mom1n)/(rot2n-rot1n);
	E3n = (mom3n-mom2n)/(rot3n-rot2n);

	Eup = E1p;
	if (E2p > Eup) Eup = E2p;
	if (E3p > Eup) Eup = E3p;

	Eun = E1n;
	if (E2n > Eun) Eun = E2n;
	if (E3n > Eun) Eun = E3n;


}

double
Trilinwp2::posEnvlpStress(double strain)
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
Trilinwp2::negEnvlpStress(double strain)
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
Trilinwp2::posEnvlpTangent(double strain)
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
Trilinwp2::negEnvlpTangent(double strain)
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
Trilinwp2::posEnvlpRotlim(double strain)
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
Trilinwp2::negEnvlpRotlim(double strain)
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
