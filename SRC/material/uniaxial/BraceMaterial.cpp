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

// Developed By: Amir Moshref; Amirkabir University of Technology (Tehran Polytechninc);a.moshref08@imperial.ac.uk
// Basic Source Written By: Chuang-Sheng Yang (Walter)

#include <BraceMaterial.h>
#include <G3Globals.h>
#include <OPS_Globals.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>
#include <Channel.h>

#include <elementAPI.h>

using namespace std;
//ErrorHandler *g3ErrorHandler;
static int numBilinMaterials = 0;


static int numCastMaterials = 0;

void *
OPS_Cast(void)
{
  if (numCastMaterials == 0) {
    numCastMaterials++;
    opserr << "Cast Fuse uniaxial material - Written by Dimitrios G. Lignos, Ph.D.\n";
  }
  
  // Pointer to a uniaxial material that will be returned
  UniaxialMaterial *theMaterial = 0;
  
  int    iData[1];
  double dData[17];
  int numData = 1;
  // Check Tag and number of Fingers
  if (OPS_GetIntInput(&numData, iData) != 0) {
    opserr << "WARNING invalid uniaxialMaterial  Cast Fuse tag" << endln;
    return 0;
  }
  
  int numRemaining = OPS_GetNumRemainingInputArgs();
  if (numRemaining != 17 || numRemaining != 13) {
    if (OPS_GetDoubleInput(&numRemaining, dData) != 0) {
      opserr << "Invalid Args want: uniaxialMaterial CastFuse tag? NLegs? bo? h? Fy? E? L? b? R0? cR1? cR2? a1? a2? a3? a4?";
      return 0;	
    }
    
    // Parsing was successful, allocate the material
    if (numRemaining == 13) {
      theMaterial = new BraceMaterial(iData[0], dData[0],
				      dData[1], dData[2], dData[3], dData[4], dData[5], 
				      dData[6], dData[7], dData[8], dData[9], dData[10],
				      dData[11], dData[12], dData[12]);
    } else {
      theMaterial = new BraceMaterial(iData[0], dData[0],
				      dData[1], dData[2], dData[3], dData[4], dData[5], 
				      dData[6], dData[7], dData[8], dData[9], dData[10],
				      dData[11], dData[12], dData[13],
				      dData[14], dData[15], dData[16]);
    }
  }

  if (theMaterial == 0) {
    opserr << "WARNING could not create uniaxialMaterial of type Cast Material\n";
    return 0;
  }
  
  return theMaterial;
}


BraceMaterial::BraceMaterial(int tag,
			     double m1p, double r1p, double m2p, double r2p, double m3p, double r3p,
			     double m1n, double r1n, double m2n, double r2n, double m3n, double r3n,
			     double px, double py, double d1, double d2, double b):
  UniaxialMaterial(tag, MAT_TAG_BraceMaterial),
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
    opserr << "BraceMaterial::BraceMaterial -- input backbone is not unique (one-to-one)\n";
    exit(-1);
  }	
  energyA = 0.5 * (rot1p*mom1p + (rot2p-rot1p)*(mom2p+mom1p) + (rot3p-rot2p)*(mom3p+mom2p) +
		   rot1n*mom1n + (rot3n-rot1n)*(mom3n+mom1n));
  
  mom2n = 0.5*(mom1n+mom3n);
  rot2n = 0.5*(rot1n+rot3n);
  
  // Set envelope slopes
  this->setEnvelope();
  
  // Initialize history variables
  this->revertToStart();
  this->revertToLastCommit();
}
BraceMaterial::BraceMaterial(int tag,
			     double m1p, double r1p, double m2p, double r2p,
			     double m1n, double r1n, double m2n, double r2n,
			     double px, double py, double d1, double d2, double b):
  UniaxialMaterial(tag, MAT_TAG_BraceMaterial),
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
    opserr << "BraceMaterial::BraceMaterial -- input backbone is not unique (one-to-one)\n";
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
BraceMaterial::BraceMaterial():
  UniaxialMaterial(0, MAT_TAG_BraceMaterial),
  pinchX(0.0), pinchY(0.0), damfc1(0.0), damfc2(0.0), beta(0.0),
  mom1p(0.0), rot1p(0.0), mom2p(0.0), rot2p(0.0), mom3p(0.0), rot3p(0.0),
  mom1n(0.0), rot1n(0.0), mom2n(0.0), rot2n(0.0), mom3n(0.0), rot3n(0.0)
{
}
BraceMaterial::~BraceMaterial()
{
  // Nothing to do
}
int
BraceMaterial::setTrialStrain(double strain, double strainRate)
{
  if (TloadIndicator == 0 && strain == 0.0)
    return 0;
  TrotMax = CrotMax;
  TrotMin = CrotMin;
  TenergyD = CenergyD;
  TrotPu = CrotPu;
  TrotNu = CrotNu;
  TmomNu = CmomNu;
  TmomPu = CmomPu;
  Tstrain = strain;
  double dStrain = Tstrain - Cstrain;
  TloadIndicator = CloadIndicator;
  Tidtime = Cidtime;
  Tidtime1 = Cidtime1;

  if (TloadIndicator == 0)
    TloadIndicator = (dStrain < 0.0) ? 2 : 1;
  if (Tstrain >= CrotMax) {
    TrotMax = Tstrain;
    Ttangent = posEnvlpTangent(Tstrain);
    Tstress = posEnvlpStress(Tstrain);
    TloadIndicator = 1;
  }
  else if (Tstrain <= CrotMin) {
    TrotMin = Tstrain;
    Ttangent = negEnvlpTangent(Tstrain);
    Tstress = negEnvlpStress(Tstrain);
    if (TrotMin <= rot1n && Tidtime == 0) {
      Tidtime = 1;
    }
    Tidtime1 = 0;
    TloadIndicator = 2;
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
BraceMaterial::getStrain(void)
{
return Tstrain;
}
double
BraceMaterial::getStress(void)
{
return Tstress;
}
double
BraceMaterial::getTangent(void)
{
return Ttangent;
}
void
BraceMaterial::positiveIncrement(double dStrain)
{
  // double kn = pow(CrotMin/rot1n,beta);
  double kn = 1.0;
  // kn = (kn < 1.0) ? 1.0 : 1.0/kn;
  // double kp = pow(CrotMax/rot1p,beta);
  double kp = 1.0;
  // kp = (kp < 1.0) ? 1.0 : 1.0/kp;
  // cout << "Hello: positive Increment" << endl;
  if (TloadIndicator == 2) {  
    TloadIndicator = 1;
    if (Cstress <= mom2p) {
      TrotNu = Cstrain;
      TmomNu = Cstress;
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
  // double rotlim = negEnvlpRotlim(CrotMin);
  // double rotrel = (rotlim > TrotNu) ? rotlim : TrotNu;
  // rotrel = TrotNu;
  // if (negEnvlpStress(CrotMin) >= 0.0)
  // rotrel = rotlim;
  double rotmp1 = TrotNu + pinchY*(TrotMax-TrotNu);
  double rotmp2 = TrotMax - (1.0-pinchY)*(maxmom-TmomNu)/(E1p*kp);
  double rotch = rotmp1 + (rotmp2-rotmp1)*pinchX;
  
  double tmpmo1;
  double tmpmo2;
  if (Tstrain < rotch) {
    Ttangent = (maxmom-TmomNu)*pinchY/(rotch-TrotNu);
    tmpmo1 = Cstress + E1p*kp*dStrain;
    tmpmo2 = TmomNu + (Tstrain-TrotNu)*Ttangent;
    if (tmpmo1 < tmpmo2) {
      Tstress = tmpmo1;
      Ttangent = E1p*kp;
      //printf("\nTstress1= %f",Tstress);
    }
    else
      Tstress = tmpmo2;
    //printf("\nTstress4= %f",Tstress);
  }
  else {
    Ttangent = (1.0-pinchY)*(maxmom-TmomNu)/(TrotMax-rotch);
    tmpmo1 = Cstress + E1p*kp*dStrain;
    tmpmo2 = TmomNu + pinchY*(maxmom-TmomNu) + (Tstrain-rotch)*Ttangent;
    if (tmpmo1 < tmpmo2) {
      Tstress = tmpmo1;
      Ttangent = E1p*kp;
      //printf("\nTstress2= %f",Tstress);
    }
    else
      Tstress = tmpmo2;
    //printf("\nTstress3= %f",Tstress);
  }
  if ( Tidtime == 1 && TrotMin > rot2n ) {
    if ( Tstress > mom2n && Tidtime1 == 0 ) {
      Tidtime1 = 1;
    }
  }
}
void
BraceMaterial::negativeIncrement(double dStrain)
{
  // double kn = pow(CrotMin/rot1n,beta);
  double kn = 1.0;
  // kn = (kn < 1.0) ? 1.0 : 1.0/kn;
  // double kp = pow(CrotMax/rot1p,beta);
  double kp = 1.0;
  // kp = (kp < 1.0) ? 1.0 : 1.0/kp;
  // cout << "Hello: negative Increment" << endl;
  // cout << "Tidtime: " << Tidtime << endl;
  
  if (TloadIndicator == 1) {
    //TloadIndicator = 2;
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
  //TrotMin = (TrotMin < rot1n) ? TrotMin : rot1n;
  //printf("\nTrotMin2= %f",TrotMin);
  if ( Tidtime == 1 && TrotMin > rot2n ) {
    double minmom = negEnvlpStress(TrotMin);
    if ( TloadIndicator == 1 ) {
      TmomPu = Cstress;
    }
    if ( Tidtime1 == 1) {
      if ( TmomPu >= mom3n ) {
	if ( Cstress > mom3n ) {
	  TrotPu = Cstrain - Cstress/(E1p*kp) + mom3n/(E1p*kp);
	}
	if ( Tstrain > TrotPu ) {
	  Ttangent = E1p*kp;
	  Tstress = Cstress + Ttangent*dStrain;
	  //printf("\nTNstress1= %f",Tstress);
	  if ( Tstress <= mom2n ) {
	    Tstress = mom2n;
	    Ttangent = E1p*1.0e-9;
	    //printf("\nTNstress2= %f",Tstress);
	  }
	}
	else {
	  Ttangent = (minmom-mom3n)/(TrotMin-TrotPu);
	  Tstress = mom3n + Ttangent*(Tstrain-TrotPu);
	  //printf("\nTNstress3= %f",Tstress);
	  if ( Tstress <= minmom ) {
	    Tstress = minmom;
	    Ttangent = Ttangent*1.0e-6;
	    //printf("\nTNstress4= %f",Tstress);
	    Tidtime1 = 0;
	  }
	}
      }
      else {
	if ( Tstrain >= TrotNu ) {
	  Ttangent = (TmomNu-Cstress)/(TrotNu-Cstrain);
	  Tstress = Cstress + Ttangent*dStrain;
	  //printf("\nTNstress5= %f",Tstress);
	  if ( Tstress <= TmomNu ) {
	    Tstress = TmomNu;
	    Ttangent = Ttangent*1.0e-6;
	    //printf("\nTNstress6= %f",Tstress);
	  }
	}
	{
	  Ttangent = (minmom-mom2n)/(TrotMin-TrotPu);
	  Tstress = mom2n + Ttangent*(Tstrain-TrotPu);
	  //printf("\nTNstress7= %f",Tstress);
	  if ( Tstress <= minmom ) {
	    Tstress = minmom;
	    Ttangent = Ttangent*1.0e-6;
	    Tidtime1 = 0;
	    //printf("\nTNstress8= %f",Tstress);
	  }
	}
      }
    }
    else {
      Ttangent = (minmom-Cstress)/(TrotMin-Cstrain);
      Tstress = Cstress + Ttangent*dStrain;
      //printf("\nTNstress9= %f",Tstress);
      if ( Tstress <= minmom ) {
	Tstress = minmom;
	Ttangent = Ttangent*1.0e-6;
	//printf("\nTNstress10= %f",Tstress);
	Tidtime1 = 0;
      }
    }
  }
  TloadIndicator = 2;
  if ( Tidtime == 0 ) {
    if ( TrotMax > rot1p ) {
      TrotMin = (TrotMin < rot1n) ? TrotMin : rot1n;
    }
    Ttangent = E1p*kp;
    Tstress = Cstress + Ttangent*dStrain;
    //printf("\nTNstress11= %f",Tstress);
    if ( Tstress <= mom1n ) {
      TrotPu = Cstrain - Cstress/(E1p*kp);
      // rot1n = rot1n + TrotPu;
      rot1n = mom1n/(E1p*kp) + Cstrain - Cstress/(E1p*kp);
      rot2n = rot1n + (mom2n-mom1n)/E2n;
      rot3n = rot2n + (mom3n-mom2n)/E3n;
      Tidtime = 1;
      //TrotMin = Tstrain;
      // cout << "Tstrain: " << Tstrain << endl;
      // cout << "rot1n: " << rot1n << endl;
      // cout << "rot2n: " << rot2n << endl;
      if (Tstrain > rot1n) {
	Tstress = mom1n;
	Ttangent = E1p*1.0e-9;
	//printf("\nTNstress12= %f",Tstress);
      }
      if (Tstrain <= rot1n && Tstrain < rot2n) {
	Tstress = mom1n + E2n*(Tstrain-rot1n);
	Ttangent = E2n;
	//printf("\nTNstress13= %f",Tstress);
	// cout << "Tstress: " << Tstress << endl;
      }
      if (Tstrain <= rot3n) {
	Tstress = mom2n + E3n*(Tstrain-rot2n);
	Ttangent = E3n;
	//printf("\nTNstress14= %f",Tstress);
	
      }
    }
  }
  if ( Tidtime >= 1 && TrotMin <= rot2n ) {
    Ttangent = E1p*kp;
    Tstress = Cstress + Ttangent*dStrain;
    //printf("\nTNstress15= %f",Tstress);
    if (Tstress <= mom3n) {
      double minmom = negEnvlpStress(TrotMin);
      Ttangent = (minmom-Cstress)/(TrotMin-Cstrain);
      Tstress = Cstress + Ttangent*(Tstrain-Cstrain);
      Tidtime = 2;
      //printf("\nTNstress16= %f",Tstress);
      //if ( Tstrain > rot3n ) {
      //Ttangent = E1p*1.0e-9;
      //Tstress = mom3n;
      //printf("\nTstress17=%f",Tstress);
      //}
      if ( Tstrain < rot3n ) {
	Ttangent = E1p*1.0e-4;
	//Tstress = mom3n;
	//Ttangent = E3n;
	Tstress = mom3n;
	//printf("\nTNstress17= %f",Tstress);
      }
    }
  }
}
int
BraceMaterial::commitState(void)
{
  CrotMax = TrotMax;
  CrotMin = TrotMin;
  CrotPu = TrotPu;
  CrotNu = TrotNu;
  CenergyD = TenergyD;
  CloadIndicator = TloadIndicator;
  Cidtime = Tidtime;
  Cidtime1 = Tidtime1;
  CmomNu = TmomNu;
  CmomPu = TmomPu;
  Cstress = Tstress;
  Cstrain = Tstrain;
  return 0;
}
int
BraceMaterial::revertToLastCommit(void)
{
  TrotMax = CrotMax;
  TrotMin = CrotMin;
  TrotPu = CrotPu;
  TrotNu = CrotNu;
  TenergyD = CenergyD;
  TloadIndicator = CloadIndicator;
  Tidtime = Cidtime;
  Tidtime1 = Cidtime1;
  TmomNu = CmomNu;
  TmomPu = CmomPu;
  Tstress = Cstress;
  Tstrain = Cstrain;
  return 0;
}
int
BraceMaterial::revertToStart(void)
{
  CrotMax = 0.0;
  CrotMin = 0.0;
  CrotPu = 0.0;
  CrotNu = 0.0;
  CenergyD = 0.0;
  CloadIndicator = 0;
  Cidtime = 0;
  Cidtime1 = 0;
  CmomNu = 0.0;
  CmomPu = 0.0;
  Cstress = 0.0;
  Cstrain = 0.0;
  Tstrain = 0;
  Tstress = 0;
  Ttangent = E1p;
  return 0;
}
UniaxialMaterial*
BraceMaterial::getCopy(void)
{
  BraceMaterial *theCopy = new BraceMaterial (this->getTag(),
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
  theCopy->Cidtime = Cidtime;
  theCopy->Cidtime1 = Cidtime1;
  theCopy->CmomNu = CmomNu;
  theCopy->CmomPu = CmomPu;
  return theCopy;
}
int
BraceMaterial::sendSelf(int commitTag, Channel &theChannel)
{
  int res = 0;
static Vector data(31);
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
 data(27) = Cidtime;
 data(28) = Cidtime1;
 data(29) = CmomNu;
 data(30) = CmomPu;
 res = theChannel.sendVector(this->getDbTag(), commitTag, data);
 if (res < 0)
   opserr << "BraceMaterial::sendSelf() - failed to send data\n";
 return res;
}
int
BraceMaterial::recvSelf(int commitTag, Channel &theChannel,
			FEM_ObjectBroker &theBroker)
{
  int res = 0;
  static Vector data(31);
  res = theChannel.recvVector(this->getDbTag(), commitTag, data);
  if (res < 0) {
    opserr << "BraceMaterial::recvSelf() - failed to receive data\n";
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
    Cidtime = int(data(27));
    Cidtime1 = int(data(28));
    CmomNu = data(29);
    CmomPu = data(30);
    // set the trial values
    TrotMax = CrotMax;
    TrotMin = CrotMin;
    TrotPu = CrotPu;
    TrotNu = CrotNu;
    TenergyD = CenergyD;
    TloadIndicator = CloadIndicator;
    Tstress = Cstress;
    Tstrain = Cstrain;
    Tidtime = Cidtime;
    Tidtime1 = Cidtime1;
    TmomNu = CmomNu;
    TmomPu = CmomPu;
  }
  return 0;
}
void
BraceMaterial::Print(OPS_Stream &s, int flag)
{
  s << "Brace Material, tag: " << this->getTag() << endln;
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
BraceMaterial::setEnvelope(void)
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
BraceMaterial::posEnvlpStress(double strain)
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
BraceMaterial::negEnvlpStress(double strain)
{
  if (strain >= 0.0)
    return 0.0;
  else if (strain >= rot1n)
    return E1n*strain;
  else if (strain >= rot2n)
    return mom1n + E2n*(strain-rot1n);
  else if (strain >= rot3n)
    return mom2n + E3n*(strain-rot2n);
  else
    return mom3n;
}
double
BraceMaterial::posEnvlpTangent(double strain)
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
BraceMaterial::negEnvlpTangent(double strain)
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
BraceMaterial::posEnvlpRotlim(double strain)
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
BraceMaterial::negEnvlpRotlim(double strain)
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









