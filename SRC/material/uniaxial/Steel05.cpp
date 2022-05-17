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
                                                                        
// $Source: /usr/local/cvs/OpenSees/SRC/material/uniaxial/Steel05.cpp,v $
                                                                      
// Written: S. A. Jalali 10/2019
// Adding Cyclic and in-cycle deterioration modes to steel02 UniaxialMaterial

#include <math.h>

#include <stdlib.h>
#include <Steel05.h>
#include <float.h>
#include <Channel.h>
#include <Information.h>
#include <Parameter.h>

#include <elementAPI.h>
#include <OPS_Globals.h>
#include "classTags.h"

static int numThisCall = 0;

void *
OPS_Steel05()
{
	if (numThisCall == 0) {
		opserr << "------ Steel05 unaxialMaterial, Written by SAJalali @ Civil Soft Science, Iran, 2019-------\n";
			numThisCall = 1;
	}
  // Pointer to a uniaxial material that will be returned
  UniaxialMaterial *theMaterial = 0;

  int    iData[1];
  double dData[16];
  int numData = 1;

  if (OPS_GetIntInput(&numData, iData) != 0) {
    opserr << "WARNING invalid uniaxialMaterial Steel05 tag" << endln;
    return 0;
  }

  numData = OPS_GetNumRemainingInputArgs();

  if (numData != 8 && numData != 11 && numData != 15 && numData != 16) {
    opserr << "Invalid Steel05 #args for: " << iData[0] << " see the syntax" << endln;
	 opserr << "-------Syntax:\n";
	 opserr << "-------UniaxialMaterial Steel05 $matTag $Fy $E $b ductilityCapacity postCapEFac gama c resFac <$R0 $cR1 $cR2> <$a1 $a2 $a3 $a4> <$sigInit>\n\n";
	 return 0;
  }
//default parameters:
	dData[15] = 0;
	dData[14] = 1;
	dData[13] = 0;
	dData[12] = 1;
	dData[11] = 0;
	dData[10] = 0.15;
	dData[9] = 0.925;
	dData[8] = 15;

    if (OPS_GetDoubleInput(&numData, dData) != 0) {
    opserr << "Invalid Steel05 #args for: " << iData[0] << " see the syntax" << endln;
	 opserr << "-------Syntax:\n";
	 opserr << "-------UniaxialMaterial Steel05 $matTag $Fy $E $b ductilityCapacity postCapEFac gama c resFac <$R0 $cR1 $cR2> <$a1 $a2 $a3 $a4> <$sigInit>\n\n";
	 return 0;
    }

    // Parsing was successful, allocate the material
    theMaterial = new Steel05(iData[0], dData[0], dData[1], dData[2], 
			      dData[3], dData[4], dData[5], dData[6], 
			      dData[7], dData[8], dData[9], dData[10],    
			      dData[11], dData[12], dData[13], dData[14], dData[15]);   


  if (theMaterial == 0) {
    opserr << "WARNING could not create uniaxialMaterial of type Steel05 Material\n";
    return 0;
  }

  return theMaterial;
}


Steel05::Steel05(int tag,
	 double _Fy, double _E0, double _b,
	 double ductC, double pcEFac, double _gama, double _c, double _resFac,
	 double _R0, double _cR1, double _cR2,
	 double _a1, double _a2, double _a3, double _a4, double sigInit) :
	 UniaxialMaterial(tag, MAT_TAG_Steel05),
	 Fy(_Fy), E0(_E0), b(_b), epsPCFac(ductC), pstcpEFac(pcEFac), gama(_gama), c(_c), resFac(_resFac),
	 R0(_R0), cR1(_cR1), cR2(_cR2), a1(_a1), a2(_a2), a3(_a3), a4(_a4), sigini(sigInit)
{
	 EnergyP = 0;
	 eP = E0;
	 epsP = sigini / E0;
	 sigP = sigini;
	 sig = 0.0;
	 eps = 0.0;
	 e = E0;

	 konP = 0;
	 epsmaxP = Fy / E0;
	 epsminP = -epsmaxP;
	 epsplP = 0.0;
	 epss0P = 0.0;
	 sigs0P = 0.0;
	 epssrP = 0.0;
	 sigsrP = 0.0;
	 FailEnerg = gama * Fy * Fy / E0;
	 ExcurEnergy = 0;
	 FydP = FydN = Fy;
}

Steel05::Steel05(void):
  UniaxialMaterial(0, MAT_TAG_Steel05)
{
	 EnergyP = 0;
	 eP = E0;
	 epsP = 0;
	 sigP = 0;
	 sig = 0.0;
	 eps = 0.0;
	 e = 0;

	 konP = 0;
	 epsmaxP = 0;
	 epsminP = 0;
	 epsplP = 0.0;
	 epss0P = 0.0;
	 sigs0P = 0.0;
	 epssrP = 0.0;
	 sigsrP = 0.0;

	 FailEnerg = gama = 0;
	 ExcurEnergy = 0;
	 FydP = FydN = 0;
}

Steel05::~Steel05(void)
{
  // Does nothing
}

UniaxialMaterial*
Steel05::getCopy(void)
{
  Steel05 *theCopy = new Steel05(this->getTag(), Fy, E0, b, epsPCFac, pstcpEFac, gama, c, resFac,
		R0, cR1, cR2, a1, a2, a3, a4, sigini);
  
  return theCopy;
}

double
Steel05::getInitialTangent(void)
{
  return E0;
}

int
Steel05::setTrialStrain(double trialStrain, double strainRate)
{
	 double Esh = b * E0;
	 const double epsy = Fy / E0;
	 eps = trialStrain + sigini / E0;
	 const double deps = eps - epsP;
	 double b2 = b;
	 const double epsPC = epsPCFac * epsy;

	 epsmax = epsmaxP;
	 epsmin = epsminP;
	 epspl = epsplP;
	 epss0 = epss0P;
	 sigs0 = sigs0P;
	 epsr = epssrP;
	 sigr = sigsrP;
	 kon = konP;

	 if (kon == 0) {
		  epsmax = epsy;
		  epsmin = -epsy;
		  if (deps < 0.0) {
				kon = 2;
				epss0 = epsmin;
				sigs0 = -FydN;
				epspl = epsmin;
		  }
		  else
		  {
				kon = 1;
				epss0 = epsmax;
				sigs0 = FydP;
				epspl = epsmax;
		  }
	 }

	 if (kon == 2) {
		  if (deps > 0.0)
		  {
				// in case of load reversal from negative to positive strain increment, 
				// update the minimum previous strain, store the last load reversal 
				// point and calculate the stress and strain (sigs0 and epss0) at the 
				// new intersection between elastic and strain hardening asymptote 
				// To include isotropic strain hardening shift the strain hardening 
				// asymptote by sigsft before calculating the intersection point 
				// Constants a3 and a4 control this stress shift on the tension side 
				kon = 1;
				epsr = epsP;
				sigr = sigP;
				if (epsP < epsmin)
					 epsmin = epsP;
				const double d1 = (epsmax - epsmin) / (2.0 * (a4 * epsy));
				const double shft = 1.0 + a3 * pow(d1, 0.8);
				epss0 = (FydP * shft - Esh * epsy * shft - sigr + E0 * epsr) / (E0 - Esh);
				sigs0 = FydP * shft + Esh * (epss0 - epsy * shft);
				epspl = epsmax;
		  }
		  else if (eps < -(epsPCFac - 1) / 2 * epsy)
		  {
				epss0 = -epsPC;
				sigs0 = -FydN - b * E0 * (epsPC - epsy);
				epsr = epsP;
				sigr = sigP;
				kon = 4;
		  }
	 }
	 else if (kon == 1) {
		  // update the maximum previous strain, store the last load reversal 
		  // point and calculate the stress and strain (sigs0 and epss0) at the 
		  // new intersection between elastic and strain hardening asymptote 
		  // To include isotropic strain hardening shift the strain hardening 
		  // asymptote by sigsft before calculating the intersection point 
		  // Constants a1 and a2 control this stress shift on compression side 
		  if (deps < 0.0)
		  {
				kon = 2;
				epsr = epsP;
				sigr = sigP;
				if (epsP > epsmax)
					 epsmax = epsP;

				const double d1 = (epsmax - epsmin) / (2.0 * (a2 * epsy));
				const double shft = 1.0 + a1 * pow(d1, 0.8);
				epss0 = (-FydN * shft + Esh * epsy * shft - sigr + E0 * epsr) / (E0 - Esh);
				sigs0 = -FydN * shft + Esh * (epss0 + epsy * shft);
				epspl = epsmin;
		  }
		  else if (eps > (epsPCFac - 1) / 2 * epsy)
		  {
				epss0 = epsPC;
				sigs0 = FydP + b * E0 * (epsPC - epsy);
				epsr = epsP;
				sigr = sigP;
				kon = 3;
		  }
	 }
	 if (kon == 4)
	 {
		  if (deps > 0.0)
		  {
				kon = 1;
				epsr = epsP;
				sigr = sigP;
				if (epsP < epsmin)
					 epsmin = epsP;
				const double d1 = (epsmax - epsmin) / (2.0 * (a4 * epsy));
				const double shft = 1.0 + a3 * pow(d1, 0.8);
				epss0 = (FydP * shft - Esh * epsy * shft - sigr + E0 * epsr) / (E0 - Esh);
				sigs0 = FydP * shft + Esh * (epss0 - epsy * shft);
				epspl = epsmax;
		  }
		  else
				b2 = pstcpEFac / b;
		 
	 }
	 if (kon == 3)
	 {
		  if (deps < 0.0)
		  {
				kon = 2;
				epsr = epsP;
				sigr = sigP;
				if (epsP > epsmax)
					 epsmax = epsP;

				const double d1 = (epsmax - epsmin) / (2.0 * (a2 * epsy));
				const double shft = 1.0 + a1 * pow(d1, 0.8);
				epss0 = (-FydN * shft + Esh * epsy * shft - sigr + E0 * epsr) / (E0 - Esh);
				sigs0 = -FydN * shft + Esh * (epss0 + epsy * shft);
				epspl = epsmin;
		  }
		  else
				b2 = pstcpEFac / b;
		 
	 }
	 double R = R0;
	 // calculate current stress sig and tangent modulus E 
	 if (cR1 != 0.0 && cR2 != 0.0)
	 {
		  const double xi = fabs((epspl - epss0) / epsy);
		  R *= (1.0 - (cR1 * xi) / (cR2 + xi));
	 }
	 const double epsrat = (eps - epsr) / (epss0 - epsr);
	 const double dum1 = 1.0 + pow(fabs(epsrat), R);
	 const double dum2 = pow(dum1, (1 / R));

	 sig = b2 * epsrat + (1.0 - b2) * epsrat / dum2;
	 sig = sig * (sigs0 - sigr) + sigr;

	 e = b2 + (1.0 - b2) / (dum1 * dum2);
	 e *= (sigs0 - sigr) / (epss0 - epsr);
	 if (kon == 3 && sig < resFac * Fy)
	 {
		  sig = resFac * Fy;
		  e = 1e-15 * E0;
	 }
	 else if (kon == 4 && sig > -resFac * Fy)
	 {
		  sig = -resFac * Fy;
		  e = 1e-15 * E0;
	 }
	 return 0;
}

double 
Steel05::getStrain(void)
{
  return eps;
}

double 
Steel05::getStress(void)
{
  return sig;
}

double 
Steel05::getTangent(void)
{
  return e;
}

int 
Steel05::commitState(void)
{
  epsminP = epsmin;
  epsmaxP = epsmax;
  epsplP = epspl;
  epss0P = epss0;
  sigs0P = sigs0;
  epssrP = epsr;
  sigsrP = sigr;

  updateDamage();
  konP = kon;
  eP = e;
  sigP = sig;
  epsP = eps;

  return 0;
}

int 
Steel05::revertToLastCommit(void)
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
  return 0;
}

int 
Steel05::revertToStart(void)
{
	 EnergyP = 0;	//by SAJalali
	 eP = E0;
	 epsP = sigini / E0;
	 sigP = sigini;
	 sig = 0.0;
	 eps = 0.0;
	 e = E0;

	 konP = 0;
	 epsmaxP = Fy / E0;
	 epsminP = -epsmaxP;
	 epsplP = 0.0;
	 epss0P = 0.0;
	 sigs0P = 0.0;
	 epssrP = 0.0;
	 sigsrP = 0.0;


	 ExcurEnergy = 0;
	 FydP = FydN = Fy;
	 return 0;
}

int 
Steel05::sendSelf(int commitTag, Channel &theChannel)
{
  static Vector data(33);
  data(0) = Fy;
  data(1) = E0;
  data(2) = b;
  data(3) = R0;
  data(4) = cR1;
  data(5) = cR2;
  data(6) = a1;
  data(7) = a2;
  data(8) = a3;
  data(9) = a4;
  data(10) = epsminP;
  data(11) = epsmaxP;
  data(12) = epsplP;
  data(13) = epss0P;
  data(14) = sigs0P;
  data(15) = epssrP;
  data(16) = sigsrP;
  data(17) = konP;  
  data(18) = epsP;  
  data(19) = sigP;  
  data(20) = eP;    
  data(21) = this->getTag();
  data(22) = sigini;
  data(23) = EnergyP;
  data(24) = epsPCFac;
  data(25) = pstcpEFac;
  data(26) = gama;
  data(27) = c;
  data(28) = resFac;
  data(29) = FydP;
  data(30) = FydN;
  data(31) = ExcurEnergy;
  data(32) = gama;
  if (theChannel.sendVector(this->getDbTag(), commitTag, data) < 0) {
    opserr << "Steel05::sendSelf() - failed to sendSelf\n";
    return -1;
  }
  return 0;
}

int 
Steel05::recvSelf(int commitTag, Channel& theChannel,
	 FEM_ObjectBroker& theBroker)
{
	 static Vector data(33);	//editted by SAJalali for energy

	 if (theChannel.recvVector(this->getDbTag(), commitTag, data) < 0) {
		  opserr << "Steel05::recvSelf() - failed to recvSelf\n";
		  return -1;
	 }

	 Fy = data(0);
	 E0 = data(1);
	 b = data(2);
	 R0 = data(3);
	 cR1 = data(4);
	 cR2 = data(5);
	 a1 = data(6);
	 a2 = data(7);
	 a3 = data(8);
	 a4 = data(9);
	 epsminP = data(10);
	 epsmaxP = data(11);
	 epsplP = data(12);
	 epss0P = data(13);
	 sigs0P = data(14);
	 epssrP = data(15);
	 sigsrP = data(16);
	 konP = int(data(17));
	 epsP = data(18);
	 sigP = data(19);
	 eP = data(20);
	 this->setTag(int(data(21)));
	 sigini = data(22);
	 EnergyP = data(23);
	 epsPCFac = data(24);
	 pstcpEFac = data(25);
	 gama = data(26);
	 c = data(27);
	 resFac = data(28);
	 FydP = data(29);
	 FydN = data(30);
	 ExcurEnergy = data(31);
	 gama = data(32);
	 FailEnerg = gama * Fy * Fy / E0;

	 e = eP;
	 sig = sigP;
	 eps = epsP;

	 return 0;
}

void 
Steel05::Print(OPS_Stream &s, int flag)
{
  if (flag == OPS_PRINT_PRINTMODEL_MATERIAL) {      
    //    s << "Steel05:(strain, stress, tangent) " << eps << " " << sig << " " << e << endln;
    s << "Steel05 tag: " << this->getTag() << endln;
    s << "  fy: " << Fy << ", ";
    s << "  E0: " << E0 << ", ";
    s << "   b: " << b << ", ";
    s << "  R0: " << R0 << ", ";
    s << " cR1: " << cR1 << ", ";
    s << " cR2: " << cR2 << ", ";    
    s << "  a1: " << a1 << ", ";
    s << "  a2: " << a2 << ", ";
    s << "  a3: " << a3 << ", ";
    s << "  a4: " << a4;    
  }
  
  if (flag == OPS_PRINT_PRINTMODEL_JSON) {
    s << "\t\t\t{";
	s << "\"name\": \"" << this->getTag() << "\", ";
	s << "\"type\": \"Steel05\", ";
	s << "\"E\": " << E0 << ", ";
	s << "\"fy\": " << Fy << ", ";
    s << "\"b\": " << b << ", ";
    s << "\"R0\": " << R0 << ", ";
    s << "\"cR1\": " << cR1 << ", ";
    s << "\"cR2\": " << cR2 << ", ";
    s << "\"a1\": " << a1 << ", ";
    s << "\"a2\": " << a2 << ", ";
    s << "\"a3\": " << a3 << ", ";
    s << "\"a4\": " << a4 << ", ";    
    s << "\"sigini\": " << sigini << "}";
  }
}

// AddingSensitivity:BEGIN ///////////////////////////////////
int
Steel05::setParameter(const char **argv, int argc, Parameter &param)
{

  if (strcmp(argv[0],"sigmaY") == 0 || strcmp(argv[0],"fy") == 0 || strcmp(argv[0],"Fy") == 0) {
    param.setValue(Fy);
    return param.addObject(1, this);
  }
  if (strcmp(argv[0],"E") == 0) {
    param.setValue(E0);
    return param.addObject(2, this);
  }
  if (strcmp(argv[0],"b") == 0) {
    param.setValue(b);
    return param.addObject(3, this);
  }
  if (strcmp(argv[0],"a1") == 0) {
    param.setValue(a1);
    return param.addObject(4, this);
  }
  if (strcmp(argv[0],"a2") == 0) {
    param.setValue(a2);
    return param.addObject(5, this);
  }
  if (strcmp(argv[0],"a3") == 0) {
    param.setValue(a3);
    return param.addObject(6, this);
  }
  if (strcmp(argv[0],"a4") == 0) {
    param.setValue(a4);
    return param.addObject(7, this);
  }

  return -1;
}



int
Steel05::updateParameter(int parameterID, Information &info)
{
  switch (parameterID) {
  case -1:
    return -1;
  case 1:
    this->Fy = info.theDouble;
    break;
  case 2:
    this->E0 = info.theDouble;
    break;
  case 3:
    this->b = info.theDouble;
    break;
  case 4:
    this->a1 = info.theDouble;
    break;
  case 5:
    this->a2 = info.theDouble;
    break;
  case 6:
    this->a3 = info.theDouble;
    break;
  case 7:
    this->a4 = info.theDouble;
    break;
  default:
    return -1;
  }
  
  return 0;
}


void Steel05::updateDamage()
{
	if (((konP == 1 || konP == 3) && sig < 0) || ((konP == 2 || konP == 4) && sig > 0))
	{
		//submit Pos damage and reset for new excursion
		double zeroSigEps = epsP - sigP/E0;
		double dE = 0.5 * sigP * (zeroSigEps - epsP);
		EnergyP += dE;
		if (EnergyP < 0) EnergyP = 0.;
		if (gama > 9999)
		{
			return;
		}
		double& Fyd = (konP == 2 || konP == 4) ? FydP : FydN;
		//double& Fyd = FydP;
		ExcurEnergy += dE;
		if (ExcurEnergy < 0) ExcurEnergy = 0.;
		double beta = pow( ExcurEnergy / ( FailEnerg - EnergyP) , c );
		if (beta > 0.999 || beta < 0)
		{
			opserr<< "\nSteel05:"<< this->getTag()<< " WARNING! Maximum Energy Absorbance Capacity Reached\n"<< endln;
			beta = 0.999;
		}
		Fyd = (1. - beta)*Fyd + beta * resFac*Fyd;
		//FydN = Fyd;
		ExcurEnergy = 0.0;
	} 
	else
	{
		double dE = 0.5 * (sig + sigP) * (eps - epsP);
		ExcurEnergy += dE;
		EnergyP += dE;
    }
}

