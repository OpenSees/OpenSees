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
 
//This file contains the definition of material Steel02Thermal, which is
// modified from Steel02 by adding functions for considering temperature
// dependent material properties

//Modified by:  Jian Zhang(j.zhang@ed.ac.uk)---------07,2010// 
//              Panagiotis Kotsovinos(P.Kotsovinos@ed.ac.uk)// 
//              Liming Jiang(liming.jiang@ed.ac.uk)



#include <stdlib.h>
#include <Steel02Thermal.h>
#include <OPS_Globals.h>
#include <float.h>
#include <Channel.h>
#include <Information.h>

#include <math.h>

#include <elementAPI.h>
#include <OPS_Globals.h>


void *
OPS_Steel02Thermal()
{
  // Pointer to a uniaxial material that will be returned
  UniaxialMaterial *theMaterial = 0;

  int    iData[1];
  double dData[12];
  int numData = 1;

  if (OPS_GetIntInput(&numData, iData) != 0) {
    opserr << "WARNING invalid uniaxialMaterial Steel02Thermal tag" << endln;
    return 0;
  }

  numData = OPS_GetNumRemainingInputArgs();

  if (numData != 3 && numData != 6 && numData != 10 && numData != 11) {
    opserr << "Invalid #args, want: uniaxialMaterial Steel02Thermal " << iData[0] << 
      " fy? E? b? <R0? cR1? cR2? <a1? a2? a3? a4?>>" << endln;
    return 0;
  }

  if (numData == 3) {
    if (OPS_GetDoubleInput(&numData, dData) != 0) {
      opserr << "Invalid arggs: uniaxialMaterial Steel02Thermal " << iData[0] << 
	" fy? E? b? <R0? cR1? cR2? <a1? a2? a3? a4?>>" << endln;
      return 0;
    }

    // Parsing was successful, allocate the material
    theMaterial = new Steel02Thermal(iData[0], dData[0], dData[1], dData[2]);    
    
  } else if (numData == 6) {
    if (OPS_GetDoubleInput(&numData, dData) != 0) {
      opserr << "Invalid arggs: uniaxialMaterial Steel02Thermal " << iData[0] << 
	" fy? E? b? <R0? cR1? cR2? <a1? a2? a3? a4?>>" << endln;
      return 0;
    }
    
    // Parsing was successful, allocate the material
    theMaterial = new Steel02Thermal(iData[0], dData[0], dData[1], dData[2], dData[3], dData[4], dData[5]);    

  } else if (numData == 10) {
    if (OPS_GetDoubleInput(&numData, dData) != 0) {
      opserr << "Invalid arggs: uniaxialMaterial Steel02Thermal " << iData[0] << 
	" fy? E? b? <R0? cR1? cR2? <a1? a2? a3? a4?>>" << endln;
      return 0;
    }

    // Parsing was successful, allocate the material
    theMaterial = new Steel02Thermal(iData[0], dData[0], dData[1], dData[2], 
				     dData[3], dData[4], dData[5], dData[6], 
				     dData[7], dData[8], dData[9]);    

  } else if (numData == 11) {
    if (OPS_GetDoubleInput(&numData, dData) != 0) {
      opserr << "Invalid arggs: uniaxialMaterial Steel02Thermal " << iData[0] << 
	" fy? E? b? <R0? cR1? cR2? <a1? a2? a3? a4?>>" << endln;
      return 0;
    }

    // Parsing was successful, allocate the material
    theMaterial = new Steel02Thermal(iData[0], dData[0], dData[1], dData[2], 
				     dData[3], dData[4], dData[5], dData[6], 
				     dData[7], dData[8], dData[9], dData[10]);    

  }   

  if (theMaterial == 0) {
    opserr << "WARNING could not create uniaxialMaterial of type Steel02Thermal Material\n";
    return 0;
  }

  return theMaterial;
}



Steel02Thermal::Steel02Thermal(int tag,
			       double _Fy, double _E0, double _b,
			       double _R0, double _cR1, double _cR2,
			       double _a1, double _a2, double _a3, double _a4, double sigInit):
  UniaxialMaterial(tag, MAT_TAG_Steel02Thermal),
  FyT(_Fy), E0T(_E0), b(_b), R0(_R0), cR1(_cR1), cR2(_cR2), a1(_a1), a2(_a2), a3(_a3), a4(_a4), 
  sigini(sigInit)
{
  
  ThermalElongation = 0; //initialize //JZ, 11/10//
  E0 = E0T;//JZ, 11/10//
  Fy = FyT;//JZ, 11/10//
  E0P= E0;  //Liming 2013
  FyP= Fy;  //Liming 2013
  FiberTP = 0; //Liming 2013

  konP = 0;
  kon = 0;
  eP = E0;
  epsP = 0.0;
  sigP = 0.0;
  sig = 0.0;
  eps = 0.0;
  e = E0;

  epsmaxP = Fy/E0;
  epsminP = -epsmaxP;
  epsplP = 0.0;
  epss0P = 0.0;
  sigs0P = 0.0;
  epssrP = 0.0;
  sigsrP = 0.0;

  if (sigini != 0.0) {
	  epsP = sigini/E0;
	  sigP = sigini;
   } 




}

Steel02Thermal::Steel02Thermal(int tag,
		 double _Fy, double _E0, double _b,
		 double _R0, double _cR1, double _cR2):
  UniaxialMaterial(tag, MAT_TAG_Steel02Thermal),
  FyT(_Fy), E0T(_E0), b(_b), R0(_R0), cR1(_cR1), cR2(_cR2), sigini(0.0)
{

	  ThermalElongation = 0; //initialize //JZ, 11/10//
	  E0 = E0T;//JZ, 11/10//
	  Fy = FyT;//JZ, 11/10//
	  E0P= E0;  //Liming 2013
      FyP= Fy;  //liming 2013
      FiberTP = 0; //Liming 2013

  konP = 0;

  // Default values for no isotropic hardening
  a1 = 0.0;
  a2 = 1.0;
  a3 = 0.0;
  a4 = 1.0;

  eP = E0;
  epsP = 0.0;
  sigP = 0.0;
  sig = 0.0;
  eps = 0.0;
  e = E0;

  epsmaxP = Fy/E0;
  epsminP = -epsmaxP;
  epsplP = 0.0;
  epss0P = 0.0;
  sigs0P = 0.0;
  epssrP = 0.0;
  sigsrP = 0.0;

  //ThermalElongation = 1e-6; //initialize
}

Steel02Thermal::Steel02Thermal(int tag, double _Fy, double _E0, double _b):
  UniaxialMaterial(tag, MAT_TAG_Steel02Thermal),
  FyT(_Fy), E0T(_E0), b(_b), sigini(0.0)
{
	  ThermalElongation = 0; //initialize //JZ, 11/10//
	  E0 = E0T;//JZ, 11/10//
	  Fy = FyT;//JZ, 11/10//
      E0P= E0;  //Liming 2013
      FyP= Fy;  //liming 2013
      FiberTP = 0; //Liming 2013
  konP = 0;

  // Default values for elastic to hardening transitions
  R0 = 15.0;
  cR1 = 0.925;
  cR2 = 0.15;

  // Default values for no isotropic hardening
  a1 = 0.0;
  a2 = 1.0;
  a3 = 0.0;
  a4 = 1.0;

  eP = E0;
  epsP = 0.0;
  sigP = 0.0;
  sig = 0.0;
  eps = 0.0;
  e = E0;

  epsmaxP = Fy/E0;
  epsminP = -epsmaxP;
  epsplP = 0.0;
  epss0P = 0.0;
  sigs0P = 0.0;
  epssrP = 0.0;
  sigsrP = 0.0;

 // ThermalElongation = 1e-6; //initialize
}

Steel02Thermal::Steel02Thermal(void):
  UniaxialMaterial(0, MAT_TAG_Steel02Thermal)
{
  konP = 0;
}

Steel02Thermal::~Steel02Thermal(void)
{
  // Does nothing
}

UniaxialMaterial*
Steel02Thermal::getCopy(void)
{
  Steel02Thermal *theCopy = new Steel02Thermal(this->getTag(), Fy, E0, b, R0, cR1, cR2, a1, a2, a3, a4, sigini);
  
  return theCopy;
}

double
Steel02Thermal::getInitialTangent(void)
{
  return E0;
}

int
Steel02Thermal::setTrialStrain(double trialStrain, double FiberTemperature, double strainRate)
{
  double Esh = b * E0;
  double epsy = Fy / E0;


  // modified C-P. Lamarche 2006
  if (sigini != 0.0) {
	  double epsini = sigini/E0;
	  eps = trialStrain+epsini;
   } else
	   eps = trialStrain;
   // modified C-P. Lamarche 2006

  double deps = eps - epsP;
  
// ---------------Initiating for New Thermal-load step, added by Liming
  double epsyP = FyP/E0P;
  if(fabs(epsmaxP-epsyP)<1e-6)
	  epsmaxP = epsy;
  
  if(fabs(epsminP+epsyP)<1e-6)
	  epsminP = -epsy;

   if(fabs(epsplP-epsyP)<1e-6)
	  epsplP = epsy;
  
  if(fabs(epsplP+epsyP)<1e-6)
	  epsplP = -epsy;

   if(fabs(epss0P-epsyP)<1e-6)
	  epss0P = epsy;
  
  if(fabs(epss0P+epsyP)<1e-6)
	  epss0P = -epsy;

   if(fabs(sigs0P-FyP)<1e-6)
	  sigs0P = Fy;
  
   if(fabs(sigs0P+FyP)<1e-6)
	  sigs0P = -Fy;
   //---------

  epsmax = epsmaxP;
  epsmin = epsminP;
  epspl  = epsplP;
  epss0  = epss0P;  
  sigs0  = sigs0P; 
  epsr   = epssrP;  
  sigr   = sigsrP;  
  kon = konP;

  if (kon == 0 || kon == 3) { // modified C-P. Lamarche 2006


    if (fabs(deps) < 10.0*DBL_EPSILON) {

      e = E0;
      sig = sigini;                // modified C-P. Lamarche 2006
      kon = 3;                     // modified C-P. Lamarche 2006 flag to impose initial stess/strain
      return 0;

    } else {

      epsmax = epsy;
      epsmin = -epsy;
      if (deps < 0.0) {
	kon = 2;
	epss0 = epsmin;
	sigs0 = -Fy;
	epspl = epsmin;
      } else {
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
  
  if (kon == 2 && FiberTemperature<FiberTP && deps > 0.0) {


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

  } else if (kon == 1 && FiberTemperature<FiberTP && deps < 0.0) {
    
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
    
    double d1 = (epsmax - epsmin) / (2.0*(a2 * epsy));
    double shft = 1.0 + a1 * pow(d1, 0.8);
    epss0 = (-Fy * shft + Esh * epsy * shft - sigr + E0 * epsr) / (E0 - Esh);
    sigs0 = -Fy * shft + Esh * (epss0 + epsy * shft);
    epspl = epsmin;
  }

  
  // calculate current stress sig and tangent modulus E 

  double xi     = fabs((epspl-epss0)/epsy);
  double R      = R0*(1.0 - (cR1*xi)/(cR2+xi));
  double epsrat = (eps-epsr)/(epss0-epsr);
  double dum1  = 1.0 + pow(fabs(epsrat),R);
  double dum2  = pow(dum1,(1/R));

  sig   = b*epsrat +(1.0-b)*epsrat/dum2;
  sig   = sig*(sigs0-sigr)+sigr;

  e = b + (1.0-b)/(dum1*dum2);
  e = e*(sigs0-sigr)/(epss0-epsr);

  FiberTP=FiberTemperature;
	  
  return 0;
}



double 
Steel02Thermal::getStrain(void)
{
  return eps;
}

double 
Steel02Thermal::getStress(void)
{
  return sig;
}

double 
Steel02Thermal::getTangent(void)
{
  return e;
}



//JZ 07/10 /////////////////////////////////////////////////////////////start
double 
Steel02Thermal::getElongTangent(double TempT, double&ET, double&Elong, double TempTmax) //PK add to include max temp
{
	// EN 1992 pt 1-2-1. Carbon steel at elevated temperatures
  if (TempT <= 80) {
		Fy = FyT;
		E0 = E0T;
  }
  else if (TempT <= 180) {
      Fy = FyT;
      E0 = E0T*(1 - (TempT - 80)*0.1/100);

  }
  else if (TempT <= 280) {
      Fy = FyT;
      E0 = E0T*(0.9 - (TempT - 180)*0.1/100);
  }
  else if (TempT <= 380) {
		Fy = FyT;
		E0 = E0T*(0.8 - (TempT - 280)*0.1/100);
  }
  else if (TempT <= 480) {
      Fy = FyT*(1 - (TempT - 380)*0.22/100);
	  E0 = E0T*(0.7 - (TempT - 380)*0.1/100);
  }
  else if (TempT <= 580) {
      Fy = FyT*(0.78 - (TempT - 480)*0.31/100);
	  E0 = E0T*(0.6 - (TempT - 480)*0.29/100);
  }
  else if (TempT <= 680) {
      Fy = FyT*(0.47 - (TempT - 580)*0.24/100);
	  E0 = E0T*(0.31 - (TempT - 580)*0.18/100);
  }
  else if (TempT <= 780) {
      Fy = FyT*(0.23 - (TempT - 680)*0.12/100);
	  E0 = E0T*(0.13 - (TempT - 680)*0.04/100);
  }
  else if (TempT <= 880) {
      Fy = FyT*(0.11 - (TempT - 780)*0.05/100);
	  E0 = E0T*(0.09 - (TempT - 780)*0.0225/100);
  }
  else if (TempT <= 980) {
      Fy = FyT*(0.06 - (TempT - 880)*0.02/100);
	  E0 = E0T*(0.0675 - (TempT - 880)*0.0225/100);
  }
  else if (TempT <= 1080) {
      Fy = FyT*(0.04 - (TempT - 980)*0.02/100);
	  E0 = E0T*(0.045 - (TempT - 980)*0.0225/100);
  }
  else if (TempT <= 1180) {
      Fy = FyT*(0.02 - (TempT - 1080)*0.02/100);
	  E0 = E0T*(0.0225 - (TempT - 1080)*0.0225/100);
  }
  else  {
      opserr << "the temperature is invalid\n"; 
  } 

  // calculation of thermal elongation of reinforcing steel. JZ
 //opserr<<TempT<<endln;
	  if (TempT <= 1) {
		  ThermalElongation = TempT * 1.2164e-5;
	  }
  else if (TempT <= 730) {
      ThermalElongation = -2.416e-4 + 1.2e-5 *(TempT+20) + 0.4e-8 *(TempT+20)*(TempT+20);
  }
  else if (TempT <= 840) {
      ThermalElongation = 11e-3;
  }
  else if (TempT <= 1180) {
      ThermalElongation = -6.2e-3 + 2e-5*(TempT+20);
  }
  else {
	  opserr << "the temperature is invalid\n";
  }
  

/*

  if (TempT <= 80){
	  E0 = E0T;
  }
E0 = E0T*(1-0.001*(TempT-80));
 if (TempT <= 380){
	Fy = FyT;
  }
 else if (TempT <= 680){
	Fy = FyT*(1-0.8/300*(TempT-380));
  }
else if (TempT <= 980){
	Fy = FyT*(0.2-0.2/300*(TempT-680));
  }
*/
  
  ET = E0;   
  Elong = ThermalElongation;
  return 0;
}
//JZ 07/10 /////////////////////////////////////////////////////////////end 

int 
Steel02Thermal::commitState(void)
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

  E0P = E0;  //Added by liming for tracking the initial stiffnes change
  FyP = Fy;  //Added by liming for tracking the yield stress change
  return 0;
}

int 
Steel02Thermal::revertToLastCommit(void)
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
   E0 = E0P;  //Added by liming for tracking the initial stiffnes change
   Fy = FyP;  //Added by liming for tracking the yield stress change
  return 0;
}

int 
Steel02Thermal::revertToStart(void)
{
  eP = E0;
  epsP = 0.0;
  sigP = 0.0;
  sig = 0.0;
  eps = 0.0;
  e = E0;  

  konP = 0;
  epsmaxP = Fy/E0;
  epsminP = -epsmaxP;
  epsplP = 0.0;
  epss0P = 0.0;
  sigs0P = 0.0;
  epssrP = 0.0;
  sigsrP = 0.0;

  if (sigini != 0.0) {
	  epsP = sigini/E0;
	  sigP = sigini;
   } 

  return 0;
}

int 
Steel02Thermal::sendSelf(int commitTag, Channel &theChannel)
{
  static Vector data(23);
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

  if (theChannel.sendVector(this->getDbTag(), commitTag, data) < 0) {
    opserr << "Steel02Thermal::sendSelf() - failed to sendSelf\n";
    return -1;
  }
  return 0;
}

int 
Steel02Thermal::recvSelf(int commitTag, Channel &theChannel, 
	     FEM_ObjectBroker &theBroker)
{
  static Vector data(23);

  if (theChannel.recvVector(this->getDbTag(), commitTag, data) < 0) {
    opserr << "Steel02Thermal::recvSelf() - failed to recvSelf\n";
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
  konP = data(17);   
  epsP = data(18);   
  sigP = data(19);   
  eP   = data(20);   
  this->setTag(data(21));
  sigini = data(22);

  e = eP;
  sig = sigP;
  eps = epsP;
  
  return 0;
}


void 
Steel02Thermal::Print(OPS_Stream &s, int flag)
{
  s << "Steel02Thermal:(strain, stress, tangent) " << eps << " " << sig << " " << e << endln;
}


//this function is no use, just for the definiation of pure virtual function.
int
Steel02Thermal::setTrialStrain(double strain, double strainRate)
{
  opserr << "Steel02Thermal::setTrialStrain(double strain, strainRate) - should never be called\n";
  return 0;
}


int 
Steel02Thermal::getVariable(const char *variable, Information &info)
{
  if (strcmp(variable,"ThermalElongation") == 0) {
    info.theDouble = ThermalElongation;    
    return 0;
  } else if (strcmp(variable,"ElongTangent") == 0) {
    Vector *theVector = info.theVector;
    if (theVector != 0) {
      double tempT, ET, Elong, TempTmax;
      tempT = (*theVector)(0);
	  ET = (*theVector)(1);
	  Elong = (*theVector)(2);
      TempTmax = (*theVector)(3);
      this->getElongTangent(tempT, ET, Elong, TempTmax);
	  (*theVector)(0) = tempT;
      (*theVector)(1) = ET;
      (*theVector)(2) = Elong;
	  (*theVector)(3) = TempTmax;
    }
    return 0;
  }

  return -1;
}
