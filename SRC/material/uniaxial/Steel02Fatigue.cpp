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
                                                                        
// $Revision: 1.5 $
// $Date: 2010-09-16 00:04:05 $
// $Source: /usr/local/cvs/OpenSees/SRC/material/uniaxial/Steel02.cpp,v $
                                                                      
// Written: fmk
// Created: 03/06
//
// Description: This file contains the class implementation of Steel02. 
// This Steel02 is based on an f2c of the FEDEAS material
// Steel02.f which is:
//-----------------------------------------------------------------------
// MENEGOTTO-PINTO STEEL MODEL WITH FILIPPOU ISOTROPIC HARDENING
//            written by MOHD YASSIN (1993)
//          adapted to FEDEAS material library
//    by E. Spacone, G. Monti and F.C. Filippou (1994)
//    Modified by Nasser A. Marafi (2018) to have strength and 
//   stiffness deterioration as per Kunnath et al. (2009)
//-----------------------------------------------------------------------

#include <math.h>

#include <stdlib.h>
#include <Steel02Fatigue.h>
#include <OPS_Globals.h>
#include <float.h>
#include <Channel.h>


#include <elementAPI.h>
#include <OPS_Globals.h>


void * OPS_ADD_RUNTIME_VPV(OPS_Steel02Fatigue)
{
  // Pointer to a uniaxial material that will be returned
  UniaxialMaterial *theMaterial = 0;

  int    iData[1];
  double dData[18];
  int numData = 1;

  if (OPS_GetIntInput(&numData, iData) != 0) {
    opserr << "WARNING invalid uniaxialMaterial Steel02Fatigue tag" << endln;
    return 0;
  }

  numData = OPS_GetNumRemainingInputArgs();

  if (numData != 9 && numData != 12 && numData != 16 && numData != 17) {
    opserr << "Invalid #args, want: uniaxialMaterial Steel02Fatigue " << iData[0] << 
      " fy? E? b? Cd? Cf? alpha? beta? minStrain? maxStrain? <R0? cR1? cR2? <a1? a2? a3? a4?>>" << endln;
    return 0;
  }

  if (numData == 9) {
    if (OPS_GetDoubleInput(&numData, dData) != 0) {
      opserr << "Invalid double: uniaxialMaterial Steel02Fatigue " << iData[0] << 
	" fy? E? b? Cd? Cf? alpha? beta? minStrain? maxStrain? <R0? cR1? cR2? <a1? a2? a3? a4?>>" << endln;
      return 0;
    }

    // Parsing was successful, allocate the material
    theMaterial = new Steel02Fatigue(iData[0], dData[0], dData[1], dData[2], dData[3], dData[4], dData[5], dData[6], dData[7], dData[8]);

  } else if (numData == 12) {
    if (OPS_GetDoubleInput(&numData, dData) != 0) {
      opserr << "Invalid int: uniaxialMaterial Steel02Fatigue " << iData[0] << 
	" fy? E? b? Cd? Cf? alpha? beta? minStrain? maxStrain? <R0? cR1? cR2? <a1? a2? a3? a4?>>" << endln;
      return 0;
    }

    // Parsing was successful, allocate the material
    theMaterial = new Steel02Fatigue(iData[0], dData[0], dData[1], dData[2], dData[3], dData[4], dData[5], dData[6], dData[7], dData[8], dData[9], dData[10], dData[11]);

  } else if (numData == 16) {
    if (OPS_GetDoubleInput(&numData, dData) != 0) {
      opserr << "Invalid arggs: uniaxialMaterial Steel02Fatigue " << iData[0] << 
	" fy? E? b? Cd? Cf? alpha? beta? minStrain? maxStrain? <R0? cR1? cR2? <a1? a2? a3? a4?>>" << endln;
      return 0;
    }

    // Parsing was successful, allocate the material
    theMaterial = new Steel02Fatigue(iData[0], dData[0], dData[1], dData[2],
			      dData[3], dData[4], dData[5], dData[6], 
			      dData[7], dData[8], dData[9], dData[10], dData[11], dData[12], dData[13], dData[14], dData[15]);

  } else if (numData == 17) {
    if (OPS_GetDoubleInput(&numData, dData) != 0) {
      opserr << "Invalid arggs: uniaxialMaterial Steel02Fatigue " << iData[0] << 
	" fy? E? b? Cd? Cf? alpha? beta? minStrain? maxStrain? <R0? cR1? cR2? <a1? a2? a3? a4?>>" << endln;
      return 0;
    }

    // Parsing was successful, allocate the material
    theMaterial = new Steel02Fatigue(iData[0], dData[0], dData[1], dData[2],
			      dData[3], dData[4], dData[5], dData[6], 
			      dData[7], dData[8], dData[9], dData[10], dData[11], dData[12], dData[13], dData[14], dData[15], dData[16]);

  }   

  if (theMaterial == 0) {
    opserr << "WARNING could not create uniaxialMaterial of type Steel02Fatigue Material\n";
    return 0;
  }

  return theMaterial;
}


Steel02Fatigue::Steel02Fatigue(int tag,
		 double _Fy, double _E0, double _b,
		 double _R0, double _cR1, double _cR2,
		 double _a1, double _a2, double _a3, double _a4,
	double _Cd, double _Cf , double _Alpha, double _Beta , double _minStrain, double _maxStrain, double sigInit ):
  UniaxialMaterial(tag, MAT_TAG_Steel02Fatigue),
  Fy(_Fy), E0(_E0), b(_b), R0(_R0), cR1(_cR1), cR2(_cR2), a1(_a1), a2(_a2), a3(_a3), a4(_a4), 
   Cd(_Cd), Cf(_Cf), Alpha(_Alpha), Beta(_Beta), minStrain(_minStrain), maxStrain(_maxStrain), sigini(sigInit)
{
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
	  epsP = sigini / E0;
	  sigP = sigini;
  }

  // Fatigue Parameters
  //Cf = Cf; // 0.16;
  //Cd = Cd; // 0.35;
  //Alpha = Alpha; //0.44;
  //Beta = Beta; //0.44;

  Fatigue_Cfailed = false;
  Fatigue_CfailedP = false;
  Fatigue_DI = 0; //Damage index
  Fatigue_X = 0; //Range in consideration
  Fatigue_Y = 0; //Previous Adjacent Range
  Fatigue_A = 0; //Peak or valley 1
  Fatigue_B = 0; //Peak or valley 2
  Fatigue_C = 0; //Peak or valley 2
  Fatigue_D = 0; //Peak or valley 4
  Fatigue_PCC = 0; /*Previous Cycle counter flag if >1 then previous 'n'
		   cycles did not flag a complete cycle */
  Fatigue_R1F = 0; //Flag for first  peak count
  Fatigue_R2F = 0; //Flag for second peak count
  Fatigue_cSlope = 0; //Current Slope
  Fatigue_PS = 0; //Previous slope
  Fatigue_EP = 0; //Previous Strain
  Fatigue_SF = 0; /*Start Flag = 0 if very first strain, (i.e. when initializing)
		  = 1 otherwise */
  Fatigue_DL = 0; //Damage if current strain was last peak.

  // added 6/9/2006
  Fatigue_SR1 = 0;
  Fatigue_NC1 = 0;
  Fatigue_SR2 = 0;
  Fatigue_NC2 = 0;
  Fatigue_SR3 = 0;
  Fatigue_NC3 = 0;

  Fatigue_Dmax = 1.0; 
  Fatigue_E0 = Cf; // 0.191; // Marafi Change 2018/01/31
  Fatigue_m = -1 * Alpha; //-0.458; // Marafi Change 2018/01/31
  Fatigue_epsmin = minStrain;
  Fatigue_epsmax = maxStrain;

  //Previous States
  Fatigue_DIP = 0; //Damage index
  Fatigue_XP = 0; //Range in consideration
  Fatigue_YP = 0; //Previous Adjacent Range
  Fatigue_AP = 0; //Peak or valley 1
  Fatigue_BP = 0; //Peak or valley 2
  Fatigue_CP = 0; //Peak or valley 2
  Fatigue_DP = 0; //Peak or valley 4
  Fatigue_PCCP = 0; /*Previous Cycle counter flag if >1 then previous 'n'
					cycles did not flag a complete cycle */
  Fatigue_R1FP = 0; //Flag for first  peak count
  Fatigue_R2FP = 0; //Flag for second peak count
  Fatigue_cSlopeP = 0; //Current Slope
  Fatigue_PSP = 0; //Previous slope
  Fatigue_EPP = 0; //Previous Strain
  Fatigue_SFP = 0; /*Start Flag = 0 if very first strain, (i.e. when initializing)
				   = 1 otherwise */
  Fatigue_DLP = 0; //Damage if current strain was last peak.

				   // added 6/9/2006
  Fatigue_SR1P = 0;
  Fatigue_NC1P = 0;
  Fatigue_SR2P = 0;
  Fatigue_NC2P = 0;
  Fatigue_SR3P = 0;
  Fatigue_NC3P = 0;

  Fatigue_DmaxP = Fatigue_Dmax;
  Fatigue_E0P = Fatigue_E0;
  Fatigue_mP = Fatigue_m;

  Zd = pow((Cf / Cd), (1. / Alpha));

  Lambda_SR = Zd * Fatigue_DI;
  Lambda_SRP = Lambda_SR;

  Fatigue_FyInitial = Fy;
  Fatigue_Fy = Fy;
  Fatigue_FyP = Fy;

}

Steel02Fatigue::Steel02Fatigue(int tag,
		 double _Fy, double _E0, double _b,
		 double _R0, double _cR1, double _cR2,
	double _Cd, double _Cf, double _Alpha, double _Beta, double _minStrain, double _maxStrain) :
  UniaxialMaterial(tag, MAT_TAG_Steel02Fatigue),
  Fy(_Fy), E0(_E0), b(_b), R0(_R0), cR1(_cR1), cR2(_cR2), Cd(_Cd), Cf(_Cf), Alpha(_Alpha), Beta(_Beta), minStrain(_minStrain), maxStrain(_maxStrain),  sigini(0.0)
{
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

  // Fatigue Parameters
  //Cf = Cf; // 0.16;
  //Cd = Cd; // 0.35;
  //Alpha = Alpha; //0.44;
  //Beta = Beta; //0.44;

  Fatigue_Cfailed = false;
  Fatigue_CfailedP = false;
  Fatigue_DI = 0; //Damage index
  Fatigue_X = 0; //Range in consideration
  Fatigue_Y = 0; //Previous Adjacent Range
  Fatigue_A = 0; //Peak or valley 1
  Fatigue_B = 0; //Peak or valley 2
  Fatigue_C = 0; //Peak or valley 2
  Fatigue_D = 0; //Peak or valley 4
  Fatigue_PCC = 0; /*Previous Cycle counter flag if >1 then previous 'n'
				   cycles did not flag a complete cycle */
  Fatigue_R1F = 0; //Flag for first  peak count
  Fatigue_R2F = 0; //Flag for second peak count
  Fatigue_cSlope = 0; //Current Slope
  Fatigue_PS = 0; //Previous slope
  Fatigue_EP = 0; //Previous Strain
  Fatigue_SF = 0; /*Start Flag = 0 if very first strain, (i.e. when initializing)
				  = 1 otherwise */
  Fatigue_DL = 0; //Damage if current strain was last peak.

				  // added 6/9/2006
  Fatigue_SR1 = 0;
  Fatigue_NC1 = 0;
  Fatigue_SR2 = 0;
  Fatigue_NC2 = 0;
  Fatigue_SR3 = 0;
  Fatigue_NC3 = 0;

  Fatigue_Dmax = 1.0;
  Fatigue_E0 = Cf; // 0.191; // Marafi Change 2018/01/31
  Fatigue_m = -1 * Alpha; //-0.458; // Marafi Change 2018/01/31
  Fatigue_epsmin = minStrain;
  Fatigue_epsmax = maxStrain;

  //Previous States
  Fatigue_DIP = 0; //Damage index
  Fatigue_XP = 0; //Range in consideration
  Fatigue_YP = 0; //Previous Adjacent Range
  Fatigue_AP = 0; //Peak or valley 1
  Fatigue_BP = 0; //Peak or valley 2
  Fatigue_CP = 0; //Peak or valley 2
  Fatigue_DP = 0; //Peak or valley 4
  Fatigue_PCCP = 0; /*Previous Cycle counter flag if >1 then previous 'n'
					cycles did not flag a complete cycle */
  Fatigue_R1FP = 0; //Flag for first  peak count
  Fatigue_R2FP = 0; //Flag for second peak count
  Fatigue_cSlopeP = 0; //Current Slope
  Fatigue_PSP = 0; //Previous slope
  Fatigue_EPP = 0; //Previous Strain
  Fatigue_SFP = 0; /*Start Flag = 0 if very first strain, (i.e. when initializing)
				   = 1 otherwise */
  Fatigue_DLP = 0; //Damage if current strain was last peak.

				   // added 6/9/2006
  Fatigue_SR1P = 0;
  Fatigue_NC1P = 0;
  Fatigue_SR2P = 0;
  Fatigue_NC2P = 0;
  Fatigue_SR3P = 0;
  Fatigue_NC3P = 0;

  Fatigue_DmaxP = Fatigue_Dmax;
  Fatigue_E0P = Fatigue_E0;
  Fatigue_mP = Fatigue_m;

  Zd = pow((Cf / Cd), (1. / Alpha));

  Lambda_SR = Zd * Fatigue_DI;
  Lambda_SRP = Lambda_SR;

  Fatigue_FyInitial = Fy;
  Fatigue_Fy = Fy;
  Fatigue_FyP = Fy;

}

Steel02Fatigue::Steel02Fatigue(int tag, double _Fy, double _E0, double _b,
	double _Cd, double _Cf, double _Alpha, double _Beta, double _minStrain, double _maxStrain) :
  UniaxialMaterial(tag, MAT_TAG_Steel02Fatigue),
  Fy(_Fy), E0(_E0), b(_b), Cd(_Cd), Cf(_Cf), Alpha(_Alpha), Beta(_Beta), minStrain(_minStrain), maxStrain(_maxStrain), sigini(0.0)
{
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

  // Fatigue Parameters
  //Cf = Cf; // 0.16;
  //Cd = Cd; // 0.35;
  //Alpha = Alpha; //0.44;
  //Beta = Beta; //0.44;

  Fatigue_Cfailed = false;
  Fatigue_CfailedP = false;
  Fatigue_DI = 0; //Damage index
  Fatigue_X = 0; //Range in consideration
  Fatigue_Y = 0; //Previous Adjacent Range
  Fatigue_A = 0; //Peak or valley 1
  Fatigue_B = 0; //Peak or valley 2
  Fatigue_C = 0; //Peak or valley 2
  Fatigue_D = 0; //Peak or valley 4
  Fatigue_PCC = 0; /*Previous Cycle counter flag if >1 then previous 'n'
				   cycles did not flag a complete cycle */
  Fatigue_R1F = 0; //Flag for first  peak count
  Fatigue_R2F = 0; //Flag for second peak count
  Fatigue_cSlope = 0; //Current Slope
  Fatigue_PS = 0; //Previous slope
  Fatigue_EP = 0; //Previous Strain
  Fatigue_SF = 0; /*Start Flag = 0 if very first strain, (i.e. when initializing)
				  = 1 otherwise */
  Fatigue_DL = 0; //Damage if current strain was last peak.

				  // added 6/9/2006
  Fatigue_SR1 = 0;
  Fatigue_NC1 = 0;
  Fatigue_SR2 = 0;
  Fatigue_NC2 = 0;
  Fatigue_SR3 = 0;
  Fatigue_NC3 = 0;

  Fatigue_Dmax = 1.0;
  Fatigue_E0 = Cf; // 0.191; // Marafi Change 2018/01/31
  Fatigue_m = -1 * Alpha; //-0.458; // Marafi Change 2018/01/31
  Fatigue_epsmin = minStrain;
  Fatigue_epsmax = maxStrain;

  //Previous States
  Fatigue_DIP = 0; //Damage index
  Fatigue_XP = 0; //Range in consideration
  Fatigue_YP = 0; //Previous Adjacent Range
  Fatigue_AP = 0; //Peak or valley 1
  Fatigue_BP = 0; //Peak or valley 2
  Fatigue_CP = 0; //Peak or valley 2
  Fatigue_DP = 0; //Peak or valley 4
  Fatigue_PCCP = 0; /*Previous Cycle counter flag if >1 then previous 'n'
					cycles did not flag a complete cycle */
  Fatigue_R1FP = 0; //Flag for first  peak count
  Fatigue_R2FP = 0; //Flag for second peak count
  Fatigue_cSlopeP = 0; //Current Slope
  Fatigue_PSP = 0; //Previous slope
  Fatigue_EPP = 0; //Previous Strain
  Fatigue_SFP = 0; /*Start Flag = 0 if very first strain, (i.e. when initializing)
				   = 1 otherwise */
  Fatigue_DLP = 0; //Damage if current strain was last peak.

				   // added 6/9/2006
  Fatigue_SR1P = 0;
  Fatigue_NC1P = 0;
  Fatigue_SR2P = 0;
  Fatigue_NC2P = 0;
  Fatigue_SR3P = 0;
  Fatigue_NC3P = 0;

  Fatigue_DmaxP = Fatigue_Dmax;
  Fatigue_E0P = Fatigue_E0;
  Fatigue_mP = Fatigue_m;

  Zd = pow((Cf / Cd), (1. / Alpha));

  Lambda_SR = Zd * Fatigue_DI;
  Lambda_SRP = Lambda_SR;

  Fatigue_FyInitial = Fy;
  Fatigue_Fy = Fy;
  Fatigue_FyP = Fy;

}

Steel02Fatigue::Steel02Fatigue(void):
  UniaxialMaterial(0, MAT_TAG_Steel02Fatigue)
{
  konP = 0;
}

Steel02Fatigue::~Steel02Fatigue(void)
{
  // Does nothing
}

static int sign(double a) {
	if (a < 0)
		return -1;
	else if (a == 0)
		return 0;
	else
		return 1;
}

int
Steel02Fatigue::commitSteel(void) {
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

	Lambda_SRP = Lambda_SR;
	Fatigue_FyP = Fy;

	Fatigue_CfailedP = Fatigue_Cfailed;

	//opserr << "FatigueMaterial: material tag " << this->getTag() << ", current Lambda_SR : " << Lambda_SR << " \n";

	return 0;
}

UniaxialMaterial*
Steel02Fatigue::getCopy(void) // Add Fatigue... Need to add additional variables later
{
	Steel02Fatigue *theCopy = new Steel02Fatigue(this->getTag(), Fy, E0, b, R0, cR1, cR2, a1, a2, a3, a4, Cd, Cf, Alpha, Beta, minStrain, maxStrain, sigini);
  
  return theCopy;
}

double
Steel02Fatigue::getInitialTangent(void) // Added Fatigue
{
  return E0;
}

int
Steel02Fatigue::setTrialStrain(double trialStrain, double strainRate)
{
  Fy = Fatigue_FyInitial * (1. - Lambda_SR);

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
  
  if (kon == 2 && deps > 0.0) {


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

  } else if (kon == 1 && deps < 0.0) {
    
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

  return 0;
}

double 
Steel02Fatigue::getStrain(void) // Added Fatigue
{
  return eps;
}

double 
Steel02Fatigue::getStress(void)  // Added Fatigue
{
 double modifier = 1.0;
 double damageloc = 1.0 - Fatigue_Dmax + Fatigue_DL;
 if (Fatigue_Cfailed) {
	 // Reduce stress to 0.0 
	 return sig * 1.0e-8;
 }
 else {
	return sig;
 }
}

double 
Steel02Fatigue::getTangent(void)  // Added Fatigue
{
 double modifier = 1.0;
 double damageloc = 1.0 - Fatigue_Dmax + Fatigue_DL;
 if (Fatigue_Cfailed)
 	// Reduce tangent to 0.0 
 	return 1.0e-8*e;
 else
	return e; 
}

int 
Steel02Fatigue::commitState(void) // Added Fatigue
{
	//if (Fatigue_Cfailed) {
	//	opserr << "FatigueMaterial: " << "True " << minStrain << " " << maxStrain << "\n";
	//}
	//else {
	//	opserr << "FatigueMaterial: " << "False " << minStrain << " " << maxStrain << "\n";
	//}

	trialStrain = eps;

	// NOTE: Do not accumulate damage if peaks are too small (e.g. if X < 1e-10)
	// may get floating point errors.  This is essentially a filter for 
	// strain cycles smaller than 1e-10 strain.  
	//  THIS FATIGE MATERIAL CODE WAS NOT INTENDED FOR HIGH CYCLE FATIGUE, 
	//  THE LINEAR ACCUMULATION OF DAMAGE IS NOT AS APPROPRIATE FOR HIGH
	//  CYCLE FATIGUE. 

	//added 6/9/2006
	// for recorder
	Fatigue_SR1 = 0;
	Fatigue_NC1 = 0;


	// No need to continue if the uniaxial material copy 
	// has already failed.
	if (Fatigue_Cfailed) {
		//commitSteel(); //return 0;
		return 0;
	}

	//opserr << "FatigueMaterial: " << trialStrain << " " << maxStrain << " " << " \n";

	//Simple check to see if we reached max strain capacities
	if (trialStrain >= maxStrain || trialStrain <= minStrain) {
		Fatigue_Cfailed = true;
		opserr << "FatigueMaterial: material tag " << this->getTag() << " failed from excessive strain\n";
		Fatigue_DI = Fatigue_Dmax;
		Fatigue_DL = Fatigue_Dmax;

		//commitSteel();  //return 0;
		return 0;
	}

	//Initialize the fatigue parameters if they have 
	// not been initialized yet
	if (Fatigue_SF == 0) {

		Fatigue_A = trialStrain;
		Fatigue_SF = 1;
		Fatigue_EP = trialStrain;
		// Initialize other params if not done so already
		Fatigue_PCC = 0;
		Fatigue_B = 0;
		Fatigue_C = 0;
		Fatigue_D = 0;

	}

	/* Now we need to determine if we are at a peak or not
	If we are, then we need to do some calcs to determine
	the amount of damage suffered. If we are not at a peak, we need to
	pretend like we are at a peak, so that we can calculate the damage
	as if it WERE a peak.
	*/

	// Determine the slope of the strain hysteresis
	if (Fatigue_EP == trialStrain) {
		Fatigue_cSlope = Fatigue_PS;         // No real slope here....
	}
	else {
		Fatigue_cSlope = trialStrain - Fatigue_EP;   // Determine Current Slope
	}


	// If we are at a peak or a valley, then check for damage
	if (sign(Fatigue_PS) != sign(Fatigue_cSlope) && sign(Fatigue_PS) != 0) {

		if (Fatigue_R1F == 0) {    // mark second peak

			Fatigue_B = Fatigue_EP;
			Fatigue_Y = fabs(Fatigue_B - Fatigue_A);
			Fatigue_R1F = 1;

		}
		else {   // start at least the third peak

				 // begin modified Rainflow cycle counting
			if (Fatigue_PCC == 1) {

				Fatigue_D = Fatigue_EP;
				Fatigue_X = fabs(Fatigue_D - Fatigue_C);

			}
			else {

				Fatigue_C = Fatigue_EP;
				Fatigue_X = fabs(Fatigue_C - Fatigue_B);

			}

			if (Fatigue_X < Fatigue_Y) {

				Fatigue_PCC = Fatigue_PCC + 1;

				if (Fatigue_PCC == 1) {
					Fatigue_Y = fabs(Fatigue_C - Fatigue_B);
				}
				else if (Fatigue_PCC == 2) {
					// Count X = |D-C| as a 1.0 cycle
					Fatigue_DI = Fatigue_DI + 1.0 / fabs(pow((Fatigue_X / Fatigue_E0), 1 / Fatigue_m));
					Lambda_SR = Zd * Fatigue_DI;
					//added 6/9/2006
					Fatigue_SR1 = Fatigue_X;
					Fatigue_NC1 = 1.0;
					// Reset parameters
					Fatigue_D = 0;
					Fatigue_C = 0;
					Fatigue_Y = fabs(Fatigue_B - Fatigue_A);
					Fatigue_PCC = 0;

				}

			}
			else {

				if (Fatigue_PCC == 1) {

					// Count Y = |C-B| as a 1.0 cycle
					Fatigue_DI = Fatigue_DI + 1.0 / fabs(pow((Fatigue_Y / Fatigue_E0), 1 / Fatigue_m));
					Lambda_SR = Zd * Fatigue_DI;
					//added 6/9/2006
					Fatigue_SR1 = Fatigue_Y;
					Fatigue_NC1 = 1.0;

					// Reset parameters
					Fatigue_B = Fatigue_D;
					Fatigue_C = 0;
					Fatigue_D = 0;
					Fatigue_Y = fabs(Fatigue_B - Fatigue_A);
					Fatigue_PCC = 0;



				}
				else {

					// Count Y = |A-B| as a 0.5 cycle
					Fatigue_DI = Fatigue_DI + 0.5 / fabs(pow((Fatigue_Y / Fatigue_E0), 1 / Fatigue_m));
					Lambda_SR = Zd * Fatigue_DI;
					//added 6/9/2006
					// For recorder
					Fatigue_SR1 = Fatigue_Y;
					Fatigue_NC1 = 0.5;

					// Reset parameters
					Fatigue_A = Fatigue_B;
					Fatigue_B = Fatigue_C;
					Fatigue_C = 0;
					Fatigue_D = 0;
					Fatigue_Y = Fatigue_X;
					Fatigue_PCC = 0;


				}

			}

		}

		// Flag failure if we have reached that point
		if (Fatigue_DI >= Fatigue_Dmax) {
			// Most likely will not fail at this point, more 
			// likely at the psuedo peak. But this step is
			// is important for accumulating damage
			Fatigue_Cfailed = true;
			opserr << "FatigueMaterial: material tag " << this->getTag() << " failed at peak\n";
			Fatigue_DL = Fatigue_DI;
		}
		else {
			Fatigue_Cfailed = false;
			Fatigue_DL = Fatigue_DI;
		}

		// Modified by Patxi 3/5/2006
		// } else {
	}
	if (Fatigue_Cfailed == false) {

		// Now check for damage, although we may not be at a peak at all.
		// Store temporary damage only as if it were the last peak: DL
		// Commit to DI only if failure occurs.
		if (Fatigue_B == 0 && Fatigue_C == 0 && Fatigue_D == 0) {

			// If we have not yet found the second peak
			Fatigue_X = fabs(trialStrain - Fatigue_A);

			if (fabs(Fatigue_X) < 1e-10) {
				Fatigue_DL = Fatigue_DI;
				// added 6/9/2006
				//values for recorder
				Fatigue_SR2 = 0.0;
				Fatigue_NC2 = 0.0;
				Fatigue_SR3 = 0.0;
				Fatigue_NC3 = 0.0;
			}
			else {
				Fatigue_DL = Fatigue_DI + 0.5 / fabs(pow((Fatigue_X / Fatigue_E0), 1 / Fatigue_m));
				// added 6/9/2006
				//values for recorder
				Fatigue_SR2 = Fatigue_X;
				Fatigue_NC2 = 0.5;
				Fatigue_SR3 = 0.0;
				Fatigue_NC3 = 0.0;
			}

		}
		else if (Fatigue_B != 0 && Fatigue_C == 0 && Fatigue_D == 0) {

			// On our way to find point C. Range Y defined, no X yet
			Fatigue_X = fabs(trialStrain - Fatigue_B);

			if (fabs(Fatigue_X) < 1e-10) {
				Fatigue_DL = Fatigue_DI;
				// added 6/9/2006
				//values for recorder
				Fatigue_SR2 = 0.0;
				Fatigue_NC2 = 0.0;
			}
			else {
				Fatigue_DL = Fatigue_DI + 0.5 / fabs(pow((Fatigue_X / Fatigue_E0), 1 / Fatigue_m));
				// added 6/9/2006
				//values for recorder
				Fatigue_SR2 = Fatigue_X;
				Fatigue_NC2 = 0.5;
			}

			if (fabs(Fatigue_Y) < 1e-10) {
				Fatigue_DL = Fatigue_DL;
				// added 6/9/2006
				//values for recorder
				Fatigue_SR3 = 0.0;
				Fatigue_NC3 = 0.0;
			}
			else {
				Fatigue_DL = Fatigue_DL + 0.5 / fabs(pow((Fatigue_Y / Fatigue_E0), 1 / Fatigue_m));
				// added 6/9/2006
				//values for recorder
				Fatigue_SR3 = Fatigue_Y;
				Fatigue_NC3 = 0.5;
			}

		}
		else if (Fatigue_B != 0 && Fatigue_C != 0 && Fatigue_D == 0) {

			// Two ranges stored, but no cycles for either stored
			// There are two scenarios that can result:
			// 1.)  |A-D| is the predominate 1/2 cycle, 
			//      and |B-C| is a small full cycle.
			//
			// 2.)  One full cylce at |D-C|, 1/2 cycle at |A-B|

			if (fabs(Fatigue_A - trialStrain)> fabs(Fatigue_A - Fatigue_B)) {

				//   count 1/2 cycle at |D-A|, and one full cycle at |B-C|.
				Fatigue_X = fabs(trialStrain - Fatigue_A);

				if (fabs(Fatigue_Y) < 1e-10) {
					Fatigue_DL = Fatigue_DI;
					// added 6/9/2006
					//values for recorder
					Fatigue_SR3 = 0.0;
					Fatigue_NC3 = 0.0;
				}
				else {
					Fatigue_DL = Fatigue_DI + 1.0 / fabs(pow((Fatigue_Y / Fatigue_E0), 1 / Fatigue_m));
					// added 6/9/2006
					//values for recorder
					Fatigue_SR3 = Fatigue_Y;
					Fatigue_NC3 = 1.0;
				}

				if (fabs(Fatigue_X) < 1e-10) {
					Fatigue_DL = Fatigue_DL;
					// added 6/9/2006
					//values for recorder
					Fatigue_SR2 = 0.0;
					Fatigue_NC2 = 0.0;
				}
				else {
					Fatigue_DL = Fatigue_DL + 0.5 / fabs(pow((Fatigue_X / Fatigue_E0), 1 / Fatigue_m));
					// added 6/9/2006
					//values for recorder
					Fatigue_SR2 = Fatigue_X;
					Fatigue_NC2 = 0.5;
				}

			}
			else {

				// One full cycle of |C-D| and 1/2 cyle of |A-B|

				if (fabs(Fatigue_C - trialStrain) < 1e-10) {
					Fatigue_DL = Fatigue_DI;
					// added 6/9/2006
					//values for recorder
					Fatigue_SR3 = 0.0;
					Fatigue_NC3 = 0.0;
				}
				else {
					Fatigue_DL = Fatigue_DI + 1.0 / fabs(pow((fabs(Fatigue_C - trialStrain) / Fatigue_E0), 1 / Fatigue_m));
					// added 6/9/2006
					//values for recorder
					Fatigue_SR3 = fabs(Fatigue_C - trialStrain);
					Fatigue_NC3 = 1.0;
				}

				if (fabs(Fatigue_A - Fatigue_B) < 1e-10) {
					Fatigue_DL = Fatigue_DL;
					// added 6/9/2006
					//values for recorder
					Fatigue_SR2 = 0.0;
					Fatigue_NC2 = 0.0;
				}
				else {
					Fatigue_DL = Fatigue_DL + 0.5 / fabs(pow((fabs(Fatigue_A - Fatigue_B) / Fatigue_E0), 1 / Fatigue_m));
					// added 6/9/2006
					//values for recorder
					Fatigue_SR2 = fabs(Fatigue_A - Fatigue_B);
					Fatigue_NC2 = 0.5;
				}
			}
		}

		// Did we fail before a peak?
		double mStress = getStress();
		if (Fatigue_DL > Fatigue_Dmax && mStress > 0.0) {
			Fatigue_DI = Fatigue_DL;
			Fatigue_Cfailed = true;
			opserr << "FatigueMaterial: material tag " << this->getTag() << " failed at pseudo peak\n";
		}
		else {
			Fatigue_Cfailed = false;
		}

	}

	Fatigue_PS = Fatigue_cSlope;        // Previous Slope
	Fatigue_EP = trialStrain;   // Keep track of previous strain

    // Check if failed at current step
	if (Fatigue_Cfailed) {
		return 0;
	}
	else {
		commitSteel();
		return 0;
	}

  return 0;
}

int 
Steel02Fatigue::revertToLastCommit(void) // Added Fatigue
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

  // Fatigue Material Previous States
  Fatigue_DI = Fatigue_DIP; //Damage index
  Fatigue_X = Fatigue_XP; //Range in consideration
  Fatigue_Y = Fatigue_YP; //Previous Adjacent Range
  Fatigue_A = Fatigue_AP; //Peak or valley 1
  Fatigue_B = Fatigue_BP; //Peak or valley 2
  Fatigue_C = Fatigue_CP; //Peak or valley 2
  Fatigue_D = Fatigue_DP; //Peak or valley 4
  Fatigue_PCC = Fatigue_PCCP; /*Previous Cycle counter flag if >1 then previous 'n'
					cycles did not flag a complete cycle */
  Fatigue_R1F = Fatigue_R1FP; //Flag for first  peak count
  Fatigue_R2F = Fatigue_R2FP; //Flag for second peak count
  Fatigue_cSlope = Fatigue_cSlopeP; //Current Slope
  Fatigue_PS = Fatigue_PSP; //Previous slope
  Fatigue_EP = Fatigue_EPP; //Previous Strain
  Fatigue_SF = Fatigue_SFP; /*Start Flag = 0 if very first strain, (i.e. when initializing)
				   = 1 otherwise */
  Fatigue_DL = Fatigue_DLP; //Damage if current strain was last peak.

				   // added 6/9/2006
  Fatigue_SR1 = Fatigue_SR1P;
  Fatigue_NC1 = Fatigue_NC1P;
  Fatigue_SR2 = Fatigue_SR2P;
  Fatigue_NC2 = Fatigue_NC2P;
  Fatigue_SR3 = Fatigue_SR3P;
  Fatigue_NC3 = Fatigue_NC3P;

  Fatigue_Dmax = Fatigue_DmaxP;
  Fatigue_E0 = Fatigue_E0P;
  Fatigue_m = Fatigue_mP;

  Lambda_SR = Lambda_SRP;

  Fy = Fatigue_FyP;

  Fatigue_Cfailed = Fatigue_CfailedP;

  return 0;
}

int 
Steel02Fatigue::revertToStart(void) // Add Fatigue
{
  eP = E0;
  epsP = 0.0;
  sigP = 0.0;
  sig = 0.0;
  eps = 0.0;
  e = E0;  

  konP = 0;
  epsmaxP = Fatigue_FyInitial/E0;
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

  // Fatigue Variables
  Fatigue_Cfailed = false;
  Fatigue_DI = 0; //Damage index
  Fatigue_X = 0; //Range in consideration
  Fatigue_Y = 0; //Previous Adjacent Range
  Fatigue_A = 0; //Peak or valley 1
  Fatigue_B = 0; //Peak or valley 2
  Fatigue_C = 0; //Peak or valley 2
  Fatigue_D = 0; //Peak or valley 4
  Fatigue_PCC = 0; /*Previous Cycle counter flag if >1 then previous 'n'
		   cycles did not flag a complete cycle */
  Fatigue_R1F = 0; //Flag for first  peak count
  Fatigue_R2F = 0; //Flag for second peak count
  Fatigue_cSlope = 0; //Current Slope
  Fatigue_PS = 0; //Previous slope
  Fatigue_EP = 0; //Previous Strain
  Fatigue_SF = 0; /*Start Flag = 0 if very first strain, (i.e. when initializing)
		  = 1 otherwise */
  Fatigue_DL = 0; //Damage if current strain was last peak.

  // added 6/9/2006
  //values for recorder
  Fatigue_SR1 = 0;
  Fatigue_NC1 = 0;
  Fatigue_SR2 = 0;
  Fatigue_NC2 = 0;
  Fatigue_SR3 = 0;
  Fatigue_NC3 = 0;

  Zd = pow((Cf / Cd), (1. / Alpha));

  Lambda_SR = Zd * Fatigue_DI;
  Lambda_SRP = Lambda_SR;

  Fatigue_Fy = Fatigue_FyInitial;
  Fatigue_FyP = Fatigue_FyInitial;

  return 0;
} 

int 
Steel02Fatigue::sendSelf(int commitTag, Channel &theChannel)
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
    opserr << "Steel02Fatigue::sendSelf() - failed to sendSelf\n";
    return -1;
  }
  return 0;
}

int 
Steel02Fatigue::recvSelf(int commitTag, Channel &theChannel,
	     FEM_ObjectBroker &theBroker)
{
  static Vector data(23);

  if (theChannel.recvVector(this->getDbTag(), commitTag, data) < 0) {
    opserr << "Steel02Fatigue::recvSelf() - failed to recvSelf\n";
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
  eP   = data(20);   
  this->setTag(int(data(21)));
  sigini = data(22);

  e = eP;
  sig = sigP;
  eps = epsP;
  
  return 0;
}

void 
Steel02Fatigue::Print(OPS_Stream &s, int flag)
{
  if (flag == OPS_PRINT_PRINTMODEL_MATERIAL) {      
    //    s << "Steel02Fatigue:(strain, stress, tangent) " << eps << " " << sig << " " << e << endln;
    s << "Steel02Fatigue tag: " << this->getTag() << endln;
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
	s << "\"type\": \"Steel02Fatigue\", ";
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
