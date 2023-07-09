/* ****************************************************************** **
**    OpenSees - Open System for Earthquake Engineering Simulation    **
**          Pacific Earthquake Engineering Research Center            **
**                                                                    **
**                                                                    **
** (C) Copyright 1999, The Regents of the University of California    **
** All Rights Reserved.                                               **
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
// $Date: 2011-04-16 00:04:05 $
// $Source: /usr/local/cvs/OpenSees/SRC/material/uniaxial/Cast.cpp,v $
                                                                      
// Written: Dimitrios Lignos
// Created: 04/2011
//
// Description: This file contains the class implementation of Cast. 
// This Cast is a material model that describes the Cast Fuse brace
// Developed by Gray et al. (2010).
// this model is based on:
// 1. Modified Menegotto-Pinto Steel Model with Filippou Isotropic Hardening
// 2. Closed Form solution for monotonic loading developed by Gray et al. (2010)
//-----------------------------------------------------------------------

#include <stdlib.h>
#include <Cast.h>
#include <Vector.h>
#include <Matrix.h>
#include <Channel.h>
#include <Information.h>
#include <Parameter.h>

#include <math.h>

#include <float.h>

#include <elementAPI.h>
#include <OPS_Globals.h>


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
	double dData[14];
	int numData = 1;
	// Check Tag and number of Fingers
	if (OPS_GetIntInput(&numData, iData) != 0) {
		opserr << "WARNING invalid uniaxialMaterial  Cast Fuse tag" << endln;
		return 0;
	}
	
	numData = OPS_GetNumRemainingInputArgs();
	if (numData < 14) {
	    opserr << "WARNING insufficient number of args want  uniaxialMaterial CastFuse tag? NLegs? bo? h? Fy? E? L? b? R0? cR1? cR2? a1? a2? a3? a4\n";
	    return 0;
	}
	if (OPS_GetDoubleInput(&numData, dData) != 0) {
		opserr << "Invalid Args want: uniaxialMaterial CastFuse tag? NLegs? bo? h? Fy? E? L? b? R0? cR1? cR2? a1? a2? a3? a4?";
		return 0;	
	}
	
	// Parsing was successful, allocate the material
	theMaterial = new Cast(iData[0], dData[0],
							   dData[1], dData[2], dData[3], dData[4], dData[5], 
							   dData[6], dData[7], dData[8], dData[9], dData[10],
						       dData[11], dData[12], dData[13]);
	
	if (theMaterial == 0) {
		opserr << "WARNING could not create uniaxialMaterial of type Cast Material\n";
		return 0;
	}
	
	return theMaterial;
}

Cast::Cast(int tag, double NLegs, double _BO, double H, double _Fy,
		 double E, double L, double _b,
		 double _R0, double _cR1, double _cR2,
		 double _a1, double _a2, double _a3, double _a4):
  UniaxialMaterial(tag, MAT_TAG_Cast),
  nLegs(NLegs), bo(_BO), h(H), fy(_Fy), eo(E), l(L), b(_b), R0(_R0), 
  cR1(_cR1), cR2(_cR2), a1(_a1), a2(_a2), a3(_a3), a4(_a4)
{
  
  kp = (1.0/6.0) * nLegs * bo * eo * pow((h/l),3.0); // Initial Stiffness of Cast Fuse
  Pp = nLegs * bo * pow(h,2.0) * fy/(4.0*l);   // Yield strength of Cast Fuse (Monotonic Backbone)

  konP = 0;
  kon = 0;
  eP = kp;
  epsP = 0.0;
  sigP = 0.0;
  sig = 0.0;
  eps = 0.0;
  e = kp;

  epsmaxP = Pp/kp;
  epsminP = -epsmaxP;
  epsplP = 0.0;
  epss0P = 0.0;
  sigs0P = 0.0;
  epssrP = 0.0;
  sigsrP = 0.0;

  // added by MG
  epsmaxrP = Pp/kp;
  epsminrP = -epsmaxrP;

}

Cast::Cast(int tag, double NLegs, double _BO, double H, double _Fy,
		   double E, double L, double _b,
		   double _R0, double _cR1, double _cR2):
  UniaxialMaterial(tag, MAT_TAG_Cast),
  nLegs(NLegs), bo(_BO), h(H), fy(_Fy), eo(E), l(L), b(_b), R0(_R0),
  cR1(_cR1), cR2(_cR2)
{

  kp = (1.0/6.0) * nLegs * bo * eo * pow((h/l),3.0); // Initial Stiffness of Cast Fuse
  Pp = nLegs * bo * pow(h,2.0) * fy/(4.0*l);         // Yield strength of Cast Fuse (Monotonic Backbone)

  konP = 0;

  // Default values for no isotropic hardening
  a1 = 0.0;
  a2 = 1.0;
  a3 = 0.0;
  a4 = 1.0;

  eP = kp;
  epsP = 0.0;
  sigP = 0.0;
  sig = 0.0;
  eps = 0.0;
  e = kp;

  epsmaxP = Pp/kp;
  epsminP = -epsmaxP;
  epsplP = 0.0;
  epss0P = 0.0;
  sigs0P = 0.0;
  epssrP = 0.0;
  sigsrP = 0.0;

  // added by MG
 double  epsmaxrP = Pp/kp;
 double  epsminrP = -epsmaxrP;
}

Cast::Cast(int tag, double NLegs, double _BO, double H,
		   double _Fy, double E, double L, double _b):
  UniaxialMaterial(tag, MAT_TAG_Cast),
  nLegs(NLegs), bo(_BO), h(H), fy(_Fy), eo(E), l(L), b(_b)
{
  kp = (1.0/6.0) * nLegs * bo * eo * pow((h/l),3.0); // Initial Stiffness of Cast Fuse
  Pp = nLegs * bo * pow(h,2.0) * fy/(4.0*l);         // Yield strength of Cast Fuse (Monotonic Backbone)

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

  eP = kp;
  epsP = 0.0;
  sigP = 0.0;
  sig = 0.0;
  eps = 0.0;
  e = kp;

  epsmaxP = Pp/kp;
  epsminP = -epsmaxP;
  epsplP = 0.0;
  epss0P = 0.0;
  sigs0P = 0.0;
  epssrP = 0.0;
  sigsrP = 0.0;

  // added by MG
double  epsmaxrP = Pp/kp;
double  epsminrP = -epsmaxrP;
}

Cast::Cast(void):
  UniaxialMaterial(0, MAT_TAG_Cast)
{
  konP = 0;
}

Cast::~Cast(void)
{
  // Does nothing
}

double
Cast::getInitialTangent(void)
{
  
  return kp;
}

int
Cast::setTrialStrain(double trialStrain, double strainRate)
{
  
  double Esh = b * kp;
  
  double epsy = Pp / kp;
	
  eps = trialStrain;

  double deps = eps - epsP;
  
  epsmax = epsmaxP;
  epsmin = epsminP;
  epspl  = epsplP;
  epss0  = epss0P;  
  sigs0  = sigs0P; 
  epsr   = epssrP;  
  sigr   = sigsrP;  
  kon = konP;

  //added by MG
double  epsminr = epsminrP;
double  epsmaxr = epsmaxrP;

  if (kon == 0 ) {
	  if (fabs(deps) < 10.0*DBL_EPSILON) {
		  
		  e = kp;
		  sig = 0;
		  return 0;
		
	  } else {
		  epsmax = epsy;
	      epsmin = -epsy;
	
		  if (deps < 0.0) {
			  kon = 2;						  
			  epss0 = epsmin;		  
			  sigs0 = -Pp;
			  epspl = epsmin;
				  
		  } else {		  
			  kon = 1;		 
			  epss0 = epsmax;		  
			  sigs0 = Pp;
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

	// MG: This variable is used to evaluate if the strain is enough since the
	// previous reversal to consider an increase in strain hardening
	double yieldcheck = (eps - epsr)/(epss0 - epsr);

    epsr = epsP;

	// MG: This ensures that in the case where post-yeilding is applied it isn't double counted
    if ((eps > 0 && sig > 0) || (eps < 0 && sig < 0)) {
	  sigr = sigP * cos(2.0*epsP/l);
  } else {
      sigr = sigP;
  }
	// commented before  
	// MG: This is redundant
    // epsmin = min(epsP, epsmin);

	// MG: This stores the min strain that was used to calculate the strain
	// hardening before this load reversal.  In the case where there
	// is a quick double back we won't consider strain hardening
	// epsminr = epsmin;

    if (epsP < epsmin)
      epsmin = epsP;

      //double d1 = (epsmax - epsmin) / (2.0*(a4 * epsy));
	  
	  // MG: recalculate d1 if the material hasn't yielded in this excursion
	  if (fabs(yieldcheck) > 1){
		  epsmaxr = epsmax;
	  }

	  double d1 = (epsmaxr - epsmin) / (2.0*(a2 * epsy));
      double shft = 1.0 + a3 * pow(d1, 0.8);
      epss0 = (Pp * shft - Esh * epsy * shft - sigr + kp * epsr) / (kp - Esh);
      sigs0 = Pp * shft + Esh * (epss0 - epsy * shft);
      epspl = epsmax;

  } else if (kon == 1 && deps < 0.0) {
    
    // update the maximum previous strain, store the last load reversal 
    // point and calculate the stress and strain (sigs0 and epss0) at the 
    // new intersection between elastic and strain hardening asymptote 
    // To include isotropic strain hardening shift the strain hardening 
    // asymptote by sigsft before calculating the intersection point 
    // Constants a1 and a2 control this stress shift on compression side 

    kon = 2;

	// MG: This variable is used to evaluate if strain the strain is enough since the
	// previous reversal to consider strain hardening
	double yieldcheck = (eps - epsr)/(epss0 - epsr);

    epsr = epsP;

	// MG: This ensures that in the case where post-yeilding is applied it isn't double counted
    if ((eps > 0 && sig > 0) || (eps < 0 && sig < 0)) {
	  sigr = sigP * cos(2.0*epsP/l);
  } else {
      sigr = sigP;
  }

	  
	// commented before
	// MG: This is redundant
    // epsmax = max(epsP, epsmax);

	// MG: This stores the max strain that was used to calculate the strain
	// hardening before this load reversal.  In the case where there
	// is a quick double back we won't consider strain hardening
	// epsmaxr = epsmax;
	  
    if (epsP > epsmax)
      epsmax = epsP;
    
    // double d1 = (epsmax - epsmin) / (2.0*(a2 * epsy));

	// MG: recalculate d1 if the material hasn't yielded in this excursion
	if (fabs(yieldcheck) > 1){
	   double epsminr = epsmin;
	} 
	double d1 = (epsmax - epsminr) / (2.0*(a2 * epsy));
    double shft = 1.0 + a1 * pow(d1, 0.8);
    epss0 = (-Pp * shft + Esh * epsy * shft - sigr + kp * epsr) / (kp - Esh);
    sigs0 = -Pp * shft + Esh * (epss0 + epsy * shft);
    epspl = epsmin;
  }

  
  // calculate current stress sig and tangent modulus E 

  double xi     = fabs((epspl-epss0)/epsy);
  double R      = R0*(1.0 - (cR1*xi)/(cR2+xi));
  double epsrat = (eps-epsr)/(epss0-epsr);
  double dum1  = 1.0 + pow(fabs(epsrat),R);
  double dum2  = pow(dum1,(1/R));

  sig   = b*epsrat +(1.0-b)*epsrat/dum2;
	
  // Added by DL to Change Paths and Pinch
  sig = sig*(sigs0-sigr);
 
  sig = sig+sigr;
  // MG: This variable is used in the tangent calc when post-yield stiffening is considered	
  double sign = 1;

  if ((eps-epsr) < 0){
	  sign = -1;
  }else {
	  sign = 1;
  }
  
  // MG: Checks which quadrant we are in and then considers post-yield stiffening when appropriate
  // this is not technically correct but it fixes the problem that occurs when there is a load
  // reversal in the elastic rebound

  if ((eps > 0 && sig > 0) || (eps < 0 && sig < 0)) {
	  sig = sig/cos(2.0*eps/l);
	  e = ((sigr - sigs0)*(b/(epsr - epss0) - (b - 1)/((epsr - epss0)*pow(pow((fabs(eps - epsr)/fabs(epsr - epss0)),R) + 1,(1/R))) + (sign*(eps - epsr)*pow(fabs(eps - epsr)/fabs(epsr - epss0),(R - 1))*(b - 1))/(fabs(epsr - epss0)*(epsr - epss0)*pow(pow(fabs(eps - epsr)/fabs(epsr - epss0),R) + 1,(1/R + 1)))))/cos((2*eps)/l) + (2*sin((2*eps)/l)*(sigr + ((b*(eps - epsr))/(epsr - epss0) - ((eps - epsr)*(b - 1))/((epsr - epss0)*pow(pow(fabs(eps - epsr)/fabs(epsr - epss0),R) + 1,(1/R))))*(sigr - sigs0)))/(l*pow(cos((2*eps)/l),2));
  } else {
      e = b + (1.0-b)/(dum1*dum2);
      e = e*(sigs0-sigr)/(epss0-epsr);
  }

//  if (fabs(eps)>fabs(epsy)) {
//	  sig = sig/cos(2.0*eps/l) + sigr;
//  }else
  //	  sig = sig + sigr;
  
 // e = b + (1.0-b)/(dum1*dum2);
 // e = e*(sigs0-sigr)/(epss0-epsr);

   

  return 0;
}

double 
Cast::getStrain(void)
{
  return eps;
}

double 
Cast::getStress(void)
{
  return sig;
}

double 
Cast::getTangent(void)
{
  return e;
}

int 
Cast::commitState(void)
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

  // Added by MG
 double epsmaxrP = epsmaxr;
 double epsminrP = epsminr;

  return 0;
}

int 
Cast::revertToLastCommit(void)
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

  // Added by MG
double  epsmaxr = epsmaxrP;
double  epsminr = epsminrP;

  return 0;
}

int 
Cast::revertToStart(void)
{
  kp = (1.0/6.0) * nLegs * bo * eo * pow((h/l),3.0); // Initial Stiffness of Cast Fuse
  Pp = nLegs * bo * pow(h,2.0) * fy/(4.0*l);   // Yield strength of Cast Fuse (Monotonic Backbone)
  
  eP = kp;
  epsP = 0.0;
  sigP = 0.0;
  sig = 0.0;
  eps = 0.0;
  e = kp;  

  konP = 0;
  epsmaxP = Pp/kp;
  epsminP = -epsmaxP;
  epsplP = 0.0;
  epss0P = 0.0;
  sigs0P = 0.0;
  epssrP = 0.0;
  sigsrP = 0.0;

  // Added by MG
double  epsmaxrP = Pp/kp;
double epsminrP = -epsmaxrP;

  return 0;
}

UniaxialMaterial*
Cast::getCopy(void)
{
	Cast *theCopy = new Cast(this->getTag(), nLegs, bo, h, fy, eo, l, b, R0, cR1, cR2, a1, a2, a3, a4);
	
	return theCopy;
}

int 
Cast::sendSelf(int commitTag, Channel &theChannel)
{
  int res = 0;
  static Vector data(30);

  data(0) = this->getTag();
	
  // Material properties
  data(1) = nLegs;
  data(2) = bo;
  data(3) = h;
  data(4) = fy;
  data(5) = eo;
  data(6) = l;
  data(7) = b;

  // State Variables from last converged state
  data(8) = R0;
  data(9) = cR1;
  data(10) = cR2;
  data(11) = a1;
  data(12) = a2;
  data(13) = a3;
  data(14) = a4;
  data(15) = epsminP;
  data(16) = epsmaxP;
  data(17) = epsplP;
  data(18) = epss0P;
  data(19) = sigs0P;
  data(20) = epssrP;
  data(21) = sigsrP;
  data(22) = konP;  
  data(23) = epsP;  
  data(24) = sigP;  
  data(25) = eP;    

  data(26) = Pp;
  data(27) = kp;

  // Added by MG
  data(28) = epsmaxrP;
  data(29) = epsminrP;

  res = theChannel.sendVector(this->getDbTag(), commitTag, data);
  if (res < 0)
	opserr << "Cast::sendSelf() - failed to send data\n";
    
  return res;
}

int 
Cast::recvSelf(int commitTag, Channel &theChannel, 
	     FEM_ObjectBroker &theBroker)
{
  int res = 0;
  static Vector data(30);
  res = theChannel.recvVector(this->getDbTag(),commitTag, data);
  
  if (res < 0) {
	  
  //if (theChannel.recvVector(this->getDbTag(), commitTag, data) < 0) {
    opserr << "Cast::recvSelf() - failed to receive data\n";
	this->setTag(0);
  }
  else {
	  this->setTag((int)data(0));

	  // Material properties
	  nLegs = data(1);
	  bo = data(2);
	  h = data(3);
	  fy = data(4);
	  eo = data(5);
	  l = data(6);
	  b = data(7); 

	  // State Variables from last converged state  
	  R0 = data(8);
	  cR1 = data(9);
	  cR2 = data(10);
	  a1 = data(11); 
	  a2 = data(12);	
	  a3 = data(13); 
	  a4 = data(14); 
	  epsminP = data(15);
	  epsmaxP = data(16);
	  epsplP = data(17); 
	  epss0P = data(18); 
	  sigs0P = data(19); 
	  epssrP = data(20);  
	  sigsrP = data(21); 
	  konP = data(22);   
	  epsP = data(23);   
  	  sigP = data(24);   
	  eP   = data(25);   
	  Pp = data(26);
	  kp = data(27);

	  // Added by MG
	  epsmaxrP = data(28);
	  epsminrP = data(29);
  }
	
  return res;
}

void 
Cast::Print(OPS_Stream &s, int flag)
{
	s << "Cast Fuse tag: " << this->getTag() << endln;
    s << "  Finger Yield Strength: " << fy << endln;	
    s << "  Finger height: " << h << endln;
    s << "  Number of fingers: " << nLegs << endln;
    s << "  Finger Length: " << l << endln;
}
