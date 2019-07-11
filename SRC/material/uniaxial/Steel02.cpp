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
//-----------------------------------------------------------------------

#include <math.h>

#include <stdlib.h>
#include <Steel02.h>
#include <OPS_Globals.h>
#include <float.h>
#include <Channel.h>


#include <elementAPI.h>
#include <OPS_Globals.h>


void *
OPS_Steel02()
{
  // Pointer to a uniaxial material that will be returned
  UniaxialMaterial *theMaterial = 0;

  int    iData[1];
  double dData[12];
  int numData = 1;

  if (OPS_GetIntInput(&numData, iData) != 0) {
    opserr << "WARNING invalid uniaxialMaterial Steel02 tag" << endln;
    return 0;
  }

  numData = OPS_GetNumRemainingInputArgs();

  if (numData != 3 && numData != 6 && numData != 10 && numData != 11) {
    opserr << "Invalid #args, want: uniaxialMaterial Steel02 " << iData[0] << 
      " fy? E? b? <R0? cR1? cR2? <a1? a2? a3? a4?>>" << endln;
    return 0;
  }

  if (numData == 3) {
    if (OPS_GetDoubleInput(&numData, dData) != 0) {
      opserr << "Invalid double: uniaxialMaterial Steel02 " << iData[0] << 
	" fy? E? b? <R0? cR1? cR2? <a1? a2? a3? a4?>>" << endln;
      return 0;
    }

    // Parsing was successful, allocate the material
    theMaterial = new Steel02(iData[0], dData[0], dData[1], dData[2]);    

  } else if (numData == 6) {
    if (OPS_GetDoubleInput(&numData, dData) != 0) {
      opserr << "Invalid int: uniaxialMaterial Steel02 " << iData[0] << 
	" fy? E? b? <R0? cR1? cR2? <a1? a2? a3? a4?>>" << endln;
      return 0;
    }

    // Parsing was successful, allocate the material
    theMaterial = new Steel02(iData[0], dData[0], dData[1], dData[2], dData[3], dData[4], dData[5]);    

  } else if (numData == 10) {
    if (OPS_GetDoubleInput(&numData, dData) != 0) {
      opserr << "Invalid arggs: uniaxialMaterial Steel02 " << iData[0] << 
	" fy? E? b? <R0? cR1? cR2? <a1? a2? a3? a4?>>" << endln;
      return 0;
    }

    // Parsing was successful, allocate the material
    theMaterial = new Steel02(iData[0], dData[0], dData[1], dData[2], 
			      dData[3], dData[4], dData[5], dData[6], 
			      dData[7], dData[8], dData[9]);    

  } else if (numData == 11) {
    if (OPS_GetDoubleInput(&numData, dData) != 0) {
      opserr << "Invalid arggs: uniaxialMaterial Steel02 " << iData[0] << 
	" fy? E? b? <R0? cR1? cR2? <a1? a2? a3? a4?>>" << endln;
      return 0;
    }

    // Parsing was successful, allocate the material
    theMaterial = new Steel02(iData[0], dData[0], dData[1], dData[2], 
			      dData[3], dData[4], dData[5], dData[6], 
			      dData[7], dData[8], dData[9], dData[10]);    

  }   

  if (theMaterial == 0) {
    opserr << "WARNING could not create uniaxialMaterial of type Steel02 Material\n";
    return 0;
  }

  return theMaterial;
}


Steel02::Steel02(int tag,
		 double _Fy, double _E0, double _b,
		 double _R0, double _cR1, double _cR2,
		 double _a1, double _a2, double _a3, double _a4, double sigInit):
  UniaxialMaterial(tag, MAT_TAG_Steel02),
  Fy(_Fy), E0(_E0), b(_b), R0(_R0), cR1(_cR1), cR2(_cR2), a1(_a1), a2(_a2), a3(_a3), a4(_a4), 
  sigini(sigInit)
{
	EnergyP = 0;	//by SAJalali
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

Steel02::Steel02(int tag,
		 double _Fy, double _E0, double _b,
		 double _R0, double _cR1, double _cR2):
  UniaxialMaterial(tag, MAT_TAG_Steel02),
  Fy(_Fy), E0(_E0), b(_b), R0(_R0), cR1(_cR1), cR2(_cR2), sigini(0.0)
{
	EnergyP = 0;	//by SAJalali
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
}

Steel02::Steel02(int tag, double _Fy, double _E0, double _b):
  UniaxialMaterial(tag, MAT_TAG_Steel02),
  Fy(_Fy), E0(_E0), b(_b), sigini(0.0)
{
	EnergyP = 0;	//by SAJalali
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
}

Steel02::Steel02(void):
  UniaxialMaterial(0, MAT_TAG_Steel02)
{
	EnergyP = 0;	//by SAJalali
	konP = 0;
}

Steel02::~Steel02(void)
{
  // Does nothing
}

UniaxialMaterial*
Steel02::getCopy(void)
{
  Steel02 *theCopy = new Steel02(this->getTag(), Fy, E0, b, R0, cR1, cR2, a1, a2, a3, a4, sigini);
  
  return theCopy;
}

double
Steel02::getInitialTangent(void)
{
  return E0;
}

int
Steel02::setTrialStrain(double trialStrain, double strainRate)
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
Steel02::getStrain(void)
{
  return eps;
}

double 
Steel02::getStress(void)
{
  return sig;
}

double 
Steel02::getTangent(void)
{
  return e;
}

int 
Steel02::commitState(void)
{
  epsminP = epsmin;
  epsmaxP = epsmax;
  epsplP = epspl;
  epss0P = epss0;
  sigs0P = sigs0;
  epssrP = epsr;
  sigsrP = sigr;
  konP = kon;

  //by SAJalali
  EnergyP += 0.5*(sig + sigP)*(eps - epsP);

  eP = e;
  sigP = sig;
  epsP = eps;

  return 0;
}

int 
Steel02::revertToLastCommit(void)
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
Steel02::revertToStart(void)
{
	EnergyP = 0;	//by SAJalali
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
Steel02::sendSelf(int commitTag, Channel &theChannel)
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
    opserr << "Steel02::sendSelf() - failed to sendSelf\n";
    return -1;
  }
  return 0;
}

int 
Steel02::recvSelf(int commitTag, Channel &theChannel, 
	     FEM_ObjectBroker &theBroker)
{
  static Vector data(23);

  if (theChannel.recvVector(this->getDbTag(), commitTag, data) < 0) {
    opserr << "Steel02::recvSelf() - failed to recvSelf\n";
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
Steel02::Print(OPS_Stream &s, int flag)
{
  if (flag == OPS_PRINT_PRINTMODEL_MATERIAL) {      
    //    s << "Steel02:(strain, stress, tangent) " << eps << " " << sig << " " << e << endln;
    s << "Steel02 tag: " << this->getTag() << endln;
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
	s << "\"type\": \"Steel02\", ";
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
