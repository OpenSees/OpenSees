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
#include <float.h>
#include <Channel.h>
#include <Information.h>
#include <Parameter.h>
#include <Matrix.h> 
#include <string.h>

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
  sigini(sigInit), parameterID(0), SHVs(0)
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
  Fy(_Fy), E0(_E0), b(_b), R0(_R0), cR1(_cR1), cR2(_cR2), sigini(0.0), parameterID(0), SHVs(0)
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
  Fy(_Fy), E0(_E0), b(_b), sigini(0.0), parameterID(0), SHVs(0)
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
  UniaxialMaterial(0, MAT_TAG_Steel02), parameterID(0), SHVs(0)
{
	EnergyP = 0;	//by SAJalali
	konP = 0;
}

Steel02::~Steel02(void)
{
    if (SHVs != 0)
        delete SHVs;
}

UniaxialMaterial*
Steel02::getCopy(void)
{
  Steel02 *theCopy = new Steel02(this->getTag(), Fy, E0, b, R0, cR1, cR2, a1, a2, a3, a4, sigini);
  
  theCopy->parameterID = parameterID;
  if (SHVs != 0)
      theCopy->SHVs = new Matrix(*SHVs);
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

  if (SHVs != 0)
      SHVs->Zero();
  return 0;
}

int 
Steel02::sendSelf(int commitTag, Channel &theChannel)
{
  static Vector data(25);
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
  data(23) = parameterID;
  data(24) = -1;
  if (SHVs != 0)
    data(24) = SHVs->noCols();

  if (theChannel.sendVector(this->getDbTag(), commitTag, data) < 0) {
    opserr << "Steel02::sendSelf() - failed to sendSelf\n";
    return -1;
  }

  if (SHVs != 0) {
    if (theChannel.sendMatrix(this->getDbTag(), commitTag, *SHVs) < 0) {
      opserr << "Steel02::sendSelf() - failed to send SHVs matrix\n";
      return -2;
	}
  }
  return 0;
}

int 
Steel02::recvSelf(int commitTag, Channel &theChannel, 
	     FEM_ObjectBroker &theBroker)
{
  static Vector data(25);

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
  parameterID = (int)data(23);

  e = eP;
  sig = sigP;
  eps = epsP;
  
  // Receive sensitivity history variables (SHVs)
  int noCols = (int)data(24);
  if (noCols > 0) {
    if (SHVs != 0)
        delete SHVs;
    SHVs = new Matrix(7, noCols); 
    if (SHVs == 0) {
        opserr << "Steel02::recvSelf() - failed to allocate SHVs matrix\n";
        return -2;
    }

    if (theChannel.recvMatrix(this->getDbTag(), commitTag, *SHVs) < 0) {
        opserr << "Steel02::recvSelf() - failed to receive SHVs matrix\n";
        return -3;
	}
  }

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

// AddingSensitivity:BEGIN ///////////////////////////////////
// Direct differentiation method (DDM) implementation
// Following DDM formulation of Barbato & Conte (2006), extended to all parameters
// Contributed by Luigi Caglio, 2026

int Steel02::activateParameter(int passedParameterID)
{
    parameterID = passedParameterID;
    return 0;
}

double Steel02::getInitialTangentSensitivity(int gradIndex)
{
    if (parameterID == 2) {
        return 1.0;  // dE0/dE0 = 1
    }
    return 0.0;
}

void Steel02::computeAsymptoticPointSensitivities(
    double& epss0Sensitivity,
    double& sigs0Sensitivity,
    double epsrSensitivity,
    double sigrSensitivity,
    int gradIndex
)
{
    // Assign parameter derivatives
    double FySensitivity = 0.0;
    double E0Sensitivity = 0.0;
    double bSensitivity = 0.0;
    double a1Sensitivity = 0.0;
    double a2Sensitivity = 0.0;
    double a3Sensitivity = 0.0;
    double a4Sensitivity = 0.0;

    if (parameterID == 1) FySensitivity = 1.0;
    else if (parameterID == 2) E0Sensitivity = 1.0;
    else if (parameterID == 3) bSensitivity = 1.0;
    else if (parameterID == 4) a1Sensitivity = 1.0;
    else if (parameterID == 5) a2Sensitivity = 1.0;
    else if (parameterID == 6) a3Sensitivity = 1.0;
    else if (parameterID == 7) a4Sensitivity = 1.0;

    // Compute derivatives of material constants
    double Esh = b * E0;
    double dEshDtheta = bSensitivity * E0 + b * E0Sensitivity;
    double fyOneMinusB = Fy * (1.0 - b);
    double dfyOneMinusBDtheta = FySensitivity * (1.0 - b) - Fy * bSensitivity;

    // Compute yield strain derivative
    double epsy = Fy / E0;
    double depsyDtheta = FySensitivity / E0 - Fy / (E0 * E0) * E0Sensitivity;

    // Determine loading direction and compute shift factor
    double sign = (kon == 1) ? 1.0 : -1.0;
    double d1, shft, dd1Dtheta, dshftDtheta;

    if (kon == 1) {  // Tension reversal
        double denomD1 = 2.0 * a4 * epsy;
        if (fabs(denomD1) < 1e-15) {
            d1 = 0.0;
            dd1Dtheta = 0.0;
        }
        else {
            d1 = (epsmax - epsmin) / denomD1;
            double d_denominator = 2.0 * epsy * a4Sensitivity + 2.0 * a4 * depsyDtheta;
            dd1Dtheta = -(epsmax - epsmin) / (denomD1 * denomD1) * d_denominator;
        }

        shft = 1.0 + a3 * pow(d1, 0.8);
        dshftDtheta = a3Sensitivity * pow(d1, 0.8) + a3 * 0.8 * pow(d1, -0.2) * dd1Dtheta;
    }
    else {  // Compression reversal (kon == 2)
        double denomD1 = 2.0 * a2 * epsy;
        if (fabs(denomD1) < 1e-15) {
            d1 = 0.0;
            dd1Dtheta = 0.0;
        }
        else {
            d1 = (epsmax - epsmin) / denomD1;
            double d_denominator = 2.0 * epsy * a2Sensitivity + 2.0 * a2 * depsyDtheta;
            dd1Dtheta = -(epsmax - epsmin) / (denomD1 * denomD1) * d_denominator;
        }

        shft = 1.0 + a1 * pow(d1, 0.8);
        dshftDtheta = a1Sensitivity * pow(d1, 0.8) + a1 * 0.8 * pow(d1, -0.2) * dd1Dtheta;
    }

    // Compute asymptotic point derivatives with shift factor
    double intercept_0 = fyOneMinusB * shft * sign;
    double dintercept0Dtheta = sign * (dfyOneMinusBDtheta * shft + fyOneMinusB * dshftDtheta);

    double sigma_intercept = sigr - E0 * epsr;
    double dsigmaIntercept = sigrSensitivity - (E0Sensitivity * epsr + E0 * epsrSensitivity);

    double numerator = intercept_0 - sigma_intercept;
    double denominator = E0 - Esh;
    double ddenomDtheta = E0Sensitivity - dEshDtheta;

    if (fabs(denominator) < 1e-15) {
        // Special case: E0 ≈ Esh (elastic-perfectly plastic limit)
        epss0Sensitivity = epsrSensitivity + sign * (depsyDtheta * shft + epsy * dshftDtheta);
        sigs0Sensitivity = sigrSensitivity + sign * (FySensitivity * shft + Fy * dshftDtheta);
    }
    else {
        // General case
        epss0Sensitivity = ((dintercept0Dtheta - dsigmaIntercept) * denominator -
            numerator * ddenomDtheta) / (denominator * denominator);
        sigs0Sensitivity = dintercept0Dtheta + epss0Sensitivity * Esh + epss0 * dEshDtheta;
    }

    // Update sensitivity history variables
    if (SHVs != 0) {
        (*SHVs)(2, gradIndex) = epsrSensitivity;
        (*SHVs)(3, gradIndex) = sigrSensitivity;
        (*SHVs)(4, gradIndex) = epss0Sensitivity;
        (*SHVs)(5, gradIndex) = sigs0Sensitivity;
    }
}

double Steel02::computeStressGradient(
    double strainSensitivity,
    double epsrSens,
    double sigrSens,
    double epss0Sens,
    double sigs0Sens,
    double epsplSens
)
{
    double gradient = 0.0;
    // Assign parameter derivatives
    double FySensitivity = 0.0;
    double E0Sensitivity = 0.0;
    double bSensitivity = 0.0;
    double R0Sensitivity = 0.0;
    double cR1Sensitivity = 0.0;
    double cR2Sensitivity = 0.0;

    if (parameterID == 1) FySensitivity = 1.0;
    else if (parameterID == 2) E0Sensitivity = 1.0;
    else if (parameterID == 3) bSensitivity = 1.0;
    else if (parameterID == 8) R0Sensitivity = 1.0;
    else if (parameterID == 9) cR1Sensitivity = 1.0;
    else if (parameterID == 10) cR2Sensitivity = 1.0;

    // Compute basic derivatives
    double epsy = Fy / E0;
    double depsyDtheta = FySensitivity / E0 - Fy / (E0 * E0) * E0Sensitivity;

    // Current values needed for derivatives
    double epsrat = (eps - epsr) / (epss0 - epsr);
    double xi = fabs((epspl - epss0) / epsy);
    double R = R0 * (1.0 - (cR1 * xi) / (cR2 + xi));

    // Derivative of xi
    double xiSign = (epspl - epss0) >= 0 ? 1.0 : -1.0;
    double dxiDtheta = xiSign * ((epsplSens - epss0Sens) * epsy -
        (epspl - epss0) * depsyDtheta) / (epsy * epsy);

    // Derivative of R
    double denomR = cR2 + xi;
    double dRDtheta = 0.0;
    if (fabs(denomR) > 1e-15) {
        double dR_dR0 = 1.0 - (cR1 * xi) / denomR;
        double dR_dcR1 = -R0 * xi / denomR;
        double dR_dcR2 = R0 * cR1 * xi / (denomR * denomR);
        double dR_dxi = -R0 * cR1 * cR2 / (denomR * denomR);
        dRDtheta = dR_dR0 * R0Sensitivity +
            dR_dcR1 * cR1Sensitivity +
            dR_dcR2 * cR2Sensitivity +
            dR_dxi * dxiDtheta;
    }

    // Derivative of normalized strain
    double denomEpsrat = epss0 - epsr;
    double numerEpsrat = eps - epsr;
    double depsratDtheta = ((strainSensitivity - epsrSens) * denomEpsrat -
        numerEpsrat * (epss0Sens - epsrSens)) / (denomEpsrat * denomEpsrat);

    // Derivative of normalized stress f_star
    double absEpsrat = fabs(epsrat);
    double term = 1.0 + pow(absEpsrat, R);
    double g = pow(term, -1.0 / R);

    // df_star/depsrat
    double dgDepsrat = 0.0;
    if (absEpsrat > 1e-15 && R > 1e-15) {
        dgDepsrat = -pow(absEpsrat, R - 1.0) * (epsrat / absEpsrat) *
            pow(term, -1.0 / R - 1.0);
    }
    double df_starDepsrat = b + (1.0 - b) * (g + epsrat * dgDepsrat);

    // df_star/dR
    double df_starDR = 0.0;
    if (absEpsrat > 1e-15 && R > 1e-15) {
        double lnTerm = log(term);
        double lnEps = log(absEpsrat);
        double dgDR = g * (lnTerm / (R * R) - pow(absEpsrat, R) * lnEps / (R * term));
        df_starDR = (1.0 - b) * epsrat * dgDR;
    }

    // df_star/db
    double df_starDb = (absEpsrat > 1e-15) ? epsrat * (1.0 - g) : 0.0;

    // Total derivative of f_star
    double df_starDtheta = df_starDepsrat * depsratDtheta +
        df_starDR * dRDtheta +
        df_starDb * bSensitivity;

    // Derivative of sigma = f_star * (sigs0 - sigr) + sigr
    double f_star_value = b * epsrat + (1.0 - b) * epsrat * g;
    gradient = df_starDtheta * (sigs0 - sigr) +
        f_star_value * (sigs0Sens - sigrSens) +
        sigrSens;

    return gradient;
}

double Steel02::getStressSensitivity(int gradIndex, bool conditional)
{ 
    double gradient = 0.0;


    // Pick up sensitivity history variables
    double CstrainSensitivity = 0.0;
    double CstressSensitivity = 0.0;
    double CepsrSensitivity = 0.0;
    double CsigrSensitivity = 0.0;
    double Cepss0Sensitivity = 0.0;
    double Csigs0Sensitivity = 0.0;
    double CepsplSensitivity = 0.0;

    if (SHVs != 0) {
        CstrainSensitivity = (*SHVs)(0, gradIndex);
        CstressSensitivity = (*SHVs)(1, gradIndex);
        CepsrSensitivity = (*SHVs)(2, gradIndex);
        CsigrSensitivity = (*SHVs)(3, gradIndex);
        Cepss0Sensitivity = (*SHVs)(4, gradIndex);
        Csigs0Sensitivity = (*SHVs)(5, gradIndex);
        CepsplSensitivity = (*SHVs)(6, gradIndex);
    }

    // Handle initial loading curve (before any reversal): Initialize asymptotic point sensitivities
    // Initial loading curve: epss0 = ±Fy/E0, sigs0 = ±Fy, epsr = sigr = 0
    bool onInitialCurve = (fabs(epsr) < 1e-15 && fabs(sigr) < 1e-15);

    if (onInitialCurve) {
        // Recompute asymptotic point sensitivities each time

        // a1-a4 have no effect on initial loading curve
        if (parameterID >= 4 && parameterID <= 7) {
            return 0.0;
        }

        double FySensitivity = 0.0;
        double E0Sensitivity = 0.0;

        if (parameterID == 1) FySensitivity = 1.0;
        else if (parameterID == 2) E0Sensitivity = 1.0;

        double epsy = Fy / E0;
        double sign = (kon == 1) ? 1.0 : -1.0;

        Cepss0Sensitivity = sign * (FySensitivity / E0 - Fy / (E0 * E0) * E0Sensitivity);
        Csigs0Sensitivity = sign * FySensitivity;
        CepsplSensitivity = Cepss0Sensitivity;
        if (SHVs != 0) {
            (*SHVs)(4, gradIndex) = Cepss0Sensitivity;
            (*SHVs)(5, gradIndex) = Csigs0Sensitivity;
            (*SHVs)(6, gradIndex) = CepsplSensitivity;
        }
    }

    //Detect reversal and update asymptotic point sensitivities
    bool reversal = ((konP == 1 && kon == 2) || (konP == 2 && kon == 1));

    if (reversal) {
        // At reversal: reversal point = previous converged state
        CepsrSensitivity = CstrainSensitivity;
        CsigrSensitivity = CstressSensitivity;
        CepsplSensitivity = CstrainSensitivity;

        // Compute asymptotic point sensitivities
        computeAsymptoticPointSensitivities(
            Cepss0Sensitivity,
            Csigs0Sensitivity,
            CepsrSensitivity,
            CsigrSensitivity,
            gradIndex
        );
        if (SHVs != 0) {
            (*SHVs)(6, gradIndex) = CepsplSensitivity;   
        }
    }

    // Compute conditional stress gradient

    gradient = computeStressGradient(
        0.0,  // Conditional: strain fixed
        CepsrSensitivity,
        CsigrSensitivity,
        Cepss0Sensitivity,
        Csigs0Sensitivity,
        CepsplSensitivity
    );

    return gradient;
}

int Steel02::commitSensitivity(double TstrainSensitivity, int gradIndex, int numGrads)
{
    if (SHVs == 0) {
        SHVs = new Matrix(7, numGrads);
        SHVs->Zero();
    }

    // Pick up sensitivity history variables
    double CstrainSensitivity = 0.0;
    double CstressSensitivity = 0.0;
    double CepsrSensitivity = 0.0;
    double CsigrSensitivity = 0.0;
    double Cepss0Sensitivity = 0.0;
    double Csigs0Sensitivity = 0.0;
    double CepsplSensitivity = 0.0;

    if (SHVs != 0) {
        CstrainSensitivity = (*SHVs)(0, gradIndex);
        CstressSensitivity = (*SHVs)(1, gradIndex);
        CepsrSensitivity = (*SHVs)(2, gradIndex);
        CsigrSensitivity = (*SHVs)(3, gradIndex);
        Cepss0Sensitivity = (*SHVs)(4, gradIndex);
        Csigs0Sensitivity = (*SHVs)(5, gradIndex);
        CepsplSensitivity = (*SHVs)(6, gradIndex);
    }

    // Handle initial loading curve (before any reversal): Initialize asymptotic point sensitivities
    // Initial loading curve: epss0 = ±Fy/E0, sigs0 = ±Fy, epsr = sigr = 0
    bool onInitialCurve = (fabs(epsr) < 1e-15 && fabs(sigr) < 1e-15);

    if (onInitialCurve) {
        // Recompute asymptotic point sensitivities each time

        if (parameterID >= 4 && parameterID <= 7) {
            // a1-a4 only matter after at least one load reversal with plastic excursion
            return 0.0; 
        }

        double FySensitivity = 0.0;
        double E0Sensitivity = 0.0;

        if (parameterID == 1) FySensitivity = 1.0;
        else if (parameterID == 2) E0Sensitivity = 1.0;

        double epsy = Fy / E0;
        double sign = (kon == 1) ? 1.0 : -1.0;

        Cepss0Sensitivity = sign * (FySensitivity / E0 - Fy / (E0 * E0) * E0Sensitivity);
        Csigs0Sensitivity = sign * FySensitivity;
        CepsplSensitivity = Cepss0Sensitivity;

        if (SHVs != 0) {
            (*SHVs)(4, gradIndex) = Cepss0Sensitivity;
            (*SHVs)(5, gradIndex) = Csigs0Sensitivity;
            (*SHVs)(6, gradIndex) = CepsplSensitivity;
        }
    }

    // Detect reversal and update
    bool reversal = ((konP == 1 && kon == 2) || (konP == 2 && kon == 1)); 

    double TepsrSensitivity = CepsrSensitivity;
    double TsigrSensitivity = CsigrSensitivity;
    double Tepss0Sensitivity = Cepss0Sensitivity;
    double Tsigs0Sensitivity = Csigs0Sensitivity;
    double TepsplSensitivity = CepsplSensitivity;

    if (reversal) {
        TepsrSensitivity = CstrainSensitivity;
        TsigrSensitivity = CstressSensitivity;
        TepsplSensitivity = CstrainSensitivity;

        int sign = (kon == 1) ? 1 : -1;
        computeAsymptoticPointSensitivities(
            Tepss0Sensitivity,
            Tsigs0Sensitivity,
            TepsrSensitivity,
            TsigrSensitivity,
            gradIndex
        );
    }

    // Compute conditional stress gradient
    double gradient = computeStressGradient(
        TstrainSensitivity,  // Total derivative
        TepsrSensitivity,
        TsigrSensitivity,
        Tepss0Sensitivity,
        Tsigs0Sensitivity,
        TepsplSensitivity
    );

    // add tangent term (only in commitSensitivity, not in getStressSensitivity)
    gradient += e * TstrainSensitivity;

    (*SHVs)(0, gradIndex) = TstrainSensitivity;
    (*SHVs)(1, gradIndex) = gradient;
    (*SHVs)(2, gradIndex) = TepsrSensitivity;
    (*SHVs)(3, gradIndex) = TsigrSensitivity;
    (*SHVs)(4, gradIndex) = Tepss0Sensitivity;
    (*SHVs)(5, gradIndex) = Tsigs0Sensitivity;
    (*SHVs)(6, gradIndex) = TepsplSensitivity;

    return 0;
}

int
Steel02::setParameter(const char **argv, int argc, Parameter &param)
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
  if (strcmp(argv[0],"R0") == 0) {
    param.setValue(R0);
    return param.addObject(8, this);
  }
  if (strcmp(argv[0],"cR1") == 0) {
    param.setValue(cR1);
    return param.addObject(9, this);
  }
  if (strcmp(argv[0],"cR2") == 0) {
    param.setValue(cR2);
    return param.addObject(10, this);
  }
  if (strcmp(argv[0],"sig0") == 0) {
    param.setValue(sigini);
    return param.addObject(11, this);
  }
	
  return -1;
}



int
Steel02::updateParameter(int parameterID, Information &info)
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
  case 8:
    this->R0 = info.theDouble;
    break;
  case 9:
    this->cR1 = info.theDouble;
    break;
  case 10:
    this->cR2 = info.theDouble;
    break;
  case 11:
    this->sigini = info.theDouble;
    break;	  
  default:
    return -1;
  }
  
  return 0;
}


