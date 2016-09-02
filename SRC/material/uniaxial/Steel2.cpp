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
// $Source: /usr/local/cvs/OpenSees/SRC/material/uniaxial/Steel2.cpp,v $
                                                                      
// Written: fmk
// Created: 03/06
//
// Description: This file contains the class implementation of Steel2. 
// This Steel2 is based on an f2c of the FEDEAS material
// Steel2.f which is:
//-----------------------------------------------------------------------
// MENEGOTTO-PINTO STEEL MODEL WITH FILIPPOU ISOTROPIC HARDENING
//            written by MOHD YASSIN (1993)
//          adapted to FEDEAS material library
//    by E. Spacone, G. Monti and F.C. Filippou (1994)
//-----------------------------------------------------------------------

#include <math.h>

#include <stdlib.h>
#include <Steel2.h>
#include <OPS_Globals.h>
#include <float.h>
#include <Channel.h>


#include <elementAPI.h>
#include <OPS_Globals.h>


void *
OPS_Steel2(void)
{
  // Pointer to a uniaxial material that will be returned
  UniaxialMaterial *theMaterial = 0;

  int    iData[1];
  double dData[12];
  int numData = 1;

  if (OPS_GetIntInput(&numData, iData) != 0) {
    opserr << "WARNING invalid uniaxialMaterial Steel2 tag" << endln;
    return 0;
  }

  numData = OPS_GetNumRemainingInputArgs();

  if (numData != 3 && numData != 6 && numData != 10 && numData != 11) {
    opserr << "Invalid #args, want: uniaxialMaterial Steel2 " << iData[0] << 
      " fy? E? b? <R0? cR1? cR2? <a1? a2? a3? a4?>>" << endln;
    return 0;
  }

  if (numData == 3) {
    if (OPS_GetDoubleInput(&numData, dData) != 0) {
      opserr << "Invalid arggs: uniaxialMaterial Steel2 " << iData[0] << 
	" fy? E? b? <R0? cR1? cR2? <a1? a2? a3? a4?>>" << endln;
      return 0;
    }

    // Parsing was successful, allocate the material
    theMaterial = new Steel2(iData[0], dData[0], dData[1], dData[2]);    

  } else if (numData == 6) {
    if (OPS_GetDoubleInput(&numData, dData) != 0) {
      opserr << "Invalid arggs: uniaxialMaterial Steel2 " << iData[0] << 
	" fy? E? b? <R0? cR1? cR2? <a1? a2? a3? a4?>>" << endln;
      return 0;
    }

    // Parsing was successful, allocate the material
    theMaterial = new Steel2(iData[0], dData[0], dData[1], dData[2], dData[3], dData[4], dData[5]);    

  } else if (numData == 10) {
    if (OPS_GetDoubleInput(&numData, dData) != 0) {
      opserr << "Invalid arggs: uniaxialMaterial Steel2 " << iData[0] << 
	" fy? E? b? <R0? cR1? cR2? <a1? a2? a3? a4?>>" << endln;
      return 0;
    }

    // Parsing was successful, allocate the material
    theMaterial = new Steel2(iData[0], dData[0], dData[1], dData[2], 
			      dData[3], dData[4], dData[5], dData[6], 
			      dData[7], dData[8], dData[9]);    

  } else if (numData == 11) {
    if (OPS_GetDoubleInput(&numData, dData) != 0) {
      opserr << "Invalid arggs: uniaxialMaterial Steel2 " << iData[0] << 
	" fy? E? b? <R0? cR1? cR2? <a1? a2? a3? a4?>>" << endln;
      return 0;
    }

    // Parsing was successful, allocate the material
    theMaterial = new Steel2(iData[0], dData[0], dData[1], dData[2], 
			      dData[3], dData[4], dData[5], dData[6], 
			      dData[7], dData[8], dData[9], dData[10]);    

  }   

  if (theMaterial == 0) {
    opserr << "WARNING could not create uniaxialMaterial of type Steel2 Material\n";
    return 0;
  }

  return theMaterial;
}

Steel2::Steel2(int tag,
	double _Fy, double _E0, double _b,
	double _R0, double _cR1, double _cR2,
	double _a1, double _a2, double _a3, double _a4, double sigInit):
UniaxialMaterial(tag, MAT_TAG_Steel2),
	Fy(_Fy), E0(_E0), b(_b), R0(_R0), cR1(_cR1), cR2(_cR2), a1(_a1), a2(_a2), a3(_a3), a4(_a4), 
	sigini(sigInit)
{
	this->revertToStart();
}

Steel2::Steel2(int tag,
	double _Fy, double _E0, double _b,
	double _R0, double _cR1, double _cR2):
UniaxialMaterial(tag, MAT_TAG_Steel2),
	Fy(_Fy), E0(_E0), b(_b), R0(_R0), cR1(_cR1), cR2(_cR2), sigini(0.0)
{
	// Default values for no isotropic hardening
	a1 = 0.0;
	a2 = 1.0;
	a3 = 0.0;
	a4 = 1.0;

	this->revertToStart();
}

Steel2::Steel2(int tag, double _Fy, double _E0, double _b):
UniaxialMaterial(tag, MAT_TAG_Steel2),
	Fy(_Fy), E0(_E0), b(_b), sigini(0.0)
{
	// Default values for elastic to hardening transitions
	R0 = 15.0;
	cR1 = 0.925;
	cR2 = 0.15;

	// Default values for no isotropic hardening
	a1 = 0.0;
	a2 = 1.0;
	a3 = 0.0;
	a4 = 1.0;

	this->revertToStart();
}

Steel2::Steel2(void):
UniaxialMaterial(0, MAT_TAG_Steel2)
{
	Fy = 0.0;
	E0 = 0.0;
	b = 0.0;
	sigini = 0.0;

	// Default values for elastic to hardening transitions
	R0 = 15.0;
	cR1 = 0.925;
	cR2 = 0.15;

	// Default values for no isotropic hardening
	a1 = 0.0;
	a2 = 1.0;
	a3 = 0.0;
	a4 = 1.0;

	this->revertToStart();
}

Steel2::~Steel2(void)
{
	// Does nothing
}

UniaxialMaterial*
	Steel2::getCopy(void)
{
	Steel2 *theCopy = new Steel2(this->getTag(), Fy, E0, b, R0, cR1, cR2, a1, a2, a3, a4, sigini);

	return theCopy;
}

double
	Steel2::getInitialTangent(void)
{
	return E0;
}

int
	Steel2::setTrialStrain(double trialStrain, double strainRate)
{
	double Esh = b * E0;
	double epsy = Fy / E0;
	double xi, R, epsrat, dum1, dum2, sigTemp;

	this->revertToLastCommit();

	if (sigini != 0.0) {
		double epsini = sigini/E0;
		eps = trialStrain+epsini;
	} else {
		eps = trialStrain;
	}

	double deps = eps - epsP;


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

				Cepss0 = epss0;
				Csigs0 = sigs0;
				Cepsr = epsr;
				Csigr = sigr;
				Cepspl = epspl;

			} else {
				kon = 1;
				epss0 = epsmax;
				sigs0 = Fy;
				epspl = epsmax;
				
				Tepsr = epsr;
				Tsigr = sigr;
				Tepss0 = epss0;
				Tsigs0 = sigs0;
				Tepspl = epspl;
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

	if (fabs(deps) < 10.0*DBL_EPSILON) {
		
		return 0;

	} else if (kon == 2 && deps > 0.0) {

		kon = 1;

		if ( (fabs(eps - epsr) <= 0.5*epsy) && (eP >= 0.99*E0) ) {	// ignore the incoming curve

			// see if we can use the old curve (calculate sigTemp at epsP)
			xi     = fabs((Tepspl-Tepss0)/epsy);
			R      = R0*(1.0 - (cR1*xi)/(cR2+xi));
			epsrat = (epsP-Tepsr)/(Tepss0-Tepsr);
			dum1  = 1.0 + pow(fabs(epsrat),R);
			dum2  = pow(dum1,(1/R));

			sigTemp   = b*epsrat +(1.0-b)*epsrat/dum2;
			sigTemp   = sigTemp*(Tsigs0-Tsigr)+Tsigr;

			if (sigP <= sigTemp) {

				// use old curve
				epsr = Tepsr;
				sigr = Tsigr;
				epss0 = Tepss0;
				sigs0 = Tsigs0;
				epspl = Tepspl;

			} else {

				// make new tension curve (don't save compression curve)
				epsr = epsP;
				sigr = sigP;
				//epsmin = min(epsP, epsmin);
				if (epsP < epsmin)
					epsmin = epsP;

				double shft = 1;
				double d1;
				if (epsmax > -epsmin) {
					d1 = (epsmax/(a4*epsy)) - 1;
				} else {
					d1 = (-epsmin/(a4*epsy)) - 1;
				}
				if (d1 > 0) {
					shft = 1.0 + a3 * d1;
				}

				epss0 = (Fy * shft - Esh * epsy * shft - sigr + E0 * epsr) / (E0 - Esh);
				sigs0 = Fy * shft + Esh * (epss0 - epsy * shft);
				epspl = epsmax;

			}

		} else {
			// save older compression curve as decent
			Cepsr = epsr;
			Csigr = sigr;
			Cepss0 = epss0;
			Csigs0 = sigs0;
			Cepspl = epspl;

			// make new tension curve
			epsr = epsP;
			sigr = sigP;

			if (epsP < epsmin)
				epsmin = epsP;

			double shft = 1;
			double d1;
			if (epsmax > -epsmin) {
				d1 = (epsmax/(a4*epsy)) - 1;
			} else {
				d1 = (-epsmin/(a4*epsy)) - 1;
			}
			if (d1 > 0) {
				shft = 1.0 + a3 * d1;
			}

			epss0 = (Fy * shft - Esh * epsy * shft - sigr + E0 * epsr) / (E0 - Esh);
			sigs0 = Fy * shft + Esh * (epss0 - epsy * shft);
			epspl = epsmax;
		}

	} else if (kon == 1 && deps < 0.0) {

		// update the maximum previous strain, store the last load reversal 
		// point and calculate the stress and strain (sigs0 and epss0) at the 
		// new intersection between elastic and strain hardening asymptote 
		// To include isotropic strain hardening shift the strain hardening 
		// asymptote by sigsft before calculating the intersection point 
		// Constants a1 and a2 control this stress shift on compression side 

		kon = 2;
		
		if ( (fabs(eps - epsr) <= 0.5*epsy) && (eP >= 0.99*E0) ) {	// ignore the incoming curve

			// see if we can use the old curve (calculate sigTemp at epsP)
			xi     = fabs((Cepspl-Cepss0)/epsy);
			R      = R0*(1.0 - (cR1*xi)/(cR2+xi));
			epsrat = (epsP-Cepsr)/(Cepss0-Cepsr);
			dum1  = 1.0 + pow(fabs(epsrat),R);
			dum2  = pow(dum1,(1/R));

			sigTemp   = b*epsrat +(1.0-b)*epsrat/dum2;
			sigTemp   = sigTemp*(Csigs0-Csigr)+Csigr;

			if (sigP >= sigTemp) {

				// use old curve
				epsr = Cepsr;
				sigr = Csigr;
				epss0 = Cepss0;
				sigs0 = Csigs0;
				epspl = Cepspl;

			} else {

				// make new compression curve (don't save tension curve)
				epsr = epsP;
				sigr = sigP;
				//      epsmax = max(epsP, epsmax);
				if (epsP > epsmax)
					epsmax = epsP;

				double shft = 1;
				double d1;
				if (epsmax > -epsmin) {
					d1 = (epsmax/(a2*epsy)) - 1;
				} else {
					d1 = (-epsmin/(a2*epsy)) - 1;
				}
				if (d1 > 0) {
					shft = 1.0 + a1 * d1;
				}

				epss0 = (-Fy * shft + Esh * epsy * shft - sigr + E0 * epsr) / (E0 - Esh);
				sigs0 = -Fy * shft + Esh * (epss0 + epsy * shft);
				epspl = epsmin;

			}

		} else {
			// save older tension curve as decent
			Tepsr = epsr;
			Tsigr = sigr;
			Tepss0 = epss0;
			Tsigs0 = sigs0;
			Tepspl = epspl;

			// make new compression curve
			epsr = epsP;
			sigr = sigP;
			//      epsmax = max(epsP, epsmax);
			if (epsP > epsmax)
				epsmax = epsP;

			double shft = 1;
			double d1;
			if (epsmax > -epsmin) {
				d1 = (epsmax/(a2*epsy)) - 1;
			} else {
				d1 = (-epsmin/(a2*epsy)) - 1;
			}
			if (d1 > 0) {
				shft = 1.0 + a1 * d1;
			}

			epss0 = (-Fy * shft + Esh * epsy * shft - sigr + E0 * epsr) / (E0 - Esh);
			sigs0 = -Fy * shft + Esh * (epss0 + epsy * shft);
			epspl = epsmin;
		}
	}


	// calculate current stress sig and tangent modulus E 

	xi     = fabs((epspl-epss0)/epsy);
	R      = R0*(1.0 - (cR1*xi)/(cR2+xi));
	epsrat = (eps-epsr)/(epss0-epsr);
	dum1  = 1.0 + pow(fabs(epsrat),R);
	dum2  = pow(dum1,(1/R));

	sig   = b*epsrat +(1.0-b)*epsrat/dum2;
	sig   = sig*(sigs0-sigr)+sigr;

	double newE = (sig - sigP)/(eps - epsP);
	if (newE > E0) {
		sig = E0*(eps - epsP) + sigP;
		e = E0;
	} else {
		e = b + (1.0-b)/(dum1*dum2);
		e = e*(sigs0-sigr)/(epss0-epsr);
	}

	return 0;
}



double 
	Steel2::getStrain(void)
{
	return eps;
}

double 
	Steel2::getStress(void)
{
	return sig;
}

double 
	Steel2::getTangent(void)
{
	return e;
}

int 
	Steel2::commitState(void)
{
	epsminP = epsmin;
	epsmaxP = epsmax;

	epsplP = epspl;
	epss0P = epss0;
	sigs0P = sigs0;
	epssrP = epsr;
	sigsrP = sigr;

	TepsplP = Tepspl;
	Tepss0P = Tepss0;
	Tsigs0P = Tsigs0;
	TepssrP = Tepsr;
	TsigsrP = Tsigr;

	CepsplP = Cepspl;
	Cepss0P = Cepss0;
	Csigs0P = Csigs0;
	CepssrP = Cepsr;
	CsigsrP = Csigr;

	konP = kon;
	eP = e;
	sigP = sig;
	epsP = eps;

	return 0;
}

int 
	Steel2::revertToLastCommit(void)
{
	epsmin = epsminP;
	epsmax = epsmaxP;

	epspl = epsplP;
	epss0 = epss0P;
	sigs0 = sigs0P;
	epsr = epssrP;
	sigr = sigsrP;

	Cepspl = CepsplP;
	Cepss0 = Cepss0P;
	Csigs0 = Csigs0P;
	Cepsr = CepssrP;
	Csigr = CsigsrP;

	Tepspl = TepsplP;
	Tepss0 = Tepss0P;
	Tsigs0 = Tsigs0P;
	Tepsr = TepssrP;
	Tsigr = TsigsrP;

	kon = konP;
	e = eP;
	sig = sigP;
	eps = epsP;

	return 0;
}

int 
	Steel2::revertToStart(void)
{	
	konP = 0;
	eP = E0;
	epsP = 0.0;
	sigP = 0.0;

	epsmaxP = Fy/E0;
	epsminP = -epsmaxP;
	epsplP = 0.0;
	epss0P = 0.0;
	sigs0P = 0.0;
	epssrP = 0.0;
	sigsrP = 0.0;

	TepsplP = epsmaxP;
	Tepss0P = epsmaxP;
	Tsigs0P = Fy;
	TepssrP = 0.0;
	TsigsrP = 0.0;

	CepsplP = epsminP;
	Cepss0P = epsminP;
	Csigs0P = -Fy;
	CepssrP = 0.0;
	CsigsrP = 0.0;

	if (sigini != 0.0) {
		epsP = sigini/E0;
		sigP = sigini;
	} 

	this->revertToLastCommit();
	return 0;
}

int 
	Steel2::sendSelf(int commitTag, Channel &theChannel)
{
	static Vector data(34);
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

	data(23) = TepsplP;
	data(24) = Tepss0P;
	data(25) = Tsigs0P;
	data(26) = TepssrP;
	data(27) = TsigsrP;

	data(28) = CepsplP;
	data(29) = Cepss0P;
	data(30) = Csigs0P;
	data(31) = CepssrP;
	data(32) = CsigsrP;

	if (theChannel.sendVector(this->getDbTag(), commitTag, data) < 0) {
		opserr << "Steel2::sendSelf() - failed to sendSelf\n";
		return -1;
	}
	return 0;
}

int 
	Steel2::recvSelf(int commitTag, Channel &theChannel, 
	FEM_ObjectBroker &theBroker)
{
	static Vector data(34);

	if (theChannel.recvVector(this->getDbTag(), commitTag, data) < 0) {
		opserr << "Steel2::recvSelf() - failed to recvSelf\n";
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

	TepsplP = data(23); 
	Tepss0P = data(24); 
	Tsigs0P = data(25); 
	TepssrP = data(26); 
	TsigsrP = data(27); 

	CepsplP = data(28); 
	Cepss0P = data(29); 
	Csigs0P = data(30); 
	CepssrP = data(31); 
	CsigsrP = data(32); 

	this->revertToLastCommit();

	return 0;
}

void 
	Steel2::Print(OPS_Stream &s, int flag)
{
	s << "Steel2:(strain, stress, tangent) " << eps << " " << sig << " " << e << endln;
}
