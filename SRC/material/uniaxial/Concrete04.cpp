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
// $Date: 2005-09-23 22:51:27 $
// $Source: /usr/local/cvs/OpenSees/SRC/material/uniaxial/Concrete04.cpp,v $
                                                                        
                                                                        
// File: ~/material/Concrete04.cpp
//
// Written: N.Mitra (nmitra@u.washington.edu) 
// Created: 09/04
// Revision: A
//
// Description: This file contains the class implementation for 
// Concrete04 based on Popovics pre and post compression curve 
// for concrete. 
//
// What: "@(#) Concrete04.C, revA"
// Revision 1. Adding in exponential tensile envelope for concrete
// Dt. 05-16-05


#include <Concrete04.h>
#include <Vector.h>
#include <Matrix.h>
#include <Channel.h>
#include <Information.h>
#include <math.h>
#include <float.h>


Concrete04::Concrete04
(int tag, double FPC, double EPSC0, double EPSCU, double EC0, double FCT, double ETU)
  :UniaxialMaterial(tag, MAT_TAG_Concrete04),
   fpc(FPC), epsc0(EPSC0), epscu(EPSCU), Ec0(EC0), fct(FCT), etu(ETU), beta(0.1),
   CminStrain(0.0), CendStrain(0.0), CcompStrain(0.0), CUtenStress(FCT),
   Cstrain(0.0), Cstress(0.0) 
{
	// Make all concrete parameters negative
	if (fpc > 0.0 || epsc0 > 0.0 || epscu > 0.0) {
		opserr << "error: negative values required for concrete stress-strain model" << endln;
	}

	if (fct < 0.0) {
		fct = 0.0;
		opserr << "warning: fct less than 0.0 so the tensile response part is being set to 0" << endln;
	}

	Ctangent = Ec0;
	CunloadSlope = Ec0;
	CUtenSlope = Ec0;

	// Set trial values
	this->revertToLastCommit();
}

Concrete04::Concrete04
(int tag, double FPC, double EPSC0, double EPSCU, double EC0, double FCT, double ETU, double BETA)
  :UniaxialMaterial(tag, MAT_TAG_Concrete04),
   fpc(FPC), epsc0(EPSC0), epscu(EPSCU), Ec0(EC0), fct(FCT), etu(ETU), beta(BETA),
   CminStrain(0.0), CendStrain(0.0), CcompStrain(0.0), CUtenStress(FCT),
   Cstrain(0.0), Cstress(0.0) 
{
	// Make all concrete parameters negative
	if (fpc > 0.0 || epsc0 > 0.0 || epscu > 0.0) {
		opserr << "error: negative values required for concrete stress-strain model" << endln;
	}

	if (fct < 0.0) {
		fct = 0.0;
		opserr << "warning: fct less than 0.0 so the tensile response part is being set to 0" << endln;
	}

	Ctangent = Ec0;
	CunloadSlope = Ec0;
	CUtenSlope = Ec0;

	// Set trial values
	this->revertToLastCommit();
}


Concrete04::Concrete04
(int tag, double FPC, double EPSC0, double EPSCU, double EC0)
  :UniaxialMaterial(tag, MAT_TAG_Concrete04),
   fpc(FPC), epsc0(EPSC0), epscu(EPSCU), Ec0(EC0), fct(0.0), etu(0.0), beta(0.0),
   CminStrain(0.0), CendStrain(0.0), CcompStrain(0.0), CUtenStress(0.0),
   Cstrain(0.0), Cstress(0.0) 
{
	// Make all concrete parameters negative
	if (fpc > 0.0 || epsc0 > 0.0 || epscu > 0.0) {
		opserr << "error: negative values required for concrete stress-strain model" << endln;
	}

	Ctangent = Ec0;
	CunloadSlope = Ec0;
	CUtenSlope = 0.0;

	// Set trial values
	this->revertToLastCommit();

}

Concrete04::Concrete04():UniaxialMaterial(0, MAT_TAG_Concrete04),
 fpc(0.0), epsc0(0.0), epscu(0.0), Ec0(0.0), fct(0.0), etu(0.0), beta(0.0),
 CminStrain(0.0), CunloadSlope(0.0), CendStrain(0.0), CcompStrain(0.0), CUtenStress(0.0),
 CUtenSlope(0.0), Cstrain(0.0), Cstress(0.0)
{
	// Set trial values
	this->revertToLastCommit();

}

Concrete04::~Concrete04()
{
   // Does nothing
}

int Concrete04::setTrialStrain (double strain, double strainRate)
{
  // Set trial strain
  Tstrain = strain;
  
  // Determine change in strain from last converged state
  double dStrain = Tstrain - Cstrain;

  if (fabs(dStrain) < DBL_EPSILON)   
    return 0;
  
  // Calculate the trial state given the change in strain
  // determineTrialState (dStrain);
  TunloadSlope = CunloadSlope;
  TcompStrain = CcompStrain;
  TUtenStress = CUtenStress;
  TUtenSlope = CUtenSlope;

//  double dpStrain = Tstrain - TcompStrain;  // June 07
  
  
  if (dStrain <= 0.0) {
	  // Material can be either in Compression-Reloading
	  // or Tension-Unloading state.
//      if (dpStrain > 0.0) { // June 07
	  if (Tstrain >= 0.0) {   
          // Material is in Tension-Unloading State

		  Ttangent = TUtenSlope;
		  Tstress = Cstress + dStrain*TUtenSlope;
		  if (Tstress < 0.0 || Tstrain == 0.0) {
			  Tstress = 0.0;
		  }
//		  if (Tstress > -DBL_EPSILON)
//			  Tstress = 0.0;

	  } else {
		  // Material is in Compression-Reloading State
            TminStrain = CminStrain;
		    TendStrain = CendStrain;
            double tempStress = Cstress + TunloadSlope*dStrain;
		    CompReload();

			if (tempStress > Tstress)
			{
				Tstress = tempStress;
				Ttangent = TunloadSlope;
			}
			if (Tstress >= 0.0) {
				Tstress = 0.0;
			}
	  }
  } else {
	  // Material can be either in Compression-Unloading
	  // or Tension-Reloading State.
 //     if (dpStrain > 0.0) { // June 07
	  if (Tstrain >= 0.0) {
          // Material is in Tension-Reloading State
          TminStrain = CminStrain;
		  TensReload();

	  } else {
		  // Material is in Compression-Unloading State
            Ttangent = TunloadSlope;
		    Tstress = Cstress + Ttangent*dStrain;
			if (Tstress > 0.0)
				Tstress = 0.0;
	  }
  }
  return 0;
}


void Concrete04::CompReload()
{
if (Tstrain <= TminStrain) {
	TminStrain = Tstrain;

	// Determine point on compression envelope
	CompEnvelope();

	// set the unload envelope in compression
	setCompUnloadEnv();
}
	else if (Tstrain <= TendStrain) {
		Ttangent = TunloadSlope;
		Tstress = Ttangent*(Tstrain-TendStrain);
	} 	else {
		Tstress = 0.0;
		Ttangent = 0.0;
	}
}

void Concrete04::CompEnvelope()
{
	if (Tstrain > epscu) {
		double Esec = fpc/epsc0;
		double r = 0.0;
		if (Esec >= Ec0) {
			r = 400.0;
		} else {
			r = Ec0/(Ec0-Esec);
		}
		double eta = Tstrain/epsc0;
		Tstress = fpc*eta*r/(r-1+pow(eta,r));
		Ttangent = fpc*r*(r-1)*(1-pow(eta,r))/(pow((r-1+pow(eta,r)),2)*epsc0);
	  } else {
		    Tstress = 0.0;
		    Ttangent = 0.0;
	  }
}

void Concrete04::setCompUnloadEnv()
{
	double tempStrain = TminStrain;

	if (tempStrain < epscu)
		tempStrain = epscu;

	double eta = tempStrain/epsc0;

	double ratio = 0.707*(eta-2.0) + 0.834; // unloading parameter as per Karsan-Jirsa

	if (eta < 2.0)
		ratio = 0.145*eta*eta + 0.13*eta;

	TendStrain = ratio*epsc0;

	double temp1 = TminStrain - TendStrain;

	double temp2 = Tstress/Ec0;

	if (temp1 > -DBL_EPSILON) {	// temp1 should always be negative
		TunloadSlope = Ec0;
	}
	else if (temp1 <= temp2) {
		TendStrain = TminStrain - temp1;
		TunloadSlope = Tstress/temp1;
		TcompStrain = 0.0; //TendStrain; // June 07
	}
	else {
		TendStrain = TminStrain - temp2;
		TunloadSlope = Ec0;
		TcompStrain = 0.0; //TendStrain; // June 07
	}	

}

void Concrete04::TensReload()
{
	// Determine point on tension envelope and 
	// also determine tension unload parameters
//	if (Tstrain > TminStrain) {
	//TcompStrain = TendStrain;
		TminStrain = Tstrain;
	    TensEnvelope();
		setTenUnload();
/*		if (Tstrain < TUtenStress/TUtenSlope) {
			Tstress = TUtenSlope*(Tstrain - TcompStrain);
			Ttangent = TUtenSlope;
		}
//	} 
/*	else if (Tstrain >= TcompStrain) {
		Ttangent = TUtenSlope;
		Tstress = Ttangent*(Tstrain-TUtenSlope);
	} 
	else {

		  Ttangent = 0.0;
		  Tstress = 0.0;
	}
*/
}

void Concrete04::TensEnvelope()
{

	if (Tstrain >= TUtenStress/TUtenSlope && Tstrain <= (etu+TcompStrain)) {
//		if (Tstrain <= (etu+TcompStrain)) {
/*		double tempStress = TUtenSlope*(Tstrain - TcompStrain);
		if (Tstrain < TUtenStress/TUtenSlope) {
			Tstress = tempStress;
			Ttangent = TUtenSlope;
		} else {
*/			double tenSlp = fct/Ec0;
			double nom = Tstrain-tenSlp;   // chnage this formulae for genaralised response
			double denom = -tenSlp+etu;
			Tstress = fct*pow(beta,nom/denom);
			Ttangent = Tstress*log(beta)/denom;
			if (nom <= 0.0) {
				Tstress = 0.0;
				Ttangent = 0.0;
			}
//		}
//		TUtenStress = Tstress;
	//	TUtenSlope = Tstress/(Tstrain - TcompStrain);
	  } 
	else if (Tstrain < TUtenStress/TUtenSlope)
	{
		Tstress = TUtenSlope*(Tstrain - TcompStrain);
		Ttangent = TUtenSlope;
	}
	else {
		    Tstress = 0.0;
		    Ttangent = 0.0;
//			TUtenStress = Tstress;
	//		TUtenSlope = 0.0;
	  }
//			TUtenStress = Tstress;
//		TUtenSlope = Tstress/(Tstrain - TcompStrain);

}

void Concrete04::setTenUnload()
{
/*	double tempStrain = TminStrain;
	if (tempStrain > (etu+TcompStrain))
		tempStrain = etu+TcompStrain;
	double tempStress = TUtenSlope*(Tstrain - TcompStrain);
		if (tempStress > TUtenStress) {
			TUtenStress = Tstress;
			TUtenSlope = Tstress/(tempStrain - TcompStrain);
		}
	*/


	if (Tstrain > TUtenStress/TUtenSlope) {
		TUtenStress = Tstress;
		TUtenSlope = Tstress/(Tstrain - TcompStrain);
	}
/*	else {
		TUtenStress = CUtenStress;
		TUtenSlope = CUtenSlope;
	}
	else {
		Tstress = TUtenSlope*(Tstrain - TcompStrain);
		Ttangent = TUtenSlope;
	}
	*/
}
double Concrete04::getStress ()
{
   return Tstress;
}

double Concrete04::getStrain ()
{
   return Tstrain;
}

double Concrete04::getTangent ()
{
   return Ttangent;
}

int Concrete04::commitState ()
{
   // History variables
   CminStrain = TminStrain;
   CunloadSlope = TunloadSlope;
   CendStrain = TendStrain;
   CcompStrain = TcompStrain;
   CUtenStress = TUtenStress;
   CUtenSlope = TUtenSlope;

   // State variables
   Cstrain = Tstrain;
   Cstress = Tstress;
   Ctangent = Ttangent;

   return 0;
}

int Concrete04::revertToLastCommit ()
{
   // Reset trial history variables to last committed state
   TminStrain = CminStrain;
   TendStrain = CendStrain;
   TunloadSlope = CunloadSlope;
   TcompStrain = CcompStrain;
   TUtenStress = CUtenStress;
   TUtenSlope = CUtenSlope;

   // Recompute trial stress and tangent
   Tstrain = Cstrain;
   Tstress = Cstress;
   Ttangent = Ctangent;

   return 0;
}

int Concrete04::revertToStart ()
{

   // History variables
   CminStrain = 0.0;
   CunloadSlope = Ec0;
   CendStrain = 0.0;
   CcompStrain = 0.0;
   CUtenStress = fct;
   CUtenSlope = Ec0;

   // State variables
   Cstrain = 0.0;
   Cstress = 0.0;
   Ctangent = Ec0;

   // Reset trial variables and state
   this->revertToLastCommit();

   return 0;
}

UniaxialMaterial* Concrete04::getCopy ()
{
   Concrete04* theCopy = new Concrete04(this->getTag(),
                                    fpc, epsc0, epscu, Ec0, fct, etu, beta);

   // Converged history variables
   theCopy->CminStrain = CminStrain;
   theCopy->CunloadSlope = CunloadSlope;
   theCopy->CendStrain = CendStrain;
   theCopy->CcompStrain = CcompStrain;
   theCopy->CUtenStress = CUtenStress;
   theCopy->CUtenSlope = CUtenSlope;

   // Converged state variables
   theCopy->Cstrain = Cstrain;
   theCopy->Cstress = Cstress;
   theCopy->Ctangent = Ctangent;

   return theCopy;
}

int Concrete04::sendSelf (int commitTag, Channel& theChannel)
{
   int res = 0;
   static Vector data(15);
   data(0) = this->getTag();

   // Material properties
   data(1) = fpc;
   data(2) = epsc0;
   data(3) = epscu;
   data(4) = Ec0;
   data(5) = fct;

   // History variables from last converged state
   data(6) = CminStrain;
   data(7) = CunloadSlope;
   data(8) = CendStrain;
   data(9) = CcompStrain;
   data(10) = CUtenStress;
   data(11) = CUtenSlope;

   // State variables from last converged state
   data(12) = Cstrain;
   data(13) = Cstress;
   data(14) = Ctangent;

   // Data is only sent after convergence, so no trial variables
   // need to be sent through data vector

   res = theChannel.sendVector(this->getDbTag(), commitTag, data);
   if (res < 0) 
      opserr << "Concrete04::sendSelf() - failed to send data\n";

   return res;
}

int Concrete04::recvSelf (int commitTag, Channel& theChannel,
                                 FEM_ObjectBroker& theBroker)
{
   int res = 0;
   static Vector data(15);
   res = theChannel.recvVector(this->getDbTag(), commitTag, data);

   if (res < 0) {
      opserr << "Concrete04::recvSelf() - failed to receive data\n";
      this->setTag(0);      
   }
   else {
      this->setTag(int(data(0)));

      // Material properties 
      fpc = data(1);
      epsc0 = data(2);
      epscu = data(3);
      Ec0 = data(4);
	  fct = data(5);

      // History variables from last converged state
      CminStrain = data(6);
      CunloadSlope = data(7);
      CendStrain = data(8);
      CcompStrain = data(9);
      CUtenStress = data(10);
      CUtenSlope = data(11);

      // State variables from last converged state
      Cstrain = data(12);
      Cstress = data(13);
      Ctangent = data(14);

      // Set trial state variables
      Tstrain = Cstrain;
      Tstress = Cstress;
      Ttangent = Ctangent;
   }

   return res;
}

void Concrete04::Print (OPS_Stream& s, int flag)
{
   s << "Concrete04, tag: " << this->getTag() << endln;
   s << "  fpc: " << fpc << endln;
   s << "  epsc0: " << epsc0 << endln;
   s << "  fct: " << fct << endln;
   s << "  epscu: " << epscu << endln;
   s << "  Ec0:  " << Ec0 << endln;
   s << "  etu:  " << etu << endln;
   s << "  beta: " << beta << endln;
}

// LOWES: add functions for variable hinge-length model
int
Concrete04::getMaterialType()
{
	return 0;
}
// LOWES: end
