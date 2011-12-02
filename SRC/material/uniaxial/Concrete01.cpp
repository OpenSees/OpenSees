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
                                                                        
// $Revision: 1.1.1.1 $
// $Date: 2000-09-15 08:23:22 $
// $Source: /usr/local/cvs/OpenSees/SRC/material/uniaxial/Concrete01.cpp,v $
                                                                        
                                                                        
// File: ~/material/Concrete01.C
//
// Written: MHS 
// Created: 06/99
// Revision: A
//
// Description: This file contains the class implementation for 
// Concrete01. 
//
// What: "@(#) Concrete01.C, revA"


#include <Concrete01.h>
#include <Vector.h>
#include <Channel.h>
#include <math.h>
#include <float.h>

Concrete01::Concrete01
(int tag, double FPC, double EPSC0, double FPCU,
 double EPSCU, double min, double max) :
   UniaxialMaterial(tag, MAT_TAG_Concrete01),
   fpc(FPC), epsc0(EPSC0), fpcu(FPCU), epscu(EPSCU), 
   epsmin(min), epsmax(max), Cstrain(0.0), Cstress(0.0), 
   CminStrain(0.0), CendStrain(0.0), Cfailed(0)
{
	// Make all concrete parameters negative
	if (fpc > 0.0)
		fpc = -fpc;

	if (epsc0 > 0.0)
		epsc0 = -epsc0;

	if (fpcu > 0.0)
		fpcu = -fpcu;

	if (epscu > 0.0)
		epscu = -epscu;

	// Initial tangent
	Ec0 = 2*fpc/epsc0;
	Ctangent = Ec0;
	CunloadSlope = Ec0;

	// Set trial values
	this->revertToLastCommit();
}

Concrete01::Concrete01():UniaxialMaterial(0, MAT_TAG_Concrete01),
 fpc(0.0), epsc0(0.0), fpcu(0.0), epscu(0.0), epsmin(0.0), epsmax(0.0),
 Cstrain(0.0), Cstress(0.0), CminStrain(0.0), CunloadSlope(0.0),
 CendStrain(0.0), Cfailed(0)
{
	// Set trial values
	this->revertToLastCommit();
}

Concrete01::~Concrete01 ()
{
   // Does nothing
}

int Concrete01::setTrialStrain (double strain, double strainRate)
{
	// Early return
	if (Cfailed)
		return 0;

	// Set trial strain
	Tstrain = strain;

	// Another early return
	if (Tstrain <= epsmin || Tstrain >= epsmax) {
		Tfailed = 1;
		Tstress = 0.0;
		Ttangent = 0.0;
		return 0;
	}

	// Determine change in strain from last converged state
	double dStrain = Tstrain - Cstrain;

	// Calculate the trial state given the change in strain
	determineTrialState (dStrain);

	return 0;
}

void Concrete01::determineTrialState (double dStrain)
{  
	TminStrain = CminStrain;
	TendStrain = CendStrain;
	TunloadSlope = CunloadSlope;

	double tempStress = Cstress + TunloadSlope*dStrain;

	// Material goes further into compression
	if (dStrain <= 0.0) {
		
		reload ();
  
		if (tempStress > Tstress) {
			Tstress = tempStress;
			Ttangent = TunloadSlope;
		}
	}

	// Material goes TOWARD tension
	else if (tempStress <= 0.0) {
		Tstress = tempStress;
		Ttangent = TunloadSlope;
	}

	// Made it into tension
	else {
		Tstress = 0.0;
		Ttangent = 0.0;
	}

}

void Concrete01::reload ()
{
	if (Tstrain <= TminStrain) {

		TminStrain = Tstrain;

		// Determine point on envelope
		envelope ();

		unload ();
	}
	else if (Tstrain <= TendStrain) {
		Ttangent = TunloadSlope;
		Tstress = Ttangent*(Tstrain-TendStrain);
	}
	else {
		Tstress = 0.0;
		Ttangent = 0.0;
	}
}

void Concrete01::envelope ()
{
	if (Tstrain > epsc0) {
		double eta = Tstrain/epsc0;
		Tstress = fpc*(2*eta-eta*eta);
		Ttangent = Ec0*(1.0-eta);
	}
	else if (Tstrain > epscu) {
		Ttangent = (fpc-fpcu)/(epsc0-epscu);
		Tstress = fpc + Ttangent*(Tstrain-epsc0);
	}
	else {
		Tstress = fpcu;
		Ttangent = 0.0;
	}
}

void Concrete01::unload ()
{
	double tempStrain = TminStrain;

	if (tempStrain < epscu)
		tempStrain = epscu;

	double eta = tempStrain/epsc0;

	double ratio = 0.707*(eta-2.0) + 0.834;

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
	}
	else {
		TendStrain = TminStrain - temp2;
		TunloadSlope = Ec0;
	}
}

double Concrete01::getStress ()
{
   return Tstress;
}

double Concrete01::getStrain ()
{
   return Tstrain;
}

double Concrete01::getTangent ()
{
   return Ttangent;
}

double Concrete01::getSecant ()
{
	if (fabs(Tstrain) > DBL_EPSILON)
		return Tstress/Tstrain;
	else
		return Ttangent;
}


int Concrete01::commitState ()
{
   // History variables
   CminStrain = TminStrain;
   CunloadSlope = TunloadSlope;
   CendStrain = TendStrain;
   Cfailed = Tfailed;

   // State variables
   Cstrain = Tstrain;
   Cstress = Tstress;
   Ctangent = Ttangent;

   return 0;
}

int Concrete01::revertToLastCommit ()
{
   // Reset trial history variables to last committed state
   TminStrain = CminStrain;
   TendStrain = CendStrain;
   TunloadSlope = CunloadSlope;
   Tfailed = Cfailed;

   // Recompute trial stress and tangent
   Tstrain = Cstrain;
   Tstress = Cstress;
   Ttangent = Ctangent;

   return 0;
}

int Concrete01::revertToStart ()
{
   // History variables
   CminStrain = 0.0;
   CunloadSlope = Ec0;
   CendStrain = 0.0;
   Cfailed = 0;

   // State variables
   Cstrain = 0.0;
   Cstress = 0.0;
   Ctangent = Ec0;

   // Reset trial variables and state
   this->revertToLastCommit();

   return 0;
}

UniaxialMaterial* Concrete01::getCopy ()
{
   Concrete01* theCopy = new Concrete01(this->getTag(),
                                    fpc, epsc0, fpcu, epscu, epsmin, epsmax);

   // Converged history variables
   theCopy->CminStrain = CminStrain;
   theCopy->CunloadSlope = CunloadSlope;
   theCopy->CendStrain = CendStrain;
   theCopy->Cfailed = Cfailed;

   // Converged state variables
   theCopy->Cstrain = Cstrain;
   theCopy->Cstress = Cstress;
   theCopy->Ctangent = Ctangent;

   return theCopy;
}

int Concrete01::sendSelf (int commitTag, Channel& theChannel)
{
   int res = 0;
   static Vector data(14);
   data(0) = this->getTag();

   // Material properties
   data(1) = fpc;
   data(2) = epsc0;
   data(3) = fpcu;
   data(4) = epscu;
   data(5) = epsmin;
   data(6) = epsmax;

   // History variables from last converged state
   data(7) = CminStrain;
   data(8) = CunloadSlope;
   data(9) = CendStrain;
   data(10) = Cfailed;

   // State variables from last converged state
   data(11) = Cstrain;
   data(12) = Cstress;
   data(13) = Ctangent;

   // Data is only sent after convergence, so no trial variables
   // need to be sent through data vector

   res = theChannel.sendVector(this->getDbTag(), commitTag, data);
   if (res < 0) 
      cerr << "Concrete01::sendSelf() - failed to send data\n";

   return res;
}

int Concrete01::recvSelf (int commitTag, Channel& theChannel,
                                 FEM_ObjectBroker& theBroker)
{
   int res = 0;
   static Vector data(14);
   res = theChannel.recvVector(this->getDbTag(), commitTag, data);

   if (res < 0) {
      cerr << "Concrete01::recvSelf() - failed to receive data\n";
      this->setTag(0);      
   }
   else {
      this->setTag(int(data(0)));

      // Material properties 
      fpc = data(1);
      epsc0 = data(2);
      fpcu = data(3);
      epscu = data(4);
      epsmin = data(5);
      epsmax = data(6);

	  Ec0 = 2*fpc/epsc0;

      // History variables from last converged state
      CminStrain = data(7);
      CunloadSlope = data(8);
      CendStrain = data(9);
      Cfailed = int(data(10));

      // State variables from last converged state
      Cstrain = data(11);
      Cstress = data(12);
	  Ctangent = data(13);

	  // Set trial state variables
	  Tstrain = Cstrain;
	  Tstress = Cstress;
	  Ttangent = Ctangent;
   }

   return res;
}

void Concrete01::Print (ostream& s, int flag)
{
   s << "Concrete01, tag: " << this->getTag() << endl;
   s << "  fpc: " << fpc << endl;
   s << "  epsc0: " << epsc0 << endl;
   s << "  fpcu: " << fpcu << endl;
   s << "  epscu: " << epscu << endl;
   if (epsmin != NEG_INF_STRAIN)
      s << "  epsmin: " << epsmin << endl;
   if (epsmax != POS_INF_STRAIN)
      s << "  epsmax: " << epsmax << endl;
}
