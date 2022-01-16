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
                                                                                                                                               
//
// Written by: Won Lee of Stanford University
// Created: 09/04

// $Revision: 1.1 $
// $Date: 2007-02-02 22:58:36 $
// $Source: /usr/local/cvs/OpenSees/SRC/material/uniaxial/ECC01.cpp,v $

// Description: This file contains the class implementation for 
// ECC01. 
//   - ECC model based on Han et al. model
//      (Han TS, Feenstra PH, Billington SL, ACI Structural Journal,
//			Nov-Dec 2003, "Simulation of Highly Ductile Fiber Reinforced
//			Cement-Based Composite Components Under Cyclic Loading")
//


#include "ECC01.h"
#include <Vector.h>
#include <Matrix.h>
#include <Channel.h>
#include <Information.h>
#include <math.h>
#include <float.h>
#include <elementAPI.h>

void * OPS_ADD_RUNTIME_VPV(OPS_ECC01)
{
    int numdata = OPS_GetNumRemainingInputArgs();
    if (numdata < 15) {
	opserr << "WARNING insufficient arguments\n";
	opserr << "Want: uniaxialMaterial ECC01 TAG? SIGT0? EPST0? SIGT1? EPST1? EPST2? SIGC0? EPSC0? EPSC1? ";
	opserr << "ALPHAT1? ALPHAT2? ALPHAC? ALPHACU? BETAT? BETAC\n";
	return 0;
    }

    int tag;
    numdata = 1;
    if (OPS_GetIntInput(&numdata,&tag) < 0) {
	return 0;
    }

    double data[14];
    numdata = 14;
    if (OPS_GetDoubleInput(&numdata,data)) {
	return 0;
    }

    UniaxialMaterial* mat = new ECC01(tag,data[0],data[1],data[2],data[3],data[4],data[5],data[6],data[7],data[8],data[9],data[10],data[11],data[12],data[13]);
    if (mat == 0) {
	opserr << "WARNING: failed to create ECC01 material\n";
	return 0;
    }

    return mat;
}

ECC01::ECC01
(int tag, double SIGT0, double EPST0, double SIGT1, double EPST1, double EPST2, double SIGC0, 
   double EPSC0, double EPSC1, double ALPHAT1, double ALPHAT2, double ALPHAC, double ALPHACU, double BETAT, double BETAC)
  :UniaxialMaterial(tag, MAT_TAG_ECC01),
   sigt0(SIGT0), epst0(EPST0), sigt1(SIGT1), epst1(EPST1), 
   epst2(EPST2), sigc0(SIGC0), epsc0(EPSC0), epsc1(EPSC1),
   alphaT1(ALPHAT1), alphaT2(ALPHAT2), alphaC(ALPHAC), alphaCU(ALPHACU), betaT(BETAT), betaC(BETAC),
   CminStrain(0.0), CmaxStrain(0.0), 
   Cstrain(0.0), Cstress(0.0), 
   Cstmp(0.0), Cetmp(0.0), Cindex(0), TmaxStrain(0.0), TminStrain(0.0), Tindex(0) 
{
	
	// Make all compressive parameters negative

	if (sigc0 > 0.0)
		sigc0 = -sigc0;

	if (epsc0 > 0.0)
		epsc0 = -epsc0;

	if (epsc1 > 0.0)
		epsc1 = -epsc1;

	// Initial tangent
	double Ec0 = sigc0/epsc0;
	Ctangent = Ec0;
	Ttangent = Ec0;

	// Set trial values
	this->revertToLastCommit();


}

ECC01::ECC01():UniaxialMaterial(0, MAT_TAG_ECC01),  
 sigt0(0.0), epst0(0.0), sigt1(0.0), epst1(0.0),    
 epst2(0.0), sigc0(0.0), epsc0(0.0), epsc1(0.0),
 alphaT1(0.0), alphaT2(0.0), alphaC(0.0), alphaCU(0.0), betaT(0.0), betaC(0.0),
 CminStrain(0.0), CmaxStrain(0.0), 
 Cstrain(0.0), Cstress(0.0),
 Cstmp(0.0), Cetmp(0.0), Cindex(0), TmaxStrain(0.0), TminStrain(0.0), Tindex(0)  
{
	// Set trial values
	this->revertToLastCommit();

}

ECC01::~ECC01 ()
{
   // Does nothing
}


int ECC01::setTrialStrain (double strain, double strainRate)
{
  double sigmax =0.0, epstul =0.0, sigmin =0.0, epscul =0.0;

  // Set trial strain
  Tstrain = strain;

  // update max and min values
  if (Tstrain > TmaxStrain) {
	  TmaxStrain = Tstrain;
  }
  if (Tstrain < TminStrain) {
	  TminStrain = Tstrain;
  }
  
  double dStrain = Tstrain - Cstrain;

  if (fabs(dStrain) < DBL_EPSILON)   
    return 0;
    
  // TENSION
  if (Tstrain > 0.0) {
	  // loading  in tension
	  if (TmaxStrain <= Tstrain) {
		  if (Tstrain <= epst0) {
			  Tindex = 1;
		  }
		  else if (Tstrain <= epst1) {
			  Tindex = 2;
		  }
		  else if (Tstrain <= epst2) {
			  Tindex = 3;
		  }
		  else {
			  Tindex = 4;
		  }
	  }
	  else {
		  // unloading/reloading in tension (hardening, first branch)
		  if (TmaxStrain <= epst0) {
			  Tindex = 1;
		  }
		  // unloading/reloading in tension (hardening, second branch)
		  else if (TmaxStrain <= epst1) {
			  // unloading (tension:hardening second branch)
			  epstul = betaT*(TmaxStrain-epst0); 
			  sigmax = sigt0 + (sigt1-sigt0)*(TmaxStrain-epst0)/(epst1-epst0);
			  if (Tstrain <= Cstrain) {
				  if (Tstrain <= epstul) {
					  Tindex = 9;
				  }
				  else {
					if (Cindex == 2) {
					  Tstmp = sigmax;
					  Tetmp = TmaxStrain;
					}
					else if (Cindex == 7) {
					  Tstmp = Cstress;
					  Tetmp = Cstrain;
					}
					Tindex = 5;
				  }
			  }
			  else {// reloading (tension:hardening second branch)
				  if (Tstrain <= epstul) {
					  Tindex = 9;
				  }
				  else {
					  if (Cindex == 5) {
						  Tstmp = Cstress;
						  Tetmp = Cstrain;
					  }
					  else if ((Cindex == 9) || (Cindex <= -1)) {
						  Tstmp = 0.0;
						  Tetmp = epstul;
					  }
					  Tindex = 7;
				  }
			  }
		  }
		  //unloading/reloading in tension (softening region)
		  else if (TmaxStrain <= epst2) {
			  epstul = betaT*(epst1-epst0);
			  sigmax = sigt1*(1.0-(TmaxStrain-epst1)/(epst2-epst1));
			  //unloading (tension:softening)
			  if (Tstrain <= Cstrain) {
				  if (Tstrain <= epstul) {
					  Tindex =9;
				  }
				  else {
					  if (Cindex == 3) {
						  Tstmp = sigmax;
						  Tetmp = TmaxStrain;
					  }
					  else if (Cindex ==8) {
						  Tstmp = Cstress;  
						  Tetmp = Cstrain;
					  }
					  Tindex = 6;
				  }
			  }
			  // reloading (tension:softening)
			  else {
				  if (Tstrain <= epstul) {
					  Tindex = 9;
				  }
				  else {
					  if (Cindex == 6) {
						  Tstmp = Cstress;  
						  Tetmp = Cstrain;   
					  }
					  else if (Cindex == 9) {
						  Tstmp = 0.0;
						  Tetmp = epstul;
					  }
					  Tindex = 8;
				  }
			  }
		  }
		  else {
			  if (Tstrain <= epst2) {
				  Tindex = 9;
			  }
			  else {
				  Tindex = 4;
			  }
		  }
	  }
		  
	  
  }

  else { // COMPRESSION
	  //Loading in compression
	  if (TminStrain >= Tstrain) {
		  if (Tstrain >= epsc0) {
			  Tindex = -1;
		  }
		  else if (Tstrain >= epsc1) {
			  Tindex = -2;
		  }
		  else {
			  Tindex = -3;
		  }
	  }
	  else {
		  //unloading/reloading in compression
		  if (TminStrain >= epsc0) {
			  // unloading/reloading in compression:pre-peak
			  Tindex = -1;
		  }
		  else if (TminStrain >= epsc1) {
			  //unloading compression:post-peak
			  epscul = betaC*(TminStrain-epsc0);
			  sigmin = sigc0*pow(((TminStrain-epsc1)/(epsc0-epsc1)),alphaCU); 
			  //                     
			  if (Tstrain >= Cstrain) {
				  if (Tstrain >= epscul) {
					  Tindex = -6;
				  }
				  else {
					  if (Cindex == -2) {
						  Tstmp = sigmin;
						  Tetmp = TminStrain;
					  }
					  else if (Cindex == -5) {
						  Tstmp = Cstress; 
						  Tetmp = Cstrain;
					  }
					  Tindex = -4;
				  }
			  }
			  // reloading compression:post-peak
			  else {
				  if (Tstrain >= epscul) {
					  Tindex = -6;
				  }
				  else {
					  if (Cindex == -4) {
						  Tstmp = Cstress;  
						  Tetmp = Cstrain;
					  }
					  else if ((Cindex == -6) || (Cindex >= 1)) {
						  Tstmp = 0.0;
						  Tetmp = epscul;
					  }
					  Tindex = -5;
				  }
			  }
		  }
		  else  {
			  if (Tstrain >= epsc1) {
				  Tindex = -6;
			  }
			  else {
				  Tindex = -3;
			  }
		  }
	  }

  }

  ECCGetStressAndStiffness (Tindex, sigmax, epstul, sigmin, epscul);

  return 0;

}

int ECC01::setTrial (double strain, double &stress, double &tangent, double strainRate)
{

  double sigmax =0.0, epstul =0.0, sigmin =0.0, epscul =0.0;

  // Set trial strain
  Tstrain = strain;
  
  if (Tstrain > TmaxStrain) {
	  TmaxStrain = Tstrain;
  }
  if (Tstrain < TminStrain) {
	  TminStrain = Tstrain;
  }
  
  double dStrain = Tstrain - Cstrain;

  if (fabs(dStrain) < DBL_EPSILON) {
    tangent = Ttangent;  
    stress = Tstress;  
    return 0;
  }
  
   
  // TENSION
  if (Tstrain > 0.0) {
	  // loading  in tension
	  if (TmaxStrain <= Tstrain) {
		  if (Tstrain <= epst0) {
			  Tindex = 1;
		  }
		  else if (Tstrain <= epst1) {
			  Tindex = 2;
		  }
		  else if (Tstrain <= epst2) {
			  Tindex = 3;
		  }
		  else {
			  Tindex = 4;
		  }
	  }
	  else {
		  // unloading/reloading in tension (hardening, first branch)
		  if (TmaxStrain <= epst0) {
			  Tindex = 1;
		  }
		  // unloading/reloading in tension (hardening, second branch)
		  else if (TmaxStrain <= epst1) {
			  // unloading (tension:hardening second branch)
			  epstul = betaT*(TmaxStrain-epst0); 
			  sigmax = sigt0 + (sigt1-sigt0)*(TmaxStrain-epst0)/(epst1-epst0);
			  if (Tstrain <= Cstrain) {
				  if (Tstrain <= epstul) {
					  Tindex = 9;
				  }
				  else {
					if (Cindex == 2) {
					  Tstmp = sigmax;
					  Tetmp = TmaxStrain;
					}
					else if (Cindex == 7) {
					  Tstmp = Cstress;
					  Tetmp = Cstrain;
					}
					Tindex = 5;
				  }
			  }
			  else {  // reloading (tension:hardening second branch)
				  if (Tstrain <= epstul) {
					  Tindex = 9;
				  }
				  else {
					  if (Cindex == 5) {
						  Tstmp = Cstress;
						  Tetmp = Cstrain;
					  }
					  else if ((Cindex == 9) || (Cindex <= -1)) {
						  Tstmp = 0.0;
						  Tetmp = epstul;
					  }
					  Tindex = 7;
				  }
			  }
		  }
		  //unloading/reloading in tension (softening region)
		  else if (TmaxStrain <= epst2) {
			  epstul = betaT*(epst1-epst0);
			  sigmax = sigt1*(1.0-(TmaxStrain-epst1)/(epst2-epst1));
			  //unloading (tension:softening)
			  if (Tstrain <= Cstrain) {
				  if (Tstrain <= epstul) {
					  Tindex =9;
				  }
				  else {
					  if (Cindex == 3) {
						  Tstmp = sigmax;
						  Tetmp = TmaxStrain;
					  }
					  else if (Cindex ==8) {
						  Tstmp = Cstress;  
						  Tetmp = Cstrain;
					  }
					  Tindex = 6;
				  }
			  }
			  // reloading (tension:softening)
			  else {
				  if (Tstrain <= epstul) {
					  Tindex = 9;
				  }
				  else {
					  if (Cindex == 6) {
						  Tstmp = Cstress;  
						  Tetmp = Cstrain;   
					  }
					  else if (Cindex == 9) {
						  Tstmp = 0.0;
						  Tetmp = epstul;
					  }
					  Tindex = 8;
				  }
			  }
		  }
		  else {
			  if (Tstrain <= epst2) {
				  Tindex = 9;
			  }
			  else {
				  Tindex = 4;
			  }
		  }
	  }
		  
  }

  else { // if it is compression
	  //Loading in compression
	  if (TminStrain >= Tstrain) {
		  if (Tstrain >= epsc0) {
			  Tindex = -1;
		  }
		  else if (Tstrain >= epsc1) {
			  Tindex = -2;
		  }
		  else {
			  Tindex = -3;
		  }
	  }
	  else {
		  //unloading/reloading in compression
		  if (TminStrain >= epsc0) {
			  // unloading/reloading in compression:pre-peak
			  Tindex = -1;
		  }
		  else if (TminStrain >= epsc1) {
			  //unloading compression:post-peak
			  epscul = betaC*(TminStrain-epsc0);                          
			  sigmin = sigc0*pow(((TminStrain-epsc1)/(epsc0-epsc1)),alphaCU); 
			  if (Tstrain >= Cstrain) {
				  if (Tstrain >= epscul) {
					  Tindex = -6;
				  }
				  else {
					  if (Cindex == -2) {
						  Tstmp = sigmin;
						  Tetmp = TminStrain;
					  }
					  else if (Cindex == -5) {
						  Tstmp = Cstress; 
						  Tetmp = Cstrain;
					  }
					  Tindex = -4;
				  }
			  }
			  // reloading compression:post-peak
			  else {
				  if (Tstrain >= epscul) {
					  Tindex = -6;
				  }
				  else {
					  if (Cindex == -4) {
						  Tstmp = Cstress;  
						  Tetmp = Cstrain;
					  }
					  else if ((Cindex == -6) || (Cindex >= 1)) {
						  Tstmp = 0.0;
						  Tetmp = epscul;
					  }
					  Tindex = -5;
				  }
			  }
		  }
		  else  {
			  if (Tstrain >= epsc1) {
				  Tindex = -6;
			  }
			  else {
				  Tindex = -3;
			  }
		  }
	  }

	  //return 0;
  }

  ECCGetStressAndStiffness (Tindex, sigmax, epstul, sigmin, epscul);
  stress = Tstress;  // 
  tangent = Ttangent;  

  return 0;

}


void ECC01::ECCGetStressAndStiffness (int index, double sigmax, double epstul, double sigmin, double epscul)
{
	// anywhere on the envelope curve
	if ((Tindex >= -3) && (Tindex <= 4)) {
		envelope ();
	}
	// tension region
	else if (Tindex == 5) {

	  if (Tetmp-epstul != 0.0) {
	    Tstress = Tstmp*pow(((Tstrain-epstul)/(Tetmp-epstul)),alphaT1);
	    Ttangent = alphaT1*Tstmp*pow(((Tstrain-epstul)/(Tetmp-epstul)),(alphaT1-1))*(1/(Tetmp-epstul));
	  }

	}
	else if (Tindex== 6) {
	  if (Tetmp-epstul != 0.0) {
		Tstress = Tstmp*pow(((Tstrain-epstul)/(Tetmp-epstul)),alphaT2);
		Ttangent = alphaT2*Tstmp*pow(((Tstrain-epstul)/(Tetmp-epstul)),(alphaT2-1))*(1/(Tetmp-epstul));
	  }
	}
	else if (Tindex== 7) {
	  if (TmaxStrain-Tetmp != 0.0) {
		Tstress = Tstmp + (sigmax-Tstmp)*(Tstrain-Tetmp)/(TmaxStrain-Tetmp);
		Ttangent = (sigmax-Tstmp)/(TmaxStrain-Tetmp);
	  }
	}
	else if (Tindex== 8) {
	  if (TmaxStrain-Tetmp != 0.0) {
		Tstress = Tstmp + (sigmax-Tstmp)*(Tstrain-Tetmp)/(TmaxStrain-Tetmp);
		Ttangent = (sigmax-Tstmp)/(TmaxStrain-Tetmp);
	  }
	}
	else if (Tindex== 9) {
		Tstress = 0.0;
		Ttangent = 0.0;
	}
	// compression region
	else if (Tindex== -4) {
	  if (Tetmp-epscul != 0.0) {
		Tstress = Tstmp*pow(((Tstrain-epscul)/(Tetmp-epscul)),alphaC);
		Ttangent = alphaC*Tstmp*pow(((Tstrain-epscul)/(Tetmp-epscul)),(alphaC-1))*(1/(Tetmp-epscul));
	  }
	}
	else if (Tindex== -5) {
	  if (TminStrain-Tetmp != 0.0) {
		Tstress = Tstmp + (sigmin-Tstmp)*(Tstrain-Tetmp)/(TminStrain-Tetmp);
		Ttangent = (sigmin-Tstmp)/(TminStrain-Tetmp);
	  }
	}
	else if (Tindex== -6) {
		Tstress = 0.0;
		Ttangent = 0.0;
	}

}


void ECC01::envelope ()
{
	double initialSlope = sigt0/epst0;
	double Ec0 = sigc0/epsc0;

	if (Tstrain > 0) { //WL: if in tension
		if (Tstrain < epst0) {
			Tstress = initialSlope*Tstrain;
			Ttangent = initialSlope;
		}
		else if (Tstrain < epst1) {
			Ttangent = (sigt1-sigt0)/(epst1-epst0);
			Tstress = sigt0 + Ttangent*(Tstrain-epst0);
		}
		else if (Tstrain < epst2) {
			Ttangent = (-sigt1)/(epst2-epst1);
			Tstress = sigt1 + Ttangent*(Tstrain-epst1);
		}
		else {
			Tstress = 0.0;
			Ttangent = 0.0;
		}
	}
	else { // WL: if in compression
		if (Tstrain > epsc0) {
			// for now hardcode in the r coefficient = 5
			Tstress = sigc0*5*(Tstrain/epsc0)*(1/(5-1+pow(Tstrain/epsc0,5)));
			//Ttangent = (1/pow(5-1+pow(Tstrain/epsc0,5),2))*((Tstrain/epsc0)*((1/epsc0)*5*pow(Tstrain/epsc0,5-1) )-(1/epsc0)*(5-1+pow(Tstrain/epsc0,5)));
			//Tstress = Ec0*Tstrain;
			Ttangent = Ec0;
		}
		else if (Tstrain > epsc1) {
			Ttangent = alphaCU*sigc0*pow(((Tstrain-epsc1)/(epsc0-epsc1)),(alphaCU-1))*(1/(epsc0-epsc1));
			Tstress = sigc0*pow(((Tstrain-epsc1)/(epsc0-epsc1)),alphaCU); 
		}
		else {
			Tstress = 0.0;
			Ttangent = 0.0;
		}
	}

}


//  the below has been modified /////////////////////
double ECC01::getStress ()
{
   return Tstress;
}

double ECC01::getStrain ()
{
   return Tstrain;
}

double ECC01::getTangent ()
{
   return Ttangent;
}

int ECC01::commitState ()
{
   // History variables
   CminStrain = TminStrain;
   CmaxStrain = TmaxStrain; 
   Cstmp = Tstmp;  
   Cetmp = Tetmp;  
   Cindex = Tindex;

   // State variables
   Cstrain = Tstrain;
   Cstress = Tstress;
   Ctangent = Ttangent;

   return 0;
}

int ECC01::revertToLastCommit ()
{
   // Reset trial history variables to last committed state
   TminStrain = CminStrain;
   TmaxStrain = CmaxStrain;
   Tstmp = Cstmp;
   Tetmp = Cetmp;
   Tindex = Cindex;

   // Recompute trial stress and tangent
   Tstrain = Cstrain;
   Tstress = Cstress;
   Ttangent = Ctangent;

   return 0;
}

int ECC01::revertToStart ()
{
	double Ec0 = sigc0/epsc0;

   // History variables
   CminStrain = 0.0;
   CmaxStrain = 0.0;
   Cstmp = 0.0;
   Cetmp = 0.0;
   Cindex = 0;

   // State variables
   Cstrain = 0.0;
   Cstress = 0.0;
   Ctangent = Ec0;

   // Reset trial variables and state
   this->revertToLastCommit();

   return 0;
}

UniaxialMaterial* ECC01::getCopy ()
{
   ECC01* theCopy = new ECC01(this->getTag(),
                                    sigt0, epst0, sigt1, epst1, epst2, sigc0, epsc0, epsc1,
									alphaT1, alphaT2, alphaC, alphaCU, betaT, betaC);

   // Converged history variables
   theCopy->CminStrain = CminStrain;
   theCopy->CmaxStrain = CmaxStrain;
   theCopy->Cstmp = Cstmp;
   theCopy->Cetmp = Cetmp;
   theCopy->Cindex = Cindex;

   // Converged state variables
   theCopy->Cstrain = Cstrain;
   theCopy->Cstress = Cstress;
   theCopy->Ctangent = Ctangent;

   return theCopy;
} 

 int ECC01::sendSelf (int commitTag, Channel& theChannel)
{
   int res = 0;
   //static Vector data(11);
   static Vector data(23);
   data(0) = this->getTag();

   // Material properties
   data(1) = sigt0;
   data(2) = epst0;
   data(3) = sigt1;
   data(4) = epst1;
   data(5) = epst2;
   data(6) = sigc0;
   data(7) = epsc0;
   data(8) = epsc1;
   data(9) = alphaT1;
   data(10) = alphaT2;
   data(11) = alphaC;
   data(12) = alphaCU;
   data(13) = betaT;
   data(14) = betaC;

   // History variables from last converged state
   data(15) = CminStrain;
   data(16) = CmaxStrain;
   data(17) = Cstmp;
   data(18) = Cetmp;
   data(19) = Cindex;

   // State variables from last converged state
   data(20) = Cstrain;
   data(21) = Cstress;
   data(22) = Ctangent;

   res = theChannel.sendVector(this->getDbTag(), commitTag, data);
   if (res < 0) 
      opserr << "ECC01::sendSelf() - failed to send data\n";

   return res;
}

int ECC01::recvSelf (int commitTag, Channel& theChannel,
                                 FEM_ObjectBroker& theBroker)
{
   int res = 0;
   static Vector data(23);
   res = theChannel.recvVector(this->getDbTag(), commitTag, data);

   if (res < 0) {
      opserr << "ECC01::recvSelf() - failed to receive data\n";
      this->setTag(0);      
   }
   else {
      this->setTag(int(data(0)));

      // Material properties 
	  sigt0 = data(1);
	  epst0 = data(2);
	  sigt1 = data(3);
	  epst1 = data(4);
	  epst2 = data(5);
	  sigc0 = data(6);
	  epsc0 = data(7);
	  epsc1 = data(8);
	  alphaT1 = data(9);
	  alphaT2 = data(10);
	  alphaC = data(11);
	  alphaCU = data(12);
	  betaT = data(13);
	  betaC = data(14);

      // History variables from last converged state
	  CminStrain = data(15);
	  CmaxStrain = data(16);
	  Cstmp = data(17);
	  Cetmp = data(18);
	  Cindex = data(19);

      // State variables from last converged state
      Cstrain = data(20);
      Cstress = data(21);
      Ctangent = data(22);

      // Set trial state variables
      Tstrain = Cstrain;
      Tstress = Cstress;
      Ttangent = Ctangent;
   }

   return res;
}

void ECC01::Print (OPS_Stream& s, int flag)
{
   s << "ECC01, tag: " << this->getTag() << endln;
   s << "  sigt0: " << sigt0 << endln;
   s << "  epst0: " << epst0 << endln;
   s << "  sigt1: " << sigt1 << endln;
   s << "  epst1: " << epst1 << endln;
   s << "  epst2: " << epst2 << endln;
   s << "  sigc0: " << sigc0 << endln;
   s << "  epsc0: " << epsc0 << endln;
   s << "  epsc1: " << epsc1 << endln;
   s << "  alphaT1: " << alphaT1 << endln;
   s << "  alphaT2: " << alphaT2 << endln;
   s << "  alphaC: " << alphaC << endln;
   s << "  alphaCU: " << alphaCU << endln;
   s << "  betaT: " << betaT << endln;
   s << "  betaC: " << betaC << endln;
}




