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
                                                                        
// $Revision: 1.11 $
// $Date: 2004-07-15 21:34:10 $
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
#include <Matrix.h>
#include <Channel.h>
#include <Information.h>
#include <math.h>
#include <float.h>

int count = 0;

Concrete01::Concrete01
(int tag, double FPC, double EPSC0, double FPCU, double EPSCU)
  :UniaxialMaterial(tag, MAT_TAG_Concrete01),
   fpc(FPC), epsc0(EPSC0), fpcu(FPCU), epscu(EPSCU), 
   CminStrain(0.0), CendStrain(0.0),
   Cstrain(0.0), Cstress(0.0) 
{
  count++;
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
	double Ec0 = 2*fpc/epsc0;
	Ctangent = Ec0;
	CunloadSlope = Ec0;
	Ttangent = Ec0;

	// Set trial values
	this->revertToLastCommit();

// AddingSensitivity:BEGIN /////////////////////////////////////
	parameterID = 0;
	SHVs = 0;
// AddingSensitivity:END //////////////////////////////////////
}

Concrete01::Concrete01():UniaxialMaterial(0, MAT_TAG_Concrete01),
 fpc(0.0), epsc0(0.0), fpcu(0.0), epscu(0.0),
 CminStrain(0.0), CunloadSlope(0.0), CendStrain(0.0),
 Cstrain(0.0), Cstress(0.0)
{
	// Set trial values
	this->revertToLastCommit();

// AddingSensitivity:BEGIN /////////////////////////////////////
	parameterID = 0;
	SHVs = 0;
// AddingSensitivity:END //////////////////////////////////////
}

Concrete01::~Concrete01 ()
{
   // Does nothing
}

int Concrete01::setTrialStrain (double strain, double strainRate)
{
  // Set trial strain
  Tstrain = strain;
  
  // check for a quick return
  if (Tstrain > 0.0) {
    Tstress = 0;
    Ttangent = 0;
    return 0;
  }
  
  // Determine change in strain from last converged state
  double dStrain = Tstrain - Cstrain;

  if (fabs(dStrain) < DBL_EPSILON)   
    return 0;
  
  // Calculate the trial state given the change in strain
  // determineTrialState (dStrain);
  TunloadSlope = CunloadSlope;
  
  double tempStress = Cstress + TunloadSlope*dStrain;
  
  // Material goes further into compression
  if (dStrain <= 0.0) {
    TminStrain = CminStrain;
    TendStrain = CendStrain;
    
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
  
  return 0;
}



int 
Concrete01::setTrial (double strain, double &stress, double &tangent, double strainRate)
{
  // Set trial strain
  Tstrain = strain;
  
  // check for a quick return
  if (Tstrain > 0.0) {
    Tstress = 0;
    Ttangent = 0;
    stress = 0;
    tangent = 0;
    return 0;
  }
  
  // Determine change in strain from last converged state
  double dStrain = Tstrain - Cstrain;
  
  if (fabs(dStrain) < DBL_EPSILON) {
    tangent = Ttangent;
    stress = Tstress;
        return 0;
  }

  // Calculate the trial state given the change in strain
  // determineTrialState (dStrain);
  TunloadSlope = CunloadSlope;
  
  double tempStress = Cstress + TunloadSlope*dStrain;
  
  // Material goes further into compression
  if (dStrain <= 0.0) {
    TminStrain = CminStrain;
    TendStrain = CendStrain;
    
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
  
  //opserr << "Concrete01::setTrial() " << strain << " " << tangent << " " << strain << endln;
  
  stress = Tstress;
  tangent =  Ttangent;
  
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
		double Ec0 = 2.0*fpc/epsc0;
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

	double Ec0 = 2.0*fpc/epsc0;

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

int Concrete01::commitState ()
{
   // History variables
   CminStrain = TminStrain;
   CunloadSlope = TunloadSlope;
   CendStrain = TendStrain;

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

   // Recompute trial stress and tangent
   Tstrain = Cstrain;
   Tstress = Cstress;
   Ttangent = Ctangent;

   return 0;
}

int Concrete01::revertToStart ()
{
	double Ec0 = 2.0*fpc/epsc0;

   // History variables
   CminStrain = 0.0;
   CunloadSlope = Ec0;
   CendStrain = 0.0;

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
                                    fpc, epsc0, fpcu, epscu);

   // Converged history variables
   theCopy->CminStrain = CminStrain;
   theCopy->CunloadSlope = CunloadSlope;
   theCopy->CendStrain = CendStrain;

   // Converged state variables
   theCopy->Cstrain = Cstrain;
   theCopy->Cstress = Cstress;
   theCopy->Ctangent = Ctangent;

   return theCopy;
}

int Concrete01::sendSelf (int commitTag, Channel& theChannel)
{
   int res = 0;
   static Vector data(11);
   data(0) = this->getTag();

   // Material properties
   data(1) = fpc;
   data(2) = epsc0;
   data(3) = fpcu;
   data(4) = epscu;

   // History variables from last converged state
   data(5) = CminStrain;
   data(6) = CunloadSlope;
   data(7) = CendStrain;

   // State variables from last converged state
   data(8) = Cstrain;
   data(9) = Cstress;
   data(10) = Ctangent;

   // Data is only sent after convergence, so no trial variables
   // need to be sent through data vector

   res = theChannel.sendVector(this->getDbTag(), commitTag, data);
   if (res < 0) 
      opserr << "Concrete01::sendSelf() - failed to send data\n";

   return res;
}

int Concrete01::recvSelf (int commitTag, Channel& theChannel,
                                 FEM_ObjectBroker& theBroker)
{
   int res = 0;
   static Vector data(11);
   res = theChannel.recvVector(this->getDbTag(), commitTag, data);

   if (res < 0) {
      opserr << "Concrete01::recvSelf() - failed to receive data\n";
      this->setTag(0);      
   }
   else {
      this->setTag(int(data(0)));

      // Material properties 
      fpc = data(1);
      epsc0 = data(2);
      fpcu = data(3);
      epscu = data(4);

      // History variables from last converged state
      CminStrain = data(5);
      CunloadSlope = data(6);
      CendStrain = data(7);

      // State variables from last converged state
      Cstrain = data(8);
      Cstress = data(9);
      Ctangent = data(10);

      // Set trial state variables
      Tstrain = Cstrain;
      Tstress = Cstress;
      Ttangent = Ctangent;
   }

   return res;
}

void Concrete01::Print (OPS_Stream& s, int flag)
{
   s << "Concrete01, tag: " << this->getTag() << endln;
   s << "  fpc: " << fpc << endln;
   s << "  epsc0: " << epsc0 << endln;
   s << "  fpcu: " << fpcu << endln;
   s << "  epscu: " << epscu << endln;
}




// AddingSensitivity:BEGIN ///////////////////////////////////
int
Concrete01::setParameter(const char **argv, int argc, Information &info)
{
	if (strcmp(argv[0],"fc") == 0) {// Compressive strength
		info.theType = DoubleType;
		return 1;
	}
	if (strcmp(argv[0],"epsco") == 0) {// Strain at compressive strength
		info.theType = DoubleType;
		return 2;
	}
	if (strcmp(argv[0],"fcu") == 0) {// Crushing strength
		info.theType = DoubleType;
		return 3;
	}
	if (strcmp(argv[0],"epscu") == 0) {// Strain at crushing strength
		info.theType = DoubleType;
		return 4;
	}
	else {
		opserr << "WARNING: Could not set parameter in Concrete01! " << endln;
		return -1;
	}
}
    
                            

int
Concrete01::updateParameter(int parameterID, Information &info)
{
	switch (parameterID) {
	case 1:
		this->fpc = info.theDouble;
		break;
	case 2:
		this->epsc0 = info.theDouble;
		break;
	case 3:
		this->fpcu = info.theDouble;
		break;
	case 4:
		this->epscu = info.theDouble;
		break;
	default:
		break;
	}
        
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
	double Ec0 = 2*fpc/epsc0;
	Ctangent = Ec0;
	CunloadSlope = Ec0;
	Ttangent = Ec0;

	return 0;
}




int
Concrete01::activateParameter(int passedParameterID)
{
	parameterID = passedParameterID;

	return 0;
}

double
Concrete01::getStressSensitivity(int gradNumber, bool conditional)
{
	// Initialize return value
	double TstressSensitivity = 0.0;
	double dktdh = 0.0;
	double TstrainSensitivity = 0.0;


	// Pick up sensitivity history variables
	double CminStrainSensitivity = 0.0;
	double CunloadSlopeSensitivity = 0.0;
	double CendStrainSensitivity = 0.0;
	double CstressSensitivity = 0.0;
	double CstrainSensitivity = 0.0;
	if (SHVs != 0) {
		CminStrainSensitivity   = (*SHVs)(0,(gradNumber-1));
		CunloadSlopeSensitivity = (*SHVs)(1,(gradNumber-1));
		CendStrainSensitivity   = (*SHVs)(2,(gradNumber-1));
		CstressSensitivity      = (*SHVs)(3,(gradNumber-1));
		CstrainSensitivity      = (*SHVs)(4,(gradNumber-1));
	}


	// Assign values to parameter derivatives (depending on what's random)
	double fpcSensitivity = 0.0;
	double epsc0Sensitivity = 0.0;
	double fpcuSensitivity = 0.0;
	double epscuSensitivity = 0.0;

	if (parameterID == 1) {
		fpcSensitivity = 1.0;
	}
	else if (parameterID == 2) {
		epsc0Sensitivity = 1.0;
	}
	else if (parameterID == 3) {
		fpcuSensitivity = 1.0;
	}
	else if (parameterID == 4) {
		epscuSensitivity = 1.0;
	}


	// Strain increment 
	double dStrain = Tstrain - Cstrain;

	// Evaluate stress sensitivity 
	if (dStrain < 0.0) {					// applying more compression to the material

		if (Tstrain < CminStrain) {			// loading along the backbone curve

			if (Tstrain > epsc0) {			//on the parabola
				
				TstressSensitivity = fpcSensitivity*(2.0*Tstrain/epsc0-(Tstrain/epsc0)*(Tstrain/epsc0))
					      + fpc*( (2.0*TstrainSensitivity*epsc0-2.0*Tstrain*epsc0Sensitivity)/(epsc0*epsc0) 
						  - 2.0*(Tstrain/epsc0)*(TstrainSensitivity*epsc0-Tstrain*epsc0Sensitivity)/(epsc0*epsc0));
				
				dktdh = 2.0*((fpcSensitivity*epsc0-fpc*epsc0Sensitivity)/(epsc0*epsc0))
					  * (1.0-Tstrain/epsc0)
					  - 2.0*(fpc/epsc0)*(TstrainSensitivity*epsc0-Tstrain*epsc0Sensitivity)
					  / (epsc0*epsc0);
			}
			else if (Tstrain > epscu) {		// on the straight inclined line
//cerr << "ON THE STRAIGHT INCLINED LINE" << endl;

				dktdh = ( (fpcSensitivity-fpcuSensitivity)
					  * (epsc0-epscu) 
					  - (fpc-fpcu)
					  * (epsc0Sensitivity-epscuSensitivity) )
					  / ((epsc0-epscu)*(epsc0-epscu));

				double kt = (fpc-fpcu)/(epsc0-epscu);

				TstressSensitivity = fpcSensitivity 
					      + dktdh*(Tstrain-epsc0)
						  + kt*(TstrainSensitivity-epsc0Sensitivity);
			}
			else {							// on the horizontal line
//cerr << "ON THE HORIZONTAL LINES" << endl;
				TstressSensitivity = fpcuSensitivity;
				dktdh = 0.0;
			
			}
		}
		else if (Tstrain < CendStrain) {	// reloading after an unloading that didn't go all the way to zero stress
//cerr << "RELOADING AFTER AN UNLOADING THAT DIDN'T GO ALL THE WAY DOWN" << endl;
			TstressSensitivity = CunloadSlopeSensitivity * (Tstrain-CendStrain)
				      + CunloadSlope * (TstrainSensitivity-CendStrainSensitivity);

			dktdh = CunloadSlopeSensitivity;
		}
		else {

			TstressSensitivity = 0.0;
			dktdh = 0.0;

		}
	}
	else if (Cstress+CunloadSlope*dStrain<0.0) {// unloading, but not all the way down to zero stress
//cerr << "UNLOADING, BUT NOT ALL THE WAY DOWN" << endl;
		TstressSensitivity = CstressSensitivity 
			               + CunloadSlopeSensitivity*dStrain
				           + CunloadSlope*(TstrainSensitivity-CstrainSensitivity);

		dktdh = CunloadSlopeSensitivity;
	}
	else {									// unloading all the way down to zero stress
//cerr << "UNLOADING ALL THE WAY DOWN" << endl;

		TstressSensitivity = 0.0;
		dktdh = 0.0;

	}

	return TstressSensitivity;
}



int
Concrete01::commitSensitivity(double TstrainSensitivity, int gradNumber, int numGrads)
{

	// Initialize unconditaional stress sensitivity
	double TstressSensitivity = 0.0;
	double dktdh = 0.0;


	// Assign values to parameter derivatives (depending on what's random)
	double fpcSensitivity = 0.0;
	double epsc0Sensitivity = 0.0;
	double fpcuSensitivity = 0.0;
	double epscuSensitivity = 0.0;

	if (parameterID == 1) {
		fpcSensitivity = 1.0;
	}
	else if (parameterID == 2) {
		epsc0Sensitivity = 1.0;
	}
	else if (parameterID == 3) {
		fpcuSensitivity = 1.0;
	}
	else if (parameterID == 4) {
		epscuSensitivity = 1.0;
	}


	// Pick up sensitivity history variables
	double CminStrainSensitivity = 0.0;
	double CunloadSlopeSensitivity = 0.0;
	double CendStrainSensitivity = 0.0;
	double CstressSensitivity = 0.0;
	double CstrainSensitivity = 0.0;
	
	if (SHVs == 0) {
		SHVs = new Matrix(5,numGrads);
		CunloadSlopeSensitivity = (2.0*fpcSensitivity*epsc0-2.0*fpc*epsc0Sensitivity) / (epsc0*epsc0);
	}
	else {
		CminStrainSensitivity   = (*SHVs)(0,(gradNumber-1));
		CunloadSlopeSensitivity = (*SHVs)(1,(gradNumber-1));
		CendStrainSensitivity   = (*SHVs)(2,(gradNumber-1));
		CstressSensitivity      = (*SHVs)(3,(gradNumber-1));
		CstrainSensitivity      = (*SHVs)(4,(gradNumber-1));
	}


	// Strain increment 
	double dStrain = Tstrain - Cstrain;

	// Evaluate stress sensitivity 
	if (dStrain < 0.0) {					// applying more compression to the material

		if (Tstrain < CminStrain) {			// loading along the backbone curve

			if (Tstrain > epsc0) {			//on the parabola
				
				TstressSensitivity = fpcSensitivity*(2.0*Tstrain/epsc0-(Tstrain/epsc0)*(Tstrain/epsc0))
					      + fpc*( (2.0*TstrainSensitivity*epsc0-2.0*Tstrain*epsc0Sensitivity)/(epsc0*epsc0) 
						  - 2.0*(Tstrain/epsc0)*(TstrainSensitivity*epsc0-Tstrain*epsc0Sensitivity)/(epsc0*epsc0));
				
				dktdh = 2.0*((fpcSensitivity*epsc0-fpc*epsc0Sensitivity)/(epsc0*epsc0))
					  * (1.0-Tstrain/epsc0)
					  - 2.0*(fpc/epsc0)*(TstrainSensitivity*epsc0-Tstrain*epsc0Sensitivity)
					  / (epsc0*epsc0);
			}
			else if (Tstrain > epscu) {		// on the straight inclined line

				dktdh = ( (fpcSensitivity-fpcuSensitivity)
					  * (epsc0-epscu) 
					  - (fpc-fpcu)
					  * (epsc0Sensitivity-epscuSensitivity) )
					  / ((epsc0-epscu)*(epsc0-epscu));

				double kt = (fpc-fpcu)/(epsc0-epscu);

				TstressSensitivity = fpcSensitivity 
					      + dktdh*(Tstrain-epsc0)
						  + kt*(TstrainSensitivity-epsc0Sensitivity);
			}
			else {							// on the horizontal line

				TstressSensitivity = fpcuSensitivity;
				dktdh = 0.0;
			
			}
		}
		else if (Tstrain < CendStrain) {	// reloading after an unloading that didn't go all the way to zero stress

			TstressSensitivity = CunloadSlopeSensitivity * (Tstrain-CendStrain)
				      + CunloadSlope * (TstrainSensitivity-CendStrainSensitivity);

			dktdh = CunloadSlopeSensitivity;
		}
		else {

			TstressSensitivity = 0.0;
			dktdh = 0.0;

		}
	}
	else if (Cstress+CunloadSlope*dStrain<0.0) {// unloading, but not all the way down to zero stress
	
		TstressSensitivity = CstressSensitivity 
			               + CunloadSlopeSensitivity*dStrain
				           + CunloadSlope*(TstrainSensitivity-CstrainSensitivity);

		dktdh = CunloadSlopeSensitivity;
	}
	else {									// unloading all the way down to zero stress

		TstressSensitivity = 0.0;
		dktdh = 0.0;

	}

	// Commit some history variables
	(*SHVs)(3,(gradNumber-1)) = TstressSensitivity;
	(*SHVs)(4,(gradNumber-1)) = TstrainSensitivity;





	// Possibly update history variables for the three ordinary history variable derivatives
	double epsTemp, epsTempSensitivity;
	double eta, etaSensitivity;
	double ratio, ratioSensitivity;
	double temp1, temp1Sensitivity;
	double temp2, temp2Sensitivity;
	double TminStrainSensitivity;
	double TunloadSlopeSensitivity;
	double TendStrainSensitivity;

	if (dStrain<0.0 && Tstrain<CminStrain) {

		TminStrainSensitivity = TstrainSensitivity;

		if (Tstrain < epscu) {

			epsTemp = epscu; 

			epsTempSensitivity = epscuSensitivity;

		}
		else {

			epsTemp = Tstrain;

			epsTempSensitivity = TstrainSensitivity;
		}

		eta = epsTemp/epsc0;

		etaSensitivity = (epsTempSensitivity*epsc0-epsTemp*epsc0Sensitivity) / (epsc0*epsc0);

		if (eta < 2.0) {

			ratio = 0.145 * eta*eta + 0.13*eta;

			ratioSensitivity = 0.29 * eta * etaSensitivity + 0.13 * etaSensitivity;

		}
		else {

			ratio = 0.707*(eta-2.0) + 0.834;

			ratioSensitivity = 0.707 * etaSensitivity;
		}

		temp1 = Tstrain - ratio * epsc0;

		temp1Sensitivity = TstrainSensitivity - ratioSensitivity * epsc0
			                                  - ratio * epsc0Sensitivity;

		temp2 = Tstress * epsc0 / (2.0*fpc); 
		
		temp2Sensitivity = (2.0*fpc*(TstressSensitivity*epsc0+Tstress*epsc0Sensitivity)
			-2.0*Tstress*epsc0*fpcSensitivity) / (4.0*fpc*fpc);

		if (temp1 == 0.0) {

			TunloadSlopeSensitivity = (2.0*fpcSensitivity*epsc0-2.0*fpc*epsc0Sensitivity) / (epsc0*epsc0);
		}
		else if (temp1 < temp2) {

			TendStrainSensitivity = TstrainSensitivity - temp1Sensitivity;

			TunloadSlopeSensitivity = (TstressSensitivity*temp1-Tstress*temp1Sensitivity) / (temp1*temp1);

		}
		else {

			TendStrainSensitivity = TstrainSensitivity - temp2Sensitivity;

			TunloadSlopeSensitivity = (2.0*fpcSensitivity*epsc0-2.0*fpc*epsc0Sensitivity) / (epsc0*epsc0);
		}
	}
	else {
		TminStrainSensitivity = CminStrainSensitivity;
		TunloadSlopeSensitivity = CunloadSlopeSensitivity;
		TendStrainSensitivity = CendStrainSensitivity;
	}



	(*SHVs)(0,(gradNumber-1)) = TminStrainSensitivity;
	(*SHVs)(1,(gradNumber-1)) = TunloadSlopeSensitivity;
	(*SHVs)(2,(gradNumber-1)) = TendStrainSensitivity;

	return 0;
}

























/* THE OLD METHODS:
double
Concrete01::getStressSensitivity(int gradNumber, bool conditional)
{
	// (This method only works for path-independent problems for now.)

	double gradient = 0.0;

	if ( parameterID == 0 ) {
		// Leave the gradient as zero if nothing is random here;
	}
	else if (Tstrain > 0.0 ) {
		gradient = 0.0;
	}
	else if (Tstrain > epsc0) {					// IN PARABOLIC AREA

		if ( parameterID == 1 ) {		// d{sigma}d{fpc}
			gradient = 2.0*Tstrain/epsc0-Tstrain*Tstrain/(epsc0*epsc0);
		}
		else if ( parameterID == 2  ) {	// d{sigma}d{epsc0}
			gradient = 2.0*fpc/(epsc0*epsc0)*(Tstrain*Tstrain/epsc0-Tstrain);
		}
		else if ( parameterID == 3  ) {	// d{sigma}d{fpcu}
			gradient = 0.0;
		}
		else if ( parameterID == 4  ) {	// d{sigma}d{epscu}
			gradient = 0.0;
		}
		else {
			gradient = 0.0;
		}
	}
	else if (Tstrain > epscu) {					// IN LINEAR AREA

		if ( parameterID == 1 ) {		// d{sigma}d{fpc}
			gradient = (epscu-Tstrain)/(epscu-epsc0);
		}
		else if ( parameterID == 2  ) {	// d{sigma}d{epsc0}
			gradient = (fpc-fpcu)*(epscu-Tstrain)/((epscu-epsc0)*(epscu-epsc0));
		}
		else if ( parameterID == 3  ) {	// d{sigma}d{fpcu}
			gradient = (Tstrain-epsc0)/(epscu-epsc0);
		}
		else if ( parameterID == 4  ) {	// d{sigma}d{epscu}
			gradient = (Tstrain-epsc0)*(fpc-fpcu)/((epsc0-epscu)*(epsc0-epscu));
		}
		else {
			gradient = 0.0;
		}
	}
	else {										// IN ZERO STIFFNESS AREA

		if ( parameterID == 1 ) {		// d{sigma}d{fpc}
			gradient = 0.0;
		}
		else if ( parameterID == 2  ) {	// d{sigma}d{epsc0}
			gradient = 0.0;
		}
		else if ( parameterID == 3  ) {	// d{sigma}d{fpcu}
			gradient = 1.0;
		}
		else if ( parameterID == 4  ) {	// d{sigma}d{epscu}
			gradient = 0.0;
		}
		else {
			gradient = 0.0;
		}
	}

	return gradient;
}



int
Concrete01::commitSensitivity(double TstrainSensitivity, int gradNumber, int numGrads)
{
	// Not treated yet. 

	return 0;
}
*/
// AddingSensitivity:END /////////////////////////////////////////////
