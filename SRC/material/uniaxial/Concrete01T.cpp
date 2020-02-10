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
                                                                        
// $Revision: 1.22 $
// $Date: 2008-08-26 16:23:47 $
// $Source: /usr/local/cvs/OpenSees/SRC/material/uniaxial/Concrete01T.cpp,v $
                                                                        
// Written: MHS 
// Created: 06/99
// Revision: A
//
// Description: This file contains the class implementation for 
// Concrete01T. 
//
// What: "@(#) Concrete01T.C, revA"

// Modified by Kristijan Kolozvari to include tension capacity up to cracking.

#include <Concrete01T.h>
#include <Vector.h>
#include <Matrix.h>
#include <Channel.h>
#include <Information.h>
#include <Parameter.h>
#include <string.h>
#include <MaterialResponse.h> // KK

#include <math.h>
#include <float.h>

#include <elementAPI.h>
#include <OPS_Globals.h>

void *
OPS_Concrete01T()
{
  // Pointer to a uniaxial material that will be returned
  UniaxialMaterial *theMaterial = 0;

  int    iData[1];
  double dData[7];
  int numData = 1;

  if (OPS_GetIntInput(&numData, iData) != 0) {
    opserr << "WARNING invalid uniaxialMaterial Concrete01T tag" << endln;
    return 0;
  }

  numData = OPS_GetNumRemainingInputArgs();

  if (numData != 4 && numData != 6 && numData != 7) {
    opserr << "Invalid #args, want: uniaxialMaterial Concrete01T " << iData[0] << "fpc? epsc0? fpcu? epscu? <ft? et0?>\n";
    return 0;
  }

  if (OPS_GetDoubleInput(&numData, dData) != 0) {
    opserr << "Invalid #args, want: uniaxialMaterial Concrete01T " << iData[0] << "fpc? epsc0? fpcu? epscu? <ft? et0?>\n";
    return 0;
  }

  if (numData == 4) {
	  // Parsing was successful, allocate the material
	  theMaterial = new Concrete01T(iData[0], dData[0], dData[1], dData[2], dData[3]);
  }

  if (numData == 6) {
	  // Parsing was successful, allocate the material
	  theMaterial = new Concrete01T(iData[0], dData[0], dData[1], dData[2], dData[3], dData[4], dData[5]);
  }

  if (numData == 7) {

	  int mon = 1;

	  // MAKE IT WORK LIKE FOR ConcreteCM !!!
	  /* numData = 1;
	  if (OPS_GetIntInput(&numData, &mon) != 0) {
		  opserr << "Invalid $mon parameter for uniaxialMaterial with tag  " << iData[0] << endln;
		  return 0;
	  }

	  if (mon != 0 && mon != 1) {
		  opserr << "Invalid $mon parameter for uniaxialMaterial with tag  " << iData[0] << endln;
		  return 0;
	  } */

	  // Parsing was successful, allocate the material
	  theMaterial = new Concrete01T(iData[0], dData[0], dData[1], dData[2], dData[3], dData[4], dData[5], mon);
  }

  if (theMaterial == 0) {
    opserr << "WARNING could not create uniaxialMaterial of type Concrete01T Material\n";
    return 0;
  }

  return theMaterial;
}

// No tensile strength
Concrete01T::Concrete01T
(int tag, double FPC, double EPSC0, double FPCU, double EPSCU)
  :UniaxialMaterial(tag, MAT_TAG_Concrete01T),
   fpc(FPC), epsc0(EPSC0), fpcu(FPCU), epscu(EPSCU), 
   ft(0.0), et(0.0), mon(0), Cracking(0),
   CminStrain(0.0), CendStrain(0.0),
   Cstrain(0.0), Cstress(0.0),
   SMALL(1e-30) // small stiffness
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
  
  if (ft < 0.0)
    ft = -ft;
  
  if (et < 0.0)
    et = -et;

  // Initial tangent
  double Ec0 = 2*fpc/epsc0;
  Ctangent = Ec0;
  CunloadSlope = Ec0;
  Ttangent = Ec0;
  
  Et = 0.0; // KK: no tension (switch)
  
  // Set trial values
  this->revertToLastCommit();
  
  // AddingSensitivity:BEGIN /////////////////////////////////////
  parameterID = 0;
  SHVs = 0;
  // AddingSensitivity:END //////////////////////////////////////
}

// With tensile strength
Concrete01T::Concrete01T
(int tag, double FPC, double EPSC0, double FPCU, double EPSCU, double FT, double ET)
  :UniaxialMaterial(tag, MAT_TAG_Concrete01T),
   fpc(FPC), epsc0(EPSC0), fpcu(FPCU), epscu(EPSCU), 
   ft(FT), et(ET), mon(0), Cracking(0),
   CminStrain(0.0), CendStrain(0.0),
   Cstrain(0.0), Cstress(0.0),
   SMALL(1e-30) // small stiffness
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

  if (ft < 0.0)
    ft = -ft;
  
  if (et < 0.0)
    et = -et;
  
  // Initial tangent
  double Ec0 = 2*fpc/epsc0;
  Ctangent = Ec0;
  CunloadSlope = Ec0;
  Ttangent = Ec0;
  
  Et = ft/et; // KK: linear elastic tension

  // Set trial values
  this->revertToLastCommit();
  
  // AddingSensitivity:BEGIN /////////////////////////////////////
  parameterID = 0;
  SHVs = 0;
  // AddingSensitivity:END //////////////////////////////////////
}

// With tensile strength
Concrete01T::Concrete01T
(int tag, double FPC, double EPSC0, double FPCU, double EPSCU, double FT, double ET, int MON)
  :UniaxialMaterial(tag, MAT_TAG_Concrete01T),
   fpc(FPC), epsc0(EPSC0), fpcu(FPCU), epscu(EPSCU), 
   ft(FT), et(ET), mon(MON), Cracking(0),
   CminStrain(0.0), CendStrain(0.0),
   Cstrain(0.0), Cstress(0.0),
   SMALL(1e-30) // small stiffness
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

  if (ft < 0.0)
    ft = -ft;
  
  if (et < 0.0)
    et = -et;
  
  // Initial tangent
  double Ec0 = 2*fpc/epsc0;
  Ctangent = Ec0;
  CunloadSlope = Ec0;
  Ttangent = Ec0;
  
  Et = ft/et; // KK: linear elastic tension

  // Set trial values
  this->revertToLastCommit();
  
  // AddingSensitivity:BEGIN /////////////////////////////////////
  parameterID = 0;
  SHVs = 0;
  // AddingSensitivity:END //////////////////////////////////////
}

Concrete01T::Concrete01T():UniaxialMaterial(0, MAT_TAG_Concrete01T),
 fpc(0.0), epsc0(0.0), fpcu(0.0), epscu(0.0),
 ft(0.0), et(0.0), mon(0), Cracking(0),
 CminStrain(0.0), CunloadSlope(0.0), CendStrain(0.0),
 Cstrain(0.0), Cstress(0.0),
 SMALL(1e-30) // small stiffness
{
  // Set trial values
  this->revertToLastCommit();
  
  // AddingSensitivity:BEGIN /////////////////////////////////////
  parameterID = 0;
  SHVs = 0;
  // AddingSensitivity:END //////////////////////////////////////
}

Concrete01T::~Concrete01T ()
{
  // Does nothing
}


int Concrete01T::setTrialStrain (double strain, double strainRate)
{
   // Reset trial history variables to last committed state
   TminStrain = CminStrain;
   TendStrain = CendStrain;
   TunloadSlope = CunloadSlope;
   Tstress = Cstress;
   Ttangent = Ctangent;
   Tstrain = Cstrain;

  // Determine change in strain from last converged state
  double dStrain = strain - Cstrain;

  if (fabs(dStrain) < DBL_EPSILON)
    return 0;

  // Set trial strain
  Tstrain = strain;
  
  // check for a quick return
  if (Tstrain > 0.0) {

	  if(Et != 0.0) { // if tension capacity defined

		  if (Cracking == 1) { // concrete has cracked
			  Tstress = 0.0;
			  Ttangent = SMALL;

		  } else {
			  if (Tstrain <= et) {
				  Tstress = Et*Tstrain;
				  Ttangent = Et;

			  } else {
				  Tstress = 0.0;
				  Ttangent = SMALL;

			  }
		  }

	  } else { // if tension capacity not defined
		  Tstress = 0.0;
		  Ttangent = SMALL;
	  }

    return 0;
  }
  
  if (mon == 1) { // monotonic compression 

	  envelope (); // go to envelope each time

	  return 0;
  }

  // Calculate the trial state given the change in strain
  // determineTrialState (dStrain);
  TunloadSlope = CunloadSlope;
  
  double tempStress = Cstress + TunloadSlope*Tstrain - TunloadSlope*Cstrain;
  
  // Material goes further into compression
  if (strain < Cstrain) {
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
	  Ttangent = SMALL;

  }
  
  return 0;
}


int 
	Concrete01T::setTrial (double strain, double &stress, double &tangent, double strainRate)
{
	// Reset trial history variables to last committed state
	TminStrain = CminStrain;
	TendStrain = CendStrain;
	TunloadSlope = CunloadSlope;
	Tstress = Cstress;
	Ttangent = Ctangent;
	Tstrain = Cstrain;

	// Determine change in strain from last converged state
	double dStrain = strain - Cstrain;

	if (fabs(dStrain) < DBL_EPSILON) {
		stress = Tstress;
		tangent = Ttangent;
		return 0;
	}

	// Set trial strain
	Tstrain = strain;

	// check for a quick return
	if (Tstrain > 0.0) {

		if(Et != 0.0) { // if tension capacity defined

			if (Cracking == 1) { // concrete has cracked
				Tstress = 0.0;
				Ttangent = SMALL;
				stress = 0.0;
				tangent = SMALL;

			} else {

				if (Tstrain <= et) {
					Tstress = Et*Tstrain;
					Ttangent = Et;
					stress = Et*Tstrain;
					tangent = Et;

				} else {
					Tstress = 0.0;
					Ttangent = SMALL;
					stress = 0.0;
					tangent = SMALL;
				}
			}

		} else { // if tension capacity not defined
			Tstress = 0.0;
			Ttangent = SMALL;
			stress = 0.0;
			tangent = SMALL;

		}

		return 0;
	}

	if (mon == 1) { // monotonic compression 

		envelope (); // go to envelope each time

		return 0;
	}

	// Calculate the trial state given the change in strain
	// determineTrialState (dStrain);
	TunloadSlope = CunloadSlope;

	double tempStress = Cstress + TunloadSlope*Tstrain - TunloadSlope*Cstrain;

	// Material goes further into compression
	if (strain <= Cstrain) {
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
		Ttangent = SMALL;

	}

	//opserr << "Concrete01T::setTrial() " << strain << " " << tangent << " " << strain << endln;

	stress = Tstress;
	tangent = Ttangent;

	return 0;
}

// NOT USED IN THIS VERSION
void Concrete01T::determineTrialState (double dStrain)
{  
  TminStrain = CminStrain;
  TendStrain = CendStrain;
  TunloadSlope = CunloadSlope;
  
  double tempStress = Cstress + TunloadSlope*dStrain;
  
  // Material goes further into compression
  if (Tstrain <= Cstrain) {
    
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
	  Ttangent = SMALL;

  }
  
}

void Concrete01T::reload ()
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
	  Ttangent = SMALL;

  }

}

void Concrete01T::envelope ()
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
    Ttangent = SMALL;
  }
}

void Concrete01T::unload ()
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

double Concrete01T::getStress ()
{
   return Tstress;
}

double Concrete01T::getStrain ()
{
   return Tstrain;
}

double Concrete01T::getTangent ()
{
   return Ttangent;
}

int Concrete01T::commitState ()
{
   // History variables
   CminStrain = TminStrain;
   CunloadSlope = TunloadSlope;
   CendStrain = TendStrain;

   // State variables
   Cstrain = Tstrain;
   Cstress = Tstress;
   Ctangent = Ttangent;

   	// Concrete has cracked
	if (Cstrain >= et && Cracking == 0 && mon == 0) {
		Cracking = 1;
	}

   return 0;
}

int Concrete01T::revertToLastCommit ()
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

int Concrete01T::revertToStart ()
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

   // Quan April 2006---
   if (SHVs !=0) {SHVs->Zero();}
   parameterID=0;

   return 0;
}

UniaxialMaterial* Concrete01T::getCopy ()
{
   Concrete01T* theCopy = new Concrete01T(this->getTag(),
                                    fpc, epsc0, fpcu, epscu, ft, et, mon);

   // Converged history variables
   theCopy->CminStrain = CminStrain;
   theCopy->CunloadSlope = CunloadSlope;
   theCopy->CendStrain = CendStrain;

   // Converged state variables
   theCopy->Cstrain = Cstrain;
   theCopy->Cstress = Cstress;
   theCopy->Ctangent = Ctangent;
   theCopy->Cracking = Cracking;

   return theCopy;
}

int Concrete01T::sendSelf (int commitTag, Channel& theChannel)
{
   int res = 0;
   static Vector data(16);
   data(0) = this->getTag();

   // Material properties
   data(1) = fpc;
   data(2) = epsc0;
   data(3) = fpcu;
   data(4) = epscu;
   data(5) = ft;
   data(6) = et;
   data(7) = mon;

   // History variables from last converged state
   data(8) = CminStrain;
   data(9) = CunloadSlope;
   data(10) = CendStrain;
   data(11) = Cracking;

   // State variables from last converged state
   data(12) = Cstrain;
   data(13) = Cstress;
   data(14) = Ctangent;

   data(15) = SMALL;

   // Data is only sent after convergence, so no trial variables
   // need to be sent through data vector

   res = theChannel.sendVector(this->getDbTag(), commitTag, data);
   if (res < 0) 
      opserr << "Concrete01T::sendSelf() - failed to send data\n";

   return res;
}

int Concrete01T::recvSelf (int commitTag, Channel& theChannel,
                                 FEM_ObjectBroker& theBroker)
{
   int res = 0;
   static Vector data(16);
   res = theChannel.recvVector(this->getDbTag(), commitTag, data);

   if (res < 0) {
      opserr << "Concrete01T::recvSelf() - failed to receive data\n";
      this->setTag(0);      
   }
   else {
      this->setTag(int(data(0)));

      // Material properties 
      fpc = data(1);
      epsc0 = data(2);
      fpcu = data(3);
      epscu = data(4);
	  ft = data(5);
	  et = data(6);
	  mon = data(7);

      // History variables from last converged state
      CminStrain = data(8);
      CunloadSlope = data(9);
      CendStrain = data(10);
	  Cracking = data(11);

      // State variables from last converged state
      Cstrain = data(12);
      Cstress = data(13);
      Ctangent = data(14);

	  SMALL = data(15);

      // Set trial state variables
      Tstrain = Cstrain;
      Tstress = Cstress;
      Ttangent = Ctangent;
   }

   return res;
}

void Concrete01T::Print (OPS_Stream& s, int flag)
{
   s << "Concrete01T, tag: " << this->getTag() << endln;
   s << "  fpc: " << fpc << endln;
   s << "  epsc0: " << epsc0 << endln;
   s << "  fpcu: " << fpcu << endln;
   s << "  epscu: " << epscu << endln;
   s << "  ft: " << ft << endln;
   s << "  et: " << et << endln;
}

// KK
Response* Concrete01T::setResponse(const char **argv, int argc,
	OPS_Stream &theOutput)
{
	Response *theResponse = 0;

	if (strcmp(argv[0],"getInputParameters") == 0) {
		Vector data1(7);
		data1.Zero();
		theResponse = new MaterialResponse(this, 100, data1);

	} else

		return this->UniaxialMaterial::setResponse(argv, argc, theOutput);

	return theResponse;
}

// KK
int Concrete01T::getResponse(int responseID, Information &matInfo)
{
	if (responseID == 100) {
		matInfo.setVector(this->getInputParameters()); 

	} else

		return this->UniaxialMaterial::getResponse(responseID, matInfo);

	return 0;
}

// KK
Vector Concrete01T::getInputParameters(void)
{
	Vector input_par(7); // size = max number of parameters (assigned + default)

	input_par.Zero();
	
	input_par(0) = this->getTag(); 
	input_par(1) = fpc; 
	input_par(2) = epsc0; 
	input_par(3) = fpcu; 
	input_par(4) = epscu;  
	input_par(5) = ft; 
	input_par(6) = et;  

	return input_par;
}

////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////

// AddingSensitivity:BEGIN ///////////////////////////////////
int
Concrete01T::setParameter(const char **argv, int argc, Parameter &param)
{

  if (strcmp(argv[0],"fc") == 0) {// Compressive strength
    param.setValue(fpc);
    return param.addObject(1, this);
  }
  else if (strcmp(argv[0],"epsco") == 0) {// Strain at compressive strength
    param.setValue(epsc0);
    return param.addObject(2, this);
  }
  else if (strcmp(argv[0],"fcu") == 0) {// Crushing strength
    param.setValue(fpcu);
    return param.addObject(3, this);
  }
  else if (strcmp(argv[0],"epscu") == 0) {// Strain at crushing strength
    param.setValue(epscu);
    return param.addObject(4, this);
  }
  
  return -1;
}          

int
Concrete01T::updateParameter(int parameterID, Information &info)
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
   	TunloadSlope = CunloadSlope;

	return 0;
}

int
Concrete01T::activateParameter(int passedParameterID)
{
	parameterID = passedParameterID;

	return 0;
}

double
Concrete01T::getStressSensitivity(int gradIndex, bool conditional)
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
		CminStrainSensitivity   = (*SHVs)(0,gradIndex);
		CunloadSlopeSensitivity = (*SHVs)(1,gradIndex);
		CendStrainSensitivity   = (*SHVs)(2,gradIndex);
		CstressSensitivity      = (*SHVs)(3,gradIndex);
		CstrainSensitivity      = (*SHVs)(4,gradIndex);
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
Concrete01T::commitSensitivity(double TstrainSensitivity, int gradIndex, int numGrads)
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
		CminStrainSensitivity   = (*SHVs)(0,gradIndex);
		CunloadSlopeSensitivity = (*SHVs)(1,gradIndex);
		CendStrainSensitivity   = (*SHVs)(2,gradIndex);
		CstressSensitivity      = (*SHVs)(3,gradIndex);
		CstrainSensitivity      = (*SHVs)(4,gradIndex);
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
	(*SHVs)(3,gradIndex) = TstressSensitivity;
	(*SHVs)(4,gradIndex) = TstrainSensitivity;





	// Possibly update history variables for the three ordinary history variable derivatives
	double epsTemp, epsTempSensitivity;
	double eta, etaSensitivity;
	double ratio, ratioSensitivity;
	double temp1, temp1Sensitivity;
	double temp2, temp2Sensitivity;
	double TminStrainSensitivity = CminStrainSensitivity;
	double TunloadSlopeSensitivity = CunloadSlopeSensitivity;
	double TendStrainSensitivity = CendStrainSensitivity;

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



	(*SHVs)(0,gradIndex) = TminStrainSensitivity;
	(*SHVs)(1,gradIndex) = TunloadSlopeSensitivity;
	(*SHVs)(2,gradIndex) = TendStrainSensitivity;

	return 0;
}

int
Concrete01T::getVariable(const char *varName, Information &theInfo)
{
  if (strcmp(varName,"ec") == 0) {
    theInfo.theDouble = epsc0;
    return 0;
  } else
    return -1;
}


/* THE OLD METHODS:
double
Concrete01T::getStressSensitivity(int gradIndex, bool conditional)
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
Concrete01T::commitSensitivity(double TstrainSensitivity, int gradIndex, int numGrads)
{
	// Not treated yet. 

	return 0;
}
*/
// AddingSensitivity:END /////////////////////////////////////////////


