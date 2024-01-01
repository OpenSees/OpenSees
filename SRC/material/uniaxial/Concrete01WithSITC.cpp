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
                                                                        
// $Revision: 1.2 $
// $Date: 2010-04-06 20:18:49 $
// $Source: /usr/local/cvs/OpenSees/SRC/material/uniaxial/Concrete01WithSITC.cpp,v $
         
// Modified by: Won Lee
// Created: 10/3/05
// Modified from Concrete01.C (see details below)
// Description: Concrete01 model modified to include SITC effect (ref. Prof. 
//	John Stanton of Univ. of Washington).  Use modified rules from his paper to include this 
//	effect (J.F. Stanton and H.D. McNiven, "The Development of a Mathematical
//	Model to Predict the Flexural Response of Reinforced Concrete Beams to Cyclic
//	Loads, Using System Identification", EERC Report Number 79/02, January 1979.
                            
// FILE BASED ON:                                   
// File: Concrete01.cpp
//
// Written: MHS 
// Created: 06/99
// Revision: A
//
// Description: This file contains the class implementation for 
// Concrete01. 
//
// What: "@(#) Concrete01.C, revA"


#include <Concrete01WithSITC.h>
#include <Vector.h>
#include <Matrix.h>
#include <Channel.h>
#include <Information.h>
#include <math.h>
#include <float.h>
#include <elementAPI.h>

//int count = 0;

void* OPS_Concrete01WithSITC()
{
    int numdata = OPS_GetNumRemainingInputArgs();
    if (numdata < 5) {
	opserr << "WARNING insufficient arguments\n";
	opserr << "Want: uniaxialMaterial Concrete01WithSITC tag? ";
	opserr << "fpc? epsc0? fpcu? epscu? <endStrainSITC?>\n";
	return 0;
    }

    int tag;
    numdata = 1;
    if (OPS_GetIntInput(&numdata,&tag) < 0) {
	opserr << "WARNING invalid tag\n";
	return 0;
    }

    double data[4];
    numdata = 4;
    if (OPS_GetDoubleInput(&numdata,data)) {
	opserr << "WARNING invalid double data\n";
	return 0;
    }
    UniaxialMaterial* mat = 0;
    
    numdata = OPS_GetNumRemainingInputArgs();
    if (numdata > 0) {
	double endStrainSITC;
	numdata = 1;
	if (OPS_GetDoubleInput(&numdata,&endStrainSITC) < 0) {
	    opserr << "WARNING invalid double data\n";
	    return 0;
	}
	mat = new Concrete01WithSITC(tag,data[0],data[1],data[2],data[3],endStrainSITC);
    } else {
	mat = new Concrete01WithSITC(tag,data[0],data[1],data[2],data[3]);
    }

    if (mat == 0) {
	opserr << "WARNING: failed to create Concrete01WithSITC material\n";
	return 0;
    }

    return mat;
}

Concrete01WithSITC::Concrete01WithSITC
(int tag, double FPC, double EPSC0, double FPCU, double EPSCU, double endStrainSITC)
  :UniaxialMaterial(tag, MAT_TAG_Concrete01WithSITC),
   fpc(FPC), epsc0(EPSC0), fpcu(FPCU), epscu(EPSCU), 
   CminStrain(0.0), CendStrain(0.0),
   Cstrain(0.0), Cstress(0.0), CmaxStrain(0.0),
   CslopeSITC(0.0), CendStrainSITC(endStrainSITC), Cindex(0), CsmallStrainIndex(0)
{
  //count++;
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
  
}

Concrete01WithSITC::Concrete01WithSITC():UniaxialMaterial(0, MAT_TAG_Concrete01WithSITC),
 fpc(0.0), epsc0(0.0), fpcu(0.0), epscu(0.0),
 CminStrain(0.0), CunloadSlope(0.0), CendStrain(0.0),
 Cstrain(0.0), Cstress(0.0), CmaxStrain(0.0),
 CslopeSITC(0.0), CendStrainSITC(0.0), Cindex(0), CsmallStrainIndex(0)

{
  // Set trial values
  this->revertToLastCommit();
  
}

Concrete01WithSITC::~Concrete01WithSITC ()
{
   // Does nothing
}

int Concrete01WithSITC::setTrialStrain (double strain, double strainRate)
{
  // Set trial strain
  Tstrain = strain;

  TminStrain = CminStrain;
  //TmaxStrain = CmaxStrain;
  Tindex = Cindex;
  
  TslopeSITC = CslopeSITC;
  TunloadSlope = CunloadSlope;

  // Determine change in strain from last converged state
  double dStrain = Tstrain - Cstrain;
  
  if (fabs(dStrain) < DBL_EPSILON) {
    return 0;
  }

  if (Tstrain < 0.0) {  // compression
    if (Tstrain <= CminStrain) { // further on envelope curve
      TminStrain = Tstrain;
      envelope();
      unload();
      Tindex = 1;
    }
    else if (Tstrain >= CendStrainSITC) {
      Tstress = 0.0;
      Ttangent = 0.0;
      Tindex = 5;
    }

    else { // anywhere in compression greater than minimum strain
      if (dStrain <= 0.0) { //loading in compression 
	if (Cindex == 2 || Cindex == 1) {   
	  Tstress = Cstress + TunloadSlope*dStrain;
	  Ttangent = TunloadSlope;
	  Tindex = 2; 
	}
	else if (Cindex == 3) { 
	  Tstress = Cstress + TslopeSITC*dStrain;
	  Ttangent = TslopeSITC;
	  Tindex = 3; 
	}
	else if (Cindex == 5) {
	  if (Tstrain <= CendStrainSITC && Cstrain >= CendStrainSITC) {
	    Ttangent = TslopeSITC;
	    Tstress = TslopeSITC*(Tstrain-CendStrainSITC);
	    Tindex = 3;
	  }
	  else if (Tstrain <= TendStrain) {
	    Ttangent = TunloadSlope;
	    Tstress = TunloadSlope*(Tstrain-TendStrain);
	    Tindex = 2;
	  }
	  else {
	    Ttangent = 0.0;
	    Tstress = 0.0;
	    Tindex = 5;
	  }
	}
	else {
	  
	  opserr << "something in compression is wrong!! Cstrain " << endln;
	}
      }
      else { // unloading in compression
	if (Cindex == 1 || Cindex == 2) { //unloading on regular branch
	  if (Tstrain >= TendStrain) {
	    Tstress = 0.0;
	    Ttangent = 0.0;
	    Tindex = 5;
	  }
	  else { 
	    Tstress = Cstress + TunloadSlope*dStrain;
	    Ttangent = TunloadSlope;
	    Tindex = 2;
	  }
	}
	else if (Cindex == 3) {
	  Tstress = Cstress + TslopeSITC*dStrain;
	  Ttangent = TslopeSITC;
	  Tindex = 3;
	  if (Tstress > 0.0) {
	    opserr << "THERE IS A PROBLEM IN UNLOADING IN COMPRESSION!!!" << endln;
	  }
	}
	else if (Cindex == 5) { // index must be 5
	  Tstress = 0.0;
	  Ttangent = 0.0;
	  Tindex = 5; // **************
	}
	else {
	  opserr << "Something is wrong in tension!!!! Cindex is " << endln;
	}
      }
    }
  }
  else { // TENSION

    Ttangent = 0.0;
    Tstress = 0.0;
    Tindex = 5;
    
  }
  
  return 0;
  
}



int 
Concrete01WithSITC::setTrial (double strain, double &stress, double &tangent, double strainRate)
{
  TminStrain = CminStrain;
  //TmaxStrain = CmaxStrain;
  Tindex = Cindex;
  
  
  // Set trial strain
  Tstrain = strain;

  TslopeSITC = CslopeSITC;
  TunloadSlope = CunloadSlope;


  // Determine change in strain from last converged state
  double dStrain = Tstrain - Cstrain;
  
  if (fabs(dStrain) < DBL_EPSILON) { 
    tangent = Ttangent;
    stress = Tstress;
    return 0;
  }
  
  if (Tstrain >=  CendStrainSITC ) {
    Ttangent = 0.0;
    Tstress = 0.0;
    Tindex = 5;
    tangent = Ttangent;
    stress = Tstress;
    return 0;
  }
  
  if (Tstrain < 0.0) {  // compression
    if (Tstrain <= CminStrain) { // further on envelope curve
      TminStrain = Tstrain;
      envelope();
      unload();
      Tindex = 1;
    }
    else if (Tstrain >= CendStrainSITC) {
      Tstress = 0.0;
      Ttangent = 0.0;
      Tindex = 5;
    }
    else { // anywhere in compression greater than minimum strain
      if (dStrain <= 0.0) { //loading in compression 
	if (Cindex == 2 || Cindex == 1) {  
	  Tstress = Cstress + TunloadSlope*dStrain;
	  Ttangent = TunloadSlope;
	  Tindex = 2;  // 
	}
	else if (Cindex == 3) { 
	  Tstress = Cstress + TslopeSITC*dStrain;
	  Ttangent = TslopeSITC;
	  Tindex = 3; 
	}
	else if (Cindex == 5) {
	  if (Tstrain <= CendStrainSITC && Cstrain >= CendStrainSITC) {
	    Ttangent = TslopeSITC;
	    Tstress = TslopeSITC*(Tstrain-CendStrainSITC);
	    Tindex = 3;
	  }
	  else if (Tstrain <= TendStrain) {
	    Ttangent = TunloadSlope;
	    Tstress = TunloadSlope*(Tstrain-TendStrain);
	    Tindex = 2;
	  }
	  else {
	    Ttangent = 0.0;
	    Tstress = 0.0;
	    Tindex = 5;
	  }
	}
	else {
	  
	  opserr << "something in compression is wrong!! Cstrain " << endln;
	}
      }
      else { // unloading in compression
	if (Cindex == 1 || Cindex == 2) { //unloading on regular branch
	  if (Tstrain >= TendStrain) {
	    Tstress = 0.0;
	    Ttangent = 0.0;
	    Tindex = 5;
	  }
	  else { 
	    Tstress = Cstress + TunloadSlope*dStrain;
	    Ttangent = TunloadSlope;
	    Tindex = 2;
	  }
	}
	else if (Cindex == 3) {
	  Tstress = Cstress + TslopeSITC*dStrain;
	  Ttangent = TslopeSITC;
	  Tindex = 3;
	  if (Tstress > 0.0) {
	    opserr << "THERE IS A PROBLEM IN UNLOADING IN COMPRESSION!!!" << endln;
	  }
	}
	else if (Cindex == 5) { // index must be 5
	  Tstress = 0.0;
	  Ttangent = 0.0;
	  Tindex = 5; // **************
	}
	else {
	  opserr << "Something is wrong in tension!!!! Cindex is " << endln;
	}
      }
    }
  }
  else { // TENSION
    
    
    if (dStrain > 0.0) { // going toward tension
      if (Cindex == 1 || Cindex == 2 || Cindex == 5 || Cindex == 0) {  
	Tstress = 0.0;
	Ttangent = 0.0;
	Tindex = 5;
      }
      else if (Cindex == 3) {
	Tstress = Cstress + TslopeSITC*dStrain;
	Ttangent = TslopeSITC;
	Tindex = 3;
      }
      else {
	opserr << " something is wrong in tension loading!!! Cindex " << endln;
      }
    } 
    else {  // going toward compression
      if (Cindex == 5) {
	if (Tstrain <= CendStrainSITC && Cstrain >= CendStrainSITC) {
	  Ttangent = TslopeSITC;
	  Tstress = TslopeSITC*(Tstrain-CendStrainSITC);
	  Tindex = 3;
	}
	else {
	  Ttangent = 0.0;
	  Tstress = 0.0;
	  Tindex = 5;
	}
      }
      else if (Cindex == 3) {
	Ttangent = TslopeSITC;
	Tstress = Cstress + TslopeSITC*dStrain;
	Tindex = 3;
      }
      else {
	opserr << "something is wrong in tension going to compression " << endln;
      }
    }
    
    
  }
  stress = Tstress;
  tangent = Ttangent;
  return 0;

}

void Concrete01WithSITC::reload ()
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

void Concrete01WithSITC::envelope ()
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

void Concrete01WithSITC::getSITCslope ()
{
  double tempStrain = Tstrain;
  double tempStress = Tstress;
  Tstrain = CminStrain;
  envelope();
  TslopeSITC = Tstress/(CminStrain-CendStrainSITC);
  Tstrain = tempStrain;
  Tstress = tempStress;
}

void Concrete01WithSITC::unload ()
{
  double tempStrain = TminStrain;
  
  if (tempStrain < epscu)
    tempStrain = epscu;
  
  double eta = tempStrain/epsc0;
  
  double ratio = 0.707*(eta-2.0) + 0.834;
  
  if (eta < 2.0)
    ratio = 0.145*eta*eta + 0.13*eta;
  
  TendStrain = ratio*epsc0;
  TslopeSITC = Tstress/(TminStrain - CendStrainSITC);
  
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

double Concrete01WithSITC::getStress ()
{
   return Tstress;
}

double Concrete01WithSITC::getStrain ()
{
   return Tstrain;
}

double Concrete01WithSITC::getTangent ()
{
   return Ttangent;
}

void Concrete01WithSITC::determineTrialState (double dStrain)

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

int Concrete01WithSITC::commitState ()
{

   // History variables
   CminStrain = TminStrain;
   CunloadSlope = TunloadSlope;
   CendStrain = TendStrain;
   CmaxStrain = TmaxStrain;
   CslopeSITC = TslopeSITC;
   Cindex = Tindex;
   CsmallStrainIndex = TsmallStrainIndex;

   // State variables
   Cstrain = Tstrain;
   Cstress = Tstress;
   Ctangent = Ttangent;

   return 0;
}

int Concrete01WithSITC::revertToLastCommit ()
{
   // Reset trial history variables to last committed state
   TminStrain = CminStrain;
   TendStrain = CendStrain;
   TunloadSlope = CunloadSlope;
   TmaxStrain = CmaxStrain;
   TslopeSITC = CslopeSITC;
   Tindex = Cindex;
   TsmallStrainIndex = CsmallStrainIndex;


   // Recompute trial stress and tangent
   Tstrain = Cstrain;
   Tstress = Cstress;
   Ttangent = Ctangent;

   return 0;
}

int Concrete01WithSITC::revertToStart ()
{
	double Ec0 = 2.0*fpc/epsc0;

   // History variables
   CminStrain = 0.0;
   CunloadSlope = Ec0;
   CendStrain = 0.0;
   CmaxStrain = 0.0;
   CslopeSITC = 0.0;
   Cindex = 0;
   CsmallStrainIndex = 0;


   // State variables
   Cstrain = 0.0;
   Cstress = 0.0;
   Ctangent = Ec0;

   // Reset trial variables and state
   this->revertToLastCommit();

   return 0;
}

UniaxialMaterial* Concrete01WithSITC::getCopy ()
{
   Concrete01WithSITC* theCopy = new Concrete01WithSITC(this->getTag(),
							fpc, epsc0, fpcu, epscu, CendStrainSITC);

   // Converged history variables
   theCopy->CminStrain = CminStrain;
   theCopy->CunloadSlope = CunloadSlope;
   theCopy->CendStrain = CendStrain;
   theCopy->CmaxStrain = CmaxStrain;
   theCopy->CslopeSITC = CslopeSITC;
   theCopy->Cindex = Cindex;
   theCopy->CsmallStrainIndex = CsmallStrainIndex;

   // Converged state variables
   theCopy->Cstrain = Cstrain;
   theCopy->Cstress = Cstress;
   theCopy->Ctangent = Ctangent;

   return theCopy;
}

int Concrete01WithSITC::sendSelf (int commitTag, Channel& theChannel)
{
   int res = 0;
   static Vector data(16);
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

   // variables added by WL
   data(11) = CmaxStrain;
   data(12) = CslopeSITC;
   data(13) = CendStrainSITC;
   data(14) = Cindex;
   data(15) = CsmallStrainIndex;

   // Data is only sent after convergence, so no trial variables
   // need to be sent through data vector

   res = theChannel.sendVector(this->getDbTag(), commitTag, data);
   if (res < 0) 
      opserr << "Concrete01WithSITC::sendSelf() - failed to send data\n";

   return res;
}

int Concrete01WithSITC::recvSelf (int commitTag, Channel& theChannel,
                                 FEM_ObjectBroker& theBroker)
{
   int res = 0;
   static Vector data(16);
   res = theChannel.recvVector(this->getDbTag(), commitTag, data);

   if (res < 0) {
      opserr << "Concrete01WithSITC::recvSelf() - failed to receive data\n";
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

      // variables added by WL
      CmaxStrain = data(11);
      CslopeSITC = data(12);
      CendStrainSITC = data(13);
      Cindex = data(14);
      CsmallStrainIndex = data(15);


      // Set trial state variables
      Tstrain = Cstrain;
      Tstress = Cstress;
      Ttangent = Ctangent;
   }

   return res;
}

void Concrete01WithSITC::Print (OPS_Stream& s, int flag)
{
	if (flag == OPS_PRINT_PRINTMODEL_MATERIAL) {
		s << "Concrete01WithSITC, tag: " << this->getTag() << endln;
		s << "  fpc: " << fpc << endln;
		s << "  epsc0: " << epsc0 << endln;
		s << "  fpcu: " << fpcu << endln;
		s << "  epscu: " << epscu << endln;
	}

	if (flag == OPS_PRINT_PRINTMODEL_JSON) {
		s << "\t\t\t{";
		s << "\"type\": \"Concrete01WithSITC\", ";
		s << "\"name\": \"" << this->getTag() << "\", ";
		s << "\"Ec\": " << 2.0*fpc/epsc0 << ", ";
		s << "\"fc\": " << fpc << ", ";
		s << "\"epsc\": " << epsc0 << ", ";
		s << "\"fcu\": " << fpcu << ", ";
		s << "\"epscu\": " << epscu << "}";
	}
}




