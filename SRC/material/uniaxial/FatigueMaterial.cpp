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
// $Date: 2003-08-14 20:23:50 $
// $Source: /usr/local/cvs/OpenSees/SRC/material/uniaxial/FatigueMaterial.cpp,v $
                                                      
// Written: Patxi
// Created: Aug 2003
//
// Description: This file contains the class definition for 
// FatigueMaterial.  FatigueMaterial wraps a UniaxialMaterial
// and imposes fatigue limits.

#include <stdlib.h>

#include <FatigueMaterial.h>
#include <ID.h>
#include <Channel.h>
#include <FEM_ObjectBroker.h>

#include <OPS_Globals.h>

FatigueMaterial::FatigueMaterial(int tag, UniaxialMaterial &material,
				   double dmax, double nf, double e0, double fe)
  :UniaxialMaterial(tag,MAT_TAG_Fatigue), theMaterial(0), 
   Tfailed(false), Cfailed(false)
{
  D   = 0; // Damage index
  X   = 0; // Range in consideration
  Y   = 0; // Previous Adjacent Range
  A   = 0; // Strain at first  cycle peak/valley
  B   = 0; // Strain at second cycle peak/valley
  R1F = 0; // Flag for first  cycle count
  R2F = 0; // Flag for second cycle count
  CS  = 0; // Current Slope
  PS  = 0; // Previous slope
  EP  = 0; // Previous Strain
  FF  = 0; // Failure Flag
  SF  = 0; // Start Flag - for initializing the very first strain
  PF  = 0; // Peak Flag --> Did we reach a peak/valley at current strain?

  Dmax = dmax;
  Nf  = nf;
  E0  = e0;
  FE  = fe;
  b   =  log10(E0) - log10(Nf)*FE; //Theoretical Intercept

  theMaterial = material.getCopy();

  if (theMaterial == 0) {
    opserr <<  "FatigueMaterial::FatigueMaterial -- failed to get copy of material\n";
    exit(-1);
  }
}

FatigueMaterial::FatigueMaterial()
  :UniaxialMaterial(0,MAT_TAG_Fatigue), theMaterial(0), 
   Tfailed(false), Cfailed(false)
{
  D   = 0; // Damage index
  X   = 0; // Range in consideration
  Y   = 0; // Previous Adjacent Range
  A   = 0; // Strain at first  cycle peak/valley
  B   = 0; // Strain at second cycle peak/valley
  R1F = 0; // Flag for first  cycle count
  R2F = 0; // Flag for second cycle count
  CS  = 0; // Current Slope
  PS  = 0; // Previous slope
  EP  = 0; // Previous Strain
  FF  = 0; // Failure Flag
  SF  = 0; // Start Flag - for initializing the very first strain
  PF  = 0; // Peak Flag --> Did we reach a peak/valley at current strain?

  Dmax = 0;
  Nf   = 0;
  E0   = 0;
  FE   = 0;
  b    =  0;
}

FatigueMaterial::~FatigueMaterial()
{
  if (theMaterial)
    delete theMaterial;
}

static int sign(double a) {
  if (a < 0)
    return -1;
  else
    return 1;
}

int 
FatigueMaterial::setTrialStrain(double strain, double strainRate)
{
  if (Cfailed)
    return 0;

  //Now check to see if we have reached failure although  
  // not at a peak or valley, assume
  // a 1/2 cycle to the current point
  if (Tfailed == false) {
    if (SF == 1) { //Do not check on the VERY first cycle i.e. SF=0
      double Ytemp = fabs(strain - B);
      double Dtemp = D +  0.5 / pow(10.0, ((log10(Ytemp)-b)/FE ) );
      if (Dtemp > Dmax) {
	D = Dtemp;
	opserr << "Fail!\n";
	Tfailed = true;
	return 0;
      }
    }
  }
  
  if (Tfailed == false) {
    if (SF == 0) {
      PS = strain;
      A  = strain;
      SF = 1;
    }
    CS = strain - EP;         // Determine Current Slope
    // Identify Peak or Valley (slope Change)
    if ((sign(PS) != sign(CS)) && (FF == 0)) {
      X =  fabs(strain - B);   // Range of current Cycle
      if (R1F == 0) {            // After Reaching First Peak
	Y = fabs(strain - A);
	B = strain;
	R1F = 1;
      }
      if ((R2F == 0) && (R1F == 1)) { // Mark the first cycle  
	A = B;
	B = strain;
	R2F = 1;
	D = 0.5 / pow(10.0, ((log10(Y)-b)/FE ) );
	PF = 1; // We reached a peak
	Y = X;
      } else if (X < Y) {
	// Do Nothing; keep counting, not a peak or 
	// valley by definition
	;
      }
      else {
	A = B;
	B = strain;
	D = D + 0.5 / pow(10.0, ((log10(Y)-b)/FE ));
	PF = 1; // We reached a peak
	Y = X;
      }
    }
    PS = CS;       // Previous Slope
    EP = strain;   // Keep track of previous strain
    
    if (D >= Dmax) {
      Tfailed = true;
      opserr << "Fail!\n";
      return  0;
    }
  }

  if(Tfailed == false) {
    Tfailed = false;
    return theMaterial->setTrialStrain(strain, strainRate);
  } else {
    return 0;
  }
  
}

double 
FatigueMaterial::getStress(void)
{
  if (Tfailed)
    return 0.0;
  else
    return theMaterial->getStress();
}

double 
FatigueMaterial::getTangent(void)
{
  if (Tfailed)
    //return 0.0;
    return 1.0e-8*theMaterial->getInitialTangent();
  else
    return theMaterial->getTangent();
}

double 
FatigueMaterial::getDampTangent(void)
{
  if (Tfailed)
    return 0.0;
  else
    return theMaterial->getDampTangent();
}



double 
FatigueMaterial::getStrain(void)
{
  return theMaterial->getStrain();
}

double 
FatigueMaterial::getStrainRate(void)
{
  return theMaterial->getStrainRate();
}

int 
FatigueMaterial::commitState(void)
{	
  Cfailed = Tfailed;

  // Check if failed at current step
  if (Tfailed) 
    return 0;
  else 
    return theMaterial->commitState();
}

int 
FatigueMaterial::revertToLastCommit(void)
{
  // Check if failed at last step
  if (Cfailed)
    return 0;
  else
    return theMaterial->revertToLastCommit();
}

int 
FatigueMaterial::revertToStart(void)
{
  Cfailed = false;
  Tfailed = false;
  
  return theMaterial->revertToStart();
}

UniaxialMaterial *
FatigueMaterial::getCopy(void)
{
  FatigueMaterial *theCopy = 
    new FatigueMaterial(this->getTag(), *theMaterial, Dmax, Nf, E0, FE);
  
  theCopy->Cfailed = Cfailed;
  theCopy->Tfailed = Tfailed;
  
  return theCopy;
}

int 
FatigueMaterial::sendSelf(int cTag, Channel &theChannel)
{
  return 0;
}

int 
FatigueMaterial::recvSelf(int cTag, Channel &theChannel, 
			 FEM_ObjectBroker &theBroker)
{
  return 0;
}

void 
FatigueMaterial::Print(OPS_Stream &s, int flag)
{
  s << "FatigueMaterial tag: " << this->getTag() << endln;
  s << "\tMaterial: " << theMaterial->getTag() << endln;
  s << "\tD: " << D << " Dmax: " << Dmax << endln;
  s << "Nf: " << Nf <<  " E0: " << E0 << " FE: " << FE << " b: " << b << endln;
}
