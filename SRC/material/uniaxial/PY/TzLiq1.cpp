/* *********************************************************************
**    Module:	TzLiq1.cpp 
**
**    Purpose:	Provide a t-z material that gets pore pressure from a
**				specified element that contains a PorousFluidSolid.
**
**    Developed by Ross W. Boulanger
**    (C) Copyright 2002, All Rights Reserved.
**
** ****************************************************************** */

// $Revision: 1.0
// $Date: 2002/2/5
// $Source: /OpenSees/SRC/material/uniaxial/TzLiq1.cpp

// Written: RWB
// Created: Feb 2002
// Revision: A
//
// Description: This file contains the class implementation for TzLiq1

#include <TzLiq1.h>
#include <NDMaterial.h>
#include <Vector.h>
#include <Channel.h>
#include <math.h>

// Controls on internal iteration between spring components
const int TZmaxIterations = 20;
const double TZtolerance = 1.0e-12;

int TzLiq1::loadStage = 0;
Vector TzLiq1::stressV3(3);

/////////////////////////////////////////////////////////////////////
//	Constructor with data

TzLiq1::TzLiq1(int tag, int classtag, int tz_type, double t_ult, double z_50,
		double dash_pot, int solid_elem1, int solid_elem2, Domain *the_Domain)
:TzSimple1(tag, classtag, tz_type, t_ult, z_50, dash_pot),
solidElem1(solid_elem1), solidElem2(solid_elem2), theDomain(the_Domain)
{
	// Initialize TzSimple variables and history variables
	//
	this->revertToStart();
    initialTangent = Tangent;
}

/////////////////////////////////////////////////////////////////////
//	Default constructor

TzLiq1::TzLiq1()
:TzSimple1(), solidElem1(0), solidElem2(0), theDomain(0)
{
}
/////////////////////////////////////////////////////////////////////
//	Default destructor
TzLiq1::~TzLiq1()
{
  // Call TzSimple1 destructor even though it does nothing.
  //
  // TzSimple1::~TzSimple1();
}

/////////////////////////////////////////////////////////////////////
int 
TzLiq1::setTrialStrain (double newz, double zRate)
{
	// Call the base class TzSimple1 to take care of the basic t-z behavior.
	//
	TzSimple1::setTrialStrain(newz, zRate);
	Tz = newz;

	// Reset mean consolidation stress if loadStage switched from 0 to 1
	//		Note: currently assuming y-axis is vertical, and that the
	//		out-of-plane normal stress equals sigma-xx.
	//
	if(lastLoadStage ==0 && loadStage ==1){

		meanConsolStress = getEffectiveStress();
		if(meanConsolStress == 0.0){
			opserr << "WARNING meanConsolStress is 0 in solid elements, ru will divide by zero";
			opserr << "TzLiq1: " << endln;
			opserr << "Adjacent solidElems: " << solidElem1 << ", " << solidElem2 << endln;
			exit(-1);
		}
	}
	lastLoadStage = loadStage;

	// Obtain the mean effective stress from the adjacent solid elements,
	//    and calculate ru for scaling of t-z base relation.
	//
	if(loadStage == 1) {
		double meanStress = getEffectiveStress();
		Tru = 1.0 - meanStress/meanConsolStress;
		if(Tru > 0.999) Tru = 0.999;
	}
	else {
		Tru = 0.0;
	}

	// Call the base class TzSimple1 to get basic t-z response,
	//
	double baseT = TzSimple1::getStress();
	double baseTangent = TzSimple1::getTangent();

	// Check: Only update Hru if not yet converged (avoiding small changes in Tru).
	//
	if(Tz !=Cz || Tt !=Ct) Hru = Tru;

	// During dilation of the soil (ru dropping), provide a stiff transition
	//   between the old and new scaled t-z relations. This avoids illogical
	//   hardening of the t-z relation (i.e., negative stiffnesses).
	//
	if (Tru < Cru){

		maxTangent = (TzSimple1::tult/TzSimple1::z50)*(1.0-Cru);

		//  If unloading, follow the old scaled t-z relation until t=0.
		//
		if(Cz>0.0 && Tz<Cz && baseT>0.0){
			Hru = Cru;
		}
		if(Cz<0.0 && Tz>Cz && baseT<0.0){
			Hru = Cru;
		}
		
		//  If above the stiff transition line (between Tru & Cru scaled surfaces)
		//
		double zref = Cz + baseT*(Cru-Tru)/maxTangent;
		if(Cz>0.0 && Tz>Cz && Tz<zref){
			Hru = 1.0 - (Ct + (Tz-Cz)*maxTangent)/baseT;
		}
		if(Cz<0.0 && Tz<Cz && Tz>zref){
			Hru = 1.0 - (Ct + (Tz-Cz)*maxTangent)/baseT;
		}
	}

	//  Now set the tangent and Tt values accordingly
	//

	Tt = baseT*(1.0-Hru);
	if(Hru==Cru || Hru==Tru){
		Tangent = (1.0-Hru)*baseTangent;
	}
	else {
		Tangent = maxTangent;
	}

	return 0;
}
/////////////////////////////////////////////////////////////////////
double 
TzLiq1::getStress(void)
{
	double dashForce = getStrainRate()*this->getDampTangent();

	// Limit the combined force to tult*(1-ru).
	//
	double tmax = (1.0-TZtolerance)*TzSimple1::tult*(1.0-Hru);
	if(fabs(Tt + dashForce) >= tmax)
		return tmax*(Tt+dashForce)/fabs(Tt+dashForce);
	else return Tt + dashForce;

}
/////////////////////////////////////////////////////////////////////
double 
TzLiq1::getTangent(void)
{
	return this->Tangent;

}
/////////////////////////////////////////////////////////////////////
double 
TzLiq1::getInitialTangent(void)
{
    return this->initialTangent;
}
/////////////////////////////////////////////////////////////////////
double 
TzLiq1::getDampTangent(void)
{
	// Call the base class TzSimple1 to get basic t-z response,
	//    and then scale by (1-ru).
	//
	double dampTangent = TzSimple1::getDampTangent();
	return dampTangent*(1.0-Hru);
}
/////////////////////////////////////////////////////////////////////
double 
TzLiq1::getStrain(void)
{
	double strain = TzSimple1::getStrain();
    return strain;
}
/////////////////////////////////////////////////////////////////////
double 
TzLiq1::getStrainRate(void)
{
    double strainrate = TzSimple1::getStrainRate();
	return strainrate;
}
/////////////////////////////////////////////////////////////////////
int 
TzLiq1::commitState(void)
{
	// Call the TzSimple1 base function to take care of details.
	//
	TzSimple1::commitState();
	Cz = Tz;
	Ct = Tt;
	Cru= Hru;
	
    return 0;
}

/////////////////////////////////////////////////////////////////////
int 
TzLiq1::revertToLastCommit(void)
{
	// reset to committed values
    //
	TzSimple1::revertToLastCommit();
	Tz = Cz;
	Tt = Ct;
	Hru= Cru;

	return 0;
}

/////////////////////////////////////////////////////////////////////
int 
TzLiq1::revertToStart(void)
{
	// Call the TzSimple1 base function to take care of most details.
	//
	TzSimple1::revertToStart();
	Tz = 0.0;
	Tt = 0.0;
	maxTangent = (TzSimple1::tult/TzSimple1::z50);

	// Excess pore pressure ratio and pointers
	//
	Tru = 0.0;
	Hru = 0.0;
	meanConsolStress = -TzSimple1::tult;
	lastLoadStage = 0;
	loadStage = 0;
	elemFlag.assign("NONE");

	// Now get all the committed variables initiated
	//
	commitState();

    return 0;
}

/////////////////////////////////////////////////////////////////////
double 
TzLiq1::getEffectiveStress(void)
{
	// Default value for meanStress
	double meanStress = meanConsolStress;
	
	// if the elemFlag has not been set yet, then set it
	//
	if(elemFlag.compare("NONE") == 0) {	//string.compare returns zero if equal

		// if theDomain pointer is nonzero, then set pointers to attached soil elements.
		//
		if(theDomain != 0)
		{	
			Element *theElement1 = theDomain->getElement(solidElem1);
			Element *theElement2 = theDomain->getElement(solidElem2);
			if (theElement1 == 0 || theElement2 == 0) {
				opserr << "WARNING solid element not found in getEffectiveStress" << endln;
				opserr << "TzLiq1: " << endln;
				opserr << "Adjacent solidElems: " << solidElem1 << ", " << solidElem2 << endln;
				exit(-1);
			}

			// Check each of the allowable element types, starting with 4NodeQuads
			theQuad1 = dynamic_cast<FourNodeQuad*>(theElement1);
			theQuad2 = dynamic_cast<FourNodeQuad*>(theElement2);
			if(theQuad1 != 0 && theQuad2 != 0) {
				elemFlag.assign("4NodeQuad");
			}

			// Check on each other type of allowable element types, only if no successful yet.
			if(elemFlag.compare("NONE") == 0) {	//string.compare returns zero if equal

				// try other elements like, 4NodeQuadUP
				elemFlag.assign("NONE");
			}
			
			// Check on acceptable soil materials
			//
			if(elemFlag.compare("4NodeQuad") == 0) {
				NDMaterial *theNDM1[4];
				NDMaterial *theNDM2[4];
				FluidSolidPorousMaterial *theFSPM1[4];
				FluidSolidPorousMaterial *theFSPM2[4];
				int dummy = 0;
				for (int i=0; i<4; i++) {
					theNDM1[i] = theQuad1->theMaterial[i];
					theNDM2[i] = theQuad2->theMaterial[i];
					theFSPM1[i] = dynamic_cast<FluidSolidPorousMaterial*>(theNDM1[i]);
					theFSPM2[i] = dynamic_cast<FluidSolidPorousMaterial*>(theNDM2[i]);
					if(theFSPM1 == 0 || theFSPM2 == 0) dummy = dummy + 1;
				}
				if(dummy == 0) elemFlag.append("-FSPM");

				// Check other acceptable soil types.
			}

			if(elemFlag.compare("NONE") == 0) {	// Never obtained a valid pointer set
				opserr << "WARNING: Adjoining solid elements did not return valid pointers";
				opserr << "TzLiq1: " << endln;
				opserr << "Adjacent solidElems: " << solidElem1 << ", " << solidElem2 << endln;
				exit(-1);
				return meanStress;
			}
	
		}
	}
	
	// Get effective stresses using appropriate pointers in elemFlag not "NONE"
	//
	if(elemFlag.compare("NONE") != 0) {

		if(elemFlag.compare("4NodeQuad-FSPM") == 0) {
			meanStress = 0.0;
			Vector *theStressVector = &stressV3;
			double excessPorePressure = 0.0;
			NDMaterial *theNDM;
			FluidSolidPorousMaterial *theFSPM;

			for(int i=0; i < 4; i++) { 
				*theStressVector = theQuad1->theMaterial[i]->getStress();
				meanStress += 2.0*(*theStressVector)[0] + (*theStressVector)[1];
				*theStressVector = theQuad2->theMaterial[i]->getStress();
				meanStress += 2.0*(*theStressVector)[0] + (*theStressVector)[1];

				theNDM = theQuad1->theMaterial[i];
				theFSPM= dynamic_cast<FluidSolidPorousMaterial*>(theNDM);
				excessPorePressure += (theFSPM->trialExcessPressure);
				theNDM = theQuad2->theMaterial[i];
				theFSPM= dynamic_cast<FluidSolidPorousMaterial*>(theNDM);
				excessPorePressure += (theFSPM->trialExcessPressure);
			}
			meanStress = meanStress/(2.0*4.0*3.0);
			excessPorePressure = excessPorePressure/(2.0*4.0);
			meanStress = meanStress - excessPorePressure;

			return meanStress;
		}

		else if(elemFlag.compare("4NodeQuadUP-FSPM") == 0) {
			
			// expect to later have UP option on quads

			return meanStress;
		}
		
		else {	// Never found a matching flag
			opserr << "WARNING: Something wrong with specification of elemFlag in getEffectiveStress";
			opserr << "TzLiq1: " << endln;
			opserr << "Adjacent solidElems: " << solidElem1 << ", " << solidElem2 << endln;
			exit(-1);
			return meanStress;
		}
	}

	return meanStress;
}

/////////////////////////////////////////////////////////////////////
int 
TzLiq1::updateParameter(int snum,Information &eleInformation)
{
	// TclUpdateMaterialStageCommand will call this routine with the
	// command:
	//
	//      updateMaterialStage - material tag -stage snum
	//
	// If snum = 0; running linear elastic for soil elements,
	//              so excess pore pressure should be zero.
	// If snum = 1; running plastic soil element behavior,
	//              so this marks the end of the "consol" gravity loading.

	if(snum !=0 && snum !=1){
		opserr << "WARNING updateMaterialStage for TzLiq1 material must be 0 or 1";
		opserr << endln;
		exit(-1);
	}
	loadStage = snum;

	return 0;
}

/////////////////////////////////////////////////////////////////////
UniaxialMaterial *
TzLiq1::getCopy(void)
{
	// Make a new instance of this class and then assign it "this" to make a copy.
	//
	TzLiq1 *clone;			// pointer to a TzLiq1 class
	clone = new TzLiq1();	// pointer gets a new instance of TzLiq1
	*clone = *this;			// the clone (dereferenced pointer) = dereferenced this.
	
	return clone;
}

/////////////////////////////////////////////////////////////////////
int 
TzLiq1::sendSelf(int cTag, Channel &theChannel)
{
	// I'm following an example by Frank here.
	//
	int res =0;

	static Vector data(16);

	TzSimple1::sendSelf(cTag, theChannel);

	data(0)  = this->getTag();
	data(1)  = Tz;
	data(2)  = Cz;
	data(3)  = Tt;
	data(4)  = Ct;
	data(5)  = Tangent;
	data(6)  = maxTangent;
	data(7)  = Tru;
	data(8)  = Cru;
	data(9)  = Hru;
	data(10) = solidElem1;
	data(11) = solidElem2;
	data(12) = meanConsolStress;
	data(13) = loadStage;
	data(14) = lastLoadStage;
	data(15) = initialTangent;
	
  res = theChannel.sendVector(this->getDbTag(), cTag, data);
  if (res < 0) 
    opserr << "TzLiq1::sendSelf() - failed to send data\n";

  return res;
}

/////////////////////////////////////////////////////////////////////
int 
TzLiq1::recvSelf(int cTag, Channel &theChannel, 
			       FEM_ObjectBroker &theBroker)
{
  int res = 0;

  static Vector data(16);
  res = theChannel.recvVector(this->getDbTag(), cTag, data);
  
  if (res < 0) {
      opserr << "TzLiq1::recvSelf() - failed to receive data\n";
      this->setTag(0);      
  }
  else {
    this->setTag((int)data(0));

	TzSimple1::recvSelf(cTag, theChannel, theBroker);

	Tz    = data(1);
	Cz    = data(2);
	Tt    = data(3);
	Ct    = data(4);
	Tangent    = data(5);
	maxTangent =data(6);
	Tru        = data(7);
	Cru        = data(8);
	Hru        = data(9);
	solidElem1        = (int)data(10);
	solidElem2        = (int)data(11);
	meanConsolStress  = data(12);
	loadStage         = (int)data(13);
	lastLoadStage     = (int)data(14);
	initialTangent    = data(15);

	// set the trial quantities
	this->revertToLastCommit();
  }

  return res;
}

/////////////////////////////////////////////////////////////////////
void 
TzLiq1::Print(OPS_Stream &s, int flag)
{
    s << "TzLiq1, tag: " << this->getTag() << endln;
    s << "  tzType: " << tzType << endln;
    s << "  tult: " << tult << endln;
    s << "  z50: " << z50 << endln;
	s << "  dashpot: " << dashpot << endln;
	s << "  solidElem1: " << solidElem1 << endln;
	s << "  solidElem2: " << solidElem2 << endln;
}

/////////////////////////////////////////////////////////////////////
