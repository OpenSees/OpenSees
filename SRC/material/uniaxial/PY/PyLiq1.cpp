/* *********************************************************************
**    Module:	PyLiq1.cpp 
**
**    Purpose:	Provide a p-y material that gets pore pressure from a
**				specified element that contains a PorousFluidSolid.
**
**    Developed by Ross W. Boulanger
**    (C) Copyright 2002, All Rights Reserved.
**
** ****************************************************************** */

// $Revision: 1.0
// $Date: 2002/5/15
// $Source: /OpenSees/SRC/material/uniaxial/PyLiq1.cpp

// Written: RWB
// Created: May 2002
// Revision: A
//
// Description: This file contains the class implementation for PyLiq1

#include <PyLiq1.h>
#include <NDMaterial.h>
#include <Vector.h>
#include <Channel.h>
#include <math.h>

// Controls on internal iteration between spring components
const int PYmaxIterations = 20;
const double PYtolerance = 1.0e-12;

int PyLiq1::loadStage = 0;
Vector PyLiq1::stressV3(3);

/////////////////////////////////////////////////////////////////////
//	Constructor with data

PyLiq1::PyLiq1(int tag, int classtag, int soil, double p_ult, double y_50,
		double dragratio, double dash_pot, double p_res, 
		int solid_elem1, int solid_elem2, Domain *the_Domain)
:PySimple1(tag, classtag, soil, p_ult, y_50, dragratio, dash_pot),
pRes(p_res), solidElem1(solid_elem1), solidElem2(solid_elem2), theDomain(the_Domain)
{
	// Initialize PySimple variables and history variables
	//
	this->revertToStart();
    initialTangent = Tangent;
}

/////////////////////////////////////////////////////////////////////
//	Default constructor

PyLiq1::PyLiq1()
:PySimple1(), pRes(0.0), solidElem1(0), solidElem2(0), theDomain(0)
{
}
/////////////////////////////////////////////////////////////////////
//	Default destructor
PyLiq1::~PyLiq1()
{
  // Call PySimple1 destructor even though it does nothing.
  //
  // PySimple1::~PySimple1();
}

/////////////////////////////////////////////////////////////////////
int 
PyLiq1::setTrialStrain (double newy, double yRate)
{
	// Call the base class PySimple1 to take care of the basic p-y behavior.
	//
	Ty = newy;
	PySimple1::setTrialStrain(Ty, yRate);

	// Reset mean consolidation stress if loadStage switched from 0 to 1
	//		Note: currently assuming y-axis is vertical, and that the
	//		out-of-plane normal stress equals sigma-xx.
	//
	if(lastLoadStage ==0 && loadStage ==1){

		meanConsolStress = getEffectiveStress();
		if(meanConsolStress == 0.0){
			opserr << "WARNING meanConsolStress is 0 in solid elements, ru will divide by zero";
			opserr << "PyLiq1: " << endln;
			opserr << "Adjacent solidElems: " << solidElem1 << ", " << solidElem2 << endln;
			exit(-1);
		}
	}
	lastLoadStage = loadStage;

	// Obtain the mean effective stress from the adjacent solid elements,
	//    and calculate ru for scaling of p-y base relation.
	//
	if(loadStage == 1) {
		double meanStress = getEffectiveStress();
		Tru = 1.0 - meanStress/meanConsolStress;
		if(Tru > 1.0-pRes/PySimple1::pult) Tru = 1.0-pRes/PySimple1::pult;
	}
	else {
		Tru = 0.0;
	}


	// Call the base class PySimple1 to get basic p-y response,
	//
	double baseP = PySimple1::getStress();
	double baseTangent = PySimple1::getTangent();

	// Check: Only update Hru if not yet converged (avoiding small changes in Tru).
	//
	Hru = Tru;
	if(Ty==Cy && Tp==Cp) { Hru = Cru;}

	// During dilation of the soil (ru dropping), provide a stiff transition
	//   between the old and new scaled p-y relations. This avoids illogical
	//   hardening of the p-y relation (i.e., negative stiffnesses).
	//
	if (Hru < Cru){

		maxTangent = (PySimple1::pult/PySimple1::y50)*(1.0-Cru);

		//  If unloading, follow the old scaled p-y relation until p=0.
		//
		if(Cy>0.0 && Ty<Cy && baseP>0.0){
			Hru = Cru;
		}
		if(Cy<0.0 && Ty>Cy && baseP<0.0){
			Hru = Cru;
		}
		
		//  If above the stiff transition line (between Tru & Cru scaled surfaces)
		//
		double yref = Cy + baseP*(Cru-Hru)/maxTangent;
		if(Cy>0.0 && Ty>Cy && Ty<yref){
			Hru = 1.0 - (Cp + (Ty-Cy)*maxTangent)/baseP;
		}
		if(Cy<0.0 && Ty<Cy && Ty>yref){
			Hru = 1.0 - (Cp + (Ty-Cy)*maxTangent)/baseP;
		}
		if(Hru > Cru) Hru = Cru;
		if(Hru < Tru) Hru = Tru;

	}

	//  Now set the tangent and Tp values accordingly
	//

	Tp = baseP*(1.0-Hru);
//	Tangent = (1.0-Hru)*baseTangent;
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
PyLiq1::getStress(void)
{
	double dashForce = getStrainRate()*this->getDampTangent();

	// Limit the combined force to pult*(1-ru).
	//
	double pmax = (1.0-PYtolerance)*PySimple1::pult*(1.0-Hru);
//	double pmax = (1.0-PYtolerance)*PySimple1::pult;
	if(fabs(Tp + dashForce) >= pmax)
		return pmax*(Tp+dashForce)/fabs(Tp+dashForce);
	else return Tp + dashForce;

}
/////////////////////////////////////////////////////////////////////
double 
PyLiq1::getTangent(void)
{
	return this->Tangent;
}
/////////////////////////////////////////////////////////////////////
double 
PyLiq1::getInitialTangent(void)
{
    return this->initialTangent;
}
/////////////////////////////////////////////////////////////////////
double 
PyLiq1::getDampTangent(void)
{
	// Call the base class PySimple1 to get basic p-y response,
	//    and then scale by (1-ru).
	//
	double dampTangent = PySimple1::getDampTangent();
	return dampTangent*(1.0-Hru);
}
/////////////////////////////////////////////////////////////////////
double 
PyLiq1::getStrain(void)
{
	double strain = PySimple1::getStrain();
    return strain;
}
/////////////////////////////////////////////////////////////////////
double 
PyLiq1::getStrainRate(void)
{
    double strainrate = PySimple1::getStrainRate();
	return strainrate;
}
/////////////////////////////////////////////////////////////////////
int 
PyLiq1::commitState(void)
{
	// Call the PySimple1 base function to take care of details.
	//
	PySimple1::commitState();
	Cy = Ty;
	Cp = Tp;
	Cru= Hru;
	
    return 0;
}

/////////////////////////////////////////////////////////////////////
int 
PyLiq1::revertToLastCommit(void)
{
	// reset to committed values
    //
	PySimple1::revertToLastCommit();
	Ty = Cy;
	Tp = Cp;
	Hru= Cru;

	return 0;
}

/////////////////////////////////////////////////////////////////////
int 
PyLiq1::revertToStart(void)
{
	// Call the PySimple1 base function to take care of most details.
	//
	PySimple1::revertToStart();
	Ty = 0.0;
	Tp = 0.0;
	maxTangent = (PySimple1::pult/PySimple1::y50);
	
	// Excess pore pressure ratio and pointers
	//
	Tru = 0.0;
	Hru = 0.0;
	meanConsolStress = -PySimple1::pult;
	lastLoadStage = 0;
	loadStage = 0;
	if(pRes <= 0.0) pRes = 0.01*PySimple1::pult;
	if(pRes > PySimple1::pult) pRes = PySimple1::pult;
	elemFlag.assign("NONE");

	// Now get all the committed variables initiated
	//
	commitState();

    return 0;
}

/////////////////////////////////////////////////////////////////////
double 
PyLiq1::getEffectiveStress(void)
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
				opserr << "PyLiq1: " << endln;
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
				opserr << "PyLiq1: " << endln;
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
			opserr << "PyLiq1: " << endln;
			opserr << "Adjacent solidElems: " << solidElem1 << ", " << solidElem2 << endln;
			exit(-1);
			return meanStress;
		}
	}

	return meanStress;
}

/////////////////////////////////////////////////////////////////////
int 
PyLiq1::updateParameter(int snum,Information &eleInformation)
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
		opserr << "WARNING updateMaterialStage for PyLiq1 material must be 0 or 1";
		opserr << endln;
		exit(-1);
	}
	loadStage = snum;

	return 0;
}

/////////////////////////////////////////////////////////////////////
UniaxialMaterial *
PyLiq1::getCopy(void)
{
	// Make a new instance of this class and then assign it "this" to make a copy.
	//
	PyLiq1 *clone;			// pointer to a PyLiq1 class
	clone = new PyLiq1();	// pointer gets a new instance of PyLiq1
	*clone = *this;			// the clone (dereferenced pointer) = dereferenced this.
	
	return clone;
}

/////////////////////////////////////////////////////////////////////
int 
PyLiq1::sendSelf(int cTag, Channel &theChannel)
{
	// I'm following an example by Frank here.
	//
	int res =0;

	static Vector data(16);

	PySimple1::sendSelf(cTag, theChannel);

	data(0)  = this->getTag();
	data(1)  = Ty;
	data(2)  = Cy;
	data(3)  = Tp;
	data(4)  = Cp;
	data(5)  = Tangent;
	data(6)  = maxTangent;
	data(7)  = Tru;
	data(8)  = Cru;
	data(9)  = Hru;
	data(10) = (int)solidElem1;
	data(11) = (int)solidElem2;
	data(12) = meanConsolStress;
	data(13) = (int)loadStage;
	data(14) = (int)lastLoadStage;
	data(15) = initialTangent;
	
  res = theChannel.sendVector(this->getDbTag(), cTag, data);
  if (res < 0) 
    opserr << "PyLiq1::sendSelf() - failed to send data\n";

  return res;
}

/////////////////////////////////////////////////////////////////////
int 
PyLiq1::recvSelf(int cTag, Channel &theChannel, 
			       FEM_ObjectBroker &theBroker)
{
  int res = 0;

  static Vector data(16);
  res = theChannel.recvVector(this->getDbTag(), cTag, data);
  
  if (res < 0) {
      opserr << "PyLiq1::recvSelf() - failed to receive data\n";
      this->setTag(0);      
  }
  else {
    this->setTag((int)data(0));

	PySimple1::recvSelf(cTag, theChannel, theBroker);

	Ty    = data(1);
	Cy    = data(2);
	Tp    = data(3);
	Cp    = data(4);
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
PyLiq1::Print(OPS_Stream &s, int flag)
{
    s << "PyLiq1, tag: " << this->getTag() << endln;
    s << "  soilType: " << soilType << endln;
    s << "  pult: " << pult << endln;
    s << "  y50: " << y50 << endln;
    s << "  drag: " << drag << endln;
	s << "  pResidual: " << pRes << endln;
	s << "  dashpot: " << dashpot << endln;
	s << "  solidElem1: " << solidElem1 << endln;
	s << "  solidElem2: " << solidElem2 << endln;
}

/////////////////////////////////////////////////////////////////////
