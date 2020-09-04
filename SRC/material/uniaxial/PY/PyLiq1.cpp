/* *********************************************************************
**    Module:	PyLiq1.cpp 
**
**    Purpose:	Provide a p-y material that gets pore pressure from a
**				specified element that contains a PorousFluidSolid.
**
** Developed by Ross W. Boulanger
**
** Copyright @ 2002 The Regents of the University of California (The Regents). All Rights Reserved.
**
** The Regents grants permission, without fee and without a written license agreement, for (a) use, 
** reproduction, modification, and distribution of this software and its documentation by educational, 
** research, and non-profit entities for noncommercial purposes only; and (b) use, reproduction and 
** modification of this software by other entities for internal purposes only. The above copyright 
** notice, this paragraph and the following three paragraphs must appear in all copies and modifications 
** of the software and/or documentation.
**
** Permission to incorporate this software into products for commercial distribution may be obtained 
** by contacting the University of California 
** Office of Technology Licensing 
** 2150 Shattuck Avenue #510, 
** Berkeley, CA 94720-1620, 
** (510) 643-7201.
**
** This software program and documentation are copyrighted by The Regents of the University of California. 
** The Regents does not warrant that the operation of the program will be uninterrupted or error-free. The 
** end-user understands that the program was developed for research purposes and is advised not to rely 
** exclusively on the program for any reason.
**
** IN NO EVENT SHALL REGENTS BE LIABLE TO ANY PARTY FOR DIRECT, INDIRECT, SPECIAL, INCIDENTAL, OR 
** CONSEQUENTIAL DAMAGES, INCLUDING LOST PROFITS, ARISING OUT OF THE USE OF THIS SOFTWARE AND ITS 
** DOCUMENTATION, EVEN IF REGENTS HAS BEEN ADVISED OF THE POSSIBILITY OF SUCH DAMAGE. REGENTS GRANTS 
** NO EXPRESS OR IMPLIED LICENSE IN ANY PATENT RIGHTS OF REGENTS BUT HAS IMPLEMENTED AN INDIVIDUAL 
** CONTRIBUTOR LICENSE AGREEMENT FOR THE OPENSEES PROJECT AT THE UNIVERSITY OF CALIFORNIA, BERKELEY 
** TO BENEFIT THE END USER.
**
** REGENTS SPECIFICALLY DISCLAIMS ANY WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES 
** OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE. THE SOFTWARE AND ACCOMPANYING DOCUMENTATION,
** IF ANY, PROVIDED HEREUNDER IS PROVIDED "AS IS". REGENTS HAS NO OBLIGATION TO PROVIDE MAINTENANCE, 
** SUPPORT, UPDATES, ENHANCEMENTS, OR MODIFICATIONS.
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
#include <Parameter.h>
#include <stdlib.h>
#include <string.h>
#include <InitialStateAnalysisWrapper.h>
#include <SSPquadUP.h>
#include <SSPquad.h>
#include <elementAPI.h>
#include <TimeSeries.h>

#include <FourNodeQuad.h>
#include <FourNodeQuadUP.h>
#include <Nine_Four_Node_QuadUP.h>
#include <FluidSolidPorousMaterial.h>
#include <PressureDependMultiYield.h>

// Controls on internal iteration between spring components
const int PYmaxIterations = 20;
const double PYtolerance = 1.0e-12;

int PyLiq1::loadStage = 0;
int PyConstructorType = 1;
Vector PyLiq1::stressV3(3);

void* OPS_PyLiq1()
{
    UniaxialMaterial* theMat = 0;
    
    int numdata = OPS_GetNumRemainingInputArgs();
    if (numdata < 9) {
	opserr << "WARNING insufficient arguments\n";
	opserr << "Want: uniaxialMaterial PyLiq1 tag? soilType? pult? y50? drag? dashpot? pRes? solidElem1? solidElem2?\n";
	opserr << "or: uniaxialMaterial PyLiq1 tag? soilType? pult? y50? drag? dashpot? -timeSeries seriesTag?\n";
	return 0;
    }

    int idata[2];
    numdata = 2;
    if (OPS_GetIntInput(&numdata, idata) < 0) {
	opserr << "WARNING invalid int inputs\n";
	return 0;
    }

    double ddata[5];
    numdata = 5;
    if (OPS_GetDoubleInput(&numdata, ddata) < 0) {
	opserr << "WARNING invalid double inputs\n";
	return 0;
    }

    const char* arg = OPS_GetString();
    Domain* theDomain = OPS_GetDomain();
    if (theDomain == 0) return 0;
    
    if (strcmp(arg, "-timeSeries") == 0) {
	int tsTag;
	numdata = 1;
	if (OPS_GetIntInput(&numdata, &tsTag) < 0) {
	    opserr << "WARNING invalid time series tag\n";
	    return 0;
	}
	TimeSeries* theSeries = OPS_getTimeSeries(tsTag);
	theMat = new PyLiq1(idata[0], MAT_TAG_PyLiq1,idata[1], ddata[0], ddata[1], ddata[2],
			    ddata[3], ddata[4],theDomain,theSeries);
	
				 
    } else {

	// back one arg
	OPS_ResetCurrentInputArg(-1);
	
	int eleTags[2];
	numdata = 2;
	if (OPS_GetIntInput(&numdata, eleTags) < 0) {
	    opserr << "WARNING invalid element tags\n";
	    return 0;
	}

	theMat = new PyLiq1(idata[0], MAT_TAG_PyLiq1,idata[1], ddata[0], ddata[1], ddata[2],
			    ddata[3], ddata[4],eleTags[0],eleTags[1],theDomain);
    }

    return theMat;
}

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
	PyConstructorType = 1;
}

PyLiq1::PyLiq1(int tag, int classtag, int soil, double p_ult, double y_50,
		double dragratio, double dash_pot, double p_res, 
		Domain *the_Domain, TimeSeries *the_Series)
:PySimple1(tag, classtag, soil, p_ult, y_50, dragratio, dash_pot),
pRes(p_res), theDomain(the_Domain), theSeries(the_Series)
{
	// Initialize PySimple variables and history variables
	//
	this->revertToStart();
    initialTangent = Tangent;
	PyConstructorType = 2;
}


/////////////////////////////////////////////////////////////////////
//	Default constructor

PyLiq1::PyLiq1()
:PySimple1(), pRes(0.0), solidElem1(0), solidElem2(0), theDomain(0), theSeries(0)
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
		if(PyConstructorType==2)
			meanConsolStress = getEffectiveStress(theSeries);
		else
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
		if(PyConstructorType==2)
			meanStress = getEffectiveStress(theSeries);
		else
			meanStress = getEffectiveStress();

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
PyLiq1::getEffectiveStress(TimeSeries *theSeries)
{
	return theSeries->getFactor(theDomain->getCurrentTime());
}
double 

PyLiq1::getEffectiveStress(void)
{

	// Default value for meanStress
	double meanStress = meanConsolStress;
	
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

		// Check that the class tags for the solid elements are either for a FourNodeQuad object, a FourNodeQuadUP object, 
		// a 9_4_QuadUP object, a SSPquadUP object, or a SSPquad object
		if(theElement1->getClassTag()!=ELE_TAG_FourNodeQuad && theElement1->getClassTag()!=ELE_TAG_FourNodeQuadUP && 
		   theElement1->getClassTag()!=ELE_TAG_Nine_Four_Node_QuadUP && theElement1->getClassTag()!=ELE_TAG_SSPquadUP && theElement1->getClassTag()!=ELE_TAG_SSPquad)
		{
			opserr << "Element: " << theElement1->getTag() << " cannot be used to read effective stress for a PyLiq1 material." << endln;
			exit(-1);
		}
		if(theElement2->getClassTag()!=ELE_TAG_FourNodeQuad && theElement2->getClassTag()!=ELE_TAG_FourNodeQuadUP && 
		   theElement2->getClassTag()!=ELE_TAG_Nine_Four_Node_QuadUP && theElement2->getClassTag()!=ELE_TAG_SSPquadUP && theElement2->getClassTag()!=ELE_TAG_SSPquad)
		{
			opserr << "Element: " << theElement2->getTag() << " cannot be used to read effective stress for a PyLiq1 material." << endln;
			exit(-1);
		}

		double excessPorePressure = 0.0;
		meanStress = 0.0;
		
		// get mean stress from element1 if it is a FourNodeQuad object
		if(theElement1->getClassTag()==ELE_TAG_FourNodeQuad)
		{
			// It's safe to cast *theElement1 onto the FourNodeQuad class because we already
			// checked the class tags.
			FourNodeQuad *theElement1 = (FourNodeQuad *)(theDomain->getElement(solidElem1));
			
			// If the element is a quad, check that the class tag for the material at each gauss point is 100 for FluidSolidPorous object
			meanStress = 0.0;
			for(int i=0;i<4;i++)
			{
				NDMaterial *NDM = theElement1->theMaterial[i];
				if(NDM->getClassTag()!=ND_TAG_FluidSolidPorousMaterial){
					opserr << "Material: " << NDM->getTag() << " cannot be used to read effective stress for a PyLiq1 material." << endln;
					exit(-1);
				}
				FluidSolidPorousMaterial *theFSPM = (FluidSolidPorousMaterial *)(NDM);
				meanStress += 1.0/8.0*(2.0/3.0*(NDM->getStress())[0] + 1.0/3.0*(NDM->getStress())[1] - theFSPM->trialExcessPressure);
			}
		}
		// get mean stress from element2 if it is a FourNodeQuad object
		if(theElement2->getClassTag()==ELE_TAG_FourNodeQuad)
		{
			// It's safe to cast *theElement1 onto the FourNodeQuad class because we already
			// checked the class tags.
			FourNodeQuad *theElement2 = (FourNodeQuad *)(theDomain->getElement(solidElem2));
			for(int i=0;i<4;i++)
			{
				NDMaterial *NDM = theElement2->theMaterial[i];
				if(NDM->getClassTag()!=ND_TAG_FluidSolidPorousMaterial){
					opserr << "Material: " << NDM->getTag() << " cannot be used to read effective stress for a PyLiq1 material." << endln;
					exit(-1);
				}
				FluidSolidPorousMaterial *theFSPM = (FluidSolidPorousMaterial *)(NDM);
				meanStress += 1.0/8.0*(2.0/3.0*(NDM->getStress())[0] + 1.0/3.0*(NDM->getStress())[1] - theFSPM->trialExcessPressure);
			}
		}

		// get mean stress from element1 if it is a FourNodeQuadUP object
		if(theElement1->getClassTag()==ELE_TAG_FourNodeQuadUP) {
			// It's safe to cast *theElement1 onto the FourNodeQuadUP class because we already checked the class tags.
			FourNodeQuadUP *theElement1 = (FourNodeQuadUP *)(theDomain->getElement(solidElem1));
			meanStress=0.0;
			
			for(int i=0;i<4;i++) {
				NDMaterial *NDM = theElement1->theMaterial[i];
				if(NDM->getClassTag()==ND_TAG_InitialStateAnalysisWrapper) {
					InitialStateAnalysisWrapper *NDM = (InitialStateAnalysisWrapper *)(theElement1->theMaterial);
					if(NDM->getMainClassTag()!=ND_TAG_PressureDependMultiYield && NDM->getMainClassTag()!=ND_TAG_PressureDependMultiYield02) {
						opserr << "Material: " << NDM->getTag() << " cannot be used to read effective stress for a PyLiq1 material." << endln;
				    	exit(-1);
					}
				} else if(NDM->getClassTag()!=ND_TAG_PressureDependMultiYield && NDM->getClassTag() !=ND_TAG_PressureDependMultiYield02){
					opserr << "Material: " << NDM->getTag() << " cannot be used to read effective stress for a PyLiq1 material." << endln;
					exit(-1);
				}
				meanStress += 1.0/8.0*(2.0/3.0*(NDM->getStress())[0] + 1.0/3.0*(NDM->getStress())[1]);
			}
		}
		// get mean stress from element 2 if it is a FourNodeQuadUP object
		if(theElement2->getClassTag()==ELE_TAG_FourNodeQuadUP) {
			// It's safe to cast *theElement2 onto the FourNodeQuadUP class because we already checked the class tags.
			FourNodeQuadUP *theElement2 = (FourNodeQuadUP *)(theDomain->getElement(solidElem2));
			for(int i=0;i<4;i++) {
				NDMaterial *NDM = theElement2->theMaterial[i];
				if(NDM->getClassTag()==ND_TAG_InitialStateAnalysisWrapper) {
					InitialStateAnalysisWrapper *NDM = (InitialStateAnalysisWrapper *)(theElement2->theMaterial);
					if(NDM->getMainClassTag()!=ND_TAG_PressureDependMultiYield && NDM->getMainClassTag()!=ND_TAG_PressureDependMultiYield02) {
						opserr << "Material: " << NDM->getTag() << " cannot be used to read effective stress for a PyLiq1 material." << endln;
				    	exit(-1);
					}
				} else if(NDM->getClassTag()!=ND_TAG_PressureDependMultiYield && NDM->getClassTag() !=ND_TAG_PressureDependMultiYield02){
					opserr << "Material: " << NDM->getTag() << " cannot be used to read effective stress for a PyLiq1 material." << endln;
					exit(-1);
				}
				meanStress += 1.0/8.0*(2.0/3.0*(NDM->getStress())[0] + 1.0/3.0*(NDM->getStress())[1]);
			}
		}

		// get mean stress from element1 if it is a 9_4_QuadUP object
		if(theElement1->getClassTag()==ELE_TAG_Nine_Four_Node_QuadUP) {
			// It's safe to cast *theElement1 onto the 9_4_QuadUP class because we already checked the class tags.
			NineFourNodeQuadUP *theElement1 = (NineFourNodeQuadUP *)(theDomain->getElement(solidElem1));
			meanStress=0.0;
			
			for(int i=0;i<9;i++) {
				NDMaterial *NDM = theElement1->theMaterial[i];
				if(NDM->getClassTag()==ND_TAG_InitialStateAnalysisWrapper) {
					InitialStateAnalysisWrapper *NDM = (InitialStateAnalysisWrapper *)(theElement1->theMaterial);
					if(NDM->getMainClassTag()!=ND_TAG_PressureDependMultiYield && NDM->getMainClassTag()!=ND_TAG_PressureDependMultiYield02) {
						opserr << "Material: " << NDM->getTag() << " cannot be used to read effective stress for a PyLiq1 material." << endln;
				    	exit(-1);
					}
				} else if(NDM->getClassTag()!=ND_TAG_PressureDependMultiYield && NDM->getClassTag() !=ND_TAG_PressureDependMultiYield02){
					opserr << "Material: " << NDM->getTag() << " cannot be used to read effective stress for a PyLiq1 material." << endln;
					exit(-1);
				}
				meanStress += 1.0/18.0*(2.0/3.0*(NDM->getStress())[0] + 1.0/3.0*(NDM->getStress())[1]);
			}
		}
		// get mean stress from element2 if it is a 9_4_QuadUP object
		if(theElement2->getClassTag()==ELE_TAG_Nine_Four_Node_QuadUP) {
            // It's safe to cast *theElement2 onto the 9_4_QuadUP class because we already checked the class tags.
			NineFourNodeQuadUP *theElement2 = (NineFourNodeQuadUP *)(theDomain->getElement(solidElem2));
			for(int i=0;i<9;i++) {
				NDMaterial *NDM = theElement2->theMaterial[i];
				if(NDM->getClassTag()==ND_TAG_InitialStateAnalysisWrapper) {
					InitialStateAnalysisWrapper *NDM = (InitialStateAnalysisWrapper *)(theElement2->theMaterial);
					if(NDM->getMainClassTag()!=ND_TAG_PressureDependMultiYield && NDM->getMainClassTag()!=ND_TAG_PressureDependMultiYield02) {
						opserr << "Material: " << NDM->getTag() << " cannot be used to read effective stress for a PyLiq1 material." << endln;
				    	exit(-1);
					}
				} else if(NDM->getClassTag()!=ND_TAG_PressureDependMultiYield && NDM->getClassTag() !=ND_TAG_PressureDependMultiYield02) {
					opserr << "Material: " << NDM->getTag() << " cannot be used to read effective stress for a PyLiq1 material." << endln;
					exit(-1);
				}
				meanStress += 1.0/18.0*(2.0/3.0*(NDM->getStress())[0] + 1.0/3.0*(NDM->getStress())[1]);
			}
		}

		// get mean stress from element1 if it is a SSPquadUP object
		if(theElement1->getClassTag()==ELE_TAG_SSPquadUP) {
			// It's safe to cast *theElement1 onto the SSPquadUP class because we already checked the class tags.
			SSPquadUP *theElement1 = (SSPquadUP *)(theDomain->getElement(solidElem1));
			meanStress = 0.0;
			
			NDMaterial *NDM = theElement1->theMaterial;
			if(NDM->getClassTag()==ND_TAG_InitialStateAnalysisWrapper) {
				InitialStateAnalysisWrapper *NDM = (InitialStateAnalysisWrapper *)(theElement1->theMaterial);
				if(NDM->getMainClassTag()!=ND_TAG_PressureDependMultiYield && NDM->getMainClassTag()!=ND_TAG_PressureDependMultiYield02) {
					opserr << "Material: " << NDM->getTag() << " cannot be used to read effective stress for a PyLiq1 material." << endln;
				    exit(-1);
				}
			} else if(NDM->getClassTag()!=ND_TAG_PressureDependMultiYield && NDM->getClassTag()!=ND_TAG_PressureDependMultiYield02) {
				opserr << "Material: " << NDM->getTag() << " cannot be used to read effective stress for a PyLiq1 material." << endln;
				exit(-1);
			}
			meanStress += 1.0/2.0*(2.0/3.0*(NDM->getStress())[0] + 1.0/3.0*(NDM->getStress())[1]);
		}
		// get mean stress from element 2 if it is a SSPquadUP object
		if(theElement2->getClassTag()==ELE_TAG_SSPquadUP) {
			// It's safe to cast *theElement1 onto the SSPquadUP class because we already checked the class tags.
			SSPquadUP *theElement2 = (SSPquadUP *)(theDomain->getElement(solidElem2));
			
			NDMaterial *NDM = theElement2->theMaterial;
			if(NDM->getClassTag()==ND_TAG_InitialStateAnalysisWrapper) {
				InitialStateAnalysisWrapper *NDM = (InitialStateAnalysisWrapper *)(theElement2->theMaterial);
				if(NDM->getMainClassTag()!=ND_TAG_PressureDependMultiYield && NDM->getMainClassTag()!=ND_TAG_PressureDependMultiYield02) {
					opserr << "Material: " << NDM->getTag() << " cannot be used to read effective stress for a PyLiq1 material." << endln;
				    exit(-1);
				}
			} else if(NDM->getClassTag()!=ND_TAG_PressureDependMultiYield && NDM->getClassTag()!=ND_TAG_PressureDependMultiYield02){
				opserr << "Material: " << NDM->getTag() << " cannot be used to read effective stress for a PyLiq1 material." << endln;
				exit(-1);
			}
			meanStress += 1.0/2.0*(2.0/3.0*(NDM->getStress())[0] + 1.0/3.0*(NDM->getStress())[1]);
		}

		// get mean stress from element1 if it is a SSPquad object
		if(theElement1->getClassTag()==ELE_TAG_SSPquad) {
			// It's safe to cast *theElement1 onto the SSPquad class because we already checked the class tags.
			SSPquad *theElement1 = (SSPquad *)(theDomain->getElement(solidElem1));
			
			// If the element is a SSPquad, check that the class tag for the material at each gauss point is FluidSolidPorous object
			meanStress = 0.0;
				
			NDMaterial *NDM = theElement1->theMaterial;
			if(NDM->getClassTag()!=ND_TAG_FluidSolidPorousMaterial){
				opserr << "Material: " << NDM->getTag() << " cannot be used to read effective stress for a PyLiq1 material." << endln;
				exit(-1);
			}
			FluidSolidPorousMaterial *theFSPM = (FluidSolidPorousMaterial *)(NDM);
			meanStress += 1.0/2.0*(2.0/3.0*(NDM->getStress())[0] + 1.0/3.0*(NDM->getStress())[1] - theFSPM->trialExcessPressure);
		}
		// get mean stress from element2 if it is a SSPquad object
		if(theElement2->getClassTag()==ELE_TAG_SSPquad) {
			// It's safe to cast *theElement2 onto the SSPquad class because we already checked the class tags.
			SSPquad *theElement2 = (SSPquad *)(theDomain->getElement(solidElem2));
			
			// If the element is a SSPquad, check that the class tag for the material at each gauss point is FluidSolidPorous object
			NDMaterial *NDM = theElement2->theMaterial;
			if(NDM->getClassTag()!=ND_TAG_FluidSolidPorousMaterial){
				opserr << "Material: " << NDM->getTag() << " cannot be used to read effective stress for a PyLiq1 material." << endln;
				exit(-1);
			}
			FluidSolidPorousMaterial *theFSPM = (FluidSolidPorousMaterial *)(NDM);
			meanStress += 1.0/2.0*(2.0/3.0*(NDM->getStress())[0] + 1.0/3.0*(NDM->getStress())[1] - theFSPM->trialExcessPressure);
		}
		
	}

	return meanStress;
}



/////////////////////////////////////////////////////////////////////

int PyLiq1::setParameter(const char **argv, int argc, Parameter &param)
{	
  if (argc >= 2)
    if (strcmp(argv[0],"updateMaterialStage") == 0 && atoi(argv[1])==this->getTag()) {
      return param.addObject(1, this);  
    }

  return -1;
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
	if(PyConstructorType==1)
	{
		data(10) = (int)solidElem1;
		data(11) = (int)solidElem2;
	}
	if(PyConstructorType==2)
	{
		data(10) = (int)theSeriesTag;
	}
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
	if(PyConstructorType==1)
	{
		solidElem1        = (int)data(10);
		solidElem2        = (int)data(11);
	}
	if(PyConstructorType==2)
	{
		theSeriesTag = (int)data(10);
	}
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
	if(PyConstructorType==1)
	{
		s << "  solidElem1: " << solidElem1 << endln;
		s << "  solidElem2: " << solidElem2 << endln;
	}
	if(PyConstructorType==2)
	{
		s << "  Time Series Tag: " << theSeries->getTag() << endln;
	}
}

/////////////////////////////////////////////////////////////////////



