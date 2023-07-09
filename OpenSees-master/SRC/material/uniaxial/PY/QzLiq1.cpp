
/*** *********************************************************************
**    Module:	QzLiq1.cpp 
**
**    Purpose:	Q-z material that incorporates liquefaction effects. It gets 
**              the pore pressure from a specified element that contains a 
**              PorousFluidSolid or from provided mean effective stress as a 
**              timeseries data.
**
**    Developed by Sumeet Kumar Sinha
**    (C) Copyright 2020, All Rights Reserved.
**
** ****************************************************************** */


// $Revision: 1.0
// $Date: 2020/6/12
// $Source: /OpenSees/SRC/material/uniaxial/QzLiq1.cpp
// Written: Sumeet Kumar Sinha
// Created: June 2020
// Revision: A
//
// Description: This file contains the class implementation for QzLiq1

#include <QzLiq1.h>
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

// Control on internal iteration between spring components
const int QZmaxIterations = 20;
const double QZtolerance = 1.0e-12;

int QzLiq1::loadStage = 0;
Vector QzLiq1::stressV3(3);
int QzConstructorType = 0;

void* OPS_QzLiq1()
{
    UniaxialMaterial* theMat = 0;
    
    int numdata = OPS_GetNumRemainingInputArgs();
    if (numdata < 8) {
	opserr << "WARNING insufficient arguments\n";
	opserr << "Want: uniaxialMaterial QzLiq1 tag? qzType? qult? z50? suction? dashpot? alpha? solidElem1? solidElem2?\n";
	opserr << "or: uniaxialMaterial QzLiq1 tag? qzType? qult? z50? suction? dashpot? alpha? -timeSeries seriesTag?\n";
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
	theMat = new QzLiq1(idata[0], idata[1], ddata[0], ddata[1], ddata[2], ddata[3], ddata[4],
			    theDomain,theSeries);
	
				 
    } else {

	// back one arg
	OPS_ResetCurrentInputArg(-1);
	
	int eleTags[2];
	numdata = 2;
	if (OPS_GetIntInput(&numdata, eleTags) < 0) {
	    opserr << "WARNING invalid element tags\n";
	    return 0;
	}

	theMat = new QzLiq1(idata[0], idata[1], ddata[0], ddata[1], ddata[2], ddata[3], ddata[4],
			    eleTags[0],eleTags[1],theDomain);
    }

    return theMat;
}

/////////////////////////////////////////////////////////////////////
//	Constructor with data

QzLiq1::QzLiq1(int tag, int qz_type, double q_ult, double z_50, double suction,
		double dash_pot, double alpha, int solid_elem1, int solid_elem2, Domain *the_Domain)
:QzSimple1(tag, qz_type, q_ult, z_50, suction, dash_pot), alpha(alpha),
solidElem1(solid_elem1), solidElem2(solid_elem2), theDomain(the_Domain)
{
	// Initialize QzSimple variables and history variables
	//
	this->revertToStart();
    initialTangent = Tangent;
	QzConstructorType = 1;
}
QzLiq1::QzLiq1(int tag, int qz_type, double q_ult, double z_50, double suction,
		double dash_pot, double alpha, Domain *the_Domain, TimeSeries *the_Series)
:QzSimple1(tag, qz_type, q_ult, z_50, suction, dash_pot),alpha(alpha),
theDomain(the_Domain), theSeries(the_Series)
{
	// Initialize QzSimple variables and history variables
	//
	this->revertToStart();
    initialTangent = Tangent;
	QzConstructorType = 2;
}
/////////////////////////////////////////////////////////////////////
//	Default constructor

QzLiq1::QzLiq1()
:QzSimple1(), solidElem1(0), solidElem2(0), theDomain(0)
{
}
/////////////////////////////////////////////////////////////////////
//	Default destructor
QzLiq1::~QzLiq1()
{
  // Call QzSimple1 destructor even though it does nothing.
  //
  // QzSimple1::~QzSimple1();
}

/////////////////////////////////////////////////////////////////////
int 
QzLiq1::setTrialStrain (double newz, double zRate)
{
	// Call the base class QzSimple1 to take care of the basic q-z behavior.
	//
	QzSimple1::setTrialStrain(newz, zRate);
	Tz = newz;

	// Reset mean consolidation stress if loadStage switched from 0 to 1
	//		Note: currently assuming y-axis is vertical, and that the
	//		out-of-plane normal stress equals sigma-xx.
	//
	if(lastLoadStage ==0 && loadStage ==1){

		if(QzConstructorType==2)
			meanConsolStress = getEffectiveStress(theSeries);
		else
			meanConsolStress = getEffectiveStress();
		if(meanConsolStress == 0.0){
			opserr << "WARNING meanConsolStress is 0 in solid elements, ru will divide by zero";
			opserr << "QzLiq1: " << endln;
			if(QzConstructorType==2)
				opserr << "Effective Stress file seriesTag: " << theSeries->getTag() << endln;
			else
				opserr << "Adjacent solidElems: " << solidElem1 << ", " << solidElem2 << endln;
			exit(-1);
		}
	}
	lastLoadStage = loadStage;

	// Obtain the mean effective stress from the adjacent solid elements,
	//    and calculate ru for scaling of q-z base relation.
	//
	double meanStress;
	if(loadStage == 1) {
		if(QzConstructorType==2)
			meanStress = getEffectiveStress(theSeries);
		else
			meanStress = getEffectiveStress();
		if(meanStress>meanConsolStress)
			meanStress=meanConsolStress;
		Tru = 1.0 - meanStress/meanConsolStress;
		if(Tru > 0.999) Tru = 0.999;
		if(Tru < 0) Tru = 0;
	}
	else {
		Tru = 0.0;
	}

	// Call the base class QzSimple1 to get basic q-z response,
	//
	double baseT = QzSimple1::getStress();
	double baseTangent = QzSimple1::getTangent();

	// Check: Only update Hru if not yet converged (avoiding small changes in Tru).
	//
	if(Tz !=Cz || Tt !=Ct) Hru = Tru;

	// During dilation of the soil (ru dropping), provide a stiff transition
	//   between the old and new scaled q-z relations. This avoids illogical
	//   hardening of the q-z relation (i.e., negative stiffnesses).
	//
	if (Tru < Cru){

		maxTangent = (QzSimple1::Qult/QzSimple1::z50)*pow(1.0-Cru,alpha);

		//  If unloading, follow the old scaled q-z relation until t=0.
		//
		if(Cz>0.0 && Tz<Cz && baseT>0.0){
			Hru = Cru;
		}
		if(Cz<0.0 && Tz>Cz && baseT<0.0){
			Hru = Cru;
		}
		
		//  If above the stiff transition line (between Tru & Cru scaled surfaces)
		
		// double zref = Cz + baseT*(Cru-Hru)/maxTangent;
		double zref = Cz + baseT*(pow(1-Hru,alpha)-pow(1-Cru,alpha))/maxTangent;
		if(Cz>0.0 && Tz>Cz && Tz<zref){
			Hru = 1.0 - pow((Ct + (Tz-Cz)*maxTangent)/baseT,1.0/alpha);
		}
		if(Cz<0.0 && Tz<Cz && Tz>zref){
			Hru = 1.0 - pow((Ct + (Tz-Cz)*maxTangent)/baseT,1.0/alpha);
		}

		if(Hru > Cru) Hru = Cru;
		if(Hru < Tru) Hru = Tru;
	}

	//  Now set the tangent and Tt values accordingly
	//

	Tt = baseT*pow(1.0-Hru,alpha);
	if(Hru==Cru || Hru==Tru){
		Tangent = pow(1.0-Hru,alpha)*baseTangent;
	}
	else {
		Tangent = maxTangent;
	}

	return 0;
}
/////////////////////////////////////////////////////////////////////
double 
QzLiq1::getStress(void)
{
	double dashForce = getStrainRate()*this->getDampTangent();

	// Limit the combined force to qult*(1-ru).
	//
	double tmax = (1.0-QZtolerance)*QzSimple1::Qult*pow(1.0-Hru,alpha);
	if(fabs(Tt + dashForce) >= tmax)
		return tmax*(Tt+dashForce)/fabs(Tt+dashForce);
	else return Tt + dashForce;

}
/////////////////////////////////////////////////////////////////////
double 
QzLiq1::getTangent(void)
{
	return this->Tangent;

}
/////////////////////////////////////////////////////////////////////
double 
QzLiq1::getInitialTangent(void)
{
    return this->initialTangent;
}
/////////////////////////////////////////////////////////////////////
double 
QzLiq1::getDampTangent(void)
{
	// Call the base class QzSimple1 to get basic q-z response,
	//    and then scale by (1-ru).
	//
	double dampTangent = QzSimple1::getDampTangent();
	return dampTangent*pow(1.0-Hru,alpha);
}
/////////////////////////////////////////////////////////////////////
double 
QzLiq1::getStrain(void)
{
	double strain = QzSimple1::getStrain();
    return strain;
}
/////////////////////////////////////////////////////////////////////
double 
QzLiq1::getStrainRate(void)
{
    double strainrate = QzSimple1::getStrainRate();
	return strainrate;
}
/////////////////////////////////////////////////////////////////////
int 
QzLiq1::commitState(void)
{
	// Call the QzSimple1 base function to take care of details.
	//
	QzSimple1::commitState();
	Cz = Tz;
	Ct = Tt;
	Cru= Hru;
	
    return 0;
}

/////////////////////////////////////////////////////////////////////
int 
QzLiq1::revertToLastCommit(void)
{
	// reset to committed values
    //
	QzSimple1::revertToLastCommit();
	Tz = Cz;
	Tt = Ct;
	Hru= Cru;

	return 0;
}

/////////////////////////////////////////////////////////////////////
int 
QzLiq1::revertToStart(void)
{
	// Call the QzSimple1 base function to take care of most details.
	//
	QzSimple1::revertToStart();
	Tz = 0.0;
	Tt = 0.0;
	maxTangent = (QzSimple1::Qult/QzSimple1::z50);

	// Excess pore pressure ratio and pointers
	//
	Tru = 0.0;
	Hru = 0.0;
	meanConsolStress = -QzSimple1::Qult;
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
QzLiq1::getEffectiveStress(TimeSeries *theSeries)
{
	return theSeries->getFactor(theDomain->getCurrentTime());
}

double 
QzLiq1::getEffectiveStress(void)
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
			opserr << "QzLiq1: " << endln;
			opserr << "Adjacent solidElems: " << solidElem1 << ", " << solidElem2 << endln;
			exit(-1);
		}

		// Check that the class tags for the solid elements are either for a FourNodeQuad object, a FourNodeQuadUP object, 
		// a 9_4_QuadUP object, a SSPquadUP object, or a SSPquad object
		if(theElement1->getClassTag()!=ELE_TAG_FourNodeQuad && theElement1->getClassTag()!=ELE_TAG_FourNodeQuadUP && 
		   theElement1->getClassTag()!=ELE_TAG_Nine_Four_Node_QuadUP && theElement1->getClassTag()!=ELE_TAG_SSPquadUP && theElement1->getClassTag()!=ELE_TAG_SSPquad)
		{
			opserr << "Element: " << theElement1->getTag() << " cannot be used to read effective stress for a QzLiq1 material." << endln;
			exit(-1);
		}
		if(theElement2->getClassTag()!=ELE_TAG_FourNodeQuad && theElement2->getClassTag()!=ELE_TAG_FourNodeQuadUP && 
		   theElement2->getClassTag()!=ELE_TAG_Nine_Four_Node_QuadUP && theElement2->getClassTag()!=ELE_TAG_SSPquadUP && theElement2->getClassTag()!=ELE_TAG_SSPquad)
		{
			opserr << "Element: " << theElement2->getTag() << " cannot be used to read effective stress for a QzLiq1 material." << endln;
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
					opserr << "Material: " << NDM->getTag() << " cannot be used to read effective stress for a QzLiq1 material." << endln;
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
					opserr << "Material: " << NDM->getTag() << " cannot be used to read effective stress for a QzLiq1 material." << endln;
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
						opserr << "Material: " << NDM->getTag() << " cannot be used to read effective stress for a QzLiq1 material." << endln;
				    	exit(-1);
					}
				} else if(NDM->getClassTag()!=ND_TAG_PressureDependMultiYield && NDM->getClassTag() !=ND_TAG_PressureDependMultiYield02){
					opserr << "Material: " << NDM->getTag() << " cannot be used to read effective stress for a QzLiq1 material." << endln;
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
						opserr << "Material: " << NDM->getTag() << " cannot be used to read effective stress for a QzLiq1 material." << endln;
				    	exit(-1);
					}
				} else if(NDM->getClassTag()!=ND_TAG_PressureDependMultiYield && NDM->getClassTag() !=ND_TAG_PressureDependMultiYield02){
					opserr << "Material: " << NDM->getTag() << " cannot be used to read effective stress for a QzLiq1 material." << endln;
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
						opserr << "Material: " << NDM->getTag() << " cannot be used to read effective stress for a QzLiq1 material." << endln;
				    	exit(-1);
					}
				} else if(NDM->getClassTag()!=ND_TAG_PressureDependMultiYield && NDM->getClassTag() !=ND_TAG_PressureDependMultiYield02){
					opserr << "Material: " << NDM->getTag() << " cannot be used to read effective stress for a QzLiq1 material." << endln;
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
						opserr << "Material: " << NDM->getTag() << " cannot be used to read effective stress for a QzLiq1 material." << endln;
				    	exit(-1);
					}
				} else if(NDM->getClassTag()!=ND_TAG_PressureDependMultiYield && NDM->getClassTag() !=ND_TAG_PressureDependMultiYield02) {
					opserr << "Material: " << NDM->getTag() << " cannot be used to read effective stress for a QzLiq1 material." << endln;
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
					opserr << "Material: " << NDM->getTag() << " cannot be used to read effective stress for a QzLiq1 material." << endln;
				    exit(-1);
				}
			} else if(NDM->getClassTag()!=ND_TAG_PressureDependMultiYield && NDM->getClassTag()!=ND_TAG_PressureDependMultiYield02) {
				opserr << "Material: " << NDM->getTag() << " cannot be used to read effective stress for a QzLiq1 material." << endln;
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
					opserr << "Material: " << NDM->getTag() << " cannot be used to read effective stress for a QzLiq1 material." << endln;
				    exit(-1);
				}
			} else if(NDM->getClassTag()!=ND_TAG_PressureDependMultiYield && NDM->getClassTag()!=ND_TAG_PressureDependMultiYield02){
				opserr << "Material: " << NDM->getTag() << " cannot be used to read effective stress for a QzLiq1 material." << endln;
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
				opserr << "Material: " << NDM->getTag() << " cannot be used to read effective stress for a QzLiq1 material." << endln;
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
				opserr << "Material: " << NDM->getTag() << " cannot be used to read effective stress for a QzLiq1 material." << endln;
				exit(-1);
			}
			FluidSolidPorousMaterial *theFSPM = (FluidSolidPorousMaterial *)(NDM);
			meanStress += 1.0/2.0*(2.0/3.0*(NDM->getStress())[0] + 1.0/3.0*(NDM->getStress())[1] - theFSPM->trialExcessPressure);
		}
		
	}

	return meanStress;
}

int QzLiq1::setParameter(const char **argv, int argc, Parameter &param)
{
  if (argc >= 2)
    if (strcmp(argv[0],"updateMaterialStage") == 0 && atoi(argv[1]) == this->getTag()) {
      return param.addObject(1, this);  
    }

  return -1;
}

/////////////////////////////////////////////////////////////////////
int 
QzLiq1::updateParameter(int snum,Information &eleInformation)
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
      opserr << "WARNING updateMaterialStage for QzLiq1 material must be 0 or 1";
      opserr << endln;
      exit(-1);
    }
    loadStage = snum;

  return 0;
}

/////////////////////////////////////////////////////////////////////
UniaxialMaterial *
QzLiq1::getCopy(void)
{
	// Make a new instance of this class and then assign it "this" to make a copy.
	//
	QzLiq1 *clone;			// pointer to a QzLiq1 class
	clone = new QzLiq1();	// pointer gets a new instance of QzLiq1
	*clone = *this;			// the clone (dereferenced pointer) = dereferenced this.
	
	return clone;
}

/////////////////////////////////////////////////////////////////////
int 
QzLiq1::sendSelf(int cTag, Channel &theChannel)
{
	// I'm following an example by Frank here.
	//
	int res =0;

	static Vector data(17);

	QzSimple1::sendSelf(cTag, theChannel);

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
	data(10)  = alpha;
	if(QzConstructorType==2)
	{
		data(11) = theSeriesTag;
		data(12) = 0.0;
	}
	if(QzConstructorType==1)
	{
		data(11) = solidElem1;
		data(12) = solidElem2;
	}
	data(13) = meanConsolStress;
	data(14) = loadStage;
	data(15) = lastLoadStage;
	data(16) = initialTangent;
	
  res = theChannel.sendVector(this->getDbTag(), cTag, data);
  if (res < 0) 
    opserr << "QzLiq1::sendSelf() - failed to send data\n";

  return res;
}

/////////////////////////////////////////////////////////////////////
int 
QzLiq1::recvSelf(int cTag, Channel &theChannel, 
			       FEM_ObjectBroker &theBroker)
{
  int res = 0;

  static Vector data(17);
  res = theChannel.recvVector(this->getDbTag(), cTag, data);
  
  if (res < 0) {
      opserr << "QzLiq1::recvSelf() - failed to receive data\n";
      this->setTag(0);      
  }
  else {
    this->setTag((int)data(0));

	QzSimple1::recvSelf(cTag, theChannel, theBroker);

	Tz    = data(1);
	Cz    = data(2);
	Tt    = data(3);
	Ct    = data(4);
	Tangent    = data(5);
	maxTangent =data(6);
	Tru        = data(7);
	Cru        = data(8);
	Hru        = data(9);
	alpha      = data(10);
	if(QzConstructorType==1)
	{
		solidElem1        = (int)data(11);
		solidElem2        = (int)data(12);
	}
	if(QzConstructorType==2)
	{
		theSeriesTag = (int)data(11);
	}
	meanConsolStress  = data(13);
	loadStage         = (int)data(14);
	lastLoadStage     = (int)data(15);
	initialTangent    = data(16);

	// set the trial quantities
	this->revertToLastCommit();
  }
 
  return res;
}

/////////////////////////////////////////////////////////////////////
void 
QzLiq1::Print(OPS_Stream &s, int flag)
{
    s << "QzLiq1, tag: " << this->getTag() << endln;
    s << "  QzType: " << QzType << endln;
    s << "  Qult: " << Qult << endln;
    s << "  z50: " << z50 << endln;
	s << "  dashpot: " << dashpot << endln;
	s << "  alpha: " << alpha << endln;
	if(QzConstructorType==1)
	{
		s << "  solidElem1: " << solidElem1 << endln;
		s << "  solidElem2: " << solidElem2 << endln;
	}
	if(QzConstructorType==2)
	{
		s << "  Time Series Tag: " << theSeries->getTag() << endln;
	}

}




