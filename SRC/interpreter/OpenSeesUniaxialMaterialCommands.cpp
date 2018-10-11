/* *****************************************************************************
Copyright (c) 2015-2017, The Regents of the University of California (Regents).
All rights reserved.

Redistribution and use in source and binary forms, with or without 
modification, are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this
   list of conditions and the following disclaimer.
2. Redistributions in binary form must reproduce the above copyright notice,
   this list of conditions and the following disclaimer in the documentation
   and/or other materials provided with the distribution.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR
ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
(INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

The views and conclusions contained in the software and documentation are those
of the authors and should not be interpreted as representing official policies,
either expressed or implied, of the FreeBSD Project.

REGENTS SPECIFICALLY DISCLAIMS ANY WARRANTIES, INCLUDING, BUT NOT LIMITED TO, 
THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
THE SOFTWARE AND ACCOMPANYING DOCUMENTATION, IF ANY, PROVIDED HEREUNDER IS 
PROVIDED "AS IS". REGENTS HAS NO OBLIGATION TO PROVIDE MAINTENANCE, SUPPORT, 
UPDATES, ENHANCEMENTS, OR MODIFICATIONS.

*************************************************************************** */


// Description: command to create uniaxial material

#include <UniaxialMaterial.h>
#include <elementAPI.h>
#include <map>
#include <LimitCurve.h>

// missing or incomplete uniaixal materials:
// Fedeas
// theMaterial = TclModelBuilder_addFedeasMaterial(clientData, interp, argc, argv);

// Drain
// theMaterial = TclModelBuilder_addDrainMaterial(clientData, interp, argc, argv);

// SNAP
// theMaterial = TclModelBuilder_addSnapMaterial(clientData, interp, argc, argv);

// Py, Tz, Qz models
// theMaterial = TclModelBuilder_addPyTzQzMaterial(clientData, interp, argc, argv, theDomain);

// LimitState
// theMaterial = Tcl_AddLimitStateMaterial(clientData, interp, argc, argv);

// UniaxialPackageCommand *matCommands = theUniaxialPackageCommands;

// material in a routine

// material class exists in a package yet to be loaded

void* OPS_ElasticMaterial();
void* OPS_ElasticPPMaterial();
void* OPS_ParallelMaterial();
void* OPS_SeriesMaterial();
void* OPS_EPPGapMaterial();
void* OPS_ENTMaterial();
void* OPS_Steel01();
void* OPS_Steel02();
void* OPS_Steel03();
void* OPS_Concrete01();
void* OPS_Steel4();
void* OPS_HystereticMaterial();
void* OPS_ReinforcingSteel();
void* OPS_Dodd_Restrepo();
void* OPS_RambergOsgoodSteel();
void* OPS_SteelMPF();
void* OPS_Concrete02();
void* OPS_Concrete04();
void* OPS_Concrete06();
void* OPS_Concrete07();
void* OPS_Concrete01WithSITC();
void* OPS_ConfinedConcrete01Material();
void* OPS_ConcreteD();
void* OPS_FRPConfinedConcrete();
void* OPS_ConcreteCM();
void* OPS_Cast();
void* OPS_ViscousDamper();
void* OPS_BilinearOilDamper();
void* OPS_Bilin();
void* OPS_Bilin02();
void* OPS_ModIMKPeakOriented();
void* OPS_ModIMKPeakOriented02();
void* OPS_ModIMKPinching();
void* OPS_ModIMKPinching02();
void* OPS_SAWSMaterial();
void* OPS_BarSlipMaterial();
void* OPS_Bond_SP01();
void* OPS_FatigueMaterial();
void* OPS_HardeningMaterial();
void* OPS_ImpactMaterial();
void* OPS_HyperbolicGapMaterial();
void* OPS_LimiStateMaterial();
void* OPS_MinMaxMaterial();
void* OPS_TensionOnlyMaterial();
void* OPS_ElasticBilin();
void* OPS_ElasticMultiLinear();
void* OPS_MultiLinear();
void* OPS_InitStrainMaterial();
void* OPS_InitStressMaterial();
void* OPS_PathIndependentMaterial();
void* OPS_Pinching4Material();
void* OPS_ECC01();
void* OPS_SelfCenteringMaterial();
void* OPS_ViscousMaterial();
void* OPS_BoucWenMaterial();
void* OPS_BWBN();
void* OPS_PySimple1();
void* OPS_TzSimple1();
void* OPS_PyLiq1();
void* OPS_TzLiq1();
void* OPS_KikuchiAikenHDR();
void* OPS_KikuchiAikenLRB();
void* OPS_AxialSp();
void* OPS_AxialSpHD();
void* OPS_PinchingLimitState();
void* OPS_CFSWSWP();
void* OPS_CFSSSWP();
void* OPS_SteelBRB();
void* OPS_SimpleFractureMaterial();
void* OPS_Maxwell();
#ifndef _NO_NEW_RESTREPO
    void* OPS_DoddRestr();
#endif
void* OPS_Steel2();
void* OPS_OriginCentered();
void* OPS_HookGap();
void* OPS_FRPConfinedConcrete02();
void* OPS_pyUCLA();
void* OPS_Steel01Thermal();
void* OPS_Steel02Thermal();
void* OPS_ConcretewBeta();
void* OPS_ConcreteSakaiKawashima();
void* OPS_Concrete02Thermal();
void* OPS_ResilienceLow();
void* OPS_ResilienceMaterialHR();
void* OPS_Elastic2();
void* OPS_Backbone();
void* OPS_ConcreteZ01Material();
void* OPS_ConcreteL01Material();
void* OPS_SteelZ01Material();
void* OPS_TendonL01Material();
void* OPS_CableMaterial();
void* OPS_ShearPanelMaterial();
void* OPS_SteelMP();
void* OPS_SmoothPSConcrete();
void* OPS_UniaxialJ2Plasticity();

namespace {

    static UniaxialMaterial *theTestingUniaxialMaterial = 0;

    struct char_cmp {
	bool operator () (const char *a,const char *b) const
	    {
		return strcmp(a,b)<0;
	    }
    };

    typedef std::map<const char *, void *(*)(void), char_cmp> OPS_ParsingFunctionMap;


    static OPS_ParsingFunctionMap uniaxialMaterialsMap;


    static int setUpUniaxialMaterials(void) {
	uniaxialMaterialsMap.insert(std::make_pair("Elastic", &OPS_ElasticMaterial));
	uniaxialMaterialsMap.insert(std::make_pair("ElasticPP", &OPS_ElasticPPMaterial));
	uniaxialMaterialsMap.insert(std::make_pair("Parallel", &OPS_ParallelMaterial));
	uniaxialMaterialsMap.insert(std::make_pair("Series", &OPS_SeriesMaterial));
	uniaxialMaterialsMap.insert(std::make_pair("ElasticPPGap", &OPS_EPPGapMaterial));
	uniaxialMaterialsMap.insert(std::make_pair("ENT", &OPS_ENTMaterial));
	uniaxialMaterialsMap.insert(std::make_pair("Steel01", &OPS_Steel01));
	uniaxialMaterialsMap.insert(std::make_pair("Steel02", &OPS_Steel02));
	uniaxialMaterialsMap.insert(std::make_pair("Steel03", &OPS_Steel03));
	uniaxialMaterialsMap.insert(std::make_pair("Concrete01", &OPS_Concrete01));
	uniaxialMaterialsMap.insert(std::make_pair("Steel4", &OPS_Steel4));
	uniaxialMaterialsMap.insert(std::make_pair("Hysteretic", &OPS_HystereticMaterial));
	uniaxialMaterialsMap.insert(std::make_pair("ReinforcingSteel", &OPS_ReinforcingSteel));
	uniaxialMaterialsMap.insert(std::make_pair("Dodd_Restrepo", &OPS_Dodd_Restrepo));
	uniaxialMaterialsMap.insert(std::make_pair("DoddRestrepo", &OPS_Dodd_Restrepo));
	uniaxialMaterialsMap.insert(std::make_pair("Restrepo", &OPS_Dodd_Restrepo));
	uniaxialMaterialsMap.insert(std::make_pair("RambergOsgoodSteel", &OPS_RambergOsgoodSteel));
	uniaxialMaterialsMap.insert(std::make_pair("RambergOsgood", &OPS_RambergOsgoodSteel));
	uniaxialMaterialsMap.insert(std::make_pair("SteelMPF", &OPS_SteelMPF));
	uniaxialMaterialsMap.insert(std::make_pair("Concrete02", &OPS_Concrete02));
	uniaxialMaterialsMap.insert(std::make_pair("Concrete04", &OPS_Concrete04));
	uniaxialMaterialsMap.insert(std::make_pair("Concrete06", &OPS_Concrete06));
	uniaxialMaterialsMap.insert(std::make_pair("Concrete07", &OPS_Concrete07));
	uniaxialMaterialsMap.insert(std::make_pair("Concrete01WithSITC", &OPS_Concrete01WithSITC));
	uniaxialMaterialsMap.insert(std::make_pair("ConfinedConcrete01", &OPS_ConfinedConcrete01Material));
	uniaxialMaterialsMap.insert(std::make_pair("ConfinedConcrete", &OPS_ConfinedConcrete01Material));
	uniaxialMaterialsMap.insert(std::make_pair("ConcreteD", &OPS_ConcreteD));
	uniaxialMaterialsMap.insert(std::make_pair("FRPConfinedConcrete", &OPS_FRPConfinedConcrete));
	uniaxialMaterialsMap.insert(std::make_pair("FRPConfinedConcrete02", &OPS_FRPConfinedConcrete02));
	uniaxialMaterialsMap.insert(std::make_pair("ConcreteCM", &OPS_ConcreteCM));
	uniaxialMaterialsMap.insert(std::make_pair("Cast", &OPS_Cast));
	uniaxialMaterialsMap.insert(std::make_pair("CastFuse", &OPS_Cast));
	uniaxialMaterialsMap.insert(std::make_pair("ViscousDamper", &OPS_ViscousDamper));
	uniaxialMaterialsMap.insert(std::make_pair("BilinearOilDamper", &OPS_BilinearOilDamper));
	uniaxialMaterialsMap.insert(std::make_pair("Bilin", &OPS_Bilin));
	uniaxialMaterialsMap.insert(std::make_pair("BilinMaterial", &OPS_Bilin));
	uniaxialMaterialsMap.insert(std::make_pair("Bilin02", &OPS_Bilin02));
	uniaxialMaterialsMap.insert(std::make_pair("ModIMKPeakOriented", &OPS_ModIMKPeakOriented));
	uniaxialMaterialsMap.insert(std::make_pair("ModIMKPeakOriented02", &OPS_ModIMKPeakOriented02));
	uniaxialMaterialsMap.insert(std::make_pair("ModIMKPinching", &OPS_ModIMKPinching));
	uniaxialMaterialsMap.insert(std::make_pair("ModIMKPinching02", &OPS_ModIMKPinching02));
	uniaxialMaterialsMap.insert(std::make_pair("SAWS", &OPS_SAWSMaterial));
	uniaxialMaterialsMap.insert(std::make_pair("SAWSMaterial", &OPS_SAWSMaterial));
	uniaxialMaterialsMap.insert(std::make_pair("BarSlip", &OPS_BarSlipMaterial));
	uniaxialMaterialsMap.insert(std::make_pair("Bond_SP01", &OPS_Bond_SP01));
	uniaxialMaterialsMap.insert(std::make_pair("Bond", &OPS_Bond_SP01));
	uniaxialMaterialsMap.insert(std::make_pair("Fatigue", &OPS_FatigueMaterial));
	uniaxialMaterialsMap.insert(std::make_pair("Hardening", &OPS_HardeningMaterial));
	uniaxialMaterialsMap.insert(std::make_pair("Impact", &OPS_ImpactMaterial));
	uniaxialMaterialsMap.insert(std::make_pair("ImpactMaterial", &OPS_ImpactMaterial));
	uniaxialMaterialsMap.insert(std::make_pair("HyperbolicGapMaterial", &OPS_HyperbolicGapMaterial));
	uniaxialMaterialsMap.insert(std::make_pair("LimitState", &OPS_LimiStateMaterial));
	uniaxialMaterialsMap.insert(std::make_pair("MinMax", &OPS_MinMaxMaterial));
	uniaxialMaterialsMap.insert(std::make_pair("MinMaxMaterial", &OPS_MinMaxMaterial));
	uniaxialMaterialsMap.insert(std::make_pair("TensionOnly", &OPS_TensionOnlyMaterial));
	uniaxialMaterialsMap.insert(std::make_pair("ElasticBilin", &OPS_ElasticBilin));
	uniaxialMaterialsMap.insert(std::make_pair("ElasticBilinear", &OPS_ElasticBilin));
	uniaxialMaterialsMap.insert(std::make_pair("ElasticMultiLinear", &OPS_ElasticMultiLinear));
	uniaxialMaterialsMap.insert(std::make_pair("MultiLinear", &OPS_MultiLinear));
	uniaxialMaterialsMap.insert(std::make_pair("InitStrainMaterial", &OPS_InitStrainMaterial));
	uniaxialMaterialsMap.insert(std::make_pair("InitStrain", &OPS_InitStrainMaterial));
	uniaxialMaterialsMap.insert(std::make_pair("InitStressMaterial", &OPS_InitStressMaterial));
	uniaxialMaterialsMap.insert(std::make_pair("InitStress", &OPS_InitStressMaterial));
	uniaxialMaterialsMap.insert(std::make_pair("PathIndependent", &OPS_PathIndependentMaterial));
	uniaxialMaterialsMap.insert(std::make_pair("Pinching4", &OPS_Pinching4Material));
	uniaxialMaterialsMap.insert(std::make_pair("ECC01", &OPS_ECC01));
	uniaxialMaterialsMap.insert(std::make_pair("SelfCentering", &OPS_SelfCenteringMaterial));
	uniaxialMaterialsMap.insert(std::make_pair("Viscous", &OPS_ViscousMaterial));
	uniaxialMaterialsMap.insert(std::make_pair("BoucWen", &OPS_BoucWenMaterial));
	uniaxialMaterialsMap.insert(std::make_pair("BWBN", &OPS_BWBN));
	uniaxialMaterialsMap.insert(std::make_pair("PySimple1", &OPS_PySimple1));
	uniaxialMaterialsMap.insert(std::make_pair("TzSimple1", &OPS_TzSimple1));
	uniaxialMaterialsMap.insert(std::make_pair("PyLiq1", &OPS_PyLiq1));
	uniaxialMaterialsMap.insert(std::make_pair("TzLiq1", &OPS_TzLiq1));
	uniaxialMaterialsMap.insert(std::make_pair("KikuchiAikenHDR", &OPS_KikuchiAikenHDR));
	uniaxialMaterialsMap.insert(std::make_pair("KikuchiAikenLRB", &OPS_KikuchiAikenLRB));
	uniaxialMaterialsMap.insert(std::make_pair("AxialSp", &OPS_AxialSp));
	uniaxialMaterialsMap.insert(std::make_pair("AxialSpHD", &OPS_AxialSpHD));
	uniaxialMaterialsMap.insert(std::make_pair("PinchingLimitStateMaterial", &OPS_PinchingLimitState));
	uniaxialMaterialsMap.insert(std::make_pair("CFSWSWP", &OPS_CFSWSWP));
	uniaxialMaterialsMap.insert(std::make_pair("CFSSSWP", &OPS_CFSSSWP));
	uniaxialMaterialsMap.insert(std::make_pair("SteelBRB", &OPS_SteelBRB));
	uniaxialMaterialsMap.insert(std::make_pair("SimpleFractureMaterial", &OPS_SimpleFractureMaterial));
	uniaxialMaterialsMap.insert(std::make_pair("SimpleFracture", &OPS_SimpleFractureMaterial));
	uniaxialMaterialsMap.insert(std::make_pair("Maxwell", &OPS_Maxwell));
	uniaxialMaterialsMap.insert(std::make_pair("MaxwellMaterial", &OPS_Maxwell));
#ifndef _NO_NEW_RESTREPO
	uniaxialMaterialsMap.insert(std::make_pair("DoddRestr", &OPS_DoddRestr));
#endif
	uniaxialMaterialsMap.insert(std::make_pair("Steel2", &OPS_Steel2));
	uniaxialMaterialsMap.insert(std::make_pair("OriginCentered", &OPS_OriginCentered));
	uniaxialMaterialsMap.insert(std::make_pair("HookGap", &OPS_HookGap));
	uniaxialMaterialsMap.insert(std::make_pair("pyUCLA", &OPS_pyUCLA));
	uniaxialMaterialsMap.insert(std::make_pair("PYUCLA", &OPS_pyUCLA));
	uniaxialMaterialsMap.insert(std::make_pair("Steel01Thermal", &OPS_Steel01Thermal));
	uniaxialMaterialsMap.insert(std::make_pair("Steel02Thermal", &OPS_Steel02Thermal));
	uniaxialMaterialsMap.insert(std::make_pair("ConcretewBeta", &OPS_ConcretewBeta));
	uniaxialMaterialsMap.insert(std::make_pair("ConcreteSakaiKawashima", &OPS_ConcreteSakaiKawashima));
	uniaxialMaterialsMap.insert(std::make_pair("Concrete02Thermal", &OPS_Concrete02Thermal));
	uniaxialMaterialsMap.insert(std::make_pair("ResilienceLow", &OPS_ResilienceLow));
	uniaxialMaterialsMap.insert(std::make_pair("ResilienceMaterialHR", &OPS_ResilienceMaterialHR));
	uniaxialMaterialsMap.insert(std::make_pair("Elastic2", &OPS_Elastic2));
	uniaxialMaterialsMap.insert(std::make_pair("Backbone", &OPS_Backbone));
	uniaxialMaterialsMap.insert(std::make_pair("ConcreteZ01Material", &OPS_ConcreteZ01Material));
	uniaxialMaterialsMap.insert(std::make_pair("ConcreteZ01", &OPS_ConcreteZ01Material));
	uniaxialMaterialsMap.insert(std::make_pair("ConcreteL01Material", &OPS_ConcreteL01Material));
	uniaxialMaterialsMap.insert(std::make_pair("ConcreteL01", &OPS_ConcreteL01Material));
	uniaxialMaterialsMap.insert(std::make_pair("SteelZ01Material", &OPS_SteelZ01Material));
	uniaxialMaterialsMap.insert(std::make_pair("SteelZ01", &OPS_SteelZ01Material));
	uniaxialMaterialsMap.insert(std::make_pair("TendonL01Material", &OPS_TendonL01Material));
	uniaxialMaterialsMap.insert(std::make_pair("TendonL01", &OPS_TendonL01Material));
	uniaxialMaterialsMap.insert(std::make_pair("Cable", &OPS_CableMaterial));
	uniaxialMaterialsMap.insert(std::make_pair("ShearPanel", &OPS_ShearPanelMaterial));
	uniaxialMaterialsMap.insert(std::make_pair("SteelMP", &OPS_SteelMP));
	uniaxialMaterialsMap.insert(std::make_pair("SmoothPSConcrete", &OPS_SmoothPSConcrete));
	uniaxialMaterialsMap.insert(std::make_pair("UniaxialJ2Plasticity", &OPS_UniaxialJ2Plasticity));

	return 0;
    }

}

int
OPS_UniaxialMaterial()
{
    static bool initDone = false;
    if (initDone == false) {
	setUpUniaxialMaterials();
	initDone = true;
    }

    if (OPS_GetNumRemainingInputArgs() < 2) {
	opserr<<"WARNING too few arguments: uniaxialMaterial type? tag? ...\n";
	return -1;
    }

    const char* matType = OPS_GetString();

    OPS_ParsingFunctionMap::const_iterator iter = uniaxialMaterialsMap.find(matType);
    if (iter == uniaxialMaterialsMap.end()) {
	opserr<<"WARNING material type " << matType << " is unknown\n";
	return -1;
    }

    UniaxialMaterial* theMaterial = (UniaxialMaterial*) (*iter->second)();
    if (theMaterial == 0) {
	return -1;
    }

    // Now add the material to the modelBuilder
    if (OPS_addUniaxialMaterial(theMaterial) == false) {
	opserr<<"ERROR could not add uniaaialMaterial.\n";
	delete theMaterial; // invoke the material objects destructor, otherwise mem leak
	return -1;
    }

    return 0;

}

int OPS_testUniaxialMaterial()
{
    if (OPS_GetNumRemainingInputArgs() != 1) {
	opserr<<"testUniaxialMaterial - You must provide a material tag.\n";
	return -1;
    }

    int tag;
    int numData = 1;
    if (OPS_GetIntInput(&numData, &tag) < 0) {
	opserr<<"invalid int value\n";
	return -1;
    }

    UniaxialMaterial* mat =OPS_getUniaxialMaterial(tag);

    if (mat == 0) {
	opserr<<"testUniaxialMaterial - Material Not Found.\n";
	return -1;
    }

    theTestingUniaxialMaterial = mat;

    return 0;
}

int OPS_setStrain()
{
    if (OPS_GetNumRemainingInputArgs() != 1) {
	opserr<<"testUniaxialMaterial - You must provide a strain value.\n";
	return -1;
    }

    UniaxialMaterial* material = theTestingUniaxialMaterial;

    if (material == 0) {
	opserr<<"setStrain WARNING no active UniaxialMaterial - use testUniaxialMaterial command.\n";
	return -1;
    }

    double strain;
    int numData = 1;
    if (OPS_GetDoubleInput(&numData, &strain) < 0) {
	opserr<<"invalid double value\n";
	return -1;
    }

    material->setTrialStrain(strain);
    material->commitState();

    return 0;
}

int OPS_getStrain()
{
    UniaxialMaterial* material = theTestingUniaxialMaterial;
    if (material == 0) {
	opserr<<"getStrain WARNING no active UniaxialMaterial - use testUniaxialMaterial command.\n";
	return -1;
    }

    double strain = material->getStrain();

    int numData = 1;

    if (OPS_SetDoubleOutput(&numData, &strain) < 0) {
	opserr<<"failed to set strain\n";
	return -1;
    }

    return 0;
}

int OPS_getStress()
{
    UniaxialMaterial* material = theTestingUniaxialMaterial;
    if (material == 0) {
	opserr<<"getStrain WARNING no active UniaxialMaterial - use testUniaxialMaterial command.\n";
	return -1;
    }

    double stress = material->getStress();

    int numData = 1;

    if (OPS_SetDoubleOutput(&numData, &stress) < 0) {
	opserr<<"failed to set stress\n";
	return -1;
    }

    return 0;
}

int OPS_getTangent()
{
    UniaxialMaterial* material = theTestingUniaxialMaterial;
    if (material == 0) {
	opserr<<"getStrain WARNING no active UniaxialMaterial - use testUniaxialMaterial command.\n";
	return -1;
    }

    double tangent = material->getTangent();

    int numData = 1;

    if (OPS_SetDoubleOutput(&numData, &tangent) < 0) {
	opserr<<"failed to set tangent\n";
	return -1;
    }

    return 0;
}

int OPS_getDampTangent()
{
    UniaxialMaterial* material = theTestingUniaxialMaterial;
    if (material == 0) {
	opserr<<"getStrain WARNING no active UniaxialMaterial - use testUniaxialMaterial command.\n";
	return -1;
    }

    double tangent = material->getDampTangent();

    int numData = 1;

    if (OPS_SetDoubleOutput(&numData, &tangent) < 0) {
	opserr<<"failed to set damp tangent\n";
	return -1;
    }

    return 0;
}

void* OPS_RotationShearCurve();
void* OPS_ThreePointCurve();
void* OPS_ShearCurve();

int OPS_LimitCurve()
{
    // Make sure there is a minimum number of arguments
    if (OPS_GetNumRemainingInputArgs() < 6)
    {
	opserr << "WARNING insufficient number of limit curve arguments\n";
	opserr << "Want: limitCurve type? tag? <specific curve args>" << endln;
	return -1;
    }

    const char* type = OPS_GetString();

    // Pointer to a limit curve that will be added to the model builder
    LimitCurve *theCurve = 0;

    if (strcmp(type, "Axial") == 0) {

	opserr << "WARNING to be implemented ...\n";
	return -1;
	
	
    } else if (strcmp(type, "RotationShearCurve") == 0) {

	void *theRSC = OPS_RotationShearCurve();
	if (theRSC != 0) {
	    theCurve = (LimitCurve *)theRSC;
	} else {
	    return -1;
	}
	  
    } else if (strcmp(type,"ThreePoint") == 0) {
	
	void* curve = OPS_RotationShearCurve();
	if (curve != 0) {
	    theCurve = (LimitCurve *)curve;
	} else {
	    return -1;
	}
	
    } else if (strcmp(type,"Shear") == 0) {

	void* curve = OPS_ShearCurve();
	if (curve != 0) {
	    theCurve = (LimitCurve *)curve;
	} else {
	    return -1;
	}
	
    } else {
	opserr << "WARNING type of limit curve is unknown\n";
	return -1;
    }

    

    // Ensure we have created the Material, out of memory if got here and no material
    if (theCurve == 0) {
	opserr << "WARNING ran out of memory creating limitCurve\n";
	return -1;
    }

    // Now add the material to the modelBuilder
    if (OPS_addLimitCurve(theCurve) == false) {
	opserr << "WARNING could not add limitCurve to the domain\n";
	delete theCurve; // invoke the material objects destructor, otherwise mem leak
	return -1;
    }

    return 0;
}
