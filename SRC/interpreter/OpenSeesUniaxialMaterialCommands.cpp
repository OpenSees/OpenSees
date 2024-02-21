/* *****************************************************************************
Copyright (c) 2015-2017, The Regents of the University of
California (Regents). All rights reserved.

Redistribution and use in source and binary forms, with or
without modification, are permitted provided that the following
conditions are met:

1. Redistributions of source code must retain the above copyright
notice, this list of conditions and the following disclaimer.
2. Redistributions in binary form must reproduce the above
copyright notice, this list of conditions and the following
disclaimer in the documentation and/or other materials provided
with the distribution.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND
CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES,
INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF
MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS
BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED
TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON
ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY
OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
POSSIBILITY OF SUCH DAMAGE.

The views and conclusions contained in the software and
documentation are those of the authors and should not be
interpreted as representing official policies, either expressed
or implied, of the FreeBSD Project.

REGENTS SPECIFICALLY DISCLAIMS ANY WARRANTIES, INCLUDING, BUT NOT
LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS
FOR A PARTICULAR PURPOSE. THE SOFTWARE AND ACCOMPANYING
DOCUMENTATION, IF ANY, PROVIDED HEREUNDER IS PROVIDED "AS IS".
REGENTS HAS NO OBLIGATION TO PROVIDE MAINTENANCE, SUPPORT,
UPDATES, ENHANCEMENTS, OR MODIFICATIONS.

***************************************************************************
*/

// Description: command to create uniaxial material

#include <HystereticBackbone.h>
#include <LimitCurve.h>
#include <StiffnessDegradation.h>
#include <StrengthDegradation.h>
#include <UniaxialMaterial.h>
#include <UnloadingRule.h>
#include <elementAPI.h>

#include <map>

// missing or incomplete uniaixal materials:
// Fedeas
// theMaterial = TclModelBuilder_addFedeasMaterial(clientData,
// interp, argc, argv);

// Drain
// theMaterial = TclModelBuilder_addDrainMaterial(clientData,
// interp, argc, argv);

// SNAP
// theMaterial = TclModelBuilder_addSnapMaterial(clientData,
// interp, argc, argv);

// Py, Tz, Qz models
// theMaterial = TclModelBuilder_addPyTzQzMaterial(clientData,
// interp, argc, argv, theDomain);

// LimitState
// theMaterial = Tcl_AddLimitStateMaterial(clientData, interp,
// argc, argv);

// UniaxialPackageCommand *matCommands =
// theUniaxialPackageCommands;

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
void* OPS_SteelFractureDI();
void* OPS_Steel02Fatigue();
void* OPS_Steel03();
void* OPS_SPSW02();
void* OPS_Concrete01();
void* OPS_Steel4();
void* OPS_HystereticMaterial();
void* OPS_HystereticSMMaterial();
void* OPS_ReinforcingSteel();
void* OPS_Dodd_Restrepo();
void* OPS_RambergOsgoodSteel();
void* OPS_SteelMPF();
void* OPS_SteelDRC();
void* OPS_Concrete02();
void* OPS_Concrete02IS();
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
void* OPS_DamperMaterial();
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
void* OPS_FlagShapeMaterial();
void* OPS_ImpactMaterial();
void* OPS_HyperbolicGapMaterial();
void* OPS_LimiStateMaterial();
void* OPS_MinMaxMaterial();
void* OPS_PenaltyMaterial();
void* OPS_TensionOnlyMaterial();
void* OPS_ElasticBilin();
void* OPS_ElasticMultiLinear();
void* OPS_ElasticPowerFunc();
void* OPS_MultiLinear();
void* OPS_ContinuumUniaxialMaterial();
void* OPS_InitStrainMaterial();
void* OPS_InitStressMaterial();
void* OPS_PathIndependentMaterial();
void* OPS_MultiplierMaterial();
void* OPS_Pinching4Material();
void* OPS_ECC01();
void* OPS_SelfCenteringMaterial();
void* OPS_ASD_SMA_3K();
void* OPS_ViscousMaterial();
void* OPS_BoucWenMaterial();
void* OPS_BoucWenOriginal();
void* OPS_BoucWenInfill();
void* OPS_BWBN();
void* OPS_PySimple1();
void* OPS_TzSimple1();
void* OPS_QzSimple1();
void* OPS_PySimple2();
void* OPS_TzSimple2();
void* OPS_QzSimple2();
void* OPS_PyLiq1();
void* OPS_TzLiq1();
void* OPS_QzLiq1();
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
void* OPS_ElasticMaterialThermal();
void* OPS_Steel01Thermal();
void* OPS_Steel02Thermal();
void* OPS_ConcretewBeta();
void* OPS_ConcreteSakaiKawashima();
void* OPS_Concrete02Thermal();
void* OPS_ResilienceLow();
void* OPS_ResilienceMaterialHR();
void* OPS_Elastic2Material();
void* OPS_BackboneMaterial();
void* OPS_ConcreteZ01Material();
void* OPS_ConcreteL01Material();
void* OPS_SteelZ01Material();
void* OPS_TendonL01Material();
void* OPS_CableMaterial();
void* OPS_ShearPanelMaterial();
void* OPS_SteelMP();
void* OPS_SmoothPSConcrete();
void* OPS_UniaxialJ2Plasticity();
void* OPS_OOHystereticMaterial();
void* OPS_UVCuniaxial();
void* OPS_GNGMaterial();
void* OPS_IMKBilin();
void* OPS_IMKPinching();
void* OPS_IMKPeakOriented();
void* OPS_SLModel();
void* OPS_SMAMaterial();
void* OPS_FRCC();
void* OPS_ConcreteZBH_original();
void* OPS_ConcreteZBH_fitted();
void* OPS_ConcreteZBH_smoothed();

void *OPS_ArctangentBackbone(void);
void *OPS_ManderBackbone(void);
void *OPS_TrilinearBackbone(void);
void *OPS_BilinearBackbone(void);
void *OPS_MultilinearBackbone(void);
void* OPS_MaterialBackbone();
void* OPS_ReeseStiffClayBelowWS();
void* OPS_ReeseStiffClayAboveWS();
void* OPS_VuggyLimestone();
void* OPS_CementedSoil();
void* OPS_WeakRock();
void* OPS_LiquefiedSand();
void* OPS_RaynorBackbone();
void* OPS_ReeseSandBackbone();
void* OPS_ReeseSoftClayBackbone();
void* OPS_CappedBackbone();
void* OPS_LinearCappedBackbone();

void* OPS_ConstantStiffnessDegradation();
void* OPS_DuctilityStiffnessDegradation();
void* OPS_EnergyStiffnessDegradation();
void* OPS_PincheiraStiffnessDegradation();

void* OPS_ConstantStrengthDegradation();
void* OPS_DuctilityStrengthDegradation();
void* OPS_EnergyStrengthDegradation();
void* OPS_ACIStrengthDegradation();
void* OPS_PetrangeliStrengthDegradation();

void* OPS_ConstantUnloadingRule();
void* OPS_TakedaUnloadingRule();
void* OPS_EnergyUnloadingRule();
void* OPS_KarsanUnloadingRule();

void* OPS_HystereticPoly();  // Salvatore Sessa 14-01-2021 Mail:
                             // salvatore.sessa2@unina.it
void* OPS_HystereticSmooth();  // Salvatore Sessa Mail:
                             // salvatore.sessa2@unina.it
void* OPS_HystereticAsym();  // Salvatore Sessa
                             // salvatore.sessa2@unina.it
void* OPS_DowelType();
void* OPS_DuctileFracture();  // Kuanshi Zhong

void* OPS_TDConcreteEXP(void);
void* OPS_TDConcrete(void);
void* OPS_TDConcreteNL(void);
void* OPS_TDConcreteMC10(void);
void* OPS_TDConcreteMC10NL(void);
void* OPS_CreepMaterial(void);

void* OPS_CoulombDamperMaterial();
void* OPS_GMG_CyclicReinforcedConcrete();

void *OPS_Hertzdamp(void);
void *OPS_JankowskiImpact(void);
void *OPS_ViscoelasticGap(void);

void *OPS_Masonry(void);
void *OPS_Trilinwp(void);
void *OPS_Trilinwp2(void);
void *OPS_Masonryt(void);

void* OPS_Ratchet(void); // Yi Xiao

namespace {

static UniaxialMaterial* theTestingUniaxialMaterial = 0;

struct char_cmp {
  bool operator()(const char* a, const char* b) const {
    return strcmp(a, b) < 0;
  }
};

typedef std::map<const char*, void* (*)(void), char_cmp>
    OPS_ParsingFunctionMap;

static OPS_ParsingFunctionMap uniaxialMaterialsMap;

static int setUpUniaxialMaterials(void) {
  uniaxialMaterialsMap.insert(
      std::make_pair("Elastic", &OPS_ElasticMaterial));
  uniaxialMaterialsMap.insert(
      std::make_pair("ElasticPP", &OPS_ElasticPPMaterial));
  uniaxialMaterialsMap.insert(
      std::make_pair("Parallel", &OPS_ParallelMaterial));
  uniaxialMaterialsMap.insert(
      std::make_pair("Series", &OPS_SeriesMaterial));
  uniaxialMaterialsMap.insert(
      std::make_pair("ElasticPPGap", &OPS_EPPGapMaterial));
  uniaxialMaterialsMap.insert(
      std::make_pair("ENT", &OPS_ENTMaterial));
  uniaxialMaterialsMap.insert(
      std::make_pair("Steel01", &OPS_Steel01));
  uniaxialMaterialsMap.insert(
      std::make_pair("Steel02", &OPS_Steel02));
  uniaxialMaterialsMap.insert(
      std::make_pair("Steel02Fatigue", &OPS_Steel02Fatigue));
  uniaxialMaterialsMap.insert(
      std::make_pair("Steel03", &OPS_Steel03));
  uniaxialMaterialsMap.insert(
      std::make_pair("SPSW02", &OPS_SPSW02));
  uniaxialMaterialsMap.insert(
      std::make_pair("Concrete01", &OPS_Concrete01));
  uniaxialMaterialsMap.insert(
      std::make_pair("Steel4", &OPS_Steel4));
  uniaxialMaterialsMap.insert(
      std::make_pair("Hysteretic", &OPS_HystereticMaterial));
  uniaxialMaterialsMap.insert(
      std::make_pair("HystereticSM", &OPS_HystereticSMMaterial));
  uniaxialMaterialsMap.insert(
      std::make_pair("ReinforcingSteel", &OPS_ReinforcingSteel));
  uniaxialMaterialsMap.insert(
      std::make_pair("Dodd_Restrepo", &OPS_Dodd_Restrepo));
  uniaxialMaterialsMap.insert(
      std::make_pair("DoddRestrepo", &OPS_Dodd_Restrepo));
  uniaxialMaterialsMap.insert(
      std::make_pair("Restrepo", &OPS_Dodd_Restrepo));
  uniaxialMaterialsMap.insert(std::make_pair(
      "RambergOsgoodSteel", &OPS_RambergOsgoodSteel));
  uniaxialMaterialsMap.insert(
      std::make_pair("RambergOsgood", &OPS_RambergOsgoodSteel));
  uniaxialMaterialsMap.insert(
      std::make_pair("SteelMPF", &OPS_SteelMPF));
  uniaxialMaterialsMap.insert(
      std::make_pair("SteelDRC", &OPS_SteelDRC));
  uniaxialMaterialsMap.insert(
      std::make_pair("Concrete02", &OPS_Concrete02));
  uniaxialMaterialsMap.insert(
      std::make_pair("Concrete02IS", &OPS_Concrete02IS));
  uniaxialMaterialsMap.insert(
      std::make_pair("Concrete04", &OPS_Concrete04));
  uniaxialMaterialsMap.insert(
      std::make_pair("Concrete06", &OPS_Concrete06));
  uniaxialMaterialsMap.insert(
      std::make_pair("Concrete07", &OPS_Concrete07));
  uniaxialMaterialsMap.insert(std::make_pair(
      "Concrete01WithSITC", &OPS_Concrete01WithSITC));
  uniaxialMaterialsMap.insert(std::make_pair(
      "ConfinedConcrete01", &OPS_ConfinedConcrete01Material));
  uniaxialMaterialsMap.insert(std::make_pair(
      "ConfinedConcrete", &OPS_ConfinedConcrete01Material));
  uniaxialMaterialsMap.insert(
      std::make_pair("ConcreteD", &OPS_ConcreteD));
  uniaxialMaterialsMap.insert(std::make_pair(
      "FRPConfinedConcrete", &OPS_FRPConfinedConcrete));
  uniaxialMaterialsMap.insert(std::make_pair(
      "FRPConfinedConcrete02", &OPS_FRPConfinedConcrete02));
  uniaxialMaterialsMap.insert(
      std::make_pair("ConcreteCM", &OPS_ConcreteCM));
  uniaxialMaterialsMap.insert(std::make_pair("Cast", &OPS_Cast));
  uniaxialMaterialsMap.insert(
      std::make_pair("CastFuse", &OPS_Cast));
  uniaxialMaterialsMap.insert(
      std::make_pair("ViscousDamper", &OPS_ViscousDamper));
  uniaxialMaterialsMap.insert(
      std::make_pair("Damper", &OPS_DamperMaterial));	
  uniaxialMaterialsMap.insert(
      std::make_pair("DamperMaterial", &OPS_DamperMaterial));	
  uniaxialMaterialsMap.insert(std::make_pair(
      "BilinearOilDamper", &OPS_BilinearOilDamper));
  uniaxialMaterialsMap.insert(
      std::make_pair("Bilin", &OPS_Bilin));
  uniaxialMaterialsMap.insert(
      std::make_pair("BilinMaterial", &OPS_Bilin));
  uniaxialMaterialsMap.insert(
      std::make_pair("Bilin02", &OPS_Bilin02));
  uniaxialMaterialsMap.insert(std::make_pair(
      "ModIMKPeakOriented", &OPS_ModIMKPeakOriented));
  uniaxialMaterialsMap.insert(std::make_pair(
      "ModIMKPeakOriented02", &OPS_ModIMKPeakOriented02));
  uniaxialMaterialsMap.insert(
      std::make_pair("ModIMKPinching", &OPS_ModIMKPinching));
  uniaxialMaterialsMap.insert(
      std::make_pair("ModIMKPinching02", &OPS_ModIMKPinching02));
  uniaxialMaterialsMap.insert(
      std::make_pair("SAWS", &OPS_SAWSMaterial));
  uniaxialMaterialsMap.insert(
      std::make_pair("SAWSMaterial", &OPS_SAWSMaterial));
  uniaxialMaterialsMap.insert(
      std::make_pair("BarSlip", &OPS_BarSlipMaterial));
  uniaxialMaterialsMap.insert(
      std::make_pair("Bond_SP01", &OPS_Bond_SP01));
  uniaxialMaterialsMap.insert(
      std::make_pair("Bond", &OPS_Bond_SP01));
  uniaxialMaterialsMap.insert(
      std::make_pair("Fatigue", &OPS_FatigueMaterial));
  uniaxialMaterialsMap.insert(
      std::make_pair("Hardening", &OPS_HardeningMaterial));
  uniaxialMaterialsMap.insert(
      std::make_pair("FlagShape", &OPS_FlagShapeMaterial));  
  uniaxialMaterialsMap.insert(
      std::make_pair("Impact", &OPS_ImpactMaterial));
  uniaxialMaterialsMap.insert(
      std::make_pair("ImpactMaterial", &OPS_ImpactMaterial));
  uniaxialMaterialsMap.insert(std::make_pair(
      "HyperbolicGapMaterial", &OPS_HyperbolicGapMaterial));
  uniaxialMaterialsMap.insert(
      std::make_pair("LimitState", &OPS_LimiStateMaterial));
  uniaxialMaterialsMap.insert(
      std::make_pair("MinMax", &OPS_MinMaxMaterial));
  uniaxialMaterialsMap.insert(
      std::make_pair("MinMaxMaterial", &OPS_MinMaxMaterial));
  uniaxialMaterialsMap.insert(
      std::make_pair("Penalty", &OPS_PenaltyMaterial));
  uniaxialMaterialsMap.insert(
      std::make_pair("TensionOnly", &OPS_TensionOnlyMaterial));  
  uniaxialMaterialsMap.insert(
      std::make_pair("ElasticBilin", &OPS_ElasticBilin));
  uniaxialMaterialsMap.insert(
      std::make_pair("ElasticBilinear", &OPS_ElasticBilin));
  uniaxialMaterialsMap.insert(std::make_pair(
      "ElasticMultiLinear", &OPS_ElasticMultiLinear));
  uniaxialMaterialsMap.insert(
      std::make_pair("ElasticPowerFunc", &OPS_ElasticPowerFunc));
  uniaxialMaterialsMap.insert(
      std::make_pair("MultiLinear", &OPS_MultiLinear));
  uniaxialMaterialsMap.insert(std::make_pair(
      "Continuum", &OPS_ContinuumUniaxialMaterial));
  uniaxialMaterialsMap.insert(std::make_pair(
      "InitStrainMaterial", &OPS_InitStrainMaterial));
  uniaxialMaterialsMap.insert(
      std::make_pair("InitStrain", &OPS_InitStrainMaterial));
  uniaxialMaterialsMap.insert(std::make_pair(
      "InitStressMaterial", &OPS_InitStressMaterial));
  uniaxialMaterialsMap.insert(
      std::make_pair("InitStress", &OPS_InitStressMaterial));
  uniaxialMaterialsMap.insert(std::make_pair(
      "PathIndependent", &OPS_PathIndependentMaterial));
  uniaxialMaterialsMap.insert(
      std::make_pair("Multiplier", &OPS_MultiplierMaterial));
  uniaxialMaterialsMap.insert(
      std::make_pair("Pinching4", &OPS_Pinching4Material));
  uniaxialMaterialsMap.insert(
      std::make_pair("ECC01", &OPS_ECC01));
  uniaxialMaterialsMap.insert(std::make_pair(
      "SelfCentering", &OPS_SelfCenteringMaterial));
  uniaxialMaterialsMap.insert(
      std::make_pair("ASD_SMA_3K", &OPS_ASD_SMA_3K));
  uniaxialMaterialsMap.insert(
      std::make_pair("Viscous", &OPS_ViscousMaterial));
  uniaxialMaterialsMap.insert(
      std::make_pair("BoucWen", &OPS_BoucWenMaterial));
  uniaxialMaterialsMap.insert(
      std::make_pair("BoucWenOriginal", &OPS_BoucWenOriginal));
  uniaxialMaterialsMap.insert(			      
      std::make_pair("BoucWenInfill", &OPS_BoucWenInfill));  
  uniaxialMaterialsMap.insert(
      std::make_pair("BWBN", &OPS_BWBN));
  uniaxialMaterialsMap.insert(
      std::make_pair("PySimple1", &OPS_PySimple1));
  uniaxialMaterialsMap.insert(
      std::make_pair("TzSimple1", &OPS_TzSimple1));
  uniaxialMaterialsMap.insert(
      std::make_pair("QzSimple1", &OPS_QzSimple1));
  uniaxialMaterialsMap.insert(
      std::make_pair("PySimple2", &OPS_PySimple2));
  uniaxialMaterialsMap.insert(
      std::make_pair("TzSimple2", &OPS_TzSimple2));
  uniaxialMaterialsMap.insert(
      std::make_pair("QzSimple2", &OPS_QzSimple2));
  uniaxialMaterialsMap.insert(
      std::make_pair("PyLiq1", &OPS_PyLiq1));
  uniaxialMaterialsMap.insert(
      std::make_pair("TzLiq1", &OPS_TzLiq1));
  uniaxialMaterialsMap.insert(
      std::make_pair("QzLiq1", &OPS_QzLiq1));
  uniaxialMaterialsMap.insert(
      std::make_pair("KikuchiAikenHDR", &OPS_KikuchiAikenHDR));
  uniaxialMaterialsMap.insert(
      std::make_pair("KikuchiAikenLRB", &OPS_KikuchiAikenLRB));
  uniaxialMaterialsMap.insert(
      std::make_pair("AxialSp", &OPS_AxialSp));
  uniaxialMaterialsMap.insert(
      std::make_pair("AxialSpHD", &OPS_AxialSpHD));
  uniaxialMaterialsMap.insert(std::make_pair(
      "PinchingLimitStateMaterial", &OPS_PinchingLimitState));
  uniaxialMaterialsMap.insert(
      std::make_pair("CFSWSWP", &OPS_CFSWSWP));
  uniaxialMaterialsMap.insert(
      std::make_pair("CFSSSWP", &OPS_CFSSSWP));
  uniaxialMaterialsMap.insert(
      std::make_pair("SteelBRB", &OPS_SteelBRB));
  uniaxialMaterialsMap.insert(std::make_pair(
      "SimpleFractureMaterial", &OPS_SimpleFractureMaterial));
  uniaxialMaterialsMap.insert(std::make_pair(
      "SimpleFracture", &OPS_SimpleFractureMaterial));
  uniaxialMaterialsMap.insert(
      std::make_pair("Maxwell", &OPS_Maxwell));
  uniaxialMaterialsMap.insert(
      std::make_pair("MaxwellMaterial", &OPS_Maxwell));
#ifndef _NO_NEW_RESTREPO
  uniaxialMaterialsMap.insert(
      std::make_pair("DoddRestr", &OPS_DoddRestr));
#endif
  uniaxialMaterialsMap.insert(
      std::make_pair("Steel2", &OPS_Steel2));
  uniaxialMaterialsMap.insert(
      std::make_pair("OriginCentered", &OPS_OriginCentered));
  uniaxialMaterialsMap.insert(
      std::make_pair("HookGap", &OPS_HookGap));
  uniaxialMaterialsMap.insert(
      std::make_pair("pyUCLA", &OPS_pyUCLA));
  uniaxialMaterialsMap.insert(
      std::make_pair("PYUCLA", &OPS_pyUCLA));
  uniaxialMaterialsMap.insert(
      std::make_pair("ElasticThermal", &OPS_ElasticMaterialThermal));
  uniaxialMaterialsMap.insert(			      
      std::make_pair("Steel01Thermal", &OPS_Steel01Thermal));
  uniaxialMaterialsMap.insert(
      std::make_pair("Steel02Thermal", &OPS_Steel02Thermal));
  uniaxialMaterialsMap.insert(
      std::make_pair("ConcretewBeta", &OPS_ConcretewBeta));
  uniaxialMaterialsMap.insert(std::make_pair(
      "ConcreteSakaiKawashima", &OPS_ConcreteSakaiKawashima));
  uniaxialMaterialsMap.insert(std::make_pair(
      "Concrete02Thermal", &OPS_Concrete02Thermal));
  uniaxialMaterialsMap.insert(
      std::make_pair("ResilienceLow", &OPS_ResilienceLow));
  uniaxialMaterialsMap.insert(std::make_pair(
      "ResilienceMaterialHR", &OPS_ResilienceMaterialHR));
  uniaxialMaterialsMap.insert(
      std::make_pair("Elastic2", &OPS_Elastic2Material));
  uniaxialMaterialsMap.insert(
      std::make_pair("Backbone", &OPS_BackboneMaterial));
  uniaxialMaterialsMap.insert(std::make_pair(
      "ConcreteZ01Material", &OPS_ConcreteZ01Material));
  uniaxialMaterialsMap.insert(
      std::make_pair("ConcreteZ01", &OPS_ConcreteZ01Material));
  uniaxialMaterialsMap.insert(std::make_pair(
      "ConcreteL01Material", &OPS_ConcreteL01Material));
  uniaxialMaterialsMap.insert(
      std::make_pair("ConcreteL01", &OPS_ConcreteL01Material));
  uniaxialMaterialsMap.insert(
      std::make_pair("SteelZ01Material", &OPS_SteelZ01Material));
  uniaxialMaterialsMap.insert(
      std::make_pair("SteelZ01", &OPS_SteelZ01Material));
  uniaxialMaterialsMap.insert(std::make_pair(
      "TendonL01Material", &OPS_TendonL01Material));
  uniaxialMaterialsMap.insert(
      std::make_pair("TendonL01", &OPS_TendonL01Material));
  uniaxialMaterialsMap.insert(
      std::make_pair("Cable", &OPS_CableMaterial));
  uniaxialMaterialsMap.insert(
      std::make_pair("ShearPanel", &OPS_ShearPanelMaterial));
  uniaxialMaterialsMap.insert(
      std::make_pair("SteelMP", &OPS_SteelMP));
  uniaxialMaterialsMap.insert(
      std::make_pair("SmoothPSConcrete", &OPS_SmoothPSConcrete));
  uniaxialMaterialsMap.insert(std::make_pair(
      "UniaxialJ2Plasticity", &OPS_UniaxialJ2Plasticity));
  uniaxialMaterialsMap.insert(
      std::make_pair("OOHysteretic", &OPS_OOHystereticMaterial));
  uniaxialMaterialsMap.insert(
      std::make_pair("UVCuniaxial", &OPS_UVCuniaxial));
  uniaxialMaterialsMap.insert(
      std::make_pair("GNG", &OPS_GNGMaterial));
  uniaxialMaterialsMap.insert(			      
      std::make_pair("SteelFractureDI", &OPS_SteelFractureDI));
  uniaxialMaterialsMap.insert(
      std::make_pair("IMKBilin", &OPS_IMKBilin));
  uniaxialMaterialsMap.insert(
      std::make_pair("IMKPinching", &OPS_IMKPinching));
  uniaxialMaterialsMap.insert(
      std::make_pair("IMKPeakOriented", &OPS_IMKPeakOriented));
  uniaxialMaterialsMap.insert(
      std::make_pair("SLModel", &OPS_SLModel));
  uniaxialMaterialsMap.insert(
      std::make_pair("SMA", &OPS_SMAMaterial));
  uniaxialMaterialsMap.insert(
      std::make_pair("FRCC", &OPS_FRCC));
  uniaxialMaterialsMap.insert(
      std::make_pair("ConcreteZBH_original", &OPS_ConcreteZBH_original));
  uniaxialMaterialsMap.insert(
      std::make_pair("ConcreteZBH_fitted", &OPS_ConcreteZBH_fitted));
  uniaxialMaterialsMap.insert(
      std::make_pair("ConcreteZBH_smoothed", &OPS_ConcreteZBH_smoothed));      
  uniaxialMaterialsMap.insert(std::make_pair(
      "HystereticPoly",
      &OPS_HystereticPoly));  // Salvatore Sessa 14-Jan-2021 Mail: salvatore.sessa2@unina.it
  uniaxialMaterialsMap.insert(std::make_pair("HystereticSmooth", &OPS_HystereticSmooth)); // Salvatore Sessa 19-Apr-2022 Mail: salvatore.sessa2@unina.it
	uniaxialMaterialsMap.insert(std::make_pair("HystereticAsym", &OPS_HystereticAsym));     // Salvatore Sessa 21-Apr-2022 Mail: salvatore.sessa2@unina.it
  uniaxialMaterialsMap.insert(
      std::make_pair("DowelType", &OPS_DowelType));
  uniaxialMaterialsMap.insert(
      std::make_pair("DuctileFracture",
                     &OPS_DuctileFracture));  // Kuanshi Zhong
  uniaxialMaterialsMap.insert(
      std::make_pair("TDConcreteEXP", &OPS_TDConcreteEXP));
  uniaxialMaterialsMap.insert(
      std::make_pair("TDConcrete", &OPS_TDConcrete));
  uniaxialMaterialsMap.insert(
      std::make_pair("TDConcreteNL", &OPS_TDConcreteNL));  
  uniaxialMaterialsMap.insert(
      std::make_pair("TDConcreteMC10", &OPS_TDConcreteMC10));
  uniaxialMaterialsMap.insert(
      std::make_pair("TDConcreteMC10NL", &OPS_TDConcreteMC10NL));
  uniaxialMaterialsMap.insert(
      std::make_pair("Creep", &OPS_CreepMaterial));  
  uniaxialMaterialsMap.insert(
      std::make_pair("CoulombDamper", &OPS_CoulombDamperMaterial));
  uniaxialMaterialsMap.insert(std::make_pair(
	  "GMG_CyclicReinforcedConcrete", &OPS_GMG_CyclicReinforcedConcrete));
  uniaxialMaterialsMap.insert(
      std::make_pair("Hertzdamp", &OPS_Hertzdamp));
  uniaxialMaterialsMap.insert(
      std::make_pair("HertzDamp", &OPS_Hertzdamp));
  uniaxialMaterialsMap.insert(
      std::make_pair("JankowskiImpact", &OPS_JankowskiImpact));
  uniaxialMaterialsMap.insert(
      std::make_pair("ViscoelasticGap", &OPS_ViscoelasticGap));
  uniaxialMaterialsMap.insert(std::make_pair("Masonry", &OPS_Masonry));
  uniaxialMaterialsMap.insert(std::make_pair("Masonryt", &OPS_Masonryt));
  uniaxialMaterialsMap.insert(std::make_pair("Trilinwp", &OPS_Trilinwp));
  uniaxialMaterialsMap.insert(std::make_pair("Trilinwp2", &OPS_Trilinwp2));
  uniaxialMaterialsMap.insert(std::make_pair("Ratchet", &OPS_Ratchet));
  
  return 0;
}

static OPS_ParsingFunctionMap hystereticBackbonesMap;

static int setUpHystereticBackbones(void) {
  hystereticBackbonesMap.insert(
      std::make_pair("Arctangent", &OPS_ArctangentBackbone));
  hystereticBackbonesMap.insert(
      std::make_pair("Bilinear", &OPS_BilinearBackbone));
  hystereticBackbonesMap.insert(
      std::make_pair("Mander", &OPS_ManderBackbone));
  hystereticBackbonesMap.insert(
      std::make_pair("Multilinear", &OPS_MultilinearBackbone));
  hystereticBackbonesMap.insert(
      std::make_pair("Trilinear", &OPS_TrilinearBackbone));
  hystereticBackbonesMap.insert(
      std::make_pair("Material", &OPS_MaterialBackbone));  
  hystereticBackbonesMap.insert(std::make_pair(
      "ReeseStiffClayBelowWS", &OPS_ReeseStiffClayBelowWS));
  hystereticBackbonesMap.insert(std::make_pair(
      "ReeseStiffClayAboveWS", &OPS_ReeseStiffClayAboveWS));
  hystereticBackbonesMap.insert(
      std::make_pair("VuggyLimestone", &OPS_VuggyLimestone));
  hystereticBackbonesMap.insert(
      std::make_pair("CementedSoil", &OPS_CementedSoil));
  hystereticBackbonesMap.insert(
      std::make_pair("WeakRock", &OPS_WeakRock));
  hystereticBackbonesMap.insert(
      std::make_pair("LiquefiedSand", &OPS_LiquefiedSand));
  hystereticBackbonesMap.insert(
      std::make_pair("Raynor", &OPS_RaynorBackbone));
  hystereticBackbonesMap.insert(
      std::make_pair("ReeseSand", &OPS_ReeseSandBackbone));
  hystereticBackbonesMap.insert(
      std::make_pair("ReeseSoftClay", &OPS_ReeseSoftClayBackbone));  
  hystereticBackbonesMap.insert(
      std::make_pair("Capped", &OPS_CappedBackbone));
  hystereticBackbonesMap.insert(
      std::make_pair("LinearCapped", &OPS_LinearCappedBackbone));
  
  return 0;
}

static OPS_ParsingFunctionMap stiffnessDegradationsMap;

static int setUpStiffnessDegradations(void) {
  stiffnessDegradationsMap.insert(std::make_pair(
      "Constant", &OPS_ConstantStiffnessDegradation));
  stiffnessDegradationsMap.insert(std::make_pair(
      "Ductility", &OPS_DuctilityStiffnessDegradation));
  stiffnessDegradationsMap.insert(
      std::make_pair("Energy", &OPS_EnergyStiffnessDegradation));
  stiffnessDegradationsMap.insert(std::make_pair(
      "Pincheira", &OPS_PincheiraStiffnessDegradation));

  return 0;
}

static OPS_ParsingFunctionMap strengthDegradationsMap;

static int setUpStrengthDegradations(void) {
  strengthDegradationsMap.insert(std::make_pair(
      "Constant", &OPS_ConstantStrengthDegradation));
  strengthDegradationsMap.insert(std::make_pair(
      "Ductility", &OPS_DuctilityStrengthDegradation));
  strengthDegradationsMap.insert(
      std::make_pair("Energy", &OPS_EnergyStrengthDegradation));
  strengthDegradationsMap.insert(
      std::make_pair("ACI", &OPS_ACIStrengthDegradation));
  strengthDegradationsMap.insert(std::make_pair(
      "Petrangeli", &OPS_PetrangeliStrengthDegradation));

  return 0;
}

static OPS_ParsingFunctionMap unloadingRulesMap;

static int setUpUnloadingRules(void) {
  unloadingRulesMap.insert(
      std::make_pair("Constant", &OPS_ConstantUnloadingRule));
  unloadingRulesMap.insert(
      std::make_pair("Ductility", &OPS_TakedaUnloadingRule));
  unloadingRulesMap.insert(
      std::make_pair("Takeda", &OPS_TakedaUnloadingRule));
  unloadingRulesMap.insert(
      std::make_pair("Energy", &OPS_EnergyUnloadingRule));
  unloadingRulesMap.insert(
      std::make_pair("Karsan", &OPS_KarsanUnloadingRule));

  return 0;
}

}  // namespace

int OPS_UniaxialMaterial() {
  static bool initDone = false;
  if (initDone == false) {
    setUpUniaxialMaterials();
    initDone = true;
  }

  if (OPS_GetNumRemainingInputArgs() < 2) {
    opserr << "WARNING too few arguments: uniaxialMaterial "
              "type? tag? ...\n";
    return -1;
  }

  const char* matType = OPS_GetString();

  OPS_ParsingFunctionMap::const_iterator iter =
      uniaxialMaterialsMap.find(matType);
  if (iter == uniaxialMaterialsMap.end()) {
    opserr << "WARNING material type " << matType
           << " is unknown\n";
    return -1;
  }

  UniaxialMaterial* theMaterial =
      (UniaxialMaterial*)(*iter->second)();
  if (theMaterial == 0) {
    return -1;
  }

  // Now add the material to the modelBuilder
  if (OPS_addUniaxialMaterial(theMaterial) == false) {
    opserr << "ERROR could not add uniaaialMaterial.\n";
    delete theMaterial;  // invoke the material objects
                         // destructor, otherwise mem leak
    return -1;
  }

  return 0;
}

int OPS_testUniaxialMaterial() {
  if (OPS_GetNumRemainingInputArgs() != 1) {
    opserr << "testUniaxialMaterial - You must provide a "
              "material tag.\n";
    return -1;
  }

  int tag;
  int numData = 1;
  if (OPS_GetIntInput(&numData, &tag) < 0) {
    opserr << "invalid int value\n";
    return -1;
  }

  UniaxialMaterial* mat = OPS_getUniaxialMaterial(tag);

  if (mat == 0) {
    opserr << "testUniaxialMaterial - Material Not Found.\n";
    return -1;
  }

  theTestingUniaxialMaterial = mat;

  return 0;
}

int OPS_setStrain() {
  if (OPS_GetNumRemainingInputArgs() < 1) {
    opserr << "testUniaxialMaterial - You must provide a strain "
              "value.\n";
    return -1;
  }

  UniaxialMaterial* material = theTestingUniaxialMaterial;

  if (material == 0) {
    opserr << "setStrain WARNING no active UniaxialMaterial - "
              "use testUniaxialMaterial command.\n";
    return -1;
  }

  double strain;
  int numData = 1;
  if (OPS_GetDoubleInput(&numData, &strain) < 0) {
    opserr << "invalid double value\n";
    return -1;
  }

  double strainRate = 0.0;
  if (OPS_GetNumRemainingInputArgs() > 0) {
      if (OPS_GetDoubleInput(&numData, &strainRate) < 0) {
          opserr << "invalid strain rate\n";
          return -1;
      }
  }

  material->setTrialStrain(strain, strainRate);
  material->commitState();

  return 0;
}

int OPS_getStrain() {
  UniaxialMaterial* material = theTestingUniaxialMaterial;
  if (material == 0) {
    opserr << "getStrain WARNING no active UniaxialMaterial - "
              "use testUniaxialMaterial command.\n";
    return -1;
  }

  double strain = material->getStrain();

  int numData = 1;

  if (OPS_SetDoubleOutput(&numData, &strain, true) < 0) {
    opserr << "failed to set strain\n";
    return -1;
  }

  return 0;
}

int OPS_getStress() {
  UniaxialMaterial* material = theTestingUniaxialMaterial;
  if (material == 0) {
    opserr << "getStrain WARNING no active UniaxialMaterial - "
              "use testUniaxialMaterial command.\n";
    return -1;
  }

  double stress = material->getStress();

  int numData = 1;

  if (OPS_SetDoubleOutput(&numData, &stress, true) < 0) {
    opserr << "failed to set stress\n";
    return -1;
  }

  return 0;
}

int OPS_getTangent() {
  UniaxialMaterial* material = theTestingUniaxialMaterial;
  if (material == 0) {
    opserr << "getStrain WARNING no active UniaxialMaterial - "
              "use testUniaxialMaterial command.\n";
    return -1;
  }

  double tangent = material->getTangent();

  int numData = 1;

  if (OPS_SetDoubleOutput(&numData, &tangent, true) < 0) {
    opserr << "failed to set tangent\n";
    return -1;
  }

  return 0;
}

int OPS_getDampTangent() {
  UniaxialMaterial* material = theTestingUniaxialMaterial;
  if (material == 0) {
    opserr << "getStrain WARNING no active UniaxialMaterial - "
              "use testUniaxialMaterial command.\n";
    return -1;
  }

  double tangent = material->getDampTangent();

  int numData = 1;

  if (OPS_SetDoubleOutput(&numData, &tangent, true) < 0) {
    opserr << "failed to set damp tangent\n";
    return -1;
  }

  return 0;
}

void* OPS_RotationShearCurve();
void* OPS_ThreePointCurve();
void* OPS_ShearCurve();
void* OPS_AxialCurve();

int OPS_LimitCurve() {
  // Make sure there is a minimum number of arguments
  if (OPS_GetNumRemainingInputArgs() < 6) {
    opserr << "WARNING insufficient number of limit curve "
              "arguments\n";
    opserr << "Want: limitCurve type? tag? <specific curve args>"
           << endln;
    return -1;
  }

  const char* type = OPS_GetString();

  // Pointer to a limit curve that will be added to the model
  // builder
  LimitCurve* theCurve = 0;

  if (strcmp(type, "Axial") == 0) {
    // opserr << "WARNING to be implemented ...\n";
    void* curve = OPS_AxialCurve();
    if (curve != 0) {
      theCurve = (LimitCurve*)curve;
    } else {
      return -1;
    }
    return -1;

    return -1;

  } else if (strcmp(type, "RotationShearCurve") == 0) {
    void* theRSC = OPS_RotationShearCurve();
    if (theRSC != 0) {
      theCurve = (LimitCurve*)theRSC;
    } else {
      return -1;
    }

  } else if (strcmp(type, "ThreePoint") == 0) {
    void* curve = OPS_ThreePointCurve();
    if (curve != 0) {
      theCurve = (LimitCurve*)curve;
    } else {
      return -1;
    }

  } else if (strcmp(type, "Shear") == 0) {
    void* curve = OPS_ShearCurve();
    if (curve != 0) {
      theCurve = (LimitCurve*)curve;
    } else {
      return -1;
    }

  } else {
    opserr << "WARNING type of limit curve is unknown\n";
    return -1;
  }

  // Ensure we have created the Material, out of memory if got
  // here and no material
  if (theCurve == 0) {
    opserr << "WARNING ran out of memory creating limitCurve\n";
    return -1;
  }

  // Now add the material to the modelBuilder
  if (OPS_addLimitCurve(theCurve) == false) {
    opserr << "WARNING could not add limitCurve to the domain\n";
    delete theCurve;  // invoke the material objects destructor,
                      // otherwise mem leak
    return -1;
  }

  return 0;
}

int OPS_hystereticBackbone() {
  static bool initDone = false;
  if (initDone == false) {
    setUpHystereticBackbones();
    initDone = true;
  }

  if (OPS_GetNumRemainingInputArgs() < 2) {
    opserr << "WARNING too few arguments: hystereticBackbone "
              "type? tag? ...\n";
    return -1;
  }

  const char* matType = OPS_GetString();

  OPS_ParsingFunctionMap::const_iterator iter =
      hystereticBackbonesMap.find(matType);
  if (iter == hystereticBackbonesMap.end()) {
    opserr << "WARNING hystereticBackbone type " << matType
           << " is unknown\n";
    return -1;
  }

  HystereticBackbone* theBackbone =
      (HystereticBackbone*)(*iter->second)();
  if (theBackbone == 0) {
    return -1;
  }

  // Now add the material to the modelBuilder
  if (OPS_addHystereticBackbone(theBackbone) == false) {
    opserr << "ERROR could not add HystereticBackbone\n";
    delete theBackbone;
    return -1;
  }

  return 0;
}

int OPS_stiffnessDegradation() {
  static bool initDone = false;
  if (initDone == false) {
    setUpStiffnessDegradations();
    initDone = true;
  }

  if (OPS_GetNumRemainingInputArgs() < 2) {
    opserr << "WARNING too few arguments: stiffnessDegradation "
              "type? tag? ...\n";
    return -1;
  }

  const char* matType = OPS_GetString();

  OPS_ParsingFunctionMap::const_iterator iter =
      stiffnessDegradationsMap.find(matType);
  if (iter == stiffnessDegradationsMap.end()) {
    opserr << "WARNING stiffnessDegradation type " << matType
           << " is unknown\n";
    return -1;
  }

  StiffnessDegradation* theBackbone =
      (StiffnessDegradation*)(*iter->second)();
  if (theBackbone == 0) {
    return -1;
  }

  // Now add the material to the modelBuilder
  if (OPS_addStiffnessDegradation(theBackbone) == false) {
    opserr << "ERROR could not add StiffnessDegradation\n";
    delete theBackbone;
    return -1;
  }

  return 0;
}

int OPS_strengthDegradation() {
  static bool initDone = false;
  if (initDone == false) {
    setUpStrengthDegradations();
    initDone = true;
  }

  if (OPS_GetNumRemainingInputArgs() < 2) {
    opserr << "WARNING too few arguments: strengthDegradation "
              "type? tag? ...\n";
    return -1;
  }

  const char* matType = OPS_GetString();

  OPS_ParsingFunctionMap::const_iterator iter =
      strengthDegradationsMap.find(matType);
  if (iter == strengthDegradationsMap.end()) {
    opserr << "WARNING strengthDegradation type " << matType
           << " is unknown\n";
    return -1;
  }

  StrengthDegradation* theBackbone =
      (StrengthDegradation*)(*iter->second)();
  if (theBackbone == 0) {
    return -1;
  }

  // Now add the material to the modelBuilder
  if (OPS_addStrengthDegradation(theBackbone) == false) {
    opserr << "ERROR could not add StrengthDegradation\n";
    delete theBackbone;
    return -1;
  }

  return 0;
}

int OPS_unloadingRule() {
  static bool initDone = false;
  if (initDone == false) {
    setUpUnloadingRules();
    initDone = true;
  }

  if (OPS_GetNumRemainingInputArgs() < 2) {
    opserr << "WARNING too few arguments: unloadingRule type? "
              "tag? ...\n";
    return -1;
  }

  const char* matType = OPS_GetString();

  OPS_ParsingFunctionMap::const_iterator iter =
      unloadingRulesMap.find(matType);
  if (iter == unloadingRulesMap.end()) {
    opserr << "WARNING unloadingRule type " << matType
           << " is unknown\n";
    return -1;
  }

  UnloadingRule* theBackbone = (UnloadingRule*)(*iter->second)();
  if (theBackbone == 0) {
    return -1;
  }

  // Now add the material to the modelBuilder
  if (OPS_addUnloadingRule(theBackbone) == false) {
    opserr << "ERROR could not add UnloadingRule\n";
    delete theBackbone;
    return -1;
  }

  return 0;
}
