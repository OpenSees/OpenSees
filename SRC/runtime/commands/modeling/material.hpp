#include <string>
#include <unordered_map>
#include <elementAPI.h>

extern OPS_Routine OPS_ASDConcrete3DMaterial;
extern OPS_Routine OPS_ReinforcedConcretePlaneStressMaterial;
extern OPS_Routine OPS_FAReinforcedConcretePlaneStressMaterial;
extern OPS_Routine OPS_FAFourSteelRCPlaneStressMaterial;
extern OPS_Routine OPS_RAFourSteelRCPlaneStressMaterial;
extern OPS_Routine OPS_PrestressedConcretePlaneStressMaterial;
extern OPS_Routine OPS_FAPrestressedConcretePlaneStressMaterial;
extern OPS_Routine OPS_FAFourSteelPCPlaneStressMaterial;
extern OPS_Routine OPS_RAFourSteelPCPlaneStressMaterial;
// extern OPS_Routine OPS_MaterialCMM;
// extern OPS_Routine OPS_NewMaterialCMM;
extern OPS_Routine OPS_NewPlasticDamageConcrete3d;
extern OPS_Routine OPS_NewPlasticDamageConcretePlaneStress;
extern OPS_Routine OPS_ElasticIsotropicMaterial;
extern OPS_Routine OPS_ElasticIsotropic3D;
extern OPS_Routine OPS_IncrementalElasticIsotropicThreeDimensional;
extern OPS_Routine OPS_ElasticOrthotropicMaterial;
extern OPS_Routine OPS_DruckerPragerMaterial;
extern OPS_Routine OPS_BoundingCamClayMaterial;
extern OPS_Routine OPS_ContactMaterial2DMaterial;
extern OPS_Routine OPS_ContactMaterial3DMaterial;
extern OPS_Routine OPS_InitialStateAnalysisWrapperMaterial;
extern OPS_Routine OPS_ManzariDafaliasMaterial;
extern OPS_Routine OPS_ManzariDafaliasMaterialRO;
extern OPS_Routine OPS_PM4SandMaterial;
extern OPS_Routine OPS_PM4SiltMaterial;
extern OPS_Routine OPS_J2CyclicBoundingSurfaceMaterial;
extern OPS_Routine OPS_CycLiqCPMaterial;
extern OPS_Routine OPS_CycLiqCPSPMaterial;
extern OPS_Routine OPS_InitStressNDMaterial;
extern OPS_Routine OPS_InitStrainNDMaterial;
extern OPS_Routine OPS_StressDensityMaterial;
extern OPS_Routine OPS_J2PlateFibreMaterial;
extern OPS_Routine OPS_PlaneStressLayeredMaterial;
extern OPS_Routine OPS_PlaneStressRebarMaterial;
extern OPS_Routine OPS_PlateFiberMaterial;
extern OPS_Routine OPS_BeamFiberMaterial;
extern OPS_Routine OPS_BeamFiberMaterial2d;
extern OPS_Routine OPS_BeamFiberMaterial2dPS;
extern OPS_Routine OPS_LinearCap;
extern OPS_Routine OPS_AcousticMedium;
extern OPS_Routine OPS_UVCmultiaxial;
extern OPS_Routine OPS_UVCplanestress;
extern OPS_Routine OPS_SAniSandMSMaterial;
extern OPS_Routine OPS_OrthotropicRotatingAngleConcreteT2DMaterial01;	// M. J. Nunez - UChile
extern OPS_Routine OPS_SmearedSteelDoubleLayerT2DMaterial01;		// M. J. Nunez - UChile

extern OPS_Routine OPS_ElasticIsotropicMaterialThermal;           // L.Jiang [SIF]
extern OPS_Routine OPS_DruckerPragerMaterialThermal;              // L.Jiang [SIF]
extern OPS_Routine OPS_PlasticDamageConcretePlaneStressThermal;   // L.Jiang [SIF]

extern OPS_Routine OPS_ElasticOrthotropicPlaneStress;
extern OPS_Routine OPS_OrthotropicMaterial;
extern OPS_Routine OPS_Series3DMaterial;
extern OPS_Routine OPS_Parallel3DMaterial;

extern OPS_Routine OPS_AllASDPlasticMaterials;

#ifdef _HAVE_Faria1998
extern OPS_Routine OPS_NewFaria1998Material;
extern OPS_Routine OPS_NewConcreteMaterial;
#endif

extern OPS_Routine OPS_FSAMMaterial; // K Kolozvari

#ifdef _HAVE_Damage2p
extern OPS_Routine OPS_Damage2p;
#endif

static std::unordered_map<std::string, OPS_Routine*> material_dispatch = {
  {"InitStressMaterial",            OPS_InitStressNDMaterial},
  {"InitStrainMaterial",            OPS_InitStrainNDMaterial},
  {"InitStrain",                    OPS_InitStrainNDMaterial},

  {"ReinforcedConcretePlaneStress", OPS_ReinforcedConcretePlaneStressMaterial},
  {"PlaneStressLayeredMaterial",    OPS_PlaneStressLayeredMaterial},
  {"PlaneStressRebarMaterial",      OPS_PlaneStressRebarMaterial},

#ifdef OPS_USE_ASDPlasticMaterials
  {"ASDPlasticMaterial",            OPS_AllASDPlasticMaterials},
#endif

  {"ASDConcrete3D",                 OPS_ASDConcrete3DMaterial},

  {"PlasticDamageConcrete",         OPS_NewPlasticDamageConcrete3d},

  {"PlasticDamageConcretePlaneStress", OPS_NewPlasticDamageConcretePlaneStress},


  {"J2PlateFibre", OPS_J2PlateFibreMaterial}, 
  {"PlateFiber",   OPS_PlateFiberMaterial},

  // Beam fiber
  {"BeamFiber",    OPS_BeamFiberMaterial},
  {"BeamFiber2d", OPS_BeamFiberMaterial2d},
  {"BeamFiber2dPS", OPS_BeamFiberMaterial2dPS},

  
  {"DruckerPragerThermal", OPS_DruckerPragerMaterialThermal},
#if 0
  {"CDPPlaneStressThermal", OPS_PlasticDamageConcretePlaneStressThermal},
#endif


#ifdef _HAVE_Faria1998
  {"Faria1998", OPS_NewFaria1998Material},  
  {"Concrete", OPS_NewConcreteMaterial},
#endif

  {"FAReinforcedConcretePlaneStress", OPS_FAReinforcedConcretePlaneStressMaterial},
  {"RAFourSteelRCPlaneStress",        OPS_RAFourSteelRCPlaneStressMaterial},
  {"FAFourSteelRCPlaneStress",        OPS_FAFourSteelRCPlaneStressMaterial},

#ifdef _HAVE_Damage2p
  {"Damage2p",                        OPS_Damage2p},
#endif

  {"PrestressedConcretePlaneStress",   OPS_PrestressedConcretePlaneStressMaterial},
  {"FAPrestressedConcretePlaneStress", OPS_FAPrestressedConcretePlaneStressMaterial},
  {"RAFourSteetPCPlaneStress",         OPS_RAFourSteelPCPlaneStressMaterial},

  {"FAFourSteelPCPlaneStress",         OPS_FAFourSteelPCPlaneStressMaterial},

  {"DruckerPrager",  OPS_DruckerPragerMaterial},

  {"TruncatedDP",    OPS_LinearCap},

  // K Kolozvari
  {"FSAM",           OPS_FSAMMaterial},

  {"AcousticMedium", OPS_AcousticMedium},

  {"UVCplanestress", OPS_UVCplanestress},

  {"UVCmultiaxial",  OPS_UVCmultiaxial},

//{"MaterialCMM",    OPS_MaterialCMM},

  {"CycLiqCP",        OPS_CycLiqCPMaterial},

  {"CycLiqCPSP",      OPS_CycLiqCPSPMaterial},

  {"BoundingCamClay", OPS_BoundingCamClayMaterial},

  {"ManzariDafalias", OPS_ManzariDafaliasMaterial},

  {"ManzariDafaliasRO", OPS_ManzariDafaliasMaterialRO},

  {"PM4Sand", OPS_PM4SandMaterial},

  {"J2CyclicBoundingSurface", OPS_J2CyclicBoundingSurfaceMaterial},

  {"PM4Silt", OPS_PM4SiltMaterial},

  {"ContactMaterial2D", OPS_ContactMaterial2DMaterial},

  {"ContactMaterial3D", OPS_ContactMaterial3DMaterial},

  {"InitialStateAnalysisWrapper", OPS_InitialStateAnalysisWrapperMaterial},

  {"stressDensity", OPS_StressDensityMaterial},

  {"ElasticIsotropic3D", OPS_ElasticIsotropic3D},

  {"ElasticIsotropic",   OPS_ElasticIsotropicMaterial},

  {"ElasticOrthotropic", OPS_ElasticOrthotropicMaterial},

  {"ElasticIsotropic3DThermal", OPS_ElasticIsotropicMaterialThermal},

  {"IncrementalElasticIsotropic3D", OPS_IncrementalElasticIsotropicThreeDimensional},

  {"OrthotropicRAConcrete",         OPS_OrthotropicRotatingAngleConcreteT2DMaterial01},
  {"SmearedSteelDoubleLayer",       OPS_SmearedSteelDoubleLayerT2DMaterial01},

  {"SAniSandMS",                    OPS_SAniSandMSMaterial},

  {"OrthotropicMaterial",           OPS_OrthotropicMaterial},
  {"ElasticOrthotropicPlaneStress", OPS_ElasticOrthotropicPlaneStress},
  {"Series3DMaterial",              OPS_Series3DMaterial},
  {"Parallel3DMaterial",            OPS_Parallel3DMaterial},
  {"Parallel3D",                    OPS_Parallel3DMaterial},
};

