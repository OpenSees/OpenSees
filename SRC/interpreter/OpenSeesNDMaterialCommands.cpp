                                                                        
// Description: command to create nD material

#include <NDMaterial.h>
#include <elementAPI.h>
#include <map>

void* OPS_ElasticIsotropicMaterial();
void* OPS_PlateFiberMaterial();
void* OPS_ReinforcedConcretePlaneStressMaterial();
void* OPS_InitStressNDMaterial();
void* OPS_J2BeamFiber2dMaterial();
void* OPS_J2BeamFiber3dMaterial();
void* OPS_J2PlateFibreMaterial();
void* OPS_FAReinforcedConcretePlaneStressMaterial();
void* OPS_RAFourSteelRCPlaneStressMaterial();
void* OPS_FAFourSteelRCPlaneStressMaterial();
void* OPS_Damage2p();
void* OPS_PrestressedConcretePlaneStressMaterial();
void* OPS_FAPrestressedConcretePlaneStressMaterial();
void* OPS_RAFourSteelPCPlaneStressMaterial();
void* OPS_FAFourSteelPCPlaneStressMaterial();
void* OPS_DruckerPragerMaterial();
void* OPS_LinearCap();
void* OPS_FSAMMaterial();
void* OPS_AcousticMedium();
void* OPS_MaterialCMM();
void* OPS_CycLiqCPMaterial();
void* OPS_CycLiqCPSPMaterial();
void* OPS_BoundingCamClayMaterial();
void* OPS_ManzariDafaliasMaterial();
void* OPS_ContactMaterial2DMaterial();
void* OPS_ContactMaterial3DMaterial();
void* OPS_InitialStateAnalysisWrapperMaterial();
void* OPS_StressDensityMaterial();
void* OPS_ElasticOrthotropicMaterial();
void* OPS_PressureDependentElastic3D();
void* OPS_J2Plasticity();
void* OPS_PlaneStressSimplifiedJ2();
void* OPS_MultiaxialCyclicPlasticity();
void* OPS_PressureIndependMultiYield();
void* OPS_MultiYieldSurfaceClay();
void* OPS_PressureDependMultiYield();
void* OPS_PressureDependMultiYield02();
void* OPS_FluidSolidPorousMaterial();
void* OPS_PlaneStress();
void* OPS_PlaneStrain();
void* OPS_CapPlasticity();
void *OPS_SimplifiedJ2();
void* OPS_PlateRebarMaterial();
void* OPS_PlateFromPlaneStressMaterial();
void* OPS_ConcreteS();
void* OPS_PlaneStressUserMaterial();
void* OPS_BeamFiberMaterial();
void* OPS_PM4Sand();

namespace {

    struct char_cmp { 
	bool operator () (const char *a,const char *b) const 
	    {
		return strcmp(a,b)<0;
	    } 
    };

    typedef std::map<const char *, void *(*)(void), char_cmp> OPS_ParsingFunctionMap;


    static OPS_ParsingFunctionMap nDMaterialsMap;

  static void* J2BeamFiber2Dor3D (void)
  {
    int NDM = OPS_GetNDM();
    if (NDM == 2)
      return OPS_J2BeamFiber2dMaterial();
    if (NDM == 3)
      return OPS_J2BeamFiber3dMaterial();
    return 0;
  }

    static int setUpNDMaterials(void)
    {
	nDMaterialsMap.insert(std::make_pair("ReinforcedConcretePlaneStress", &OPS_ReinforcedConcretePlaneStressMaterial));
	nDMaterialsMap.insert(std::make_pair("ReinforceConcretePlaneStress", &OPS_ReinforcedConcretePlaneStressMaterial));
	nDMaterialsMap.insert(std::make_pair("InitStressNDMaterial", &OPS_InitStressNDMaterial));
	nDMaterialsMap.insert(std::make_pair("InitStressND", &OPS_InitStressNDMaterial));
	nDMaterialsMap.insert(std::make_pair("J2BeamFiber", &J2BeamFiber2Dor3D));
	nDMaterialsMap.insert(std::make_pair("J2PlateFibre", &OPS_J2PlateFibreMaterial));
	nDMaterialsMap.insert(std::make_pair("FAReinforcedConcretePlaneStress", &OPS_FAReinforcedConcretePlaneStressMaterial));
	nDMaterialsMap.insert(std::make_pair("FAReinforceConcretePlaneStress", &OPS_FAReinforcedConcretePlaneStressMaterial));
	nDMaterialsMap.insert(std::make_pair("RAFourSteelRCPlaneStress", &OPS_RAFourSteelRCPlaneStressMaterial));
	nDMaterialsMap.insert(std::make_pair("FAFourSteelRCPlaneStress", &OPS_FAFourSteelRCPlaneStressMaterial));
#ifdef _HAVE_Damage2P
	nDMaterialsMap.insert(std::make_pair("Damage2p", &OPS_Damage2p));
#endif
	nDMaterialsMap.insert(std::make_pair("PrestressedConcretePlaneStress", &OPS_PrestressedConcretePlaneStressMaterial));
	nDMaterialsMap.insert(std::make_pair("FAPrestressedConcretePlaneStress", &OPS_FAPrestressedConcretePlaneStressMaterial));
	nDMaterialsMap.insert(std::make_pair("RAFourSteetPCPlaneStress", &OPS_RAFourSteelPCPlaneStressMaterial));
	nDMaterialsMap.insert(std::make_pair("FAFourSteelPCPlaneStress", &OPS_FAFourSteelPCPlaneStressMaterial));
	nDMaterialsMap.insert(std::make_pair("DruckerPrager", &OPS_DruckerPragerMaterial));
	nDMaterialsMap.insert(std::make_pair("TruncatedDP", &OPS_LinearCap));
	nDMaterialsMap.insert(std::make_pair("FSAM", &OPS_FSAMMaterial));
	nDMaterialsMap.insert(std::make_pair("AcousticMedium", &OPS_AcousticMedium));
	nDMaterialsMap.insert(std::make_pair("MaterialCMM", &OPS_MaterialCMM));
	nDMaterialsMap.insert(std::make_pair("CycLiqCP", &OPS_CycLiqCPMaterial));
	nDMaterialsMap.insert(std::make_pair("CycLiqCPSP", &OPS_CycLiqCPSPMaterial));
	nDMaterialsMap.insert(std::make_pair("BoundingCamClay", &OPS_BoundingCamClayMaterial));
	nDMaterialsMap.insert(std::make_pair("ManzariDafalias", &OPS_ManzariDafaliasMaterial));
	nDMaterialsMap.insert(std::make_pair("ContactMaterial2D", &OPS_ContactMaterial2DMaterial));
	nDMaterialsMap.insert(std::make_pair("ContactMaterial3D", &OPS_ContactMaterial3DMaterial));
	nDMaterialsMap.insert(std::make_pair("InitialStateAnalysisWrapper", &OPS_InitialStateAnalysisWrapperMaterial));
	nDMaterialsMap.insert(std::make_pair("StressDensityModel", &OPS_StressDensityMaterial));
	nDMaterialsMap.insert(std::make_pair("ElasticIsotropic", &OPS_ElasticIsotropicMaterial));
	nDMaterialsMap.insert(std::make_pair("ElasticIsotropic3D", &OPS_ElasticIsotropicMaterial));
	nDMaterialsMap.insert(std::make_pair("ElasticOrthotropic3D", &OPS_ElasticOrthotropicMaterial));
	nDMaterialsMap.insert(std::make_pair("ElasticOrthotropic", &OPS_ElasticOrthotropicMaterial));
	nDMaterialsMap.insert(std::make_pair("PressureDependentElastic3D", &OPS_PressureDependentElastic3D));
	nDMaterialsMap.insert(std::make_pair("J2Plasticity", &OPS_J2Plasticity));
	nDMaterialsMap.insert(std::make_pair("J2", &OPS_J2Plasticity));
	nDMaterialsMap.insert(std::make_pair("PlaneStressSimplifiedJ2", &OPS_PlaneStressSimplifiedJ2));
	nDMaterialsMap.insert(std::make_pair("MultiaxialCyclicPlasticity", &OPS_MultiaxialCyclicPlasticity));
	nDMaterialsMap.insert(std::make_pair("MCP", &OPS_MultiaxialCyclicPlasticity));
	nDMaterialsMap.insert(std::make_pair("PressureIndependMultiYield", &OPS_PressureIndependMultiYield));
	nDMaterialsMap.insert(std::make_pair("MultiYieldSurfaceClay", &OPS_MultiYieldSurfaceClay));
	nDMaterialsMap.insert(std::make_pair("PressureDependMultiYield", &OPS_PressureDependMultiYield));
	nDMaterialsMap.insert(std::make_pair("PressureDependMultiYield02", &OPS_PressureDependMultiYield02));
	nDMaterialsMap.insert(std::make_pair("FluidSolidPorous", &OPS_FluidSolidPorousMaterial));
	nDMaterialsMap.insert(std::make_pair("PlaneStressMaterial", &OPS_PlaneStress));
	nDMaterialsMap.insert(std::make_pair("PlaneStress", &OPS_PlaneStress));
	nDMaterialsMap.insert(std::make_pair("PlaneStrainMaterial", &OPS_PlaneStrain));
	nDMaterialsMap.insert(std::make_pair("PlaneStrain", &OPS_PlaneStrain));
	nDMaterialsMap.insert(std::make_pair("PlateFiber", &OPS_PlateFiberMaterial));
	nDMaterialsMap.insert(std::make_pair("PlateFiberMaterial", &OPS_PlateFiberMaterial));
	nDMaterialsMap.insert(std::make_pair("CapPlasticity", &OPS_CapPlasticity));
	nDMaterialsMap.insert(std::make_pair("Simplified3DJ2", &OPS_SimplifiedJ2));
	nDMaterialsMap.insert(std::make_pair("3DJ2", &OPS_SimplifiedJ2));
	nDMaterialsMap.insert(std::make_pair("PlateRebarMaterial", &OPS_PlateRebarMaterial));
	nDMaterialsMap.insert(std::make_pair("PlateRebar", &OPS_PlateRebarMaterial));
	nDMaterialsMap.insert(std::make_pair("PlateFromPlaneStressMaterial", &OPS_PlateFromPlaneStressMaterial));
	nDMaterialsMap.insert(std::make_pair("PlateFromPlaneStress", &OPS_PlateFromPlaneStressMaterial));
	nDMaterialsMap.insert(std::make_pair("ConcreteS", &OPS_ConcreteS));
	nDMaterialsMap.insert(std::make_pair("PlaneStressUserMaterial", &OPS_PlaneStressUserMaterial));
	nDMaterialsMap.insert(std::make_pair("BeamFiberMaterial", &OPS_BeamFiberMaterial));
	nDMaterialsMap.insert(std::make_pair("BeamFiber", &OPS_BeamFiberMaterial));
	nDMaterialsMap.insert(std::make_pair("PM4Sand", &OPS_BeamFiberMaterial));
	
	return 0;
    }
}

int
OPS_NDMaterial()
{
    static bool initDone = false;
    if (initDone == false) {
	setUpNDMaterials();
	initDone = true;
    }

    if (OPS_GetNumRemainingInputArgs() < 2) {
	opserr<<"WARNING too few arguments: nDMaterial type? tag? ...\n";
	return -1;
    }

    const char* matType = OPS_GetString();
    
    OPS_ParsingFunctionMap::const_iterator iter = nDMaterialsMap.find(matType);
    if (iter == nDMaterialsMap.end()) {
	opserr<<"WARNING material type " << matType << " is unknown\n";
	return -1;
    }

    NDMaterial* theMaterial = (NDMaterial*) (*iter->second)();
    if (theMaterial == 0) {
	return -1;
    }

    // Now add the material to the modelBuilder
    if (OPS_addNDMaterial(theMaterial) == false) {
	opserr<<"ERROR could not add NDMaterial.\n";
	delete theMaterial;
	return -1;
    }

    return 0;

}
