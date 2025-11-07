//===----------------------------------------------------------------------===//
//
//                                   xara
//                              https://xara.so
//
//===----------------------------------------------------------------------===//
//
// Copyright (c) 2025, OpenSees/Xara Developers
// All rights reserved.  No warranty, explicit or implicit, is provided.
//
// This source code is licensed under the BSD 2-Clause License.
// See LICENSE file or https://opensource.org/licenses/BSD-2-Clause
//
//===----------------------------------------------------------------------===//
//
//
//
//                 stress         |          strain
//          1   2   3   4   5   6 |  
//         xx  yy  zz  xy  yz  xz |  11  22  33  12  23  31
//    PSn:  1   2       3         |   1   2   0   3   0   0
//    PSe:  1   2   0   3   -   ? |   1   2  [1]  3  [2] [3]
//    PF :  1   2   -   3   4   5 |   1   2  [1]  3   4   5
//    BF :  1   0   0   2   0   3 |   1  [1] [2]  2  [3]  3
//    AS?:  1   2   3   4   -   
//
// strains ordered  00, 11, 22, 01  
//            i.e.  11, 22, 33, 12 
//
//            strain(0) =   eps_00
//            strain(1) =   eps_11
//            strain(2) =   eps_22
//            strain(3) = 2*eps_01
//
//  same ordering for stresses but no 2 
//
//
// strains ordered : eps11, eps22, eps33, 2*eps12, 2*eps23, 2*eps31 
// NDmaterial  strain order       = 11, 22, 33, 12, 23, 31 
// PlaneStress strain order       = 11, 22, 12, 33, 23, 31
// BeamFiber   strain order       = 11, 12, 31, 22, 33, 23
// PlateFiber strain order        = 11, 22, 12, 23, 31, 33
//                                   0   1   2   3   4   5

// Platefiber: 22, 33, 13, and 23 are condensed out.

// PlateFiberMaterial strain order =  11, 22, 12, 23, 31, 33

//
//      0  1  2  3  4  5
// ND : 11 22 33   12   23   31
// PS : 11 22 12   33   23   31 
// PF : 11 22 12   23   31   33 | 
// BF : 11 12 13 | 22   33   23 
// AS : 11 22 33 12
//
#include <tcl.h>
#include <string>
#include <assert.h>
#include <Parsing.h>
#include <unordered_map>
#include <elementAPI.h>
#include <runtimeAPI.h>
#include <NDMaterial.h>

extern Tcl_CmdProc TclCommand_addPlaneWrapper;
extern Tcl_CmdProc TclCommand_newJ2Material;
extern Tcl_CmdProc TclCommand_newJ2Simplified;
extern Tcl_CmdProc TclCommand_newPlasticMaterial;
extern Tcl_CmdProc TclCommand_newConcreteMaterial;
extern Tcl_CmdProc TclCommand_newElasticMaterial;
extern Tcl_CmdProc TclCommand_addWrappingMaterial;
extern Tcl_CmdProc TclCommand_newPlateRebar;
extern Tcl_CmdProc TclCommand_newPlateFiber;

extern OPS_Routine OPS_ElasticOrthotropicPlaneStress;
extern OPS_Routine OPS_OrthotropicMaterial;
extern OPS_Routine OPS_Series3DMaterial;
extern OPS_Routine OPS_Parallel3DMaterial;
extern OPS_Routine OPS_J2PlateFibreMaterial;
extern OPS_Routine OPS_J2CyclicBoundingSurfaceMaterial;
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
extern OPS_Routine OPS_NewPlasticDamageConcretePlaneStress;
extern OPS_Routine OPS_ElasticIsotropicMaterial;
extern OPS_Routine OPS_IncrementalElasticIsotropicThreeDimensional;
extern OPS_Routine OPS_ElasticOrthotropicMaterial;
extern OPS_Routine OPS_BoundingCamClayMaterial;
extern OPS_Routine OPS_ContactMaterial2DMaterial;
extern OPS_Routine OPS_ContactMaterial3DMaterial;
extern OPS_Routine OPS_InitialStateAnalysisWrapperMaterial;
extern OPS_Routine OPS_ManzariDafaliasMaterial;
extern OPS_Routine OPS_ManzariDafaliasMaterialRO;
extern OPS_Routine OPS_PM4SandMaterial;
extern OPS_Routine OPS_PM4SiltMaterial;
extern OPS_Routine OPS_LinearElasticGGmaxMaterial;
extern OPS_Routine OPS_CycLiqCPMaterial;
extern OPS_Routine OPS_CycLiqCPSPMaterial;
extern OPS_Routine OPS_InitStressNDMaterial;
extern OPS_Routine OPS_StressDensityMaterial;
extern OPS_Routine OPS_PlaneStressLayeredMaterial;
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

extern OPS_Routine OPS_AllASDPlasticMaterials;

#ifdef _HAVE_Faria1998
extern OPS_Routine OPS_NewFaria1998Material;
extern OPS_Routine OPS_NewConcreteMaterial;
#endif

extern OPS_Routine OPS_FSAMMaterial; // K Kolozvari

#ifdef _HAVE_Damage2p
extern OPS_Routine OPS_Damage2p;
#endif


extern "C"
int OPS_ResetInputNoBuilder(ClientData, Tcl_Interp *, int cArg,
                            int mArg, TCL_Char ** const argv, Domain *);


template <OPS_Routine fn> static int
dispatch(ClientData clientData, Tcl_Interp* interp, int argc, G3_Char** const argv)
{
  BasicModelBuilder *builder = static_cast<BasicModelBuilder*>(clientData);

  OPS_ResetInputNoBuilder(clientData, interp, 2, argc, argv, 0);

  G3_Runtime *rt = G3_getRuntime(interp);
  NDMaterial* theMaterial = (NDMaterial*)fn( rt, argc, argv );
  if (theMaterial == nullptr) {
    return TCL_ERROR;
  }

  if (builder->addTaggedObject<NDMaterial>(*theMaterial) != TCL_OK) {
    opserr << OpenSees::PromptValueError << "Failed to add material to the model builder.\n";
    delete theMaterial;
    return TCL_ERROR;
  }
  return TCL_OK;
}

template <int (*fn)(ClientData clientData, Tcl_Interp* interp, int, G3_Char** const)> 
static int
dispatch(ClientData clientData, Tcl_Interp* interp, int argc, G3_Char** const argv)
{
  assert(clientData != nullptr);
  return fn( clientData, interp, argc, argv );
}

namespace OpenSees {
static std::unordered_map<std::string, Tcl_CmdProc*> MaterialLibrary = {
//
// Elastic 
//
// Isotropic
  {"ElasticIsotropic3D",               dispatch<TclCommand_newElasticMaterial>},
  {"ElasticIsotropic",                 dispatch<TclCommand_newElasticMaterial>},
  {"ElasticIsotropic3DThermal",        dispatch<OPS_ElasticIsotropicMaterialThermal>},
// Orthotropic
  {"ElasticOrthotropic",               dispatch<OPS_ElasticOrthotropicMaterial>},
  {"ElasticOrthotropicPlaneStress",    dispatch<OPS_ElasticOrthotropicPlaneStress>},
//
// Plasticity
//
  {"J2",                               dispatch<TclCommand_newPlasticMaterial>},
  {"J2Plasticity",                     dispatch<TclCommand_newPlasticMaterial>},
  {"J2N",                              dispatch<TclCommand_newPlasticMaterial>},
  {"J2L",                              dispatch<TclCommand_newPlasticMaterial>},
  {"J2Thermal",                        dispatch<TclCommand_newPlasticMaterial>},
  {"J2PlasticityThermal",              dispatch<TclCommand_newPlasticMaterial>},
  {"J2BeamFiber",                      dispatch<TclCommand_newPlasticMaterial>},

  {"SimplifiedJ2",                     dispatch<TclCommand_newPlasticMaterial>},
  {"J2Simplified",                     dispatch<TclCommand_newPlasticMaterial>},
  {"Simplified3DJ2",                   dispatch<TclCommand_newPlasticMaterial>},
  {"3DJ2",                             dispatch<TclCommand_newPlasticMaterial>},
  {"PlaneStressSimplifiedJ2",          dispatch<TclCommand_newPlasticMaterial>},
  {"DruckerPrager",                    dispatch<TclCommand_newPlasticMaterial>},

  {"UVCplanestress",                   dispatch<OPS_UVCplanestress       > },
  {"UVCmultiaxial",                    dispatch<OPS_UVCmultiaxial        > },
  {"J2PlateFibre",                     dispatch<OPS_J2PlateFibreMaterial >},
//
  {"ManzariDafalias",                  dispatch<OPS_ManzariDafaliasMaterial>},
  {"ManzariDafaliasRO",                dispatch<OPS_ManzariDafaliasMaterialRO>},


  {"DruckerPragerThermal",             dispatch<OPS_DruckerPragerMaterialThermal> },
  {"TruncatedDP",                      dispatch<OPS_LinearCap     > },
  {"FSAM",                             dispatch<OPS_FSAMMaterial  > },
  {"AcousticMedium",                   dispatch<OPS_AcousticMedium> },
  {"CycLiqCP",                         dispatch<OPS_CycLiqCPMaterial>},
  {"CycLiqCPSP",                       dispatch<OPS_CycLiqCPSPMaterial>},
  {"BoundingCamClay",                  dispatch<OPS_BoundingCamClayMaterial>},
//
// Wrapper
//
  {"InitStrainMaterial",               dispatch<TclCommand_addWrappingMaterial>},
  {"InitStrain",                       dispatch<TclCommand_addWrappingMaterial>},
  {"InitialStrain",                    dispatch<TclCommand_addWrappingMaterial>},
  {"InitStressMaterial",               dispatch<OPS_InitStressNDMaterial>},
  {"OrthotropicMaterial",              dispatch<OPS_OrthotropicMaterial>},
  {"Series3DMaterial",                 dispatch<OPS_Series3DMaterial>},
  {"Parallel3DMaterial",               dispatch<OPS_Parallel3DMaterial>},
  {"Parallel3D",                       dispatch<OPS_Parallel3DMaterial>},
  // Beam fiber (             22, 33, and 23 == 0)
  {"BeamFiber",                        dispatch<TclCommand_newPlateFiber>},
  {"BeamFiber2d",                      dispatch<TclCommand_newPlateFiber>},
  {"BeamFiber2dPS",                    dispatch<TclCommand_newPlateFiber>},
  // Plane 
  {"PlaneStressMaterial",              dispatch<TclCommand_addPlaneWrapper>},
  {"PlaneStress",                      dispatch<TclCommand_addPlaneWrapper>},
  {"PlaneStrainMaterial",              dispatch<TclCommand_addPlaneWrapper>},
  {"PlaneStrain",                      dispatch<TclCommand_addPlaneWrapper>},
  {"PlaneStressRebarMaterial",         dispatch<TclCommand_newPlateRebar>},
  // Plate  (constrain stress 33 == 13 == 23 == 0) 
  {"PlateRebarMaterial",               dispatch<TclCommand_newPlateRebar>},
  {"PlateRebar",                       dispatch<TclCommand_newPlateRebar>},
  {"PlateFiberMaterial",               dispatch<TclCommand_newPlateFiber>},
  {"PlateFiber",                       dispatch<TclCommand_newPlateFiber>},
//
// Other
//
  {"ReinforcedConcretePlaneStress",    dispatch<OPS_ReinforcedConcretePlaneStressMaterial>},
  {"PlaneStressLayeredMaterial",       dispatch<OPS_PlaneStressLayeredMaterial>},
  {"ASDConcrete3D",                    dispatch<OPS_ASDConcrete3DMaterial>},
  {"PlasticDamageConcrete",            dispatch<TclCommand_newConcreteMaterial>},
  {"FariaPlasticDamage",               dispatch<TclCommand_newConcreteMaterial>},
  {"PlasticDamageConcretePlaneStress", dispatch<OPS_NewPlasticDamageConcretePlaneStress>},
};

static std::unordered_map<std::string, OPS_Routine*> OldMaterialCommands = {
#ifdef OPS_USE_ASDPlasticMaterials
  {"ASDPlasticMaterial",            OPS_AllASDPlasticMaterials},
#endif
#ifdef _HAVE_Faria1998
  {"Faria1998", OPS_NewFaria1998Material},  
  {"Concrete", OPS_NewConcreteMaterial},
#endif
#ifdef _HAVE_Damage2p
  {"Damage2p",                        OPS_Damage2p},
#endif

  {"FAReinforcedConcretePlaneStress", OPS_FAReinforcedConcretePlaneStressMaterial},
  {"RAFourSteelRCPlaneStress",        OPS_RAFourSteelRCPlaneStressMaterial},
  {"FAFourSteelRCPlaneStress",        OPS_FAFourSteelRCPlaneStressMaterial},


  {"PrestressedConcretePlaneStress",   OPS_PrestressedConcretePlaneStressMaterial},
  {"FAPrestressedConcretePlaneStress", OPS_FAPrestressedConcretePlaneStressMaterial},
  {"RAFourSteetPCPlaneStress",         OPS_RAFourSteelPCPlaneStressMaterial},
  {"FAFourSteelPCPlaneStress",         OPS_FAFourSteelPCPlaneStressMaterial},

//{"MaterialCMM",    OPS_MaterialCMM},

  {"PM4Sand",                       OPS_PM4SandMaterial},
  {"J2CyclicBoundingSurface",       OPS_J2CyclicBoundingSurfaceMaterial},
  {"PM4Silt",                       OPS_PM4SiltMaterial},
  {"LinearElasticGGmax",            OPS_LinearElasticGGmaxMaterial},
  {"ContactMaterial2D",             OPS_ContactMaterial2DMaterial},
  {"ContactMaterial3D",             OPS_ContactMaterial3DMaterial},
  {"InitialStateAnalysisWrapper",   OPS_InitialStateAnalysisWrapperMaterial},
  {"stressDensity",                 OPS_StressDensityMaterial},
  {"IncrementalElasticIsotropic3D", OPS_IncrementalElasticIsotropicThreeDimensional},
  {"OrthotropicRAConcrete",         OPS_OrthotropicRotatingAngleConcreteT2DMaterial01},
  {"SmearedSteelDoubleLayer",       OPS_SmearedSteelDoubleLayerT2DMaterial01},
  {"SAniSandMS",                    OPS_SAniSandMSMaterial},
};
} // namespace OpenSees
