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
// Description: This file contains the function invoked when the user invokes
// the uniaxialMaterial command in the interpreter.
//
// Written: fmk, MHS, cmp
// Created: 07/99
//
#include <tcl.h>
#include <assert.h>
#include <string>
#include <Parsing.h>
#include <unordered_map>
#include <runtimeAPI.h>
#include <UniaxialMaterial.h>

// Viscous
extern OPS_Routine OPS_Maxwell;
extern OPS_Routine OPS_ViscousDamper;
extern OPS_Routine OPS_ViscousMaterial;
extern OPS_Routine OPS_ViscoelasticGap;
extern OPS_Routine OPS_DamperMaterial;

extern OPS_Routine OPS_ASD_SMA_3K;
extern OPS_Routine OPS_ASDConcrete1DMaterial;
extern OPS_Routine OPS_APDFMD;
extern OPS_Routine OPS_APDMD;
extern OPS_Routine OPS_APDVFD;
extern OPS_Routine OPS_Bilin02;
extern OPS_Routine OPS_Bilin;
extern OPS_Routine OPS_BilinearOilDamper;
extern OPS_Routine OPS_Bond_SP01;
extern OPS_Routine OPS_CFSSSWP;
extern OPS_Routine OPS_CFSWSWP;
extern OPS_Routine OPS_CableMaterial;
extern OPS_Routine OPS_Cast;
extern OPS_Routine OPS_CreepMaterial;
// Concrete
extern OPS_Routine OPS_Concrete02;
extern OPS_Routine OPS_Concrete02IS;
extern OPS_Routine OPS_Concrete02Thermal;
extern OPS_Routine OPS_ConcreteCM;
extern OPS_Routine OPS_ConcreteD;
extern OPS_Routine OPS_ConcreteECThermal;
extern OPS_Routine OPS_ConcreteL01Material;
extern OPS_Routine OPS_ConcreteSakaiKawashima;
extern OPS_Routine OPS_ConcreteZ01Material;
extern OPS_Routine OPS_ConfinedConcrete01Material;

extern OPS_Routine OPS_DoddRestr;
extern OPS_Routine OPS_Dodd_Restrepo;
extern OPS_Routine OPS_EPPGapMaterial;
extern OPS_Routine OPS_ElasticBilin;
extern OPS_Routine OPS_ElasticMaterialThermal;
extern OPS_Routine OPS_ElasticMultiLinear;
extern OPS_Routine OPS_ElasticPPMaterial;
extern OPS_Routine OPS_ElasticPowerFunc;
extern OPS_Routine OPS_FRPConfinedConcrete02;
extern OPS_Routine OPS_FRPConfinedConcrete02;
extern OPS_Routine OPS_FRPConfinedConcrete;
extern OPS_Routine OPS_GMG_CyclicReinforcedConcrete;
extern OPS_Routine OPS_FRCC;
extern OPS_Routine OPS_GNGMaterial;
extern OPS_Routine OPS_HoehlerStanton;
extern OPS_Routine OPS_HookGap;
extern OPS_Routine OPS_HystereticMaterial;
extern OPS_Routine OPS_HystereticPoly;
extern OPS_Routine OPS_HystereticAsym;
extern OPS_Routine OPS_HystereticSmooth;
extern OPS_Routine OPS_HystereticSMMaterial;
extern OPS_Routine OPS_IMKBilin;
extern OPS_Routine OPS_IMKPeakOriented;
extern OPS_Routine OPS_IMKPinching;
extern OPS_Routine OPS_JankowskiImpact;
extern OPS_Routine OPS_ImpactMaterial;
extern OPS_Routine OPS_Masonry;
extern OPS_Routine OPS_Masonryt;
extern OPS_Routine OPS_MinMaxMaterial;
extern OPS_Routine OPS_ModIMKPeakOriented02;
extern OPS_Routine OPS_ModIMKPeakOriented;
extern OPS_Routine OPS_ModIMKPinching;
extern OPS_Routine OPS_ModIMKPinching02;
extern OPS_Routine OPS_MultiLinear;
extern OPS_Routine OPS_OriginCentered;
extern OPS_Routine OPS_PinchingLimitState;
extern OPS_Routine OPS_PinchingLimitStateMaterial;
extern OPS_Routine OPS_PySimple3;
extern OPS_Routine OPS_RambergOsgoodSteel;
extern OPS_Routine OPS_Ratchet;
extern OPS_Routine OPS_ReinforcingSteel;
extern OPS_Routine OPS_ResilienceLow;
extern OPS_Routine OPS_ResilienceMaterialHR;
extern OPS_Routine OPS_SAWSMaterial;
extern OPS_Routine OPS_SLModel;
extern OPS_Routine OPS_SMAMaterial;
extern OPS_Routine OPS_SPSW02;           // SAJalali
extern OPS_Routine OPS_SimpleFractureMaterial;

extern OPS_Routine OPS_StainlessECThermal;
extern OPS_Routine OPS_SteelBRB;
extern OPS_Routine OPS_SteelECThermal;
extern OPS_Routine OPS_SteelFractureDI; // galvisf
extern OPS_Routine OPS_SteelMPF;
extern OPS_Routine OPS_SteelZ01Material;
extern OPS_Routine OPS_Steel02Fatigue;
extern OPS_Routine OPS_Steel4;

extern OPS_Routine OPS_TDConcrete;       // ntosic
extern OPS_Routine OPS_TDConcreteEXP;    // ntosic
extern OPS_Routine OPS_TDConcreteMC10;   // ntosic
extern OPS_Routine OPS_TDConcreteMC10NL; // ntosic
extern OPS_Routine OPS_TendonL01Material;
extern OPS_Routine OPS_Trilinwp2;
extern OPS_Routine OPS_Trilinwp;
extern OPS_Routine OPS_UVCuniaxial;
extern OPS_Routine OPS_pyUCLA;
extern void *OPS_ConcretewBeta();

extern Tcl_CmdProc Create_OOHystereticMaterial;



#if 0
const char** DeprecatedUniaxialMaterials {
  {"Bilin02", "This material is superceded by \"IMKBilin\" and \"HystereticSM\""},
  {"CFSSSWP", ""},
  {"CFSWSWP", ""},
  {"APDVFD",  ""},
  {"APDMD",   ""},
  {"APDFMD",  ""},
};
#endif

typedef UniaxialMaterial*(G3_TclUniaxialPackage)(ClientData, Tcl_Interp *, int, TCL_Char ** const);
G3_TclUniaxialPackage TclBasicBuilder_addFedeasMaterial;
G3_TclUniaxialPackage TclBasicBuilder_addSnapMaterial;
G3_TclUniaxialPackage TclBasicBuilder_addDrainMaterial;
std::unordered_map<std::string, G3_TclUniaxialPackage *> tcl_uniaxial_package_table {
  {"DRAIN",              TclBasicBuilder_addDrainMaterial },
  {"SNAP",               TclBasicBuilder_addSnapMaterial  },
  {"snap",               TclBasicBuilder_addSnapMaterial  },
// #if defined(_STEEL2) || defined(OPSDEF_UNIAXIAL_FEDEAS)
//{"FEDEAS",             TclBasicBuilder_addFedeasMaterial},
// #endif
};


// Elastic
extern Tcl_CmdProc TclCommand_newElasticUniaxialMaterial;
// Plastic
extern Tcl_CmdProc TclCommand_newPlasticMaterial;
extern Tcl_CmdProc TclCommand_newUniaxialJ2Plasticity;
// Fedeas
extern Tcl_CmdProc TclCommand_newFedeasSteel;
extern Tcl_CmdProc TclCommand_newFedeasUniaxialDamage;
extern Tcl_CmdProc TclCommand_newFedeasConcrete;
// Wrapper
extern Tcl_CmdProc TclCommand_addWrappingMaterial;
extern Tcl_CmdProc TclCommand_newFatigueMaterial;
extern Tcl_CmdProc TclCommand_newParallelMaterial;
extern OPS_Routine OPS_SeriesMaterial;
// Bearing
Tcl_CmdProc TclCommand_AxialSp;
Tcl_CmdProc TclCommand_AxialSpHD;
Tcl_CmdProc TclCommand_KikuchiAikenHDR;
Tcl_CmdProc TclCommand_KikuchiAikenLRB;
// Concrete
Tcl_CmdProc TclCommand_newUniaxialConcrete04;
Tcl_CmdProc TclCommand_newUniaxialConcrete06;
Tcl_CmdProc TclCommand_newUniaxialConcrete07;
// Bouc-Wen
extern Tcl_CmdProc TclCommand_newBoucWen;
extern Tcl_CmdProc TclCommand_newUniaxialBoucWen;
extern Tcl_CmdProc TclCommand_newBoucWenMG;
extern OPS_Routine OPS_DegradingPinchedBW;
// Abutment
Tcl_CmdProc TclCommand_HyperbolicGapMaterial;
// Other
static Tcl_CmdProc TclDispatch_newUniaxialPinching4;
extern Tcl_CmdProc TclDispatch_LegacyUniaxials;


typedef UniaxialMaterial* (TclDispatch_UniaxialMaterial)(G3_Runtime*, int, TCL_Char ** const);
TclDispatch_UniaxialMaterial TclCommand_ReinforcingSteel;


template <OPS_Routine fn> static int
dispatch(ClientData clientData, Tcl_Interp* interp, int argc, G3_Char** const argv)
{
  BasicModelBuilder *builder = static_cast<BasicModelBuilder*>(clientData);
  G3_Runtime *rt = G3_getRuntime(interp);
  UniaxialMaterial* theMaterial = (UniaxialMaterial*)fn( rt, argc, argv );

  if (builder->addTaggedObject<UniaxialMaterial>(*theMaterial) != TCL_OK) {
    opserr << OpenSees::PromptValueError << "Failed to add UniaxialMaterial to the model builder.\n";
    delete theMaterial;
    return TCL_ERROR;
  }
  return TCL_OK;
}

template <UniaxialMaterial*(*fn)(G3_Runtime*, int, TCL_Char** const)> static int
dispatch(ClientData clientData, Tcl_Interp* interp, int argc, TCL_Char** const argv)
{
  assert(clientData != nullptr);
  BasicModelBuilder *builder = static_cast<BasicModelBuilder*>(clientData);
  G3_Runtime *rt = G3_getRuntime(interp);
  UniaxialMaterial* theMaterial = fn( rt, argc, argv );

  if (builder->addTaggedObject<UniaxialMaterial>(*theMaterial) != TCL_OK) {
    opserr << OpenSees::PromptValueError << "Could not add uniaxialMaterial to the model builder.\n";
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

std::unordered_map<std::string, Tcl_CmdProc*> uniaxial_dispatch {

    {"Elastic",                  dispatch<TclCommand_newElasticUniaxialMaterial>},
//
// Plasticity
//
    {"ElasticPP",                dispatch<OPS_ElasticPPMaterial>              },
    {"UniaxialJ2Plasticity",     dispatch<TclCommand_newUniaxialJ2Plasticity> },
    {"UVCuniaxial",              dispatch<OPS_UVCuniaxial>                    },
    {"Hardening",                dispatch<TclCommand_newPlasticMaterial>      },
    {"Hardening2",               dispatch<TclCommand_newPlasticMaterial>      },
//
// Bouc
//
    {"BWBF",                   dispatch<TclCommand_newBoucWen>         },
    {"BWBN",                   dispatch<TclCommand_newBoucWen>         },
    {"BoucWenOriginal",        dispatch<TclCommand_newBoucWen>         },
    {"BoucWen",                dispatch<TclCommand_newUniaxialBoucWen> },
    {"BoucWenMG",              dispatch<TclCommand_newBoucWenMG>       },
    {"DegradingPinchedBW",     dispatch<OPS_DegradingPinchedBW>        },
//
// Steel
//
    {"Steel",                  dispatch<TclCommand_newPlasticMaterial>},
    {"Steel01",                dispatch<TclCommand_newFedeasSteel>     },
    {"Steel02",                dispatch<TclCommand_newFedeasSteel>     },
    {"Steel2",                 dispatch<TclCommand_newFedeasSteel>     },
    {"Steel4",                 dispatch<OPS_Steel4>                    },
    {"RambergOsgood",          dispatch<OPS_RambergOsgoodSteel>        },
    {"RambergOsgoodSteel",     dispatch<OPS_RambergOsgoodSteel>        },
    {"ReinforcingSteel",       dispatch<OPS_ReinforcingSteel>          },
//  {"ReinforcingSteel",       dispatch<TclCommand_ReinforcingSteel>   },
    {"SteelBRB",               dispatch<OPS_SteelBRB>                  },
    {"SteelFractureDI",        dispatch<OPS_SteelFractureDI>           },
    {"Steel02Fatigue",         dispatch<OPS_Steel02Fatigue>            },

    {"Dodd_Restrepo",          dispatch<OPS_Dodd_Restrepo>             },
    {"DoddRestrepo" ,          dispatch<OPS_Dodd_Restrepo>             },
    {"Restrepo",               dispatch<OPS_Dodd_Restrepo>             },

#if !defined(_NO_NEW_RESTREPO)
    {"DoddRestr",              dispatch<OPS_DoddRestr>                 },
#endif

    {"SteelZ01Material",       dispatch<OPS_SteelZ01Material>          },
    {"SteelZ01",               dispatch<OPS_SteelZ01Material>          },

// Concretes
    {"Concrete01",             dispatch<TclCommand_newFedeasConcrete>  },
    {"Concrete02",             dispatch<TclCommand_newFedeasConcrete>     },
    {"Concrete04",             dispatch<TclCommand_newUniaxialConcrete04> },
    {"Concrete06",             dispatch<TclCommand_newUniaxialConcrete06> },
    {"Concrete07",             dispatch<TclCommand_newUniaxialConcrete07> },

    {"Concrete02IS",           dispatch<OPS_Concrete02IS>              },
    {"ConcreteCM",             dispatch<OPS_ConcreteCM>                },
    {"ConfinedConcrete01",     dispatch<OPS_ConfinedConcrete01Material>},
    {"ConfinedConcrete",       dispatch<OPS_ConfinedConcrete01Material>},

    {"ConcreteZ01Material",    dispatch<OPS_ConcreteZ01Material>       },
    {"ConcreteZ01",            dispatch<OPS_ConcreteZ01Material>       },


#if 0
    { "ConcretewBeta",       dispatch<OPS_ConcretewBeta>    }
#endif
//
// Viscous
//
    {"Maxwell",                dispatch<OPS_Maxwell>                   },
    {"MaxwellMaterial",        dispatch<OPS_Maxwell>                   },
    {"ViscousDamper",          dispatch<OPS_ViscousDamper>             },
    {"DamperMaterial",         dispatch<OPS_DamperMaterial>            },
    {"BilinearOilDamper",      dispatch<OPS_BilinearOilDamper>         },
//
// Multilinear
//
    {"BilinMaterial",           dispatch<OPS_Bilin>                     },
    {"Bilin",                   dispatch<OPS_Bilin>                     },
    {"MultiLinear",             dispatch<OPS_MultiLinear>               },
    {"IMKBilin",                dispatch<OPS_IMKBilin>                  },
    {"IMKPeakOriented",         dispatch<OPS_IMKPeakOriented>           },
    {"IMKPinching",             dispatch<OPS_IMKPinching>               },
    {"ModIMKPinching",          dispatch<OPS_ModIMKPinching>            },
    {"ModIMKPinching02",        dispatch<OPS_ModIMKPinching02>          },
    {"ModIMKPeakOriented",      dispatch<OPS_ModIMKPeakOriented>        },
    {"ModIMKPeakOriented02",    dispatch<OPS_ModIMKPeakOriented02>      },
    {"Bilin02",                 dispatch<OPS_Bilin02>                   },

// Piles
    {"PySimple3",               dispatch<OPS_PySimple3>                 },

//
// Wrappers
//
    {"InitStrainMaterial",     dispatch<TclCommand_addWrappingMaterial>},
    {"InitialStrain",          dispatch<TclCommand_addWrappingMaterial>},
    {"InitStrain",             dispatch<TclCommand_addWrappingMaterial>},
    {"InitStressMaterial",     dispatch<TclCommand_addWrappingMaterial>},
    {"InitialStress",          dispatch<TclCommand_addWrappingMaterial>},
    {"InitStress",             dispatch<TclCommand_addWrappingMaterial>},
    {"ContinuumUniaxial",      dispatch<TclCommand_addWrappingMaterial>},
    {"MinMaxMaterial",         dispatch<OPS_MinMaxMaterial>            },
    {"MinMax",                 dispatch<OPS_MinMaxMaterial>            },
    {"Series",                 dispatch<OPS_SeriesMaterial>            },
    {"Parallel",               dispatch<TclCommand_newParallelMaterial>},
    {"Ratchet",                dispatch<OPS_Ratchet>                   },
    {"Fatigue",                dispatch<TclCommand_newFatigueMaterial> },
  
// Other

    {"GNG",                  dispatch<OPS_GNGMaterial>                 },
    {"Bond_SP01",            dispatch<OPS_Bond_SP01>                 },
    {"Bond",                 dispatch<OPS_Bond_SP01>                 },
    {"APDFMD",               dispatch<OPS_APDFMD> },
    {"APDMD",                dispatch<OPS_APDMD> },
    {"APDVFD",               dispatch<OPS_APDVFD> },

    {"FedeasUniaxialDamage", dispatch<TclCommand_newFedeasUniaxialDamage>  },
    {"KikuchiAikenHDR",      dispatch<TclCommand_KikuchiAikenHDR>       },
    {"KikuchiAikenLRB",      dispatch<TclCommand_KikuchiAikenLRB>       },

    {"AxialSp",              dispatch<TclCommand_AxialSp>               },
    {"AxialSpHD",            dispatch<TclCommand_AxialSpHD>             },

/*
  {"PlateBearingConnectionThermal",  OPS_PlateBearingConnectionThermal},
  {"PinchingLimitStateMaterial",     OPS_PinchingLimitState           },
*/

// Other
    {"ElasticBilin",           dispatch<OPS_ElasticBilin>              },
    {"ElasticBilinear",        dispatch<OPS_ElasticBilin>              },

    {"ImpactMaterial",         dispatch<OPS_ImpactMaterial>            },
    {"Impact",                 dispatch<OPS_ImpactMaterial>            },

    {"SimpleFractureMaterial", dispatch<OPS_SimpleFractureMaterial>    },
    {"SimpleFracture",         dispatch<OPS_SimpleFractureMaterial>    },
//
    {"Cast",                   dispatch<OPS_Cast>                      },
    {"CastFuse",               dispatch<OPS_Cast>                      },

    {"ElasticMultiLinear",     dispatch<OPS_ElasticMultiLinear>        },
    {"ElasticPowerFunc",       dispatch<OPS_ElasticPowerFunc>          },

/* 
    {"HoehlerStanton",         dispatch<OPS_HoehlerStanton>            },
*/

    {"SLModel",                dispatch<OPS_SLModel>                   },

    {"OriginCentered",         dispatch<OPS_OriginCentered>            },

    {"HookGap",                dispatch<OPS_HookGap>                   },

    {"HyperbolicGapMaterial",  dispatch<TclCommand_HyperbolicGapMaterial>},

    {"FRPConfinedConcrete02",  dispatch<OPS_FRPConfinedConcrete02>     },
    {"FRCC",                   dispatch<OPS_FRCC>                      },
    {"GMG_CyclicReinforcedConcrete", dispatch<OPS_GMG_CyclicReinforcedConcrete>},

    {"PinchingLimitState",     dispatch<OPS_PinchingLimitState>        },

    {"pyUCLA",                 dispatch<OPS_pyUCLA>                    },
    {"PYUCLA",                 dispatch<OPS_pyUCLA>                    },


    {"JankowskiImpact",        dispatch<OPS_JankowskiImpact>           },

// Thermal
    {"Steel01Thermal",         dispatch<TclCommand_newFedeasSteel>     },
    {"Steel02Thermal",         dispatch<TclCommand_newFedeasSteel>     },

    {"SteelECThermal",         dispatch<OPS_SteelECThermal>            },

    {"StainlessECThermal",     dispatch<OPS_StainlessECThermal>        },

    {"ElasticThermal",         dispatch<OPS_ElasticMaterialThermal>    },

    {"ConcreteECThermal",      dispatch<OPS_ConcreteECThermal>         },

    {"Concrete02Thermal",      dispatch<OPS_Concrete02Thermal>         },
//
    {"ConcreteD",              dispatch<OPS_ConcreteD>                 },

    {"ConcreteSakaiKawashima", dispatch<OPS_ConcreteSakaiKawashima>    },

    {"SteelMPF",               dispatch<OPS_SteelMPF>                  },

    {"ResilienceLow",          dispatch<OPS_ResilienceLow>             },

    {"ResilienceMaterialHR",   dispatch<OPS_ResilienceMaterialHR>      },

    {"CFSWSWP",                dispatch<OPS_CFSWSWP>                   },
    {"CFSSSWP",                dispatch<OPS_CFSSSWP>                   },

    {"FRPConfinedConcrete",    dispatch<OPS_FRPConfinedConcrete>       },

    {"Masonry",                dispatch<OPS_Masonry>                   },

    {"Trilinwp",               dispatch<OPS_Trilinwp>                  },

    {"Trilinwp2",              dispatch<OPS_Trilinwp2>                 },

    {"Masonryt",               dispatch<OPS_Masonryt>                  },

    
    {"Hysteretic",             dispatch<OPS_HystereticMaterial>        },
    {"HystereticAsym",         dispatch<OPS_HystereticAsym>            },
    {"HystereticSmooth",       dispatch<OPS_HystereticSmooth>          },
    {"HystereticSMMaterial",   dispatch<OPS_HystereticSMMaterial>      },

    {"ElasticPPGap",           dispatch<OPS_EPPGapMaterial>            },


    {"OOHysteretic",           dispatch<Create_OOHystereticMaterial>   },

    {"Viscous",                dispatch<OPS_ViscousMaterial>           },
    {"ViscoelasticGap",        dispatch<OPS_ViscoelasticGap>           },

    {"SAWSMaterial",           dispatch<OPS_SAWSMaterial>              },
    {"SAWS",                   dispatch<OPS_SAWSMaterial>              },


    {"ConcreteL01Material",    dispatch<OPS_ConcreteL01Material>       },
    {"ConcreteL01",            dispatch<OPS_ConcreteL01Material>       },

    {"Creep",                  dispatch<OPS_CreepMaterial>             },

    {"TendonL01Material",      dispatch<OPS_TendonL01Material>         },
    {"TendonL01",              dispatch<OPS_TendonL01Material>         },

    {"Cable",                  dispatch<OPS_CableMaterial>             },

    {"SMA",                    dispatch<OPS_SMAMaterial>               },

    {"ASD_SMA_3K",             dispatch<OPS_ASD_SMA_3K>                },

    {"ASDConcrete1D",          dispatch<OPS_ASDConcrete1DMaterial>     },

    {"HystereticPoly",         dispatch<OPS_HystereticPoly>            },

    {"SPSW02",                 dispatch<OPS_SPSW02>                    },

    {"TDConcreteEXP",          dispatch<OPS_TDConcreteEXP>             },
    {"TDConcrete",             dispatch<OPS_TDConcrete>                },
    {"TDConcreteMC10",         dispatch<OPS_TDConcreteMC10>            },
    {"TDConcreteMC10NL",       dispatch<OPS_TDConcreteMC10NL>          },

    {"Pinching4",              TclDispatch_newUniaxialPinching4        },

// Legacy
    {"Elastic2",               TclDispatch_LegacyUniaxials             },
    {"ENT",                    TclDispatch_LegacyUniaxials             },
    {"BarSlip",                TclDispatch_LegacyUniaxials             },
    {"ShearPanel",             TclDispatch_LegacyUniaxials             },
    {"Concrete01WithSITC",     TclDispatch_LegacyUniaxials             },
};

