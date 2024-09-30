/* ****************************************************************** **
**    Opensee - Open System for Earthquake Engineering Simulation    **
**          Pacific Earthquake Engineering Research Center            **
** ****************************************************************** */
//
// Description: This file contains the function invoked when the user invokes
// the uniaxialMaterial command in the interpreter.
//
// Written: fmk, MHS
// Created: 07/99
//
//
#include <unordered_map>
#include <tcl.h>
#include <runtimeAPI.h>
#include <elementAPI.h>
#include <BasicModelBuilder.h>
#include <Elastic2Material.h>   // ZHY
#include <HardeningMaterial.h>  // MHS
#include <HardeningMaterial2.h> // MHS
#include <Steel03.h>            // KM
#include <Concrete01WithSITC.h> // Won Lee
#include <ECC01.h>              // Won Lee
#include <Concrete04.h>
#include <Concrete05.h>
#include <Concrete06.h>              // LMS
#include <Concrete07.h>              // JDW
#include <HystereticBackbone.h>      // MHS
#include <EPPGapMaterial.h>          // Mackie
#include <PathIndependentMaterial.h> // MHS
#include <BackboneMaterial.h>        // MHS
#include <FatigueMaterial.h>         // Patxi
#include <ENTMaterial.h>             // MHS
#include <Pinching4Material.h>       // NM
#include <ShearPanelMaterial.h>      // NM
#include <BarSlipMaterial.h>         // NM
#include <Bond_SP01.h>               // JZ

#include <SteelMP.h>               //Quan & Michele
#include <SmoothPSConcrete.h>      //Quan & Michele
#include <SelfCenteringMaterial.h> //JAE

#include <SMAMaterial.h> // Davide Fugazza
#include <Masonry.h>
#include <Trilinwp.h>
#include <Trilinwp2.h>
#include <Masonryt.h>

#include <Vector.h>
#include <string.h>

#include <UniaxialJ2Plasticity.h> // Quan

extern OPS_Routine OPS_ElasticPPMaterial;
extern OPS_Routine OPS_EPPGapMaterial;
extern OPS_Routine OPS_HardeningMaterial;
extern OPS_Routine OPS_HystereticMaterial;
extern OPS_Routine OPS_CableMaterial;
extern OPS_Routine OPS_Bilin;
extern OPS_Routine OPS_Bilin02;
extern OPS_Routine OPS_FRPConfinedConcrete02;
extern OPS_Routine OPS_SteelFractureDI; // galvisf
extern OPS_Routine OPS_RambergOsgoodSteel;
extern OPS_Routine OPS_ReinforcingSteel;
extern OPS_Routine OPS_PinchingLimitStateMaterial;
extern OPS_Routine OPS_SAWSMaterial;
extern OPS_Routine OPS_ConcreteZ01Material;
extern OPS_Routine OPS_ConcreteL01Material;
extern OPS_Routine OPS_SteelZ01Material;
extern OPS_Routine OPS_TendonL01Material;
extern OPS_Routine OPS_ConfinedConcrete01Material;
extern OPS_Routine OPS_SimpleFractureMaterial;
extern OPS_Routine OPS_HoehlerStanton;
extern OPS_Routine OPS_BilinearOilDamper;
extern OPS_Routine OPS_Cast;
extern OPS_Routine OPS_DoddRestr;
extern OPS_Routine OPS_ElasticMultiLinear;
extern OPS_Routine OPS_HookGap;
extern OPS_Routine OPS_HyperbolicGapMaterial;
extern OPS_Routine OPS_FRPConfinedConcrete;
extern OPS_Routine OPS_FRPConfinedConcrete02;
extern OPS_Routine OPS_UVCuniaxial;
extern OPS_Routine OPS_Steel01Thermal;
extern OPS_Routine OPS_Concrete02Thermal;
extern OPS_Routine OPS_StainlessECThermal;     // L.Jiang [SIF]
extern OPS_Routine OPS_SteelECThermal;         // L.Jiang [SIF]
extern OPS_Routine OPS_ConcreteECThermal;      // L.Jiang [SIF]
extern OPS_Routine OPS_ElasticMaterialThermal; // L.Jiang[SIF]
// extern OPS_Routine OPS_PlateBearingConnectionThermal;
extern OPS_Routine OPS_ASD_SMA_3K; // Luca Aceto
extern OPS_Routine OPS_IMKPeakOriented;
extern OPS_Routine OPS_IMKBilin;
extern OPS_Routine OPS_IMKPinching;
extern OPS_Routine OPS_ModIMKPeakOriented;
extern OPS_Routine OPS_ModIMKPeakOriented02;
extern OPS_Routine OPS_ModIMKPinching;
extern OPS_Routine OPS_ModIMKPinching02;
extern void *OPS_ConcretewBeta(void);
extern OPS_Routine OPS_ConcreteD;
extern OPS_Routine OPS_PinchingLimitState;
extern OPS_Routine OPS_OriginCentered;
extern OPS_Routine OPS_Steel2;
extern OPS_Routine OPS_ConcreteSakaiKawashima;
extern OPS_Routine OPS_ResilienceMaterialHR;
extern OPS_Routine OPS_CFSSSWP;
extern OPS_Routine OPS_CFSWSWP;
extern OPS_Routine OPS_ResilienceLow;
extern OPS_Routine OPS_ViscousMaterial;
extern OPS_Routine OPS_SteelMPF;   // K Kolozvari
extern OPS_Routine OPS_Bond_SP01;  // K Kolozvari
extern OPS_Routine OPS_OOHystereticMaterial;
extern OPS_Routine OPS_ElasticPowerFunc;
extern OPS_Routine OPS_UVCuniaxial;
extern OPS_Routine OPS_DegradingPinchedBW;
extern OPS_Routine OPS_SMAMaterial;
extern OPS_Routine OPS_HystereticPoly; // Salvatore Sessa 14-Jan-2021 Mail: salvatore.sessa2@unina.it
extern OPS_Routine OPS_Masonry;
extern OPS_Routine OPS_Trilinwp;
extern OPS_Routine OPS_Trilinwp2;
extern OPS_Routine OPS_Masonryt;

// extern int TclCommand_ConfinedConcrete02(ClientData clientData, Tcl_Interp
// *interp, int argc, 					 TCL_Char ** const argv, TclBasicBuilder *theTclBuilder);

extern UniaxialMaterial *Tcl_AddLimitStateMaterial(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char **arg);


#include <packages.h>

typedef struct uniaxialPackageCommand {
  char *funcName;
  void *(*funcPtr)();
  struct uniaxialPackageCommand *next;
} UniaxialPackageCommand;

static UniaxialPackageCommand *theUniaxialPackageCommands = NULL;

static void
printCommand(int argc, TCL_Char ** const argv)
{
  opserr << "Input command: ";
  for (int i = 0; i < argc; ++i)
    opserr << argv[i] << " ";
  opserr << endln;
}

// external functions

Tcl_CmdProc TclCommand_KikuchiAikenHDR;
Tcl_CmdProc TclCommand_KikuchiAikenLRB;

UniaxialMaterial *TclBasicBuilder_addDrainMaterial(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char ** const argv);

UniaxialMaterial *TclBasicBuilder_addSnapMaterial(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char ** const argv);

UniaxialMaterial *TclBasicBuilder_addPyTzQzMaterial(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char ** const argv,
                                                    Domain *theDomain);

UniaxialMaterial *TclBasicBuilder_FRPCnfinedConcrete(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char ** const argv,
                                                     Domain *theDomain);

UniaxialMaterial *TclBasicBuilder_addDegradingMaterial(ClientData, Tcl_Interp *, int, TCL_Char **);

extern "C" int OPS_ResetInputNoBuilder(ClientData clientData, Tcl_Interp *interp, int cArg, int mArg, TCL_Char ** const argv, Domain *domain);

int
TclBasicBuilderUniaxialMaterialCommand(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char ** const argv, Domain *theDomain)
{

  G3_Runtime *rt = G3_getRuntime(interp);
  BasicModelBuilder *builder = static_cast<BasicModelBuilder*>(clientData);

  // Make sure there is a minimum number of arguments
  if (argc < 3) {
    opserr << "WARNING insufficient number of uniaxial material arguments\n";
    opserr << "Want: uniaxialMaterial type? tag? <specific material args>"
           << endln;
    return TCL_ERROR;
  }

    OPS_ResetInputNoBuilder(clientData, interp, 2, argc, argv, theDomain);

  // Pointer to a uniaxial material that will be added to the model builder
  UniaxialMaterial *theMaterial = nullptr;

  // Check argv[2] for uniaxial material type
  if (strcmp(argv[1], "Elastic") == 0) {
      return TCL_ERROR;


#if 0
  } else if (strcmp(argv[1],"HoehlerStanton") == 0) {
    void *theMat = OPS_HoehlerStanton(rt, argc, argv);
    if (theMat != 0)
      theMaterial = (UniaxialMaterial *)theMat;
    else
      return TCL_ERROR;
#endif

  } else if ((strcmp(argv[1], "BilinearOilDamper") == 0)) {
    void *theMat = OPS_BilinearOilDamper(rt, argc, argv);
    if (theMat != 0)
      theMaterial = (UniaxialMaterial *)theMat;
    else
      return TCL_ERROR;

  } else if ((strcmp(argv[1], "Bond_SP01") == 0) ||
             (strcmp(argv[1], "Bond") == 0)) {
    void *theMat = OPS_Bond_SP01(rt, argc, argv);
    if (theMat != 0)
      theMaterial = (UniaxialMaterial *)theMat;
    else
      return TCL_ERROR;

#ifndef _NO_NEW_RESTREPO
  } else if ((strcmp(argv[1], "DoddRestr") == 0)) {

    void *theMat = OPS_DoddRestr(rt, argc, argv);
    if (theMat != 0)
      theMaterial = (UniaxialMaterial *)theMat;
    else
      return TCL_ERROR;
#endif

#if 0
  } else if (strcmp(argv[1],"HoehlerStanton") == 0) {
    void *theMat = OPS_HoehlerStanton(rt, argc, argv);
    if (theMat != 0)
      theMaterial = (UniaxialMaterial *)theMat;
    else
      return TCL_ERROR;
#endif


  } else if (strcmp(argv[1], "ReinforcingSteel") == 0) {
    void *theMat = OPS_ReinforcingSteel(rt, argc, argv);
    if (theMat != 0)
      theMaterial = (UniaxialMaterial *)theMat;
    else
      return TCL_ERROR;

  } else if (strcmp(argv[1], "FRPConfinedConcrete02") == 0) {
    void *theMat = OPS_FRPConfinedConcrete02(rt, argc, argv);
    if (theMat != 0)
      theMaterial = (UniaxialMaterial *)theMat;
    else
      return TCL_ERROR;

  } else if ((strcmp(argv[1], "PinchingLimitState") == 0) ||
             (strcmp(argv[1], "PinchingLimitStateMaterial") == 0)) {
    void *theMat = OPS_PinchingLimitState(rt, argc, argv);
    if (theMat != 0)
      theMaterial = (UniaxialMaterial *)theMat;
    else
      return TCL_ERROR;

  } else if (strcmp(argv[1], "ModIMKPinching") == 0) {
    void *theMat = OPS_ModIMKPinching(rt, argc, argv);
    if (theMat != 0)
      theMaterial = (UniaxialMaterial *)theMat;
    else
      return TCL_ERROR;

  } else if (strcmp(argv[1], "ModIMKPinching02") == 0) {
    void *theMat = OPS_ModIMKPinching02(rt, argc, argv);
    if (theMat != 0)
      theMaterial = (UniaxialMaterial *)theMat;
    else
      return TCL_ERROR;

  } else if (strcmp(argv[1], "ModIMKPeakOriented02") == 0) {
    void *theMat = OPS_ModIMKPeakOriented02(rt, argc, argv);
    if (theMat != 0)
      theMaterial = (UniaxialMaterial *)theMat;
    else
      return TCL_ERROR;

  } else if (strcmp(argv[1], "PlateBearingConnectionThermal") == 0) {
    // void *theMat = OPS_PlateBearingConnectionThermal(rt, argc, argv);
    void *theMat = 0;
    if (theMat != 0)
      theMaterial = (UniaxialMaterial *)theMat;
    else
      return TCL_ERROR;

  } else if (strcmp(argv[1], "Steel01Thermal") == 0) {
    void *theMat = OPS_Steel01Thermal(rt, argc, argv);
    if (theMat != 0)
      theMaterial = (UniaxialMaterial *)theMat;
    else
      return TCL_ERROR;

  } else if (strcmp(argv[1], "ConcretewBeta") == 0) {
    void *theMat = OPS_ConcretewBeta();
    if (theMat != 0)
      theMaterial = (UniaxialMaterial *)theMat;
    else
      return TCL_ERROR;



  } else if (strcmp(argv[1], "Concrete02Thermal") == 0) {
    void *theMat = OPS_Concrete02Thermal(rt, argc, argv);
    if (theMat != 0)
      theMaterial = (UniaxialMaterial *)theMat;
    else
      return TCL_ERROR;

  } else if (strcmp(argv[1], "Elastic2") == 0) {
    if (argc < 4 || argc > 5) {
      opserr << "WARNING invalid number of arguments\n";
      printCommand(argc, argv);
      opserr << "Want: uniaxialMaterial Elastic tag? E? <eta?>" << endln;
      return TCL_ERROR;
    }

    int tag;
    double E;
    double eta = 0.0;

    if (Tcl_GetInt(interp, argv[2], &tag) != TCL_OK) {
      opserr << "WARNING invalid uniaxialMaterial Elastic tag" << endln;
      return TCL_ERROR;
    }

    if (Tcl_GetDouble(interp, argv[3], &E) != TCL_OK) {
      opserr << "WARNING invalid E\n";
      opserr << "uniaxiaMaterial Elastic: " << tag << endln;
      return TCL_ERROR;
    }

    if (argc == 5) {
      if (Tcl_GetDouble(interp, argv[4], &eta) != TCL_OK) {
        opserr << "WARNING invalid eta\n";
        opserr << "uniaxialMaterial Elastic: " << tag << endln;
        return TCL_ERROR;
      }
    }

    // Parsing was successful, allocate the material
    theMaterial = new Elastic2Material(tag, E, eta);

  } else if (strcmp(argv[1], "ENT") == 0) {
    if (argc < 4) {
      opserr << "WARNING invalid number of arguments\n";
      printCommand(argc, argv);
      opserr << "Want: uniaxialMaterial ENT tag? E?" << endln;
      return TCL_ERROR;
    }

    int tag;
    double E;

    if (Tcl_GetInt(interp, argv[2], &tag) != TCL_OK) {
      opserr << "WARNING invalid uniaxialMaterial ENT tag" << endln;
      return TCL_ERROR;
    }

    if (Tcl_GetDouble(interp, argv[3], &E) != TCL_OK) {
      opserr << "WARNING invalid E\n";
      opserr << "uniaxiaMaterial ENT: " << tag << endln;
      return TCL_ERROR;
    }

    // Parsing was successful, allocate the material
    theMaterial = new ENTMaterial(tag, E);

  }

  else if (strcmp(argv[1], "ElasticPP") == 0) {
    void *theMat = OPS_ElasticPPMaterial(rt, argc, argv);
    if (theMat != 0)
      theMaterial = (UniaxialMaterial *)theMat;
    else
      return TCL_ERROR;
  }

  else if (strcmp(argv[1], "Steel03") == 0) {
    // Check that there is the minimum number of arguments
    if (argc < 9) {
      opserr << "WARNING insufficient arguments\n";
      printCommand(argc, argv);
      opserr << "Want: uniaxialMaterial Steel03 tag? fy? E0? b? r? cR1 cR2?";
      opserr << " <a1? a2? a3? a4?>" << endln;
      return TCL_ERROR;
    }

    int tag;

    if (Tcl_GetInt(interp, argv[2], &tag) != TCL_OK) {
      opserr << "WARNING invalid uniaxialMaterial Steel03 tag" << endln;
      return TCL_ERROR;
    }

    // Read required Steel01 material parameters
    double fy, E, b, r, cR1, cR2;

    if (Tcl_GetDouble(interp, argv[3], &fy) != TCL_OK) {
      opserr << "WARNING invalid fy\n";
      opserr << "uniaxialMaterial Steel03: " << tag << endln;
      return TCL_ERROR;
    }

    if (Tcl_GetDouble(interp, argv[4], &E) != TCL_OK) {
      opserr << "WARNING invalid E0\n";
      opserr << "uniaxialMaterial Steel03: " << tag << endln;
      return TCL_ERROR;
    }

    if (Tcl_GetDouble(interp, argv[5], &b) != TCL_OK) {
      opserr << "WARNING invalid b\n";
      opserr << "uniaxialMaterial Steel03: " << tag << endln;
      return TCL_ERROR;
    }

    if (Tcl_GetDouble(interp, argv[6], &r) != TCL_OK) {
      opserr << "WARNING invalid r\n";
      opserr << "uniaxialMaterial Steel03: " << tag << endln;
      return TCL_ERROR;
    }

    if (Tcl_GetDouble(interp, argv[7], &cR1) != TCL_OK) {
      opserr << "WARNING invalid cR1\n";
      opserr << "uniaxialMaterial Steel03: " << tag << endln;
      return TCL_ERROR;
    }

    if (Tcl_GetDouble(interp, argv[8], &cR2) != TCL_OK) {
      opserr << "WARNING invalid cR2\n";
      opserr << "uniaxialMaterial Steel03: " << tag << endln;
      return TCL_ERROR;
    }

    // Read optional Steel01 material parameters
    double a1, a2, a3, a4;
    if (argc > 9) {
      if (argc < 13) {
        opserr << "WARNING insufficient number of hardening parameters\n";
        opserr << "uniaxialMaterial Steel03: " << tag << endln;
        return TCL_ERROR;
      }

      if (Tcl_GetDouble(interp, argv[9], &a1) != TCL_OK) {
        opserr << "WARNING invalid a1\n";
        opserr << "uniaxialMaterial Steel03: " << tag << endln;
        return TCL_ERROR;
      }

      if (Tcl_GetDouble(interp, argv[10], &a2) != TCL_OK) {
        opserr << "WARNING invalid a2\n";
        opserr << "uniaxialMaterial Steel03: " << tag << endln;
        return TCL_ERROR;
      }

      if (Tcl_GetDouble(interp, argv[11], &a3) != TCL_OK) {
        opserr << "WARNING invalid a3\n";
        opserr << "uniaxialMaterial Steel03: " << tag << endln;
        return TCL_ERROR;
      }

      if (Tcl_GetDouble(interp, argv[12], &a4) != TCL_OK) {
        opserr << "WARNING invalid a4\n";
        opserr << "uniaxialMaterial Steel03: " << tag << endln;
        return TCL_ERROR;
      }

      // Parsing was successful, allocate the material
      theMaterial = new Steel03(tag, fy, E, b, r, cR1, cR2, a1, a2, a3, a4);
    } else
      // Parsing was successful, allocate the material
      theMaterial = new Steel03(tag, fy, E, b, r, cR1, cR2);

  }

  else if (strcmp(argv[1], "Hysteretic") == 0) {

    void *theMat = OPS_HystereticMaterial(rt, argc, argv);
    if (theMat != 0)
      theMaterial = (UniaxialMaterial *)theMat;
    else
      return TCL_ERROR;

  }

  else if (strcmp(argv[1], "OOHysteretic") == 0) {

    void *theMat = OPS_OOHystereticMaterial(rt, argc, argv);
    if (theMat != 0)
      theMaterial = (UniaxialMaterial *)theMat;
    else
      return TCL_ERROR;

  }

  else if (strcmp(argv[1], "Concrete04") == 0) {
    //        opserr << argc << endln;
    if (argc != 10 && argc != 9 && argc != 7) {
      opserr << "WARNING insufficient arguments\n";
      printCommand(argc, argv);
      opserr << "Want: uniaxialMaterial Concrete04 tag? fpc? epsc0? epscu? "
                "Ec0? <ft? etu? <beta?> >"
             << endln;
      return TCL_ERROR;
    }

    int tag;

    if (Tcl_GetInt(interp, argv[2], &tag) != TCL_OK) {
      opserr << "WARNING invalid uniaxialMaterial Concrete04 tag" << endln;
      return TCL_ERROR;
    }

    // Read required Concrete04 material parameters
    double fpc, epsc0, ft, epscu, Ec0, etu, beta;

    if (Tcl_GetDouble(interp, argv[3], &fpc) != TCL_OK) {
      opserr << "WARNING invalid fpc\n";
      opserr << "Concrete04 material: " << tag << endln;
      return TCL_ERROR;
    }

    if (Tcl_GetDouble(interp, argv[4], &epsc0) != TCL_OK) {
      opserr << "WARNING invalid epsc0\n";
      opserr << "Concrete04 material: " << tag << endln;
      return TCL_ERROR;
    }

    if (Tcl_GetDouble(interp, argv[5], &epscu) != TCL_OK) {
      opserr << "WARNING invalid epscu\n";
      opserr << "Concrete04 material: " << tag << endln;
      return TCL_ERROR;
    }

    if (Tcl_GetDouble(interp, argv[6], &Ec0) != TCL_OK) {
      opserr << "WARNING invalid Ec0\n";
      opserr << "Concrete04 material: " << tag << endln;
      return TCL_ERROR;
    }
    if (argc == 9 || argc == 10) {
      if (Tcl_GetDouble(interp, argv[7], &ft) != TCL_OK) {
        opserr << "WARNING invalid ft\n";
        opserr << "Concrete04 material: " << tag << endln;
        return TCL_ERROR;
      }
      if (Tcl_GetDouble(interp, argv[8], &etu) != TCL_OK) {
        opserr << "WARNING invalid etu\n";
        opserr << "Concrete04 material: " << tag << endln;
        return TCL_ERROR;
      }
    }
    if (argc == 10) {
      if (Tcl_GetDouble(interp, argv[9], &beta) != TCL_OK) {
        opserr << "WARNING invalid beta\n";
        opserr << "Concrete04 material: " << tag << endln;
        return TCL_ERROR;
      }
    }

    // Parsing was successful, allocate the material
    if (argc == 10) {
      theMaterial = new Concrete04(tag, fpc, epsc0, epscu, Ec0, ft, etu, beta);
    } else if (argc == 9) {
      theMaterial = new Concrete04(tag, fpc, epsc0, epscu, Ec0, ft, etu);
    } else if (argc == 7) {
      theMaterial = new Concrete04(tag, fpc, epsc0, epscu, Ec0);
    }
  }

  else if (strcmp(argv[1], "Concrete06") == 0) {
    if (argc < 12) {
      opserr << "WARNING insufficient arguments\n";
      printCommand(argc, argv);
      opserr << "Want: uniaxialMaterial Concrete06 tag? fc? eo? r? k? alphaC? "
                "fcr? ecr? b? alphaT?"
             << endln;
      return TCL_ERROR;
    }

    int tag;

    if (Tcl_GetInt(interp, argv[2], &tag) != TCL_OK) {
      opserr << "WARNING invalid uniaxialMaterial Concrete06 tag" << endln;
      return TCL_ERROR;
    }

    // Read required Concrete01 material parameters
    double fc, eo, r, k, fcr, ecr, b, alphaC, alphaT;

    if (Tcl_GetDouble(interp, argv[3], &fc) != TCL_OK) {
      opserr << "WARNING invalid fc\n";
      opserr << "Concrete06 material: " << tag << endln;
      return TCL_ERROR;
    }

    if (Tcl_GetDouble(interp, argv[4], &eo) != TCL_OK) {
      opserr << "WARNING invalid eo\n";
      opserr << "Concrete06 material: " << tag << endln;
      return TCL_ERROR;
    }

    if (Tcl_GetDouble(interp, argv[5], &r) != TCL_OK) {
      opserr << "WARNING invalid r\n";
      opserr << "Concrete06 material: " << tag << endln;
      return TCL_ERROR;
    }

    if (Tcl_GetDouble(interp, argv[6], &k) != TCL_OK) {
      opserr << "WARNING invalid k\n";
      opserr << "Concrete06 material: " << tag << endln;
      return TCL_ERROR;
    }

    if (Tcl_GetDouble(interp, argv[7], &alphaC) != TCL_OK) {
      opserr << "WARNING invalid alphaC\n";
      opserr << "Concrete06 material: " << tag << endln;
      return TCL_ERROR;
    }

    if (Tcl_GetDouble(interp, argv[8], &fcr) != TCL_OK) {
      opserr << "WARNING invalid fcr\n";
      opserr << "Concrete06 material: " << tag << endln;
      return TCL_ERROR;
    }

    if (Tcl_GetDouble(interp, argv[9], &ecr) != TCL_OK) {
      opserr << "WARNING invalid ecr\n";
      opserr << "Concrete06 material: " << tag << endln;
      return TCL_ERROR;
    }

    if (Tcl_GetDouble(interp, argv[10], &b) != TCL_OK) {
      opserr << "WARNING invalid b\n";
      opserr << "Concrete06 material: " << tag << endln;
      return TCL_ERROR;
    }

    if (Tcl_GetDouble(interp, argv[11], &alphaT) != TCL_OK) {
      opserr << "WARNING invalid alphaT\n";
      opserr << "Concrete06 material: " << tag << endln;
      return TCL_ERROR;
    }

    // Parsing was successful, allocate the material
    theMaterial =
        new Concrete06(tag, fc, eo, r, k, alphaC, fcr, ecr, b, alphaT);
  }

  else if (strcmp(argv[1], "Concrete07") == 0) {
    // Check to see if there are enough arquements
    if (argc < 11) {
      opserr << "WARNING: Insufficient arguments\n";
      printCommand(argc, argv);
      opserr << "Want: uniaxialMaterial Concrete07 tag? fpc? epsc0? Ec? fpt? "
                "epst0? xcrp? xcrn? r?\n";
      return TCL_ERROR;
    }

    int tag;

    if (Tcl_GetInt(interp, argv[2], &tag) != TCL_OK) {
      opserr << "WARNING: Invalid uniaxial Concrete07 tag\n";
      return TCL_ERROR;
    }

    // Read in the faluves required for the model
    double fpc, epsc0, Ec, fpt, epst0, xcrp, xcrn, r;

    if (Tcl_GetDouble(interp, argv[3], &fpc) != TCL_OK) {
      opserr << "WARNING: Invalid peak compression stress\n";
      opserr << "uniaxialMaterial Concrete07: " << tag << endln;
      return TCL_ERROR;
    }

    if (Tcl_GetDouble(interp, argv[4], &epsc0) != TCL_OK) {
      opserr << "WARNING: Invalid peak compression strain\n";
      opserr << "uniaxialMaterial Concrete07: " << tag << endln;
      return TCL_ERROR;
    }

    if (Tcl_GetDouble(interp, argv[5], &Ec) != TCL_OK) {
      opserr << "WARNING: Invalid Young's Modulus\n";
      opserr << "uniaxialMaterial Concrete07: " << tag << endln;
      return TCL_ERROR;
    }

    if (Tcl_GetDouble(interp, argv[6], &fpt) != TCL_OK) {
      opserr << "WARNING: Invalid peak tension stress\n";
      opserr << "uniaxialMaterial Concrete07: " << tag << endln;
      return TCL_ERROR;
    }

    if (Tcl_GetDouble(interp, argv[7], &epst0) != TCL_OK) {
      opserr << "WARNING: Invalid peak tension strain\n";
      opserr << "uniaxialMaterial Concrete07: " << tag << endln;
      return TCL_ERROR;
    }

    if (Tcl_GetDouble(interp, argv[8], &xcrp) != TCL_OK) {
      opserr << "WARNING: Invalid critical nondimensional strain in tension\n";
      opserr << "uniaxialMaterial Concrete07: " << tag << endln;
      return TCL_ERROR;
    }

    if (Tcl_GetDouble(interp, argv[9], &xcrn) != TCL_OK) {
      opserr
          << "WARNING: Invalid critical nondimensional strain in compression\n";
      opserr << "uniaxialMaterial Concrete07: " << tag << endln;
      return TCL_ERROR;
    }

    if (Tcl_GetDouble(interp, argv[10], &r) != TCL_OK) {
      opserr << "WARNING: Invalid value for r\n";
      opserr << "uniaxialMaterial Concrete07: " << tag << endln;
    }

    //		opserr << "fpc: " << fpc << endln << "epsc0: " << epsc0 << endln <<
    //"Ec: " << Ec << endln; 		opserr << "fpt: " << fpt << endln << "epst0: " <<
    //epst0 << endln << "xcrp: " << xcrp << endln; 		opserr << "xcrn: " << xcrn <<
    //endln << "r: " << r << endln;

    // Parsing was successful, allocate the material

    theMaterial =
        new Concrete07(tag, fpc, epsc0, Ec, fpt, epst0, xcrp, xcrn, r);
  }

  else if (strcmp(argv[1], "PathIndependent") == 0) {
    if (argc < 4) {
      opserr << "WARNING insufficient arguments\n";
      printCommand(argc, argv);
      opserr << "Want: uniaxialMaterial PathIndependent tag? matTag?" << endln;
      return TCL_ERROR;
    }

    int tag, matTag;

    if (Tcl_GetInt(interp, argv[2], &tag) != TCL_OK) {
      opserr << "WARNING invalid tag\n";
      opserr << "PathIndependent material: " << tag << endln;
      return TCL_ERROR;
    }

    if (Tcl_GetInt(interp, argv[3], &matTag) != TCL_OK) {
      opserr << "WARNING invalid matTag\n";
      opserr << "PathIndependent material: " << tag << endln;
      return TCL_ERROR;
    }

    UniaxialMaterial *material = builder->getTypedObject<UniaxialMaterial>(matTag);
    if (material == nullptr)
      return TCL_ERROR;

    theMaterial = new PathIndependentMaterial(tag, *material);
  }

  else if (strcmp(argv[1], "Backbone") == 0) {
    if (argc < 4) {
      opserr << "WARNING insufficient arguments\n";
      printCommand(argc, argv);
      opserr << "Want: uniaxialMaterial Backbone tag? bbTag?" << endln;
      return TCL_ERROR;
    }

    int tag, bbTag;

    if (Tcl_GetInt(interp, argv[2], &tag) != TCL_OK) {
      opserr << "WARNING invalid tag\n";
      opserr << "Backbone material: " << tag << endln;
      return TCL_ERROR;
    }

    if (Tcl_GetInt(interp, argv[3], &bbTag) != TCL_OK) {
      opserr << "WARNING invalid bTag\n";
      opserr << "Backbone material: " << tag << endln;
      return TCL_ERROR;
    }

    HystereticBackbone *backbone = OPS_getHystereticBackbone(bbTag);

    if (backbone == 0) {
      opserr << "WARNING backbone does not exist\n";
      opserr << "backbone: " << bbTag;
      opserr << "\nuniaxialMaterial Backbone: " << tag << endln;
      return TCL_ERROR;
    }

    theMaterial = new BackboneMaterial(tag, *backbone);
  }

  else if (strcmp(argv[1], "Fatigue") == 0) {
    if (argc < 4) {
      opserr << "WARNING insufficient arguments\n";
      printCommand(argc, argv);
      opserr << "Want: uniaxialMaterial Fatigue tag? matTag?";
      opserr << " <-D_max dmax?> <-e0 e0?> <-m m?>" << endln;
      opserr << " <-min min?> <-max max?>" << endln;
      return TCL_ERROR;
    }

    int tag, matTag;

    if (Tcl_GetInt(interp, argv[2], &tag) != TCL_OK) {
      opserr << "WARNING invalid uniaxialMaterial Fatigue tag" << endln;
      return TCL_ERROR;
    }

    if (Tcl_GetInt(interp, argv[3], &matTag) != TCL_OK) {
      opserr << "WARNING invalid component tag\n";
      opserr << "uniaxialMaterial Fatigue: " << tag << endln;
      return TCL_ERROR;
    }

    double Dmax = 1.0;
    double E0 = 0.191;
    double m = -0.458;
    double epsmin = NEG_INF_STRAIN;
    double epsmax = POS_INF_STRAIN;

    for (int j = 4; j < argc; j++) {
      if (strcmp(argv[j], "-Dmax") == 0) {
        if ((j + 1 >= argc) ||
            (Tcl_GetDouble(interp, argv[j + 1], &Dmax) != TCL_OK)) {
          opserr << "WARNING invalid -Dmax";
          opserr << "uniaxialMaterial Fatigue: " << tag << endln;
          return TCL_ERROR;
        }
      } else if (strcmp(argv[j], "-E0") == 0) {
        if ((j + 1 >= argc) ||
            (Tcl_GetDouble(interp, argv[j + 1], &E0) != TCL_OK)) {
          opserr << "WARNING invalid -E0";
          opserr << "uniaxialMaterial Fatigue: " << tag << endln;
          return TCL_ERROR;
        }
      } else if (strcmp(argv[j], "-m") == 0) {
        if ((j + 1 >= argc) ||
            (Tcl_GetDouble(interp, argv[j + 1], &m) != TCL_OK)) {
          opserr << "WARNING invalid -m";
          opserr << "uniaxialMaterial Fatigue: " << tag << endln;
          return TCL_ERROR;
        }
      } else if (strcmp(argv[j], "-min") == 0) {
        if ((j + 1 >= argc) ||
            (Tcl_GetDouble(interp, argv[j + 1], &epsmin) != TCL_OK)) {
          opserr << "WARNING invalid -min ";
          opserr << "uniaxialMaterial Fatigue: " << tag << endln;
          return TCL_ERROR;
        }
      } else if (strcmp(argv[j], "-max") == 0) {
        if ((j + 1 >= argc) ||
            (Tcl_GetDouble(interp, argv[j + 1], &epsmax) != TCL_OK)) {
          opserr << "WARNING invalid -max";
          opserr << "uniaxialMaterial Fatigue: " << tag << endln;
          return TCL_ERROR;
        }
      }
      j++;
    }

    UniaxialMaterial *theMat = builder->getTypedObject<UniaxialMaterial>(matTag);
    if (theMat == nullptr)
      return TCL_ERROR;

    // Parsing was successful, allocate the material
    theMaterial =
        new FatigueMaterial(tag, *theMat, Dmax, E0, m, epsmin, epsmax);

  }

  else if ((strcmp(argv[1], "SAWSMaterial") == 0) ||
           (strcmp(argv[1], "SAWS") == 0)) {

    void *theMat = OPS_SAWSMaterial(rt, argc, argv);
    if (theMat != 0)
      theMaterial = (UniaxialMaterial *)theMat;
    else
      return TCL_ERROR;

  }

  else if ((strcmp(argv[1], "BilinMaterial") == 0) ||
           (strcmp(argv[1], "Bilin") == 0)) {

    void *theMat = OPS_Bilin(rt, argc, argv);
    if (theMat != 0)
      theMaterial = (UniaxialMaterial *)theMat;
    else
      return TCL_ERROR;
  }

  else if ((strcmp(argv[1], "ConcreteZ01Material") == 0) ||
           (strcmp(argv[1], "ConcreteZ01") == 0)) {

    void *theMat = OPS_ConcreteZ01Material(rt, argc, argv);
    if (theMat != 0)
      theMaterial = (UniaxialMaterial *)theMat;
    else
      return TCL_ERROR;
  }

  else if ((strcmp(argv[1], "ConcreteL01Material") == 0) ||
           (strcmp(argv[1], "ConcreteL01") == 0)) {

    void *theMat = OPS_ConcreteL01Material(rt, argc, argv);
    if (theMat != 0)
      theMaterial = (UniaxialMaterial *)theMat;
    else
      return TCL_ERROR;
  }

  else if ((strcmp(argv[1], "SteelZ01Material") == 0) ||
           (strcmp(argv[1], "SteelZ01") == 0)) {

    void *theMat = OPS_SteelZ01Material(rt, argc, argv);
    if (theMat != 0)
      theMaterial = (UniaxialMaterial *)theMat;
    else
      return TCL_ERROR;
  }

  else if ((strcmp(argv[1], "TendonL01Material") == 0) ||
           (strcmp(argv[1], "TendonL01") == 0)) {
    void *theMat = OPS_TendonL01Material(rt, argc, argv);
    if (theMat != 0)
      theMaterial = (UniaxialMaterial *)theMat;
    else
      return TCL_ERROR;
  }


  else if (strcmp(argv[1], "Pinching4") == 0) {
    if (argc != 42 && argc != 31) {
      opserr << "WARNING insufficient arguments\n";
      printCommand(argc, argv);
      opserr << "Want: uniaxialMaterial Pinching4 tag? stress1p? strain1p? "
                "stress2p? strain2p? stress3p? strain3p? stress4p? strain4p? "
             << "\n<stress1n? strain1n? stress2n? strain2n? stress3n? "
                "strain3n? stress4n? strain4n?> rDispP? rForceP? uForceP? "
             << "\n<rDispN? rForceN? uForceN?> gammaK1? gammaK2? gammaK3? "
                "gammaK4? gammaKLimit? gammaD1? gammaD2? gammaD3? gammaD4? "
             << "\ngammaDLimit? gammaF1? gammaF2? gammaF3? gammaF4? "
                "gammaFLimit? gammaE? CycleOrEnergyDamage? ";
      return TCL_ERROR;
    }

    int tag, tDmg;
    double stress1p, stress2p, stress3p, stress4p;
    double strain1p, strain2p, strain3p, strain4p;
    double stress1n, stress2n, stress3n, stress4n;
    double strain1n, strain2n, strain3n, strain4n;
    double rDispP, rForceP, uForceP, rDispN, rForceN, uForceN;
    double gammaK1, gammaK2, gammaK3, gammaK4, gammaKLimit;
    double gammaD1, gammaD2, gammaD3, gammaD4, gammaDLimit;
    double gammaF1, gammaF2, gammaF3, gammaF4, gammaFLimit;
    double gammaE;

    int i = 2;

    if (Tcl_GetInt(interp, argv[i++], &tag) != TCL_OK) {
      opserr << "WARNING invalid uniaxialMaterial Pinching4 tag" << endln;
      return TCL_ERROR;
    }

    if (Tcl_GetDouble(interp, argv[i++], &stress1p) != TCL_OK) {
      opserr << "WARNING invalid stress1p\n";
      opserr << "Pinching4 material: " << tag << endln;
      return TCL_ERROR;
    }

    if (Tcl_GetDouble(interp, argv[i++], &strain1p) != TCL_OK) {
      opserr << "WARNING invalid strain1p\n";
      opserr << "Pinching4 material: " << tag << endln;
      return TCL_ERROR;
    }

    if (Tcl_GetDouble(interp, argv[i++], &stress2p) != TCL_OK) {
      opserr << "WARNING invalid stress2p\n";
      opserr << "Pinching4 material: " << tag << endln;
      return TCL_ERROR;
    }

    if (Tcl_GetDouble(interp, argv[i++], &strain2p) != TCL_OK) {
      opserr << "WARNING invalid strain2p\n";
      opserr << "Pinching4 material: " << tag << endln;
      return TCL_ERROR;
    }

    if (Tcl_GetDouble(interp, argv[i++], &stress3p) != TCL_OK) {
      opserr << "WARNING invalid stress3p\n";
      opserr << "Pinching4 material: " << tag << endln;
      return TCL_ERROR;
    }

    if (Tcl_GetDouble(interp, argv[i++], &strain3p) != TCL_OK) {
      opserr << "WARNING invalid strain3p\n";
      opserr << "Pinching4 material: " << tag << endln;
      return TCL_ERROR;
    }

    if (Tcl_GetDouble(interp, argv[i++], &stress4p) != TCL_OK) {
      opserr << "WARNING invalid stress4p\n";
      opserr << "Pinching4 material: " << tag << endln;
      return TCL_ERROR;
    }

    if (Tcl_GetDouble(interp, argv[i++], &strain4p) != TCL_OK) {
      opserr << "WARNING invalid strain4p\n";
      opserr << "Pinching4 material: " << tag << endln;
      return TCL_ERROR;
    }

    if (argc == 42) {
      if (Tcl_GetDouble(interp, argv[i++], &stress1n) != TCL_OK) {
        opserr << "WARNING invalid stress1n\n";
        opserr << "Pinching4 material: " << tag << endln;
        return TCL_ERROR;
      }

      if (Tcl_GetDouble(interp, argv[i++], &strain1n) != TCL_OK) {
        opserr << "WARNING invalid strain1n\n";
        opserr << "Pinching4 material: " << tag << endln;
        return TCL_ERROR;
      }

      if (Tcl_GetDouble(interp, argv[i++], &stress2n) != TCL_OK) {
        opserr << "WARNING invalid stress2n\n";
        opserr << "Pinching4 material: " << tag << endln;
        return TCL_ERROR;
      }

      if (Tcl_GetDouble(interp, argv[i++], &strain2n) != TCL_OK) {
        opserr << "WARNING invalid strain2n\n";
        opserr << "Pinching4 material: " << tag << endln;
        return TCL_ERROR;
      }

      if (Tcl_GetDouble(interp, argv[i++], &stress3n) != TCL_OK) {
        opserr << "WARNING invalid stress3n\n";
        opserr << "Pinching4 material: " << tag << endln;
        return TCL_ERROR;
      }

      if (Tcl_GetDouble(interp, argv[i++], &strain3n) != TCL_OK) {
        opserr << "WARNING invalid strain3n\n";
        opserr << "Pinching4 material: " << tag << endln;
        return TCL_ERROR;
      }

      if (Tcl_GetDouble(interp, argv[i++], &stress4n) != TCL_OK) {
        opserr << "WARNING invalid stress4n\n";
        opserr << "Pinching4 material: " << tag << endln;
        return TCL_ERROR;
      }

      if (Tcl_GetDouble(interp, argv[i++], &strain4n) != TCL_OK) {
        opserr << "WARNING invalid strain4n\n";
        opserr << "Pinching4 material: " << tag << endln;
        return TCL_ERROR;
      }
    }

    if (Tcl_GetDouble(interp, argv[i++], &rDispP) != TCL_OK) {
      opserr << "WARNING invalid rDispP\n";
      opserr << "Pinching4 material: " << tag << endln;
      return TCL_ERROR;
    }

    if (Tcl_GetDouble(interp, argv[i++], &rForceP) != TCL_OK) {
      opserr << "WARNING invalid rForceP\n";
      opserr << "Pinching4 material: " << tag << endln;
      return TCL_ERROR;
    }

    if (Tcl_GetDouble(interp, argv[i++], &uForceP) != TCL_OK) {
      opserr << "WARNING invalid uForceP\n";
      opserr << "Pinching4 material: " << tag << endln;
      return TCL_ERROR;
    }

    if (argc == 42) {
      if (Tcl_GetDouble(interp, argv[i++], &rDispN) != TCL_OK) {
        opserr << "WARNING invalid rDispN\n";
        opserr << "Pinching4 material: " << tag << endln;
        return TCL_ERROR;
      }

      if (Tcl_GetDouble(interp, argv[i++], &rForceN) != TCL_OK) {
        opserr << "WARNING invalid rForceN\n";
        opserr << "Pinching4 material: " << tag << endln;
        return TCL_ERROR;
      }

      if (Tcl_GetDouble(interp, argv[i++], &uForceN) != TCL_OK) {
        opserr << "WARNING invalid uForceN\n";
        opserr << "Pinching4 material: " << tag << endln;
        return TCL_ERROR;
      }
    }

    if (Tcl_GetDouble(interp, argv[i++], &gammaK1) != TCL_OK) {
      opserr << "WARNING invalid gammaK1\n";
      opserr << "Pinching4 material: " << tag << endln;
      return TCL_ERROR;
    }
    if (Tcl_GetDouble(interp, argv[i++], &gammaK2) != TCL_OK) {
      opserr << "WARNING invalid gammaK2\n";
      opserr << "Pinching4 material: " << tag << endln;
      return TCL_ERROR;
    }
    if (Tcl_GetDouble(interp, argv[i++], &gammaK3) != TCL_OK) {
      opserr << "WARNING invalid gammaK3\n";
      opserr << "Pinching4 material: " << tag << endln;
      return TCL_ERROR;
    }
    if (Tcl_GetDouble(interp, argv[i++], &gammaK4) != TCL_OK) {
      opserr << "WARNING invalid gammaK4\n";
      opserr << "Pinching4 material: " << tag << endln;
      return TCL_ERROR;
    }
    if (Tcl_GetDouble(interp, argv[i++], &gammaKLimit) != TCL_OK) {
      opserr << "WARNING invalid gammaKLimit\n";
      opserr << "Pinching4 material: " << tag << endln;
      return TCL_ERROR;
    }
    if (Tcl_GetDouble(interp, argv[i++], &gammaD1) != TCL_OK) {
      opserr << "WARNING invalid gammaD1\n";
      opserr << "Pinching4 material: " << tag << endln;
      return TCL_ERROR;
    }
    if (Tcl_GetDouble(interp, argv[i++], &gammaD2) != TCL_OK) {
      opserr << "WARNING invalid gammaD2\n";
      opserr << "Pinching4 material: " << tag << endln;
      return TCL_ERROR;
    }
    if (Tcl_GetDouble(interp, argv[i++], &gammaD3) != TCL_OK) {
      opserr << "WARNING invalid gammaD3\n";
      opserr << "Pinching4 material: " << tag << endln;
      return TCL_ERROR;
    }
    if (Tcl_GetDouble(interp, argv[i++], &gammaD4) != TCL_OK) {
      opserr << "WARNING invalid gammaD4\n";
      opserr << "Pinching4 material: " << tag << endln;
      return TCL_ERROR;
    }
    if (Tcl_GetDouble(interp, argv[i++], &gammaDLimit) != TCL_OK) {
      opserr << "WARNING invalid gammaDLimit\n";
      opserr << "Pinching4 material: " << tag << endln;
      return TCL_ERROR;
    }
    if (Tcl_GetDouble(interp, argv[i++], &gammaF1) != TCL_OK) {
      opserr << "WARNING invalid gammaF1\n";
      opserr << "Pinching4 material: " << tag << endln;
      return TCL_ERROR;
    }
    if (Tcl_GetDouble(interp, argv[i++], &gammaF2) != TCL_OK) {
      opserr << "WARNING invalid gammaF2\n";
      opserr << "Pinching4 material: " << tag << endln;
      return TCL_ERROR;
    }
    if (Tcl_GetDouble(interp, argv[i++], &gammaF3) != TCL_OK) {
      opserr << "WARNING invalid gammaF3\n";
      opserr << "Pinching4 material: " << tag << endln;
      return TCL_ERROR;
    }
    if (Tcl_GetDouble(interp, argv[i++], &gammaF4) != TCL_OK) {
      opserr << "WARNING invalid gammaF4\n";
      opserr << "Pinching4 material: " << tag << endln;
      return TCL_ERROR;
    }
    if (Tcl_GetDouble(interp, argv[i++], &gammaFLimit) != TCL_OK) {
      opserr << "WARNING invalid gammaFLimit\n";
      opserr << "Pinching4 material: " << tag << endln;
      return TCL_ERROR;
    }

    if (Tcl_GetDouble(interp, argv[i++], &gammaE) != TCL_OK) {
      opserr << "WARNING invalid gammaE\n";
      opserr << "Pinching4 material: " << tag << endln;
      return TCL_ERROR;
    }

    int y;
    y = i;

    if ((strcmp(argv[y], "cycle") == 0) || (strcmp(argv[y], "Cycle") == 0) ||
        (strcmp(argv[y], "DamageCycle") == 0) ||
        (strcmp(argv[y], "damageCycle") == 0)) {
      tDmg = 1;
    } else if ((strcmp(argv[y], "energy") == 0) ||
               (strcmp(argv[y], "Energy") == 0) ||
               (strcmp(argv[y], "DamageEnergy") == 0) ||
               (strcmp(argv[y], "damageEnergy") == 0)) {
      tDmg = 0;
    } else {
      opserr << "WARNING invalid type of damage calculation specified\n";
      opserr << "Pinching4 material: " << tag << endln;
      return TCL_ERROR;
    }

    // allocate the pinching material
    if (argc == 42) {
      theMaterial = new Pinching4Material(
          tag, stress1p, strain1p, stress2p, strain2p, stress3p, strain3p,
          stress4p, strain4p, stress1n, strain1n, stress2n, strain2n, stress3n,
          strain3n, stress4n, strain4n, rDispP, rForceP, uForceP, rDispN,
          rForceN, uForceN, gammaK1, gammaK2, gammaK3, gammaK4, gammaKLimit,
          gammaD1, gammaD2, gammaD3, gammaD4, gammaDLimit, gammaF1, gammaF2,
          gammaF3, gammaF4, gammaFLimit, gammaE, tDmg);
    }
    if (argc == 31) {
      theMaterial = new Pinching4Material(
          tag, stress1p, strain1p, stress2p, strain2p, stress3p, strain3p,
          stress4p, strain4p, rDispP, rForceP, uForceP, gammaK1, gammaK2,
          gammaK3, gammaK4, gammaKLimit, gammaD1, gammaD2, gammaD3, gammaD4,
          gammaDLimit, gammaF1, gammaF2, gammaF3, gammaF4, gammaFLimit, gammaE,
          tDmg);
    }
  }

  else if (strcmp(argv[1], "BarSlip") == 0) {
    if (argc != 17 && argc != 15) {
      opserr << "WARNING insufficient arguments\n";
      printCommand(argc, argv);
      opserr << "Want: uniaxialMaterial BarSlip tag? fc? fy? Es? fu? Eh? db? "
                "ld? nb? width? depth? bsflag? type? <damage? unit?>"
             << endln;
      return TCL_ERROR;
    }

    int tag, nb, bsf, typ = 2, dmg = 0, unt = 1;
    double fc, fy, Es, fu, Eh, ld, width, depth, db;

    int argStart = 2;

    if (Tcl_GetInt(interp, argv[argStart++], &tag) != TCL_OK) {
      opserr << "WARNING invalid tag\n";
      opserr << "BarSlip: " << tag << endln;
      return TCL_ERROR;
    }
    if (Tcl_GetDouble(interp, argv[argStart++], &fc) != TCL_OK) {
      opserr << "WARNING invalid fc\n";
      opserr << "BarSlip: " << tag << endln;
      return TCL_ERROR;
    }
    if (Tcl_GetDouble(interp, argv[argStart++], &fy) != TCL_OK) {
      opserr << "WARNING invalid fy\n";
      opserr << "BarSlip: " << tag << endln;
      return TCL_ERROR;
    }
    if (Tcl_GetDouble(interp, argv[argStart++], &Es) != TCL_OK) {
      opserr << "WARNING invalid Es\n";
      opserr << "BarSlip: " << tag << endln;
      return TCL_ERROR;
    }
    if (Tcl_GetDouble(interp, argv[argStart++], &fu) != TCL_OK) {
      opserr << "WARNING invalid fu\n";
      opserr << "BarSlip: " << tag << endln;
      return TCL_ERROR;
    }
    if (Tcl_GetDouble(interp, argv[argStart++], &Eh) != TCL_OK) {
      opserr << "WARNING invalid Eh\n";
      opserr << "BarSlip: " << tag << endln;
      return TCL_ERROR;
    }
    if (Tcl_GetDouble(interp, argv[argStart++], &db) != TCL_OK) {
      opserr << "WARNING invalid db\n";
      opserr << "BarSlip: " << tag << endln;
      return TCL_ERROR;
    }
    if (Tcl_GetDouble(interp, argv[argStart++], &ld) != TCL_OK) {
      opserr << "WARNING invalid ld\n";
      opserr << "BarSlip: " << tag << endln;
      return TCL_ERROR;
    }
    if (Tcl_GetInt(interp, argv[argStart++], &nb) != TCL_OK) {
      opserr << "WARNING invalid nbars\n";
      opserr << "BarSlip: " << tag << endln;
      return TCL_ERROR;
    }
    if (Tcl_GetDouble(interp, argv[argStart++], &width) != TCL_OK) {
      opserr << "WARNING invalid width\n";
      opserr << "BarSlip: " << tag << endln;
      return TCL_ERROR;
    }
    if (Tcl_GetDouble(interp, argv[argStart++], &depth) != TCL_OK) {
      opserr << "WARNING invalid depth\n";
      opserr << "BarSlip: " << tag << endln;
      return TCL_ERROR;
    }

    int y;
    y = argStart;

    if ((strcmp(argv[y], "strong") == 0) || (strcmp(argv[y], "Strong") == 0) ||
        (strcmp(argv[y], "weak") == 0) || (strcmp(argv[y], "Weak") == 0)) {
      if ((strcmp(argv[y], "strong") == 0) ||
          (strcmp(argv[y], "Strong") == 0)) {
        bsf = 0;
      }

      if ((strcmp(argv[y], "weak") == 0) || (strcmp(argv[y], "Weak") == 0)) {
        bsf = 1;
      }
    } else {
      opserr << "WARNING invalid bond strength specified\n";
      opserr << "BarSlip: " << tag << endln;
      return TCL_ERROR;
    }
    y++;

    if ((strcmp(argv[y], "beamtop") == 0) ||
        (strcmp(argv[y], "beamTop") == 0) ||
        (strcmp(argv[y], "beambot") == 0) ||
        (strcmp(argv[y], "beamBot") == 0) ||
        (strcmp(argv[y], "beambottom") == 0) ||
        (strcmp(argv[y], "beamBottom") == 0) ||
        (strcmp(argv[y], "beam") == 0) || (strcmp(argv[y], "Beam") == 0) ||
        (strcmp(argv[y], "Column") == 0) || (strcmp(argv[y], "column") == 0)) {
      if ((strcmp(argv[y], "beamtop") == 0) ||
          (strcmp(argv[y], "beamTop") == 0) || (strcmp(argv[y], "beam") == 0) ||
          (strcmp(argv[y], "Beam") == 0)) {
        typ = 0;
      }

      if ((strcmp(argv[y], "beambot") == 0) ||
          (strcmp(argv[y], "beamBot") == 0) ||
          (strcmp(argv[y], "beambottom") == 0) ||
          (strcmp(argv[y], "beamBottom") == 0)) {
        typ = 1;
      }

      if ((strcmp(argv[y], "column") == 0) ||
          (strcmp(argv[y], "Column") == 0)) {
        typ = 2;
      }
    } else {
      opserr << "WARNING invalid location of bar specified\n";
      opserr << "BarSlip: " << tag << endln;
      return TCL_ERROR;
    }
    if (argc == 17) {
      y++;

      if ((strcmp(argv[y], "damage1") == 0) ||
          (strcmp(argv[y], "Damage1") == 0) ||
          (strcmp(argv[y], "damage2") == 0) ||
          (strcmp(argv[y], "Damage2") == 0) ||
          (strcmp(argv[y], "nodamage") == 0) ||
          (strcmp(argv[y], "Nodamage") == 0) ||
          (strcmp(argv[y], "NoDamage") == 0) ||
          (strcmp(argv[y], "noDamage") == 0)) {
        if ((strcmp(argv[y], "damage1") == 0) ||
            (strcmp(argv[y], "Damage1") == 0)) {
          dmg = 1;
        } else if ((strcmp(argv[y], "damage2") == 0) ||
                   (strcmp(argv[y], "Damage2") == 0)) {
          dmg = 2;
        } else if ((strcmp(argv[y], "nodamage") == 0) ||
                   (strcmp(argv[y], "Nodamage") == 0) ||
                   (strcmp(argv[y], "NoDamage") == 0) ||
                   (strcmp(argv[y], "noDamage") == 0)) {
          dmg = 0;
        }

      } else {
        opserr << "WARNING invalid damage specified\n";
        opserr << "BarSlip: " << tag << endln;
        return TCL_ERROR;
      }

      y++;

      if ((strcmp(argv[y], "mpa") == 0) || (strcmp(argv[y], "MPa") == 0) ||
          (strcmp(argv[y], "mPa") == 0) || (strcmp(argv[y], "Mpa") == 0) ||
          (strcmp(argv[y], "psi") == 0) || (strcmp(argv[y], "Psi") == 0) ||
          (strcmp(argv[y], "PSI") == 0) || (strcmp(argv[y], "Pa") == 0) ||
          (strcmp(argv[y], "pa") == 0) || (strcmp(argv[y], "psf") == 0) ||
          (strcmp(argv[y], "Psf") == 0) || (strcmp(argv[y], "PSF") == 0) ||
          (strcmp(argv[y], "ksi") == 0) || (strcmp(argv[y], "Ksi") == 0) ||
          (strcmp(argv[y], "KSI") == 0) || (strcmp(argv[y], "ksf") == 0) ||
          (strcmp(argv[y], "Ksf") == 0) || (strcmp(argv[y], "KSF") == 0)) {
        if ((strcmp(argv[y], "mpa") == 0) || (strcmp(argv[y], "MPa") == 0) ||
            (strcmp(argv[y], "mPa") == 0) || (strcmp(argv[y], "Mpa") == 0)) {
          unt = 1;
        } else if ((strcmp(argv[y], "psi") == 0) ||
                   (strcmp(argv[y], "Psi") == 0) ||
                   (strcmp(argv[y], "PSI") == 0)) {
          unt = 2;
        } else if ((strcmp(argv[y], "Pa") == 0) ||
                   (strcmp(argv[y], "pa") == 0)) {
          unt = 3;
        } else if ((strcmp(argv[y], "psf") == 0) ||
                   (strcmp(argv[y], "Psf") == 0) ||
                   (strcmp(argv[y], "PSF") == 0)) {
          unt = 4;
        } else if ((strcmp(argv[y], "ksi") == 0) ||
                   (strcmp(argv[y], "Ksi") == 0) ||
                   (strcmp(argv[y], "KSI") == 0)) {
          unt = 5;
        } else if ((strcmp(argv[y], "ksf") == 0) ||
                   (strcmp(argv[y], "Ksf") == 0) ||
                   (strcmp(argv[y], "KSF") == 0)) {
          unt = 6;
        }
      } else {
        opserr << "WARNING invalid unit specified\n";
        opserr << "BarSlip: " << tag << endln;
        return TCL_ERROR;
      }
    }

    // allocate the material
    if (argc == 15) {
      theMaterial = new BarSlipMaterial(tag, fc, fy, Es, fu, Eh, db, ld, nb,
                                        width, depth, bsf, typ);
    }

    if (argc == 17) {
      theMaterial = new BarSlipMaterial(tag, fc, fy, Es, fu, Eh, db, ld, nb,
                                        width, depth, bsf, typ, dmg, unt);
    }

  }

  else if (strcmp(argv[1], "ShearPanel") == 0) {
    if (argc != 42 && argc != 31) {
      opserr << "WARNING insufficient arguments\n";
      printCommand(argc, argv);
      opserr << "Want: uniaxialMaterial ShearPanel tag? stress1p? strain1p? "
                "stress2p? strain2p? stress3p? strain3p? stress4p? strain4p? "
             << "\n<stress1n? strain1n? stress2n? strain2n? stress3n? "
                "strain3n? stress4n? strain4n?> rDispP? rForceP? uForceP? "
             << "\n<rDispN? rForceN? uForceN?> gammaK1? gammaK2? gammaK3? "
                "gammaK4? gammaKLimit? gammaD1? gammaD2? gammaD3? gammaD4? "
             << "\ngammaDLimit? gammaF1? gammaF2? gammaF3? gammaF4? "
                "gammaFLimit? gammaE? YieldStress? ";
      return TCL_ERROR;
    }

    int tag;
    double stress1p, stress2p, stress3p, stress4p;
    double strain1p, strain2p, strain3p, strain4p;
    double stress1n, stress2n, stress3n, stress4n;
    double strain1n, strain2n, strain3n, strain4n;
    double rDispP, rForceP, uForceP, rDispN, rForceN, uForceN;
    double gammaK1, gammaK2, gammaK3, gammaK4, gammaKLimit;
    double gammaD1, gammaD2, gammaD3, gammaD4, gammaDLimit;
    double gammaF1, gammaF2, gammaF3, gammaF4, gammaFLimit;
    double gammaE, yStr;

    int i = 2;

    if (Tcl_GetInt(interp, argv[i++], &tag) != TCL_OK) {
      opserr << "WARNING invalid uniaxialMaterial ShearPanel tag" << endln;
      return TCL_ERROR;
    }

    if (Tcl_GetDouble(interp, argv[i++], &stress1p) != TCL_OK) {
      opserr << "WARNING invalid stress1p\n";
      opserr << "ShearPanel material: " << tag << endln;
      return TCL_ERROR;
    }

    if (Tcl_GetDouble(interp, argv[i++], &strain1p) != TCL_OK) {
      opserr << "WARNING invalid strain1p\n";
      opserr << "ShearPanel material: " << tag << endln;
      return TCL_ERROR;
    }

    if (Tcl_GetDouble(interp, argv[i++], &stress2p) != TCL_OK) {
      opserr << "WARNING invalid stress2p\n";
      opserr << "ShearPanel material: " << tag << endln;
      return TCL_ERROR;
    }

    if (Tcl_GetDouble(interp, argv[i++], &strain2p) != TCL_OK) {
      opserr << "WARNING invalid strain2p\n";
      opserr << "ShearPanel material: " << tag << endln;
      return TCL_ERROR;
    }

    if (Tcl_GetDouble(interp, argv[i++], &stress3p) != TCL_OK) {
      opserr << "WARNING invalid stress3p\n";
      opserr << "ShearPanel material: " << tag << endln;
      return TCL_ERROR;
    }

    if (Tcl_GetDouble(interp, argv[i++], &strain3p) != TCL_OK) {
      opserr << "WARNING invalid strain3p\n";
      opserr << "ShearPanel material: " << tag << endln;
      return TCL_ERROR;
    }

    if (Tcl_GetDouble(interp, argv[i++], &stress4p) != TCL_OK) {
      opserr << "WARNING invalid stress4p\n";
      opserr << "ShearPanel material: " << tag << endln;
      return TCL_ERROR;
    }

    if (Tcl_GetDouble(interp, argv[i++], &strain4p) != TCL_OK) {
      opserr << "WARNING invalid strain4p\n";
      opserr << "ShearPanel material: " << tag << endln;
      return TCL_ERROR;
    }

    if (argc == 42) {
      if (Tcl_GetDouble(interp, argv[i++], &stress1n) != TCL_OK) {
        opserr << "WARNING invalid stress1n\n";
        opserr << "ShearPanel material: " << tag << endln;
        return TCL_ERROR;
      }

      if (Tcl_GetDouble(interp, argv[i++], &strain1n) != TCL_OK) {
        opserr << "WARNING invalid strain1n\n";
        opserr << "ShearPanel material: " << tag << endln;
        return TCL_ERROR;
      }

      if (Tcl_GetDouble(interp, argv[i++], &stress2n) != TCL_OK) {
        opserr << "WARNING invalid stress2n\n";
        opserr << "ShearPanel material: " << tag << endln;
        return TCL_ERROR;
      }

      if (Tcl_GetDouble(interp, argv[i++], &strain2n) != TCL_OK) {
        opserr << "WARNING invalid strain2n\n";
        opserr << "ShearPanel material: " << tag << endln;
        return TCL_ERROR;
      }

      if (Tcl_GetDouble(interp, argv[i++], &stress3n) != TCL_OK) {
        opserr << "WARNING invalid stress3n\n";
        opserr << "ShearPanel material: " << tag << endln;
        return TCL_ERROR;
      }

      if (Tcl_GetDouble(interp, argv[i++], &strain3n) != TCL_OK) {
        opserr << "WARNING invalid strain3n\n";
        opserr << "ShearPanel material: " << tag << endln;
        return TCL_ERROR;
      }

      if (Tcl_GetDouble(interp, argv[i++], &stress4n) != TCL_OK) {
        opserr << "WARNING invalid stress4n\n";
        opserr << "ShearPanel material: " << tag << endln;
        return TCL_ERROR;
      }

      if (Tcl_GetDouble(interp, argv[i++], &strain4n) != TCL_OK) {
        opserr << "WARNING invalid strain4n\n";
        opserr << "ShearPanel material: " << tag << endln;
        return TCL_ERROR;
      }
    }

    if (Tcl_GetDouble(interp, argv[i++], &rDispP) != TCL_OK) {
      opserr << "WARNING invalid rDispP\n";
      opserr << "ShearPanel material: " << tag << endln;
      return TCL_ERROR;
    }

    if (Tcl_GetDouble(interp, argv[i++], &rForceP) != TCL_OK) {
      opserr << "WARNING invalid rForceP\n";
      opserr << "ShearPanel material: " << tag << endln;
      return TCL_ERROR;
    }

    if (Tcl_GetDouble(interp, argv[i++], &uForceP) != TCL_OK) {
      opserr << "WARNING invalid uForceP\n";
      opserr << "ShearPanel material: " << tag << endln;
      return TCL_ERROR;
    }

    if (argc == 42) {
      if (Tcl_GetDouble(interp, argv[i++], &rDispN) != TCL_OK) {
        opserr << "WARNING invalid rDispN\n";
        opserr << "ShearPanel material: " << tag << endln;
        return TCL_ERROR;
      }

      if (Tcl_GetDouble(interp, argv[i++], &rForceN) != TCL_OK) {
        opserr << "WARNING invalid rForceN\n";
        opserr << "ShearPanel material: " << tag << endln;
        return TCL_ERROR;
      }

      if (Tcl_GetDouble(interp, argv[i++], &uForceN) != TCL_OK) {
        opserr << "WARNING invalid uForceN\n";
        opserr << "ShearPanel material: " << tag << endln;
        return TCL_ERROR;
      }
    }

    if (Tcl_GetDouble(interp, argv[i++], &gammaK1) != TCL_OK) {
      opserr << "WARNING invalid gammaK1\n";
      opserr << "ShearPanel material: " << tag << endln;
      return TCL_ERROR;
    }
    if (Tcl_GetDouble(interp, argv[i++], &gammaK2) != TCL_OK) {
      opserr << "WARNING invalid gammaK2\n";
      opserr << "ShearPanel material: " << tag << endln;
      return TCL_ERROR;
    }
    if (Tcl_GetDouble(interp, argv[i++], &gammaK3) != TCL_OK) {
      opserr << "WARNING invalid gammaK3\n";
      opserr << "ShearPanel material: " << tag << endln;
      return TCL_ERROR;
    }
    if (Tcl_GetDouble(interp, argv[i++], &gammaK4) != TCL_OK) {
      opserr << "WARNING invalid gammaK4\n";
      opserr << "ShearPanel material: " << tag << endln;
      return TCL_ERROR;
    }
    if (Tcl_GetDouble(interp, argv[i++], &gammaKLimit) != TCL_OK) {
      opserr << "WARNING invalid gammaKLimit\n";
      opserr << "ShearPanel material: " << tag << endln;
      return TCL_ERROR;
    }
    if (Tcl_GetDouble(interp, argv[i++], &gammaD1) != TCL_OK) {
      opserr << "WARNING invalid gammaD1\n";
      opserr << "ShearPanel material: " << tag << endln;
      return TCL_ERROR;
    }
    if (Tcl_GetDouble(interp, argv[i++], &gammaD2) != TCL_OK) {
      opserr << "WARNING invalid gammaD2\n";
      opserr << "ShearPanel material: " << tag << endln;
      return TCL_ERROR;
    }
    if (Tcl_GetDouble(interp, argv[i++], &gammaD3) != TCL_OK) {
      opserr << "WARNING invalid gammaD3\n";
      opserr << "ShearPanel material: " << tag << endln;
      return TCL_ERROR;
    }
    if (Tcl_GetDouble(interp, argv[i++], &gammaD4) != TCL_OK) {
      opserr << "WARNING invalid gammaD4\n";
      opserr << "ShearPanel material: " << tag << endln;
      return TCL_ERROR;
    }
    if (Tcl_GetDouble(interp, argv[i++], &gammaDLimit) != TCL_OK) {
      opserr << "WARNING invalid gammaDLimit\n";
      opserr << "ShearPanel material: " << tag << endln;
      return TCL_ERROR;
    }
    if (Tcl_GetDouble(interp, argv[i++], &gammaF1) != TCL_OK) {
      opserr << "WARNING invalid gammaF1\n";
      opserr << "ShearPanel material: " << tag << endln;
      return TCL_ERROR;
    }
    if (Tcl_GetDouble(interp, argv[i++], &gammaF2) != TCL_OK) {
      opserr << "WARNING invalid gammaF2\n";
      opserr << "ShearPanel material: " << tag << endln;
      return TCL_ERROR;
    }
    if (Tcl_GetDouble(interp, argv[i++], &gammaF3) != TCL_OK) {
      opserr << "WARNING invalid gammaF3\n";
      opserr << "ShearPanel material: " << tag << endln;
      return TCL_ERROR;
    }
    if (Tcl_GetDouble(interp, argv[i++], &gammaF4) != TCL_OK) {
      opserr << "WARNING invalid gammaF4\n";
      opserr << "ShearPanel material: " << tag << endln;
      return TCL_ERROR;
    }
    if (Tcl_GetDouble(interp, argv[i++], &gammaFLimit) != TCL_OK) {
      opserr << "WARNING invalid gammaFLimit\n";
      opserr << "ShearPanel material: " << tag << endln;
      return TCL_ERROR;
    }

    if (Tcl_GetDouble(interp, argv[i++], &gammaE) != TCL_OK) {
      opserr << "WARNING invalid gammaE\n";
      opserr << "ShearPanel material: " << tag << endln;
      return TCL_ERROR;
    }

    if (Tcl_GetDouble(interp, argv[i++], &yStr) != TCL_OK) {
      opserr << "WARNING invalid yield stress\n";
      opserr << "ShearPanel material: " << tag << endln;
      return TCL_ERROR;
    }

    // allocate the pinching material
    if (argc == 42) {
      theMaterial = new ShearPanelMaterial(
          tag, stress1p, strain1p, stress2p, strain2p, stress3p, strain3p,
          stress4p, strain4p, stress1n, strain1n, stress2n, strain2n, stress3n,
          strain3n, stress4n, strain4n, rDispP, rForceP, uForceP, rDispN,
          rForceN, uForceN, gammaK1, gammaK2, gammaK3, gammaK4, gammaKLimit,
          gammaD1, gammaD2, gammaD3, gammaD4, gammaDLimit, gammaF1, gammaF2,
          gammaF3, gammaF4, gammaFLimit, gammaE, yStr);
    }
    if (argc == 31) {
      theMaterial = new ShearPanelMaterial(
          tag, stress1p, strain1p, stress2p, strain2p, stress3p, strain3p,
          stress4p, strain4p, rDispP, rForceP, uForceP, gammaK1, gammaK2,
          gammaK3, gammaK4, gammaKLimit, gammaD1, gammaD2, gammaD3, gammaD4,
          gammaDLimit, gammaF1, gammaF2, gammaF3, gammaF4, gammaFLimit, gammaE,
          yStr);
    }
  }

  else if (strcmp(argv[1], "Concrete01WithSITC") == 0) {
    if (argc < 7) {
      opserr << "WARNING insufficient arguments\n";
      printCommand(argc, argv);
      opserr << "Want: uniaxialMaterial Concrete01 tag? fpc? epsc0? fpcu? "
                "epscu? <endStrainSITC?>"
             << endln;
      return TCL_ERROR;
    }

    int tag;

    if (Tcl_GetInt(interp, argv[2], &tag) != TCL_OK) {
      opserr << "WARNING invalid uniaxialMaterial Concrete01 tag" << endln;
      return TCL_ERROR;
    }

    // Read required Concrete01 material parameters
    double fpc, epsc0, fpcu, epscu;

    if (Tcl_GetDouble(interp, argv[3], &fpc) != TCL_OK) {
      opserr << "WARNING invalid fpc\n";
      opserr << "Concrete01 material: " << tag << endln;
      return TCL_ERROR;
    }

    if (Tcl_GetDouble(interp, argv[4], &epsc0) != TCL_OK) {
      opserr << "WARNING invalid epsc0\n";
      opserr << "Concrete01 material: " << tag << endln;
      return TCL_ERROR;
    }

    if (Tcl_GetDouble(interp, argv[5], &fpcu) != TCL_OK) {
      opserr << "WARNING invalid fpcu\n";
      opserr << "Concrete01 material: " << tag << endln;
      return TCL_ERROR;
    }

    if (Tcl_GetDouble(interp, argv[6], &epscu) != TCL_OK) {
      opserr << "WARNING invalid epscu\n";
      opserr << "Concrete01 material: " << tag << endln;
      return TCL_ERROR;
    }

    if (argc == 7)
      // Parsing was successful, allocate the material
      theMaterial = new Concrete01WithSITC(tag, fpc, epsc0, fpcu, epscu);
    else {
      double endStrainSITC;
      if (Tcl_GetDouble(interp, argv[7], &endStrainSITC) != TCL_OK) {
        opserr << "WARNING invalid epscu\n";
        opserr << "Concrete01 material: " << tag << endln;
        return TCL_ERROR;
      }
      theMaterial =
          new Concrete01WithSITC(tag, fpc, epsc0, fpcu, epscu, endStrainSITC);
    }
  }

  else if (strcmp(argv[1], "ECC01") == 0) {
    if (argc < 16) {
      opserr << "WARNING insufficient arguments\n";
      printCommand(argc, argv);
      opserr << "Want: uniaxialMaterial ECC01 TAG? SIGT0? EPST0? SIGT1? EPST1? "
                "EPST2? SIGC0? EPSC0? EPSC1? ";
      opserr << "ALPHAT1? ALPHAT2? ALPHAC? ALPHACU? BETAT? BETAC\n";
      return TCL_ERROR;
    }

    int tag;
    double SIGT0, EPST0, SIGT1, EPST1, EPST2, SIGC0, EPSC0, EPSC1, ALPHAT1,
        ALPHAT2, ALPHAC, ALPHACU, BETAT, BETAC;

    if (Tcl_GetInt(interp, argv[2], &tag) != TCL_OK) {
      opserr << "WARNING invalid uniaxialMaterial ECC01 tag" << endln;
      return TCL_ERROR;
    }

    if (Tcl_GetDouble(interp, argv[3], &SIGT0) != TCL_OK) {
      opserr << "WARNING invalid SIGTO\n";
      opserr << "ECC01 material: " << tag << endln;
      return TCL_ERROR;
    }

    if (Tcl_GetDouble(interp, argv[4], &EPST0) != TCL_OK) {
      opserr << "WARNING invalid EPSTO\n";
      opserr << "ECC01 material: " << tag << endln;
      return TCL_ERROR;
    }

    if (Tcl_GetDouble(interp, argv[5], &SIGT1) != TCL_OK) {
      opserr << "WARNING invalid SIGT1\n";
      opserr << "ECC01 material: " << tag << endln;
      return TCL_ERROR;
    }

    if (Tcl_GetDouble(interp, argv[6], &EPST1) != TCL_OK) {
      opserr << "WARNING invalid EPST1\n";
      opserr << "ECC01 material: " << tag << endln;
      return TCL_ERROR;
    }

    if (Tcl_GetDouble(interp, argv[7], &EPST2) != TCL_OK) {
      opserr << "WARNING invalid epscu\n";
      opserr << "ECC01 material: " << tag << endln;
      return TCL_ERROR;
    }
    if (Tcl_GetDouble(interp, argv[8], &SIGC0) != TCL_OK) {
      opserr << "WARNING invalid epscu\n";
      opserr << "ECC01 material: " << tag << endln;
      return TCL_ERROR;
    }
    if (Tcl_GetDouble(interp, argv[9], &EPSC0) != TCL_OK) {
      opserr << "WARNING invalid epscu\n";
      opserr << "ECC01 material: " << tag << endln;
      return TCL_ERROR;
    }
    if (Tcl_GetDouble(interp, argv[10], &EPSC1) != TCL_OK) {
      opserr << "WARNING invalid epscu\n";
      opserr << "ECC01 material: " << tag << endln;
      return TCL_ERROR;
    }
    if (Tcl_GetDouble(interp, argv[11], &ALPHAT1) != TCL_OK) {
      opserr << "WARNING invalid epscu\n";
      opserr << "ECC01 material: " << tag << endln;
      return TCL_ERROR;
    }
    if (Tcl_GetDouble(interp, argv[12], &ALPHAT2) != TCL_OK) {
      opserr << "WARNING invalid epscu\n";
      opserr << "ECC01 material: " << tag << endln;
      return TCL_ERROR;
    }
    if (Tcl_GetDouble(interp, argv[13], &ALPHAC) != TCL_OK) {
      opserr << "WARNING invalid epscu\n";
      opserr << "ECC01 material: " << tag << endln;
      return TCL_ERROR;
    }
    if (Tcl_GetDouble(interp, argv[14], &ALPHACU) != TCL_OK) {
      opserr << "WARNING invalid epscu\n";
      opserr << "ECC01 material: " << tag << endln;
      return TCL_ERROR;
    }
    if (Tcl_GetDouble(interp, argv[15], &BETAT) != TCL_OK) {
      opserr << "WARNING invalid epscu\n";
      opserr << "ECC01 material: " << tag << endln;
      return TCL_ERROR;
    }
    if (Tcl_GetDouble(interp, argv[16], &BETAC) != TCL_OK) {
      opserr << "WARNING invalid epscu\n";
      opserr << "ECC01 material: " << tag << endln;
      return TCL_ERROR;
    }

    theMaterial =
        new ECC01(tag, SIGT0, EPST0, SIGT1, EPST1, EPST2, SIGC0, EPSC0, EPSC1,
                  ALPHAT1, ALPHAT2, ALPHAC, ALPHACU, BETAT, BETAC);
  }

  else if (strcmp(argv[1], "ASD_SMA_3K") == 0) {
    void *theMat = OPS_ASD_SMA_3K(rt, argc, argv);
    if (theMat != 0)
      theMaterial = (UniaxialMaterial *)theMat;
    else
      return TCL_ERROR;
  }

  else if (strcmp(argv[1], "SelfCentering") == 0) {
    if (argc < 7) {
      opserr << "WARNING insufficient arguments\n";
      printCommand(argc, argv);
      opserr << "Want: uniaxialMaterial SelfCentering tag? k1? k2? ActF? beta? "
                "<SlipDef? BearDef? rBear?>"
             << endln;
      return TCL_ERROR;
    }

    int tag;
    double k1, k2, ActF, beta, rBear, SlipDef, BearDef;

    if (Tcl_GetInt(interp, argv[2], &tag) != TCL_OK) {
      opserr << "WARNING invalid uniaxialMaterial SelfCentering tag" << endln;
      return TCL_ERROR;
    }

    if (Tcl_GetDouble(interp, argv[3], &k1) != TCL_OK) {
      opserr << "WARNING invalid k1\n";
      opserr << "uniaxialMaterial SelfCentering: " << tag << endln;
      return TCL_ERROR;
    }

    if (Tcl_GetDouble(interp, argv[4], &k2) != TCL_OK) {
      opserr << "WARNING invalid k2\n";
      opserr << "uniaxialMaterial SelfCentering: " << tag << endln;
      return TCL_ERROR;
    }

    if (Tcl_GetDouble(interp, argv[5], &ActF) != TCL_OK) {
      opserr << "WARNING invalid ActF\n";
      opserr << "uniaxialMaterial SelfCentering: " << tag << endln;
      return TCL_ERROR;
    }

    if (Tcl_GetDouble(interp, argv[6], &beta) != TCL_OK) {
      opserr << "WARNING invalid beta\n";
      opserr << "uniaxialMaterial SelfCentering: " << tag << endln;
      return TCL_ERROR;
    }

    if (argc == 8) {
      if (Tcl_GetDouble(interp, argv[7], &SlipDef) != TCL_OK) {
        opserr << "WARNING invalid SlipDef\n";
        opserr << "uniaxialMaterial SelfCentering: " << tag << endln;
        return TCL_ERROR;
      }
      // Parsing was successful, allocate the material
      theMaterial =
          new SelfCenteringMaterial(tag, k1, k2, ActF, beta, SlipDef, 0, 0);
    }

    else if (argc > 8) {
      if (Tcl_GetDouble(interp, argv[7], &SlipDef) != TCL_OK) {
        opserr << "WARNING invalid SlipDef\n";
        opserr << "uniaxialMaterial SelfCentering: " << tag << endln;
        return TCL_ERROR;
      }
      if (Tcl_GetDouble(interp, argv[8], &BearDef) != TCL_OK) {
        opserr << "WARNING invalid BearDef\n";
        opserr << "uniaxialMaterial SelfCentering: " << tag << endln;
        return TCL_ERROR;
      }
      if (Tcl_GetDouble(interp, argv[9], &rBear) != TCL_OK) {
        opserr << "WARNING invalid rBear\n";
        opserr << "uniaxialMaterial SelfCentering: " << tag << endln;
        return TCL_ERROR;
      }
      // Parsing was successful, allocate the material
      theMaterial = new SelfCenteringMaterial(tag, k1, k2, ActF, beta, SlipDef,
                                              BearDef, rBear);
    }

    else {
      // Parsing was successful, allocate the material
      theMaterial = new SelfCenteringMaterial(tag, k1, k2, ActF, beta, 0, 0, 0);
    }
  }

  else if (strcmp(argv[1], "SteelMP") == 0) {
    // Check that there is the minimum number of arguments
    if (argc < 4) {
      opserr << "WARNING insufficient arguments\n";
      printCommand(argc, argv);
      opserr << "Want: uniaxialMaterial SteelMP tag? fy? E0? b? ";
      opserr << " <coeffR1?  coeffR2? a1? a2?>" << endln;
      return TCL_ERROR;
    }

    int tag;

    if (Tcl_GetInt(interp, argv[2], &tag) != TCL_OK) {
      opserr << "WARNING invalid uniaxialMaterial SteelMP tag" << endln;
      return TCL_ERROR;
    }

    // Read required Steel01 material parameters
    double fy, E, b;

    if (Tcl_GetDouble(interp, argv[3], &fy) != TCL_OK) {
      opserr << "WARNING invalid fy\n";
      opserr << "uniaxialMaterial SteelMP: " << tag << endln;
      return TCL_ERROR;
    }

    if (Tcl_GetDouble(interp, argv[4], &E) != TCL_OK) {
      opserr << "WARNING invalid E0\n";
      opserr << "uniaxialMaterial SteelMP: " << tag << endln;
      return TCL_ERROR;
    }

    if (Tcl_GetDouble(interp, argv[5], &b) != TCL_OK) {
      opserr << "WARNING invalid b\n";
      opserr << "uniaxialMaterial SteelMP: " << tag << endln;
      return TCL_ERROR;
    }

    if (argc < 5) {
      opserr << "WARNING insufficient number of hardening parameters\n";
      opserr << "uniaxialMaterial Steel03: " << tag << endln;
      return TCL_ERROR;
    }

    // Read optional Steel01 material parameters
    double r, coeffR1, coeffR2, a1, a2;
    r = 20.0;
    coeffR1 = 18.5;
    coeffR2 = .15;
    a1 = 0;
    a2 = 0;

    if (argc > 6) {
      if (Tcl_GetDouble(interp, argv[6], &r) != TCL_OK) {
        opserr << "WARNING invalid r\n";
        opserr << "uniaxialMaterial SteelMP: " << tag << endln;
        return TCL_ERROR;
      }

      if (Tcl_GetDouble(interp, argv[7], &coeffR1) != TCL_OK) {
        opserr << "WARNING invalid CR1\n";
        opserr << "uniaxialMaterial SteelMP: " << tag << endln;
        return TCL_ERROR;
      }

      if (Tcl_GetDouble(interp, argv[8], &coeffR2) != TCL_OK) {
        opserr << "WARNING invalid CR2\n";
        opserr << "uniaxialMaterial SteelMP: " << tag << endln;
        return TCL_ERROR;
      }

      if (Tcl_GetDouble(interp, argv[9], &a1) != TCL_OK) {
        opserr << "WARNING invalid a1\n";
        opserr << "uniaxialMaterial SteelMP: " << tag << endln;
        return TCL_ERROR;
      }

      if (Tcl_GetDouble(interp, argv[10], &a2) != TCL_OK) {
        opserr << "WARNING invalid a2\n";
        opserr << "uniaxialMaterial SteelMP: " << tag << endln;
        return TCL_ERROR;
      }
    } // if

    theMaterial = new SteelMP(tag, fy, E, b, r, coeffR1, coeffR2, a1, a2);
  }

  else if (strcmp(argv[1], "SmoothPSConcrete") == 0) {
    if (argc < 6 || argc > 9) {
      opserr << "WARNING invalid number of arguments\n";
      printCommand(argc, argv);
      opserr << "Want: uniaxialMaterial SmoothPSConcrete tag? fc? fu? Ec? "
                "<eps0?> <epsu?> <eta?>"
             << endln;
      return TCL_ERROR;
    }

    int tag;
    double fu, Ec, fc;
    double eps0 = 0.002;
    double epsu = 0.005;
    double eta = 0.2;

    if (Tcl_GetInt(interp, argv[2], &tag) != TCL_OK) {
      opserr << "WARNING invalid uniaxialMaterial SmoothPSConcrete tag"
             << endln;
      return TCL_ERROR;
    }

    if (Tcl_GetDouble(interp, argv[3], &fc) != TCL_OK) {
      opserr << "WARNING invalid fc\n";
      opserr << "uniaxiaMaterial SmoothPSConcrete: " << tag << endln;
      return TCL_ERROR;
    }

    if (Tcl_GetDouble(interp, argv[4], &fu) != TCL_OK) {
      opserr << "WARNING invalid fu\n";
      opserr << "uniaxiaMaterial SmoothPSConcrete: " << tag << endln;
      return TCL_ERROR;
    }

    if (Tcl_GetDouble(interp, argv[5], &Ec) != TCL_OK) {
      opserr << "WARNING invalid Ec\n";
      opserr << "uniaxiaMaterial SmoothPSConcrete: " << tag << endln;
      return TCL_ERROR;
    }

    if (argc >= 7)
      if (Tcl_GetDouble(interp, argv[6], &eps0) != TCL_OK) {
        opserr << "WARNING invalid eps0\n";
        opserr << "uniaxialMaterial SmoothPSConcrete: " << tag << endln;
        return TCL_ERROR;
      }

    if (argc >= 8)
      if (Tcl_GetDouble(interp, argv[7], &epsu) != TCL_OK) {
        opserr << "WARNING invalid epsu\n";
        opserr << "uniaxialMaterial SmoothPSConcrete: " << tag << endln;
        return TCL_ERROR;
      }

    if (argc >= 9)
      if (Tcl_GetDouble(interp, argv[8], &eta) != TCL_OK) {
        opserr << "WARNING invalid eta\n";
        opserr << "uniaxialMaterial SmoothPSConcrete: " << tag << endln;
        return TCL_ERROR;
      }

    // Parsing was successful, allocate the material
    theMaterial = new SmoothPSConcrete(tag, fc, fu, Ec, eps0, epsu, eta);
  }


  else if (strcmp(argv[1], "HystereticPoly") ==
      0) { // BEGIN Salvatore Sessa 14-Jan-2021 Mail: salvatore.sessa2@unina.it
    void *theMat = OPS_HystereticPoly(rt, argc, argv);
    if (theMat != 0)
      theMaterial = (UniaxialMaterial *)theMat;
    else
      return TCL_ERROR;
  } // END Salvatore Sessa 14-Jan-2021 Mail: salvatore.sessa2@unina.it
    // Fedeas
#if defined(_STEEL2) || defined(OPSDEF_UNIAXIAL_FEDEAS)
  if (theMaterial == 0)
    theMaterial =
        TclBasicBuilder_addFedeasMaterial(clientData, interp, argc, argv);
#endif
  // Drain
  if (theMaterial == 0)
    theMaterial =
        TclBasicBuilder_addDrainMaterial(clientData, interp, argc, argv);

  // SNAP
  if (theMaterial == 0)
    theMaterial =
        TclBasicBuilder_addSnapMaterial(clientData, interp, argc, argv);

  // Py, Tz, Qz models
  if (theMaterial == 0)
    theMaterial = TclBasicBuilder_addPyTzQzMaterial(clientData, interp, argc,
                                                    argv, theDomain);

  // LimitState
  // if (theMaterial == 0)
  //   theMaterial = Tcl_AddLimitStateMaterial(clientData, interp, argc, argv);

#if defined(OPSDEF_DAMAGE_FEDEAS)
  if (theMaterial == 0)
    theMaterial =
        TclBasicBuilder_addDegradingMaterial(clientData, interp, argc, argv);
#endif

  if (theMaterial == 0) {
    //
    // maybe element in a class package already loaded
    //  loop through linked list of loaded functions comparing names & if find
    //  call it
    //

    UniaxialPackageCommand *matCommands = theUniaxialPackageCommands;
    bool found = false;
    while (matCommands != NULL && found == false) {
      if (strcmp(argv[1], matCommands->funcName) == 0) {
        theMaterial = (UniaxialMaterial *)(*(matCommands->funcPtr))();
        found = true;
        ;
      } else
        matCommands = matCommands->next;
    }
  }


  //
  // if still here the element command does not exist
  //
  if (theMaterial == nullptr) {
    opserr << "WARNING could not create uniaxialMaterial " << argv[1] << endln;
    return TCL_ERROR;
  }

  // Now add the material to the modelBuilder
  if (builder->addTaggedObject<UniaxialMaterial>(*theMaterial) != TCL_OK) {
    delete theMaterial; // invoke the material objects destructor, otherwise mem
    return TCL_ERROR;
  }

  return TCL_OK;
}
