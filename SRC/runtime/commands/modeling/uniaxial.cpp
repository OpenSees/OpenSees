//===----------------------------------------------------------------------===//
//
//        OpenSees - Open System for Earthquake Engineering Simulation
//
//===----------------------------------------------------------------------===//
//
// Description: This file contains the function invoked when the user invokes
// the uniaxialMaterial command in the interpreter.
//
// Written: fmk, MHS, cmp
// Created: 07/99
//

#include <G3_Logging.h>
#include <iostream>
#include <BasicModelBuilder.h>
#include "uniaxial.hpp"
#include <packages.h>

#include <BackboneMaterial.h>        // MHS
#include <Concrete01WithSITC.h>      // Won Lee
#include <ECC01.h>                   // Won Lee
#include <ENTMaterial.h>             // MHS
#include <EPPGapMaterial.h>          // Mackie
#include <Elastic2Material.h>        // ZHY
#include <FatigueMaterial.h>         // Patxi
#include <HystereticBackbone.h>      // MHS
#include <PathIndependentMaterial.h> // MHS
#include <ShearPanelMaterial.h>      // NM
#include <Steel03.h>                 // KM

#include <SelfCenteringMaterial.h> //JAE
#include <SmoothPSConcrete.h>      //Quan & Michele
#include <SteelBRB.h>              //Quan & Michele
#include <SteelMP.h>               //Quan & Michele
#include <SMAMaterial.h> // Davide Fugazza

#include <string.h>
#include <assert.h>

#include <UniaxialJ2Plasticity.h> // Quan

class G3_Runtime;

extern "C" int OPS_ResetInputNoBuilder(ClientData clientData,
                                       Tcl_Interp *interp, int cArg, int mArg,
                                       TCL_Char ** const argv, Domain *domain);

// extern void *OPS_PlateBearingConnectionThermal(G3_Runtime*);

// extern int TclCommand_ConfinedConcrete02(ClientData clientData, Tcl_Interp
// *interp, int argc, 					 TCL_Char ** const argv, TclBasicBuilder
// *theTclBuilder);

// extern UniaxialMaterial *Tcl_AddLimitStateMaterial(ClientData clientData,
//                                                    Tcl_Interp *interp, int argc,
//                                                    TCL_Char **arg);

#if 0
extern UniaxialMaterial *
Tcl_addWrapperUniaxialMaterial(matObj *, ClientData clientData,
                               Tcl_Interp *interp, int argc, TCL_Char ** const argv);
#endif
typedef struct uniaxialPackageCommand {
  char *funcName;
  void *(*funcPtr)();
  struct uniaxialPackageCommand *next;
} UniaxialPackageCommand;

static UniaxialPackageCommand *theUniaxialPackageCommands = NULL;

static void printCommand(int argc, TCL_Char ** const argv) {
  opserr << "Input command: ";
  for (int i = 0; i < argc; ++i)
    opserr << argv[i] << " ";
  opserr << "\n";
}

//
// external functions
//
UniaxialMaterial *TclBasicBuilder_addPyTzQzMaterial(ClientData clientData,
                                                    Tcl_Interp *interp,
                                                    int argc, TCL_Char ** const argv,
                                                    Domain *theDomain);

UniaxialMaterial *TclBasicBuilder_FRPCnfinedConcrete(ClientData clientData,
                                                     Tcl_Interp *interp,
                                                     int argc, TCL_Char ** const argv,
                                                     Domain *theDomain);

UniaxialMaterial *TclBasicBuilder_addDegradingMaterial(ClientData, Tcl_Interp *,
                                                       int, TCL_Char **);

int
TclCommand_addUniaxialMaterial(ClientData clientData, Tcl_Interp *interp,
                                  int argc, TCL_Char ** const argv)
{

  assert(clientData != nullptr);
  BasicModelBuilder *builder = static_cast<BasicModelBuilder*>(clientData);
  Domain *theDomain = builder->getDomain();

  // Make sure there is a minimum number of arguments
  if (argc < 3) {
    opserr << G3_ERROR_PROMPT 
           << "insufficient number of uniaxial material arguments\n"
           << "Want: uniaxialMaterial type? tag? <specific material args>\n";
    return TCL_ERROR;
  }

  OPS_ResetInputNoBuilder(clientData, interp, 2, argc, argv, theDomain);


  auto tcl_cmd = uniaxial_dispatch.find(std::string(argv[1]));
  if (tcl_cmd != uniaxial_dispatch.end())
    return (*tcl_cmd->second)(clientData, interp, argc, &argv[0]);


  // Pointer to a uniaxial material that will be added to the model builder
  UniaxialMaterial *theMaterial = nullptr;

  //
  // Check for packages
  //
  char *mat_name;
  if ((mat_name = strstr((char *)argv[1], "::"))) {
    // TODO: clean this up!!!!!!!!!!!!!!
    char **new_argv = new char*[argc];
    for (int i=0; i<argc; ++i)
      new_argv[i] = (char*)argv[i];
    new_argv[1] = mat_name+2;
    char pack_name[40];
    int i = 0;
    while (argv[1][i] != ':')
      pack_name[i] = argv[1][i], i++;

    pack_name[i] = '\0';
    theMaterial = (*tcl_uniaxial_package_table[pack_name])
                  (clientData,interp,argc,(const char**)new_argv);
    delete[] new_argv;
  }

  // Fedeas
#if defined(_STEEL2) || defined(OPSDEF_UNIAXIAL_FEDEAS)
  if (theMaterial == nullptr)
    theMaterial =
        TclBasicBuilder_addFedeasMaterial(clientData, interp, argc, argv);
#endif



  // Py, Tz, Qz models
  if (theMaterial == nullptr)
    theMaterial = TclBasicBuilder_addPyTzQzMaterial(clientData, interp, argc, argv, theDomain);

  // // LimitState
  // if (theMaterial == 0)
  //   theMaterial = Tcl_AddLimitStateMaterial(clientData, interp, argc, argv);

  if (theMaterial == nullptr) {
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
  // check to see if element is a procedure
  //   the proc may already have been loaded from a package or may exist in a
  //   package yet to be loaded
  //
  if (theMaterial == nullptr) {
#if 0
    //
    // maybe material in a routine
    //
    char *matType = new char[strlen(argv[1]) + 1];
    strcpy(matType, argv[1]);
    matObj *matObject = OPS_GetMaterialType(matType, strlen(matType));

    delete[] matType;

    if (matObject != 0) {

      theMaterial = Tcl_addWrapperUniaxialMaterial(matObject, clientData, interp, argc, argv);

      if (theMaterial == 0)
        delete matObject;
    }
#endif
  }

  //
  // maybe material class exists in a package yet to be loaded
  //

  if (theMaterial == nullptr) {

    void *libHandle;
    void *(*funcPtr)();

    int matNameLength = strlen(argv[1]);
    char *tclFuncName = new char[matNameLength + 12];
    strcpy(tclFuncName, "OPS_");
    strcpy(&tclFuncName[4], argv[1]);
    int res =
        getLibraryFunction(argv[1], tclFuncName, &libHandle, (void **)&funcPtr);

    delete[] tclFuncName;

    if (res == 0) {
      //
      // add loaded function to list of functions
      //
      char *matName = new char[matNameLength + 1];
      strcpy(matName, argv[1]);
      UniaxialPackageCommand *theMatCommand = new UniaxialPackageCommand;
      theMatCommand->funcPtr = funcPtr;
      theMatCommand->funcName = matName;
      theMatCommand->next = theUniaxialPackageCommands;
      theUniaxialPackageCommands = theMatCommand;

      theMaterial = (UniaxialMaterial *)(*funcPtr)();
    }
  }

  //
  // if still here the element command does not exist
  //
  if (theMaterial == nullptr) {
    opserr << G3_ERROR_PROMPT << "Could not create uniaxialMaterial " << argv[1] << "\n";
    return TCL_ERROR;
  }

  // Now add the material to the modelBuilder
  if (builder->addTaggedObject<UniaxialMaterial>(*theMaterial) != TCL_OK) {
    opserr << G3_ERROR_PROMPT << "Could not add uniaxialMaterial to the model builder.\n";
    delete theMaterial;
    return TCL_ERROR;
  }

  return TCL_OK;
}


int
TclCommand_newFatigueMaterial(ClientData clientData, Tcl_Interp* interp, int argc, TCL_Char ** const argv)
{
  assert(clientData != nullptr);
  BasicModelBuilder *builder = static_cast<BasicModelBuilder*>(clientData);

  if (argc < 4) {
    opserr << G3_ERROR_PROMPT << "insufficient arguments\n";
    printCommand(argc, argv);
    opserr << "Want: uniaxialMaterial Fatigue tag? matTag?";
    opserr << " <-D_max dmax?> <-e0 e0?> <-m m?>" << "\n";
    opserr << " <-min min?> <-max max?>" << "\n";
    return TCL_ERROR;
  }

  int tag, matTag;

  if (Tcl_GetInt(interp, argv[2], &tag) != TCL_OK) {
    opserr << G3_ERROR_PROMPT << "invalid uniaxialMaterial Fatigue tag" << "\n";
    return TCL_ERROR;
  }

  if (Tcl_GetInt(interp, argv[3], &matTag) != TCL_OK) {
    opserr << G3_ERROR_PROMPT << "invalid component tag\n";
    opserr << "uniaxialMaterial Fatigue: " << tag << "\n";
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
        opserr << G3_ERROR_PROMPT << "invalid -Dmax";
        opserr << "uniaxialMaterial Fatigue: " << tag << "\n";
        return TCL_ERROR;
      }
    } else if (strcmp(argv[j], "-E0") == 0) {
      if ((j + 1 >= argc) ||
          (Tcl_GetDouble(interp, argv[j + 1], &E0) != TCL_OK)) {
        opserr << G3_ERROR_PROMPT << "invalid -E0";
        opserr << "uniaxialMaterial Fatigue: " << tag << "\n";
        return TCL_ERROR;
      }
    } else if (strcmp(argv[j], "-m") == 0) {
      if ((j + 1 >= argc) ||
          (Tcl_GetDouble(interp, argv[j + 1], &m) != TCL_OK)) {
        opserr << G3_ERROR_PROMPT << "invalid -m";
        opserr << "uniaxialMaterial Fatigue: " << tag << "\n";
        return TCL_ERROR;
      }
    } else if (strcmp(argv[j], "-min") == 0) {
      if ((j + 1 >= argc) ||
          (Tcl_GetDouble(interp, argv[j + 1], &epsmin) != TCL_OK)) {
        opserr << G3_ERROR_PROMPT << "invalid -min ";
        opserr << "uniaxialMaterial Fatigue: " << tag << "\n";
        return TCL_ERROR;
      }
    } else if (strcmp(argv[j], "-max") == 0) {
      if ((j + 1 >= argc) ||
          (Tcl_GetDouble(interp, argv[j + 1], &epsmax) != TCL_OK)) {
        opserr << G3_ERROR_PROMPT << "invalid -max";
        opserr << "uniaxialMaterial Fatigue: " << tag << "\n";
        return TCL_ERROR;
      }
    }
    j++;
  }

  UniaxialMaterial *theMat = builder->getTypedObject<UniaxialMaterial>(matTag);

  if (theMat == nullptr) {
    opserr << G3_ERROR_PROMPT << "component material does not exist\n";
    opserr << "Component material: " << matTag;
    opserr << "\nuniaxialMaterial Fatigue: " << tag << "\n";
    return TCL_ERROR;
  }

  // Parsing was successful, allocate the material
  UniaxialMaterial *theMaterial =
      new FatigueMaterial(tag, *theMat, Dmax, E0, m, epsmin, epsmax);

  if (builder->addTaggedObject<UniaxialMaterial>(*theMaterial) != TCL_OK) {
    delete theMaterial;
    return TCL_ERROR;
  }
  return TCL_OK;

}


static int
TclDispatch_LegacyUniaxials(ClientData clientData, Tcl_Interp* interp, int argc, TCL_Char**const argv)
{
  assert(clientData != nullptr);
  BasicModelBuilder *builder = static_cast<BasicModelBuilder*>(clientData);
  UniaxialMaterial  *theMaterial = nullptr;

  if (strcmp(argv[1], "Elastic2") == 0) {
    if (argc < 4 || argc > 5) {
      opserr << G3_ERROR_PROMPT << "invalid number of arguments\n";
      printCommand(argc, argv);
      opserr << "Want: uniaxialMaterial Elastic tag? E? <eta?>" << "\n";
      return TCL_ERROR;
    }

    int tag;
    double E;
    double eta = 0.0;

    if (Tcl_GetInt(interp, argv[2], &tag) != TCL_OK) {
      opserr << G3_ERROR_PROMPT << "invalid uniaxialMaterial Elastic tag" << "\n";
      return TCL_ERROR;
    }

    if (Tcl_GetDouble(interp, argv[3], &E) != TCL_OK) {
      opserr << G3_ERROR_PROMPT << "invalid E\n";
      opserr << "uniaxiaMaterial Elastic: " << tag << "\n";
      return TCL_ERROR;
    }

    if (argc == 5) {
      if (Tcl_GetDouble(interp, argv[4], &eta) != TCL_OK) {
        opserr << G3_ERROR_PROMPT << "invalid eta\n";
        opserr << "uniaxialMaterial Elastic: " << tag << "\n";
        return TCL_ERROR;
      }
    }

    // Parsing was successful, allocate the material
    theMaterial = new Elastic2Material(tag, E, eta);

  } else if (strcmp(argv[1], "ENT") == 0) {
    if (argc < 4) {
      opserr << G3_ERROR_PROMPT << "invalid number of arguments\n";
      printCommand(argc, argv);
      opserr << "Want: uniaxialMaterial ENT tag? E?" << "\n";
      return TCL_ERROR;
    }

    int tag;
    double E;

    if (Tcl_GetInt(interp, argv[2], &tag) != TCL_OK) {
      opserr << G3_ERROR_PROMPT << "invalid uniaxialMaterial ENT tag" << "\n";
      return TCL_ERROR;
    }

    if (Tcl_GetDouble(interp, argv[3], &E) != TCL_OK) {
      opserr << G3_ERROR_PROMPT << "invalid E\n";
      opserr << "uniaxiaMaterial ENT: " << tag << "\n";
      return TCL_ERROR;
    }

    // Parsing was successful, allocate the material
    theMaterial = new ENTMaterial(tag, E);

  }

  else if (strcmp(argv[1], "BarSlip") == 0) {
    if (argc != 17 && argc != 15) {
      opserr << G3_ERROR_PROMPT << "insufficient arguments\n";
      printCommand(argc, argv);
      opserr << "Want: uniaxialMaterial BarSlip tag? fc? fy? Es? fu? Eh? db? "
                "ld? nb? width? depth? bsflag? type? <damage? unit?>"
             << "\n";
      return TCL_ERROR;
    }

    int tag, nb, bsf, typ, dmg=1, unt=1;
    double fc, fy, Es, fu, Eh, ld, width, depth, db;

    int argStart = 2;

    if (Tcl_GetInt(interp, argv[argStart++], &tag) != TCL_OK) {
      opserr << G3_ERROR_PROMPT << "invalid tag\n";
      return TCL_ERROR;
    }
    if (Tcl_GetDouble(interp, argv[argStart++], &fc) != TCL_OK) {
      opserr << G3_ERROR_PROMPT << "invalid fc\n";
      return TCL_ERROR;
    }
    if (Tcl_GetDouble(interp, argv[argStart++], &fy) != TCL_OK) {
      opserr << G3_ERROR_PROMPT << "invalid fy\n";
      return TCL_ERROR;
    }
    if (Tcl_GetDouble(interp, argv[argStart++], &Es) != TCL_OK) {
      opserr << G3_ERROR_PROMPT << "invalid Es\n";
      return TCL_ERROR;
    }
    if (Tcl_GetDouble(interp, argv[argStart++], &fu) != TCL_OK) {
      opserr << G3_ERROR_PROMPT << "invalid fu\n";
      return TCL_ERROR;
    }
    if (Tcl_GetDouble(interp, argv[argStart++], &Eh) != TCL_OK) {
      opserr << G3_ERROR_PROMPT << "invalid Eh\n";
      return TCL_ERROR;
    }
    if (Tcl_GetDouble(interp, argv[argStart++], &db) != TCL_OK) {
      opserr << G3_ERROR_PROMPT << "invalid db\n";
      return TCL_ERROR;
    }
    if (Tcl_GetDouble(interp, argv[argStart++], &ld) != TCL_OK) {
      opserr << G3_ERROR_PROMPT << "invalid ld\n";
      return TCL_ERROR;
    }
    if (Tcl_GetInt(interp, argv[argStart++], &nb) != TCL_OK) {
      opserr << G3_ERROR_PROMPT << "invalid nbars\n";
      return TCL_ERROR;
    }
    if (Tcl_GetDouble(interp, argv[argStart++], &width) != TCL_OK) {
      opserr << G3_ERROR_PROMPT << "invalid width\n";
      return TCL_ERROR;
    }
    if (Tcl_GetDouble(interp, argv[argStart++], &depth) != TCL_OK) {
      opserr << G3_ERROR_PROMPT << "invalid depth\n";
      return TCL_ERROR;
    }

    int y;
    y = argStart;

    if ((strcmp(argv[y], "strong") == 0) ||
        (strcmp(argv[y], "Strong") == 0) || (strcmp(argv[y], "weak") == 0) ||
        (strcmp(argv[y], "Weak") == 0)) {
      if ((strcmp(argv[y], "strong") == 0) ||
          (strcmp(argv[y], "Strong") == 0)) {
        bsf = 0;
      }

      if ((strcmp(argv[y], "weak") == 0) || (strcmp(argv[y], "Weak") == 0)) {
        bsf = 1;
      }
    } else {
      opserr << G3_ERROR_PROMPT << "invalid bond strength specified\n";
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
        (strcmp(argv[y], "Column") == 0) ||
        (strcmp(argv[y], "column") == 0)) {
      if ((strcmp(argv[y], "beamtop") == 0) ||
          (strcmp(argv[y], "beamTop") == 0) ||
          (strcmp(argv[y], "beam") == 0) || (strcmp(argv[y], "Beam") == 0)) {
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
      opserr << G3_ERROR_PROMPT << "invalid location of bar specified\n";
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
        opserr << G3_ERROR_PROMPT << "invalid damage specified\n";
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
        opserr << G3_ERROR_PROMPT << "invalid unit specified\n";
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
      opserr << G3_ERROR_PROMPT << "insufficient arguments\n";
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
      opserr << G3_ERROR_PROMPT << "invalid uniaxialMaterial ShearPanel tag" << "\n";
      return TCL_ERROR;
    }

    if (Tcl_GetDouble(interp, argv[i++], &stress1p) != TCL_OK) {
      opserr << G3_ERROR_PROMPT << "invalid stress1p\n";
      opserr << "ShearPanel material: " << tag << "\n";
      return TCL_ERROR;
    }

    if (Tcl_GetDouble(interp, argv[i++], &strain1p) != TCL_OK) {
      opserr << G3_ERROR_PROMPT << "invalid strain1p\n";
      opserr << "ShearPanel material: " << tag << "\n";
      return TCL_ERROR;
    }

    if (Tcl_GetDouble(interp, argv[i++], &stress2p) != TCL_OK) {
      opserr << G3_ERROR_PROMPT << "invalid stress2p\n";
      opserr << "ShearPanel material: " << tag << "\n";
      return TCL_ERROR;
    }

    if (Tcl_GetDouble(interp, argv[i++], &strain2p) != TCL_OK) {
      opserr << G3_ERROR_PROMPT << "invalid strain2p\n";
      opserr << "ShearPanel material: " << tag << "\n";
      return TCL_ERROR;
    }

    if (Tcl_GetDouble(interp, argv[i++], &stress3p) != TCL_OK) {
      opserr << G3_ERROR_PROMPT << "invalid stress3p\n";
      opserr << "ShearPanel material: " << tag << "\n";
      return TCL_ERROR;
    }

    if (Tcl_GetDouble(interp, argv[i++], &strain3p) != TCL_OK) {
      opserr << G3_ERROR_PROMPT << "invalid strain3p\n";
      opserr << "ShearPanel material: " << tag << "\n";
      return TCL_ERROR;
    }

    if (Tcl_GetDouble(interp, argv[i++], &stress4p) != TCL_OK) {
      opserr << G3_ERROR_PROMPT << "invalid stress4p\n";
      opserr << "ShearPanel material: " << tag << "\n";
      return TCL_ERROR;
    }

    if (Tcl_GetDouble(interp, argv[i++], &strain4p) != TCL_OK) {
      opserr << G3_ERROR_PROMPT << "invalid strain4p\n";
      opserr << "ShearPanel material: " << tag << "\n";
      return TCL_ERROR;
    }

    if (argc == 42) {
      if (Tcl_GetDouble(interp, argv[i++], &stress1n) != TCL_OK) {
        opserr << G3_ERROR_PROMPT << "invalid stress1n\n";
        opserr << "ShearPanel material: " << tag << "\n";
        return TCL_ERROR;
      }

      if (Tcl_GetDouble(interp, argv[i++], &strain1n) != TCL_OK) {
        opserr << G3_ERROR_PROMPT << "invalid strain1n\n";
        opserr << "ShearPanel material: " << tag << "\n";
        return TCL_ERROR;
      }

      if (Tcl_GetDouble(interp, argv[i++], &stress2n) != TCL_OK) {
        opserr << G3_ERROR_PROMPT << "invalid stress2n\n";
        opserr << "ShearPanel material: " << tag << "\n";
        return TCL_ERROR;
      }

      if (Tcl_GetDouble(interp, argv[i++], &strain2n) != TCL_OK) {
        opserr << G3_ERROR_PROMPT << "invalid strain2n\n";
        opserr << "ShearPanel material: " << tag << "\n";
        return TCL_ERROR;
      }

      if (Tcl_GetDouble(interp, argv[i++], &stress3n) != TCL_OK) {
        opserr << G3_ERROR_PROMPT << "invalid stress3n\n";
        opserr << "ShearPanel material: " << tag << "\n";
        return TCL_ERROR;
      }

      if (Tcl_GetDouble(interp, argv[i++], &strain3n) != TCL_OK) {
        opserr << G3_ERROR_PROMPT << "invalid strain3n\n";
        opserr << "ShearPanel material: " << tag << "\n";
        return TCL_ERROR;
      }

      if (Tcl_GetDouble(interp, argv[i++], &stress4n) != TCL_OK) {
        opserr << G3_ERROR_PROMPT << "invalid stress4n\n";
        opserr << "ShearPanel material: " << tag << "\n";
        return TCL_ERROR;
      }

      if (Tcl_GetDouble(interp, argv[i++], &strain4n) != TCL_OK) {
        opserr << G3_ERROR_PROMPT << "invalid strain4n\n";
        opserr << "ShearPanel material: " << tag << "\n";
        return TCL_ERROR;
      }
    }

    if (Tcl_GetDouble(interp, argv[i++], &rDispP) != TCL_OK) {
      opserr << G3_ERROR_PROMPT << "invalid rDispP\n";
      opserr << "ShearPanel material: " << tag << "\n";
      return TCL_ERROR;
    }

    if (Tcl_GetDouble(interp, argv[i++], &rForceP) != TCL_OK) {
      opserr << G3_ERROR_PROMPT << "invalid rForceP\n";
      opserr << "ShearPanel material: " << tag << "\n";
      return TCL_ERROR;
    }

    if (Tcl_GetDouble(interp, argv[i++], &uForceP) != TCL_OK) {
      opserr << G3_ERROR_PROMPT << "invalid uForceP\n";
      opserr << "ShearPanel material: " << tag << "\n";
      return TCL_ERROR;
    }

    if (argc == 42) {
      if (Tcl_GetDouble(interp, argv[i++], &rDispN) != TCL_OK) {
        opserr << G3_ERROR_PROMPT << "invalid rDispN\n";
        opserr << "ShearPanel material: " << tag << "\n";
        return TCL_ERROR;
      }

      if (Tcl_GetDouble(interp, argv[i++], &rForceN) != TCL_OK) {
        opserr << G3_ERROR_PROMPT << "invalid rForceN\n";
        opserr << "ShearPanel material: " << tag << "\n";
        return TCL_ERROR;
      }

      if (Tcl_GetDouble(interp, argv[i++], &uForceN) != TCL_OK) {
        opserr << G3_ERROR_PROMPT << "invalid uForceN\n";
        opserr << "ShearPanel material: " << tag << "\n";
        return TCL_ERROR;
      }
    }

    if (Tcl_GetDouble(interp, argv[i++], &gammaK1) != TCL_OK) {
      opserr << G3_ERROR_PROMPT << "invalid gammaK1\n";
      opserr << "ShearPanel material: " << tag << "\n";
      return TCL_ERROR;
    }
    if (Tcl_GetDouble(interp, argv[i++], &gammaK2) != TCL_OK) {
      opserr << G3_ERROR_PROMPT << "invalid gammaK2\n";
      opserr << "ShearPanel material: " << tag << "\n";
      return TCL_ERROR;
    }
    if (Tcl_GetDouble(interp, argv[i++], &gammaK3) != TCL_OK) {
      opserr << G3_ERROR_PROMPT << "invalid gammaK3\n";
      opserr << "ShearPanel material: " << tag << "\n";
      return TCL_ERROR;
    }
    if (Tcl_GetDouble(interp, argv[i++], &gammaK4) != TCL_OK) {
      opserr << G3_ERROR_PROMPT << "invalid gammaK4\n";
      opserr << "ShearPanel material: " << tag << "\n";
      return TCL_ERROR;
    }
    if (Tcl_GetDouble(interp, argv[i++], &gammaKLimit) != TCL_OK) {
      opserr << G3_ERROR_PROMPT << "invalid gammaKLimit\n";
      opserr << "ShearPanel material: " << tag << "\n";
      return TCL_ERROR;
    }
    if (Tcl_GetDouble(interp, argv[i++], &gammaD1) != TCL_OK) {
      opserr << G3_ERROR_PROMPT << "invalid gammaD1\n";
      opserr << "ShearPanel material: " << tag << "\n";
      return TCL_ERROR;
    }
    if (Tcl_GetDouble(interp, argv[i++], &gammaD2) != TCL_OK) {
      opserr << G3_ERROR_PROMPT << "invalid gammaD2\n";
      opserr << "ShearPanel material: " << tag << "\n";
      return TCL_ERROR;
    }
    if (Tcl_GetDouble(interp, argv[i++], &gammaD3) != TCL_OK) {
      opserr << G3_ERROR_PROMPT << "invalid gammaD3\n";
      opserr << "ShearPanel material: " << tag << "\n";
      return TCL_ERROR;
    }
    if (Tcl_GetDouble(interp, argv[i++], &gammaD4) != TCL_OK) {
      opserr << G3_ERROR_PROMPT << "invalid gammaD4\n";
      opserr << "ShearPanel material: " << tag << "\n";
      return TCL_ERROR;
    }
    if (Tcl_GetDouble(interp, argv[i++], &gammaDLimit) != TCL_OK) {
      opserr << G3_ERROR_PROMPT << "invalid gammaDLimit\n";
      opserr << "ShearPanel material: " << tag << "\n";
      return TCL_ERROR;
    }
    if (Tcl_GetDouble(interp, argv[i++], &gammaF1) != TCL_OK) {
      opserr << G3_ERROR_PROMPT << "invalid gammaF1\n";
      opserr << "ShearPanel material: " << tag << "\n";
      return TCL_ERROR;
    }
    if (Tcl_GetDouble(interp, argv[i++], &gammaF2) != TCL_OK) {
      opserr << G3_ERROR_PROMPT << "invalid gammaF2\n";
      opserr << "ShearPanel material: " << tag << "\n";
      return TCL_ERROR;
    }
    if (Tcl_GetDouble(interp, argv[i++], &gammaF3) != TCL_OK) {
      opserr << G3_ERROR_PROMPT << "invalid gammaF3\n";
      opserr << "ShearPanel material: " << tag << "\n";
      return TCL_ERROR;
    }
    if (Tcl_GetDouble(interp, argv[i++], &gammaF4) != TCL_OK) {
      opserr << G3_ERROR_PROMPT << "invalid gammaF4\n";
      opserr << "ShearPanel material: " << tag << "\n";
      return TCL_ERROR;
    }
    if (Tcl_GetDouble(interp, argv[i++], &gammaFLimit) != TCL_OK) {
      opserr << G3_ERROR_PROMPT << "invalid gammaFLimit\n";
      opserr << "ShearPanel material: " << tag << "\n";
      return TCL_ERROR;
    }

    if (Tcl_GetDouble(interp, argv[i++], &gammaE) != TCL_OK) {
      opserr << G3_ERROR_PROMPT << "invalid gammaE\n";
      opserr << "ShearPanel material: " << tag << "\n";
      return TCL_ERROR;
    }

    if (Tcl_GetDouble(interp, argv[i++], &yStr) != TCL_OK) {
      opserr << G3_ERROR_PROMPT << "invalid yield stress\n";
      opserr << "ShearPanel material: " << tag << "\n";
      return TCL_ERROR;
    }

    // allocate the pinching material
    if (argc == 42) {
      theMaterial = new ShearPanelMaterial(
          tag, stress1p, strain1p, stress2p, strain2p, stress3p, strain3p,
          stress4p, strain4p, stress1n, strain1n, stress2n, strain2n,
          stress3n, strain3n, stress4n, strain4n, rDispP, rForceP, uForceP,
          rDispN, rForceN, uForceN, gammaK1, gammaK2, gammaK3, gammaK4,
          gammaKLimit, gammaD1, gammaD2, gammaD3, gammaD4, gammaDLimit,
          gammaF1, gammaF2, gammaF3, gammaF4, gammaFLimit, gammaE, yStr);
    }
    if (argc == 31) {
      theMaterial = new ShearPanelMaterial(
          tag, stress1p, strain1p, stress2p, strain2p, stress3p, strain3p,
          stress4p, strain4p, rDispP, rForceP, uForceP, gammaK1, gammaK2,
          gammaK3, gammaK4, gammaKLimit, gammaD1, gammaD2, gammaD3, gammaD4,
          gammaDLimit, gammaF1, gammaF2, gammaF3, gammaF4, gammaFLimit,
          gammaE, yStr);
    }
  }

  else if (strcmp(argv[1], "Concrete01WithSITC") == 0) {
    if (argc < 7) {
      opserr << G3_ERROR_PROMPT << "insufficient arguments\n";
      printCommand(argc, argv);
      opserr << "Want: uniaxialMaterial Concrete01 tag? fpc? epsc0? fpcu? "
                "epscu? <endStrainSITC?>"
             << "\n";
      return TCL_ERROR;
    }

    int tag;

    if (Tcl_GetInt(interp, argv[2], &tag) != TCL_OK) {
      opserr << G3_ERROR_PROMPT << "invalid uniaxialMaterial Concrete01 tag" << "\n";
      return TCL_ERROR;
    }

    // Read required Concrete01 material parameters
    double fpc, epsc0, fpcu, epscu;

    if (Tcl_GetDouble(interp, argv[3], &fpc) != TCL_OK) {
      opserr << G3_ERROR_PROMPT << "invalid fpc\n";
      opserr << "Concrete01 material: " << tag << "\n";
      return TCL_ERROR;
    }

    if (Tcl_GetDouble(interp, argv[4], &epsc0) != TCL_OK) {
      opserr << G3_ERROR_PROMPT << "invalid epsc0\n";
      opserr << "Concrete01 material: " << tag << "\n";
      return TCL_ERROR;
    }

    if (Tcl_GetDouble(interp, argv[5], &fpcu) != TCL_OK) {
      opserr << G3_ERROR_PROMPT << "invalid fpcu\n";
      opserr << "Concrete01 material: " << tag << "\n";
      return TCL_ERROR;
    }

    if (Tcl_GetDouble(interp, argv[6], &epscu) != TCL_OK) {
      opserr << G3_ERROR_PROMPT << "invalid epscu\n";
      opserr << "Concrete01 material: " << tag << "\n";
      return TCL_ERROR;
    }

    if (argc == 7)
      // Parsing was successful, allocate the material
      theMaterial = new Concrete01WithSITC(tag, fpc, epsc0, fpcu, epscu);
    else {
      double endStrainSITC;
      if (Tcl_GetDouble(interp, argv[7], &endStrainSITC) != TCL_OK) {
        opserr << G3_ERROR_PROMPT << "invalid epscu\n";
        opserr << "Concrete01 material: " << tag << "\n";
        return TCL_ERROR;
      }
      theMaterial =
          new Concrete01WithSITC(tag, fpc, epsc0, fpcu, epscu, endStrainSITC);
    }
  }

  if (theMaterial == nullptr)
    return TCL_ERROR;

  if (builder->addTaggedObject<UniaxialMaterial>(*theMaterial) != TCL_OK) {
    delete theMaterial;
    return TCL_ERROR;
  }

  return TCL_OK;
}

#include <UniaxialJ2Plasticity.h>
static int
TclCommand_newUniaxialJ2Plasticity(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char** const argv)
{
  // ----- 1D J2 Plasticity ----
    if (argc < 7) {
      opserr << "WARNING invalid number of arguments\n";
      printCommand(argc, argv);
      opserr << "Want: uniaxialMaterial UniaxialJ2Plasticity tag? E? sigmaY? "
                "Hkin? <Hiso?>"
             << endln;
      return TCL_ERROR;
    }

    int tag;
    double E, sigmaY, Hkin, Hiso;
    Hiso = 0.0;

    if (Tcl_GetInt(interp, argv[2], &tag) != TCL_OK) {
      opserr << "WARNING invalid uniaxialMaterial UniaxialJ2Plasticity tag"
             << endln;
      return TCL_ERROR;
    }

    if (Tcl_GetDouble(interp, argv[3], &E) != TCL_OK) {
      opserr << "WARNING invalid E\n";
      opserr << "uniaxiaMaterial UniaxialJ2Plasticity: " << tag << endln;
      return TCL_ERROR;
    }

    if (Tcl_GetDouble(interp, argv[4], &sigmaY) != TCL_OK) {
      opserr << "WARNING invalid sigmaY\n";
      opserr << "uniaxiaMaterial UniaxialJ2Plasticity: " << tag << endln;
      return TCL_ERROR;
    }

    if (Tcl_GetDouble(interp, argv[5], &Hkin) != TCL_OK) {
      opserr << "WARNING invalid Hkin\n";
      opserr << "uniaxiaMaterial SmoothPSConcrete: " << tag << endln;
      return TCL_ERROR;
    }

    if (argc >= 7)
      if (Tcl_GetDouble(interp, argv[6], &Hiso) != TCL_OK) {
        opserr << "WARNING invalid Hiso\n";
        opserr << "uniaxialMaterial UniaxialJ2Plasticity: " << tag << endln;
        return TCL_ERROR;
      }

    // Parsing was successful, allocate the material
    UniaxialMaterial* theMaterial = new UniaxialJ2Plasticity(tag, E, sigmaY, Hkin, Hiso);

   assert(clientData != nullptr);
   BasicModelBuilder *builder = static_cast<BasicModelBuilder*>(clientData);
   builder->addTaggedObject<UniaxialMaterial>(*theMaterial);
   return TCL_OK;

}


#include <Pinching4Material.h>       // NM
static int
TclDispatch_newUniaxialPinching4(ClientData clientData, Tcl_Interp* interp, int argc, TCL_Char ** const argv)
{
   assert(clientData != nullptr);
   BasicModelBuilder *builder = static_cast<BasicModelBuilder*>(clientData);
   UniaxialMaterial* theMaterial = nullptr;

   if (strcmp(argv[1], "Pinching4") == 0) {
      if (argc != 42 && argc != 31) {
        opserr << G3_ERROR_PROMPT << "insufficient arguments\n";
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
        opserr << G3_ERROR_PROMPT << "invalid uniaxialMaterial Pinching4 tag" << "\n";
        return TCL_ERROR;
      }

      if (Tcl_GetDouble(interp, argv[i++], &stress1p) != TCL_OK) {
        opserr << G3_ERROR_PROMPT << "invalid stress1p\n";
        opserr << "Pinching4 material: " << tag << "\n";
        return TCL_ERROR;
      }

      if (Tcl_GetDouble(interp, argv[i++], &strain1p) != TCL_OK) {
        opserr << G3_ERROR_PROMPT << "invalid strain1p\n";
        opserr << "Pinching4 material: " << tag << "\n";
        return TCL_ERROR;
      }

      if (Tcl_GetDouble(interp, argv[i++], &stress2p) != TCL_OK) {
        opserr << G3_ERROR_PROMPT << "invalid stress2p\n";
        opserr << "Pinching4 material: " << tag << "\n";
        return TCL_ERROR;
      }

      if (Tcl_GetDouble(interp, argv[i++], &strain2p) != TCL_OK) {
        opserr << G3_ERROR_PROMPT << "invalid strain2p\n";
        opserr << "Pinching4 material: " << tag << "\n";
        return TCL_ERROR;
      }

      if (Tcl_GetDouble(interp, argv[i++], &stress3p) != TCL_OK) {
        opserr << G3_ERROR_PROMPT << "invalid stress3p\n";
        opserr << "Pinching4 material: " << tag << "\n";
        return TCL_ERROR;
      }

      if (Tcl_GetDouble(interp, argv[i++], &strain3p) != TCL_OK) {
        opserr << G3_ERROR_PROMPT << "invalid strain3p\n";
        opserr << "Pinching4 material: " << tag << "\n";
        return TCL_ERROR;
      }

      if (Tcl_GetDouble(interp, argv[i++], &stress4p) != TCL_OK) {
        opserr << G3_ERROR_PROMPT << "invalid stress4p\n";
        opserr << "Pinching4 material: " << tag << "\n";
        return TCL_ERROR;
      }

      if (Tcl_GetDouble(interp, argv[i++], &strain4p) != TCL_OK) {
        opserr << G3_ERROR_PROMPT << "invalid strain4p\n";
        opserr << "Pinching4 material: " << tag << "\n";
        return TCL_ERROR;
      }

      if (argc == 42) {
        if (Tcl_GetDouble(interp, argv[i++], &stress1n) != TCL_OK) {
          opserr << G3_ERROR_PROMPT << "invalid stress1n\n";
          opserr << "Pinching4 material: " << tag << "\n";
          return TCL_ERROR;
        }

        if (Tcl_GetDouble(interp, argv[i++], &strain1n) != TCL_OK) {
          opserr << G3_ERROR_PROMPT << "invalid strain1n\n";
          opserr << "Pinching4 material: " << tag << "\n";
          return TCL_ERROR;
        }

        if (Tcl_GetDouble(interp, argv[i++], &stress2n) != TCL_OK) {
          opserr << G3_ERROR_PROMPT << "invalid stress2n\n";
          opserr << "Pinching4 material: " << tag << "\n";
          return TCL_ERROR;
        }

        if (Tcl_GetDouble(interp, argv[i++], &strain2n) != TCL_OK) {
          opserr << G3_ERROR_PROMPT << "invalid strain2n\n";
          opserr << "Pinching4 material: " << tag << "\n";
          return TCL_ERROR;
        }

        if (Tcl_GetDouble(interp, argv[i++], &stress3n) != TCL_OK) {
          opserr << G3_ERROR_PROMPT << "invalid stress3n\n";
          opserr << "Pinching4 material: " << tag << "\n";
          return TCL_ERROR;
        }

        if (Tcl_GetDouble(interp, argv[i++], &strain3n) != TCL_OK) {
          opserr << G3_ERROR_PROMPT << "invalid strain3n\n";
          opserr << "Pinching4 material: " << tag << "\n";
          return TCL_ERROR;
        }

        if (Tcl_GetDouble(interp, argv[i++], &stress4n) != TCL_OK) {
          opserr << G3_ERROR_PROMPT << "invalid stress4n\n";
          opserr << "Pinching4 material: " << tag << "\n";
          return TCL_ERROR;
        }

        if (Tcl_GetDouble(interp, argv[i++], &strain4n) != TCL_OK) {
          opserr << G3_ERROR_PROMPT << "invalid strain4n\n";
          opserr << "Pinching4 material: " << tag << "\n";
          return TCL_ERROR;
        }
      }

      if (Tcl_GetDouble(interp, argv[i++], &rDispP) != TCL_OK) {
        opserr << G3_ERROR_PROMPT << "invalid rDispP\n";
        opserr << "Pinching4 material: " << tag << "\n";
        return TCL_ERROR;
      }

      if (Tcl_GetDouble(interp, argv[i++], &rForceP) != TCL_OK) {
        opserr << G3_ERROR_PROMPT << "invalid rForceP\n";
        opserr << "Pinching4 material: " << tag << "\n";
        return TCL_ERROR;
      }

      if (Tcl_GetDouble(interp, argv[i++], &uForceP) != TCL_OK) {
        opserr << G3_ERROR_PROMPT << "invalid uForceP\n";
        opserr << "Pinching4 material: " << tag << "\n";
        return TCL_ERROR;
      }

      if (argc == 42) {
        if (Tcl_GetDouble(interp, argv[i++], &rDispN) != TCL_OK) {
          opserr << G3_ERROR_PROMPT << "invalid rDispN\n";
          opserr << "Pinching4 material: " << tag << "\n";
          return TCL_ERROR;
        }

        if (Tcl_GetDouble(interp, argv[i++], &rForceN) != TCL_OK) {
          opserr << G3_ERROR_PROMPT << "invalid rForceN\n";
          opserr << "Pinching4 material: " << tag << "\n";
          return TCL_ERROR;
        }

        if (Tcl_GetDouble(interp, argv[i++], &uForceN) != TCL_OK) {
          opserr << G3_ERROR_PROMPT << "invalid uForceN\n";
          opserr << "Pinching4 material: " << tag << "\n";
          return TCL_ERROR;
        }
      }

      if (Tcl_GetDouble(interp, argv[i++], &gammaK1) != TCL_OK) {
        opserr << G3_ERROR_PROMPT << "invalid gammaK1\n";
        opserr << "Pinching4 material: " << tag << "\n";
        return TCL_ERROR;
      }
      if (Tcl_GetDouble(interp, argv[i++], &gammaK2) != TCL_OK) {
        opserr << G3_ERROR_PROMPT << "invalid gammaK2\n";
        opserr << "Pinching4 material: " << tag << "\n";
        return TCL_ERROR;
      }
      if (Tcl_GetDouble(interp, argv[i++], &gammaK3) != TCL_OK) {
        opserr << G3_ERROR_PROMPT << "invalid gammaK3\n";
        opserr << "Pinching4 material: " << tag << "\n";
        return TCL_ERROR;
      }
      if (Tcl_GetDouble(interp, argv[i++], &gammaK4) != TCL_OK) {
        opserr << G3_ERROR_PROMPT << "invalid gammaK4\n";
        opserr << "Pinching4 material: " << tag << "\n";
        return TCL_ERROR;
      }
      if (Tcl_GetDouble(interp, argv[i++], &gammaKLimit) != TCL_OK) {
        opserr << "WARNING invalid gammaKLimit\n";
        opserr << "Pinching4 material: " << tag << "\n";
        return TCL_ERROR;
      }
      if (Tcl_GetDouble(interp, argv[i++], &gammaD1) != TCL_OK) {
        opserr << "WARNING invalid gammaD1\n";
        opserr << "Pinching4 material: " << tag << "\n";
        return TCL_ERROR;
      }
      if (Tcl_GetDouble(interp, argv[i++], &gammaD2) != TCL_OK) {
        opserr << "WARNING invalid gammaD2\n";
        opserr << "Pinching4 material: " << tag << "\n";
        return TCL_ERROR;
      }
      if (Tcl_GetDouble(interp, argv[i++], &gammaD3) != TCL_OK) {
        opserr << "WARNING invalid gammaD3\n";
        opserr << "Pinching4 material: " << tag << "\n";
        return TCL_ERROR;
      }
      if (Tcl_GetDouble(interp, argv[i++], &gammaD4) != TCL_OK) {
        opserr << "WARNING invalid gammaD4\n";
        opserr << "Pinching4 material: " << tag << "\n";
        return TCL_ERROR;
      }
      if (Tcl_GetDouble(interp, argv[i++], &gammaDLimit) != TCL_OK) {
        opserr << "WARNING invalid gammaDLimit\n";
        opserr << "Pinching4 material: " << tag << "\n";
        return TCL_ERROR;
      }
      if (Tcl_GetDouble(interp, argv[i++], &gammaF1) != TCL_OK) {
        opserr << "WARNING invalid gammaF1\n";
        opserr << "Pinching4 material: " << tag << "\n";
        return TCL_ERROR;
      }
      if (Tcl_GetDouble(interp, argv[i++], &gammaF2) != TCL_OK) {
        opserr << "WARNING invalid gammaF2\n";
        opserr << "Pinching4 material: " << tag << "\n";
        return TCL_ERROR;
      }
      if (Tcl_GetDouble(interp, argv[i++], &gammaF3) != TCL_OK) {
        opserr << "WARNING invalid gammaF3\n";
        opserr << "Pinching4 material: " << tag << "\n";
        return TCL_ERROR;
      }
      if (Tcl_GetDouble(interp, argv[i++], &gammaF4) != TCL_OK) {
        opserr << "WARNING invalid gammaF4\n";
        opserr << "Pinching4 material: " << tag << "\n";
        return TCL_ERROR;
      }
      if (Tcl_GetDouble(interp, argv[i++], &gammaFLimit) != TCL_OK) {
        opserr << "WARNING invalid gammaFLimit\n";
        opserr << "Pinching4 material: " << tag << "\n";
        return TCL_ERROR;
      }

      if (Tcl_GetDouble(interp, argv[i++], &gammaE) != TCL_OK) {
        opserr << "WARNING invalid gammaE\n";
        opserr << "Pinching4 material: " << tag << "\n";
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
        opserr << "Pinching4 material: " << tag << "\n";
        return TCL_ERROR;
      }

      // allocate the pinching material
      if (argc == 42) {
        theMaterial = new Pinching4Material(
            tag, stress1p, strain1p, stress2p, strain2p, stress3p, strain3p,
            stress4p, strain4p, stress1n, strain1n, stress2n, strain2n,
            stress3n, strain3n, stress4n, strain4n, rDispP, rForceP, uForceP,
            rDispN, rForceN, uForceN, gammaK1, gammaK2, gammaK3, gammaK4,
            gammaKLimit, gammaD1, gammaD2, gammaD3, gammaD4, gammaDLimit,
            gammaF1, gammaF2, gammaF3, gammaF4, gammaFLimit, gammaE, tDmg);
      }
      if (argc == 31) {
        theMaterial = new Pinching4Material(
            tag, stress1p, strain1p, stress2p, strain2p, stress3p, strain3p,
            stress4p, strain4p, rDispP, rForceP, uForceP, gammaK1, gammaK2,
            gammaK3, gammaK4, gammaKLimit, gammaD1, gammaD2, gammaD3, gammaD4,
            gammaDLimit, gammaF1, gammaF2, gammaF3, gammaF4, gammaFLimit,
            gammaE, tDmg);
      }
  }
  if (!theMaterial || (builder->addTaggedObject<UniaxialMaterial>(*theMaterial) != TCL_OK)) {
    delete theMaterial;
    return TCL_ERROR;
  }
  return TCL_OK;
}


#if 0
   else if (strcmp(argv[1], "Backbone") == 0) {
      if (argc < 4) {
        opserr << G3_ERROR_PROMPT << "insufficient arguments\n";
        printCommand(argc, argv);
        opserr << "Want: uniaxialMaterial Backbone tag? bbTag?" << "\n";
        return TCL_ERROR;
      }

      int tag, bbTag;

      if (Tcl_GetInt(interp, argv[2], &tag) != TCL_OK) {
        opserr << G3_ERROR_PROMPT << "invalid tag\n";
        opserr << "Backbone material: " << tag << "\n";
        return TCL_ERROR;
      }

      if (Tcl_GetInt(interp, argv[3], &bbTag) != TCL_OK) {
        opserr << G3_ERROR_PROMPT << "invalid bTag\n";
        opserr << "Backbone material: " << tag << "\n";
        return TCL_ERROR;
      }

      HystereticBackbone *backbone = OPS_getHystereticBackbone(bbTag);

      if (backbone == 0) {
        opserr << G3_ERROR_PROMPT << "backbone does not exist\n";
        opserr << "backbone: " << bbTag;
        opserr << "\nuniaxialMaterial Backbone: " << tag << "\n";
        return TCL_ERROR;
      }

      theMaterial = new BackboneMaterial(tag, *backbone);
    }
#endif
