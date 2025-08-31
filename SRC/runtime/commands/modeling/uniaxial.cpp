//===----------------------------------------------------------------------===//
//
//                                   xara
//
//===----------------------------------------------------------------------===//
//                              https://xara.so
//===----------------------------------------------------------------------===//
//
// Description: This file contains the function invoked when the user invokes
// the uniaxialMaterial command in the interpreter.
//
// Written: fmk, MHS, cmp
// Created: 07/99
//
#include <Logging.h>
#include <iostream>
#include <BasicModelBuilder.h>
#include "uniaxial.hpp"
#include <packages.h>

#include <BackboneMaterial.h>        // MHS
#include <ECC01.h>                   // Won Lee
#include <EPPGapMaterial.h>          // Mackie
#include <HystereticBackbone.h>      // MHS
#include <PathIndependentMaterial.h> // MHS
#include <Steel03.h>                 // KM

#include <SelfCenteringMaterial.h> //JAE
#include <SmoothPSConcrete.h>      //Quan & Michele
#include <SteelBRB.h>              //Quan & Michele
#include <SteelMP.h>               //Quan & Michele
#include <SMAMaterial.h> // Davide Fugazza

#include <string.h>
#include <assert.h>


class G3_Runtime;

extern "C" int OPS_ResetInputNoBuilder(ClientData clientData,
                                       Tcl_Interp *interp, int cArg, int mArg,
                                       TCL_Char ** const argv, Domain *domain);



typedef struct uniaxialPackageCommand {
  char *funcName;
  void *(*funcPtr)();
  struct uniaxialPackageCommand *next;
} UniaxialPackageCommand;

static UniaxialPackageCommand *theUniaxialPackageCommands = NULL;

//
// external functions
//
UniaxialMaterial *TclBasicBuilder_addPyTzQzMaterial(ClientData,
                                                    Tcl_Interp *,
                                                    int argc, TCL_Char ** const argv,
                                                    Domain *);

UniaxialMaterial *TclBasicBuilder_FRPCnfinedConcrete(ClientData,
                                                     Tcl_Interp *,
                                                     int argc, TCL_Char ** const argv,
                                                     Domain *);

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
    opserr << OpenSees::PromptValueError 
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
    opserr << OpenSees::PromptValueError << "Could not create uniaxialMaterial " << argv[1] << "\n";
    return TCL_ERROR;
  }

  // Now add the material to the modelBuilder
  if (builder->addTaggedObject<UniaxialMaterial>(*theMaterial) != TCL_OK) {
    opserr << OpenSees::PromptValueError << "Could not add uniaxialMaterial to the model builder.\n";
    delete theMaterial;
    return TCL_ERROR;
  }

  return TCL_OK;
}




#include <Pinching4Material.h>  // NM
static int
TclDispatch_newUniaxialPinching4(ClientData clientData, Tcl_Interp* interp, int argc, TCL_Char ** const argv)
{
   assert(clientData != nullptr);
   BasicModelBuilder *builder = static_cast<BasicModelBuilder*>(clientData);
   UniaxialMaterial* theMaterial = nullptr;

   if (strcmp(argv[1], "Pinching4") == 0) {
      if (argc != 42 && argc != 31) {
        opserr << OpenSees::PromptValueError << "insufficient arguments\n";
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
        opserr << OpenSees::PromptValueError << "invalid uniaxialMaterial Pinching4 tag" << "\n";
        return TCL_ERROR;
      }

      if (Tcl_GetDouble(interp, argv[i++], &stress1p) != TCL_OK) {
        opserr << OpenSees::PromptValueError << "invalid stress1p\n";
        opserr << "Pinching4 material: " << tag << "\n";
        return TCL_ERROR;
      }

      if (Tcl_GetDouble(interp, argv[i++], &strain1p) != TCL_OK) {
        opserr << OpenSees::PromptValueError << "invalid strain1p\n";
        opserr << "Pinching4 material: " << tag << "\n";
        return TCL_ERROR;
      }

      if (Tcl_GetDouble(interp, argv[i++], &stress2p) != TCL_OK) {
        opserr << OpenSees::PromptValueError << "invalid stress2p\n";
        opserr << "Pinching4 material: " << tag << "\n";
        return TCL_ERROR;
      }

      if (Tcl_GetDouble(interp, argv[i++], &strain2p) != TCL_OK) {
        opserr << OpenSees::PromptValueError << "invalid strain2p\n";
        opserr << "Pinching4 material: " << tag << "\n";
        return TCL_ERROR;
      }

      if (Tcl_GetDouble(interp, argv[i++], &stress3p) != TCL_OK) {
        opserr << OpenSees::PromptValueError << "invalid stress3p\n";
        opserr << "Pinching4 material: " << tag << "\n";
        return TCL_ERROR;
      }

      if (Tcl_GetDouble(interp, argv[i++], &strain3p) != TCL_OK) {
        opserr << OpenSees::PromptValueError << "invalid strain3p\n";
        opserr << "Pinching4 material: " << tag << "\n";
        return TCL_ERROR;
      }

      if (Tcl_GetDouble(interp, argv[i++], &stress4p) != TCL_OK) {
        opserr << OpenSees::PromptValueError << "invalid stress4p\n";
        opserr << "Pinching4 material: " << tag << "\n";
        return TCL_ERROR;
      }

      if (Tcl_GetDouble(interp, argv[i++], &strain4p) != TCL_OK) {
        opserr << OpenSees::PromptValueError << "invalid strain4p\n";
        opserr << "Pinching4 material: " << tag << "\n";
        return TCL_ERROR;
      }

      if (argc == 42) {
        if (Tcl_GetDouble(interp, argv[i++], &stress1n) != TCL_OK) {
          opserr << OpenSees::PromptValueError << "invalid stress1n\n";
          opserr << "Pinching4 material: " << tag << "\n";
          return TCL_ERROR;
        }

        if (Tcl_GetDouble(interp, argv[i++], &strain1n) != TCL_OK) {
          opserr << OpenSees::PromptValueError << "invalid strain1n\n";
          opserr << "Pinching4 material: " << tag << "\n";
          return TCL_ERROR;
        }

        if (Tcl_GetDouble(interp, argv[i++], &stress2n) != TCL_OK) {
          opserr << OpenSees::PromptValueError << "invalid stress2n\n";
          opserr << "Pinching4 material: " << tag << "\n";
          return TCL_ERROR;
        }

        if (Tcl_GetDouble(interp, argv[i++], &strain2n) != TCL_OK) {
          opserr << OpenSees::PromptValueError << "invalid strain2n\n";
          opserr << "Pinching4 material: " << tag << "\n";
          return TCL_ERROR;
        }

        if (Tcl_GetDouble(interp, argv[i++], &stress3n) != TCL_OK) {
          opserr << OpenSees::PromptValueError << "invalid stress3n\n";
          opserr << "Pinching4 material: " << tag << "\n";
          return TCL_ERROR;
        }

        if (Tcl_GetDouble(interp, argv[i++], &strain3n) != TCL_OK) {
          opserr << OpenSees::PromptValueError << "invalid strain3n\n";
          opserr << "Pinching4 material: " << tag << "\n";
          return TCL_ERROR;
        }

        if (Tcl_GetDouble(interp, argv[i++], &stress4n) != TCL_OK) {
          opserr << OpenSees::PromptValueError << "invalid stress4n\n";
          opserr << "Pinching4 material: " << tag << "\n";
          return TCL_ERROR;
        }

        if (Tcl_GetDouble(interp, argv[i++], &strain4n) != TCL_OK) {
          opserr << OpenSees::PromptValueError << "invalid strain4n\n";
          opserr << "Pinching4 material: " << tag << "\n";
          return TCL_ERROR;
        }
      }

      if (Tcl_GetDouble(interp, argv[i++], &rDispP) != TCL_OK) {
        opserr << OpenSees::PromptValueError << "invalid rDispP\n";
        opserr << "Pinching4 material: " << tag << "\n";
        return TCL_ERROR;
      }

      if (Tcl_GetDouble(interp, argv[i++], &rForceP) != TCL_OK) {
        opserr << OpenSees::PromptValueError << "invalid rForceP\n";
        opserr << "Pinching4 material: " << tag << "\n";
        return TCL_ERROR;
      }

      if (Tcl_GetDouble(interp, argv[i++], &uForceP) != TCL_OK) {
        opserr << OpenSees::PromptValueError << "invalid uForceP\n";
        opserr << "Pinching4 material: " << tag << "\n";
        return TCL_ERROR;
      }

      if (argc == 42) {
        if (Tcl_GetDouble(interp, argv[i++], &rDispN) != TCL_OK) {
          opserr << OpenSees::PromptValueError << "invalid rDispN\n";
          opserr << "Pinching4 material: " << tag << "\n";
          return TCL_ERROR;
        }

        if (Tcl_GetDouble(interp, argv[i++], &rForceN) != TCL_OK) {
          opserr << OpenSees::PromptValueError << "invalid rForceN\n";
          opserr << "Pinching4 material: " << tag << "\n";
          return TCL_ERROR;
        }

        if (Tcl_GetDouble(interp, argv[i++], &uForceN) != TCL_OK) {
          opserr << OpenSees::PromptValueError << "invalid uForceN\n";
          opserr << "Pinching4 material: " << tag << "\n";
          return TCL_ERROR;
        }
      }

      if (Tcl_GetDouble(interp, argv[i++], &gammaK1) != TCL_OK) {
        opserr << OpenSees::PromptValueError << "invalid gammaK1\n";
        opserr << "Pinching4 material: " << tag << "\n";
        return TCL_ERROR;
      }
      if (Tcl_GetDouble(interp, argv[i++], &gammaK2) != TCL_OK) {
        opserr << OpenSees::PromptValueError << "invalid gammaK2\n";
        opserr << "Pinching4 material: " << tag << "\n";
        return TCL_ERROR;
      }
      if (Tcl_GetDouble(interp, argv[i++], &gammaK3) != TCL_OK) {
        opserr << OpenSees::PromptValueError << "invalid gammaK3\n";
        opserr << "Pinching4 material: " << tag << "\n";
        return TCL_ERROR;
      }
      if (Tcl_GetDouble(interp, argv[i++], &gammaK4) != TCL_OK) {
        opserr << OpenSees::PromptValueError << "invalid gammaK4\n";
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
