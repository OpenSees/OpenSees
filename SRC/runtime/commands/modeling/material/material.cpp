/* ****************************************************************** **
**    OpenSees - Open System for Earthquake Engineering Simulation    **
**          Pacific Earthquake Engineering Research Center            **
**                                                                    **
** Developed by:                                                      **
**   Frank McKenna (fmckenna@ce.berkeley.edu)                         **
**   Gregory L. Fenves (fenves@ce.berkeley.edu)                       **
**   Filip C. Filippou (filippou@ce.berkeley.edu)                     **
**                                                                    **
** With a lot of additions by                                         **
**   Boris Jeremic    (jeremic@ucdavis.edu)                           **
**   Zaohui Yang      (zhyang@ucdavis.edu)                            **
**   Zhao Cheng       (zcheng@ucdavis.edu)                            **
**                                                                    **
** ****************************************************************** */
//
#include <tcl.h>
#include <string.h>
#include <Logging.h>
#include <Parsing.h>
#include <BasicModelBuilder.h>

#include <PressureDependentElastic3D.h>

#include <PlaneStressMaterial.h>
#include <PlateFiberMaterial.h>
#include <BeamFiberMaterial.h>

#include <PressureIndependMultiYield.h>
#include <PressureDependMultiYield.h>
#include <PressureDependMultiYield02.h>
#include <FluidSolidPorousMaterial.h>
// #include <isotropy.h>

// #include <Template3Dep.h>
// #include <NewTemplate3Dep.h>
// #include <FiniteDeformationElastic3D.h>
// #include <FiniteDeformationEP3D.h>


Tcl_CmdProc TclCommand_newPlasticMaterial;
Tcl_CmdProc TclCommand_newElasticMaterial;
// Tcl_CmdProc TclCommand_newIsotropicMaterial;



int
TclCommand_addMaterial(ClientData clientData, Tcl_Interp* interp, 
                        Tcl_Size argc, TCL_Char** const argv)
{
  static
  std::unordered_map<std::string, Tcl_CmdProc*> MaterialLibrary = {
    {"ElasticIsotropic",          TclCommand_newElasticMaterial},
    {"Elastic",                   TclCommand_newElasticMaterial},
    {"Isotropic",                 TclCommand_newElasticMaterial},
    {"J2",                        TclCommand_newPlasticMaterial},
    {"J2Simplified",              TclCommand_newPlasticMaterial},
    {"J2BeamFiber",               TclCommand_newPlasticMaterial},
  };


  if (argc < 2) {
    opserr << OpenSees::PromptValueError
           << "missing argument type"
           << "\n";
    return TCL_ERROR;
  }


  auto cmd = MaterialLibrary.find(std::string(argv[1]));
  if (cmd != MaterialLibrary.end())
    return (*cmd->second)(clientData, interp, argc, &argv[0]);


#if 0
  // Check argv[1] for ND material type

  // Pointer to an ND material that will be added to the model builder
  NDMaterial* theMaterial = nullptr;

  if (strcmp(argv[1], "PressureDependentElastic3D") == 0) {
    //Jul. 07, 2001 Boris Jeremic & ZHaohui Yang jeremic|zhyang@ucdavis.edu
    // Pressure dependent elastic material
    if (argc < 6) {
      opserr << "WARNING insufficient arguments\n";
      opserr << "Want: nDMaterial PressureDependentElastic3D tag? E? v? rho?" << "\n";
      return TCL_ERROR;
    }

    int tag     = 0;
    double E    = 0.0;
    double v    = 0.0;
    double rho  = 0.0;
    double expp = 0.0;
    double prp  = 0.0;
    double pop  = 0.0;

    if (Tcl_GetInt(interp, argv[2], &tag) != TCL_OK) {
      opserr << "WARNING invalid PressureDependentElastic3D tag" << "\n";
      return TCL_ERROR;
    }

    if (Tcl_GetDouble(interp, argv[3], &E) != TCL_OK) {
      opserr << "WARNING invalid E\n";
      opserr << "nDMaterial PressureDependentElastic3D: E" << tag << "\n";
      return TCL_ERROR;
    }

    if (Tcl_GetDouble(interp, argv[4], &v) != TCL_OK) {
      opserr << "WARNING invalid v\n";
      opserr << "nDMaterial PressureDependentElastic3D: v" << tag << "\n";
      return TCL_ERROR;
    }

    if (Tcl_GetDouble(interp, argv[5], &rho) != TCL_OK) {
      opserr << "WARNING invalid v\n";
      opserr << "nDMaterial PressureDependentElastic3D: rho" << tag << "\n";
      return TCL_ERROR;
    }

    if (argc == 6) {
      theMaterial = new PressureDependentElastic3D(tag, E, v, rho);
    } else if (argc == 7) {
      //get the exponent of the pressure sensitive elastic material)
      if (Tcl_GetDouble(interp, argv[6], &expp) != TCL_OK) {
        opserr << "WARNING invalid v\n";
        opserr << "nDMaterial PressureDependentElastic3D: " << tag << "\n";
        return TCL_ERROR;
      }
      theMaterial = new PressureDependentElastic3D(tag, E, v, rho, expp);
    } else if (argc == 8) {
      //get the exponent pressure of the pressure sensitive elastic material)
      if (Tcl_GetDouble(interp, argv[6], &expp) != TCL_OK) {
        opserr << "WARNING invalid v\n";
        opserr << "nDMaterial PressureDependentElastic3D: expp" << tag << "\n";
        return TCL_ERROR;
      }
      //get the reference pressure of the pressure sensitive elastic material)
      if (Tcl_GetDouble(interp, argv[7], &prp) != TCL_OK) {
        opserr << "WARNING invalid v\n";
        opserr << "nDMaterial PressureDependentElastic3D: prp " << tag << "\n";
        return TCL_ERROR;
      }
      theMaterial = new PressureDependentElastic3D(tag, E, v, rho, expp, prp);
    } else if (argc >= 9) {
      //get the exponent of the pressure sensitive elastic material)
      if (Tcl_GetDouble(interp, argv[6], &expp) != TCL_OK) {
        opserr << "WARNING invalid v\n";
        opserr << "nDMaterial PressureDependentElastic3D: expp" << tag << "\n";
        return TCL_ERROR;
      }
      //get the reference pressure of the pressure sensitive elastic material)
      if (Tcl_GetDouble(interp, argv[7], &prp) != TCL_OK) {
        opserr << "WARNING invalid v\n";
        opserr << "nDMaterial PressureDependentElastic3D: prp" << tag << "\n";
        return TCL_ERROR;
      }
      //get the cutoff pressure po of the pressure sensitive elastic material)
      if (Tcl_GetDouble(interp, argv[8], &pop) != TCL_OK) {
        opserr << "WARNING invalid v\n";
        opserr << "nDMaterial PressureDependentElastic3D: pop" << tag << "\n";
        return TCL_ERROR;
      }
      theMaterial = new PressureDependentElastic3D(tag, E, v, rho, expp, prp, pop);
    }

  }


  else if ((strcmp(argv[1], "MultiaxialCyclicPlasticity") == 0) || 
           (strcmp(argv[1], "MCP") == 0)) {

    //
    //  MultiAxialCyclicPlasticity Model   by Gang Wang
    //
    //  nDMaterial MultiaxialCyclicPlasticity $tag, $rho, $K, $G,
    //      $Su , $Ho , $h, $m, $beta, $KCoeff
    // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    if (argc < 12) {
      opserr << "WARNING insufficient arguments\n";
      opserr << "Want: nDMaterial MultiaxialCyclicPlasticity tag? rho? K? G? Su? Ho? h? m? beta? "
                "KCoeff? <eta?>"
             << "\n";
      return TCL_ERROR;
    }

    int tag;
    double K, G, rho, Su, Ho, h, m, beta, Kcoeff;
    double eta = 0.0;

    if (Tcl_GetInt(interp, argv[2], &tag) != TCL_OK) {
      opserr << "WARNING invalid MultiaxialCyclicPlasticity tag" << "\n";
      return TCL_ERROR;
    }

    if (Tcl_GetDouble(interp, argv[3], &rho) != TCL_OK) {
      opserr << "WARNING invalid rho\n";
      opserr << "nDMaterial MultiaxialCyclicPlasticity: " << tag << "\n";
      return TCL_ERROR;
    }

    if (Tcl_GetDouble(interp, argv[4], &K) != TCL_OK) {
      opserr << "WARNING invalid K\n";
      opserr << "nDMaterial MultiaxialCyclicPlasticity: " << tag << "\n";
      return TCL_ERROR;
    }

    if (Tcl_GetDouble(interp, argv[5], &G) != TCL_OK) {
      opserr << "WARNING invalid G\n";
      opserr << "nDMaterial MultiaxialCyclicPlasticity: " << tag << "\n";
      return TCL_ERROR;
    }


    if (Tcl_GetDouble(interp, argv[6], &Su) != TCL_OK) {
      opserr << "WARNING invalid alpha1\n";
      opserr << "nDMaterial MultiaxialCyclicPlasticity: " << tag << "\n";
      return TCL_ERROR;
    }

    if (Tcl_GetDouble(interp, argv[7], &Ho) != TCL_OK) {
      opserr << "WARNING invalid Ho\n";
      opserr << "nDMaterial MultiaxialCyclicPlasticity: " << tag << "\n";
      return TCL_ERROR;
    }

    if (Tcl_GetDouble(interp, argv[8], &h) != TCL_OK) {
      opserr << "WARNING invalid h\n";
      opserr << "nDMaterial MultiaxialCyclicPlasticity: " << tag << "\n";
      return TCL_ERROR;
    }

    if (Tcl_GetDouble(interp, argv[9], &m) != TCL_OK) {
      opserr << "WARNING invalid m\n";
      opserr << "nDMaterial MultiaxialCyclicPlasticity: " << tag << "\n";
      return TCL_ERROR;
    }

    if (Tcl_GetDouble(interp, argv[10], &beta) != TCL_OK) {
      opserr << "WARNING invalid beta\n";
      opserr << "nDMaterial MultiaxialCyclicPlasticity: " << tag << "\n";
      return TCL_ERROR;
    }
    if (Tcl_GetDouble(interp, argv[11], &Kcoeff) != TCL_OK) {
      opserr << "WARNING invalid Kcoeff\n";
      opserr << "nDMaterial MultiaxialCyclicPlasticity: " << tag << "\n";
      return TCL_ERROR;
    }


    if (argc > 12 && Tcl_GetDouble(interp, argv[12], &eta) != TCL_OK) {
      opserr << "WARNING invalid eta\n";
      opserr << "nDMaterial MultiaxialCyclicPlasticity: " << tag << "\n";
      return TCL_ERROR;
    }

    theMaterial =
        new MultiaxialCyclicPlasticity(tag, 0, rho, K, G, Su, Ho, h, m, beta, Kcoeff, eta);
  }


  // Pressure Independend Multi-yield, by ZHY
  else if (strcmp(argv[1], "PressureIndependMultiYield") == 0) {
    const int numParam = 6;
    const int totParam = 10;
    int tag;
    double param[totParam];
    param[6] = 0.0;
    param[7] = 100.;
    param[8] = 0.0;
    param[9] = 20;

    char* arg[] = {"nd",
                   "rho",
                   "refShearModul",
                   "refBulkModul",
                   "cohesi",
                   "peakShearStra",
                   "frictionAng (=0)",
                   "refPress (=100)",
                   "pressDependCoe (=0.0)",
                   "numberOfYieldSurf (=20)"};
    if (argc < (3 + numParam)) {
      opserr << "WARNING insufficient arguments\n";
      opserr << "Want: nDMaterial PressureIndependMultiYield tag? " << arg[0];
      opserr << "? "
             << "\n";
      opserr << arg[1] << "? " << arg[2] << "? " << arg[3] << "? "
             << "\n";
      opserr << arg[4] << "? " << arg[5] << "? " << arg[6] << "? "
             << "\n";
      opserr << arg[7] << "? " << arg[8] << "? " << arg[9] << "? " << "\n";
      return TCL_ERROR;
    }

    if (Tcl_GetInt(interp, argv[2], &tag) != TCL_OK) {
      opserr << "WARNING invalid PressureIndependMultiYield tag" << "\n";
      return TCL_ERROR;
    }

    for (int i = 3; (i < argc && i < 13); i++)
      if (Tcl_GetDouble(interp, argv[i], &param[i - 3]) != TCL_OK) {
        opserr << "WARNING invalid " << arg[i - 3] << "\n";
        opserr << "nDMaterial PressureIndependMultiYield: " << tag << "\n";
        return TCL_ERROR;
      }

    static double* gredu = 0;
    // user defined yield surfaces
    if (param[9] < 0 && param[9] > -40) {
      param[9] = -int(param[9]);
      gredu    = new double[int(2 * param[9])];
      for (int i = 0; i < 2 * param[9]; i++)
        if (Tcl_GetDouble(interp, argv[i + 13], &gredu[i]) != TCL_OK) {
          opserr << "WARNING invalid " << arg[i - 3] << "\n";
          opserr << "nDMaterial PressureIndependMultiYield: " << tag << "\n";
          return TCL_ERROR;
        }
    }

    PressureIndependMultiYield* temp =
        new PressureIndependMultiYield(tag, param[0], param[1], param[2], param[3], param[4],
                                       param[5], param[6], param[7], param[8], param[9], gredu);
    theMaterial = temp;

    if (gredu != 0) {
      delete[] gredu;
      gredu = 0;
    }
  }

  // Pressure Dependend Multi-yield, by ZHY
  else if (strcmp(argv[1], "PressureDependMultiYield") == 0) {
    const int numParam = 15;
    const int totParam = 24;
    int tag;
    double param[totParam];
    param[15] = 20;
    param[16] = 0.6;
    param[17] = 0.9;
    param[18] = 0.02;
    param[19] = 0.7;
    param[20] = 101.;
    param[21] = .3;
    param[22] = 0.;
    param[23] = 1.;

    char* arg[] = {"nd",
                   "rho",
                   "refShearModul",
                   "refBulkModul",
                   "frictionAng",
                   "peakShearStra",
                   "refPress",
                   "pressDependCoe",
                   "phaseTransformAngle",
                   "contractionParam1",
                   "dilationParam1",
                   "dilationParam2",
                   "liquefactionParam1",
                   "liquefactionParam2",
                   "liquefactionParam4",
                   "numberOfYieldSurf (=20)",
                   "e (=0.6)",
                   "volLimit1 (=0.9)",
                   "volLimit2 (=0.02)",
                   "volLimit3 (=0.7)",
                   "Atmospheric pressure (=101)",
                   "cohesi (=.5)",
                   "Hv (=0)",
                   "Pv (=1.)"};
    if (argc < (3 + numParam)) {
      opserr << "WARNING insufficient arguments\n";
      opserr << "Want: nDMaterial PressureDependMultiYield tag? " << arg[0];
      opserr << "? "
             << "\n";
      opserr << arg[1] << "? " << arg[2] << "? " << arg[3] << "? "
             << "\n";
      opserr << arg[4] << "? " << arg[5] << "? " << arg[6] << "? "
             << "\n";
      opserr << arg[7] << "? " << arg[8] << "? " << arg[9] << "? "
             << "\n";
      opserr << arg[10] << "? " << arg[11] << "? " << arg[12] << "? "
             << "\n";
      opserr << arg[13] << "? " << arg[14] << "? " << arg[15] << "? "
             << "\n";
      opserr << arg[16] << "? " << arg[17] << "? " << arg[18] << "? "
             << "\n";
      opserr << arg[19] << "? " << arg[20] << "? " << arg[21] << "? " << "\n";
      return TCL_ERROR;
    }

    if (Tcl_GetInt(interp, argv[2], &tag) != TCL_OK) {
      opserr << "WARNING invalid PressureDependMultiYield tag" << "\n";
      return TCL_ERROR;
    }

    for (int i = 3; (i < argc && i < 19); i++)
      if (Tcl_GetDouble(interp, argv[i], &param[i - 3]) != TCL_OK) {
        opserr << "WARNING invalid " << arg[i - 3] << "\n";
        opserr << "nDMaterial PressureDependMultiYield: " << tag << "\n";
        return TCL_ERROR;
      }

    static double* gredu = 0;
    // user defined yield surfaces
    if (param[15] < 0 && param[15] > -40) {
      param[15] = -int(param[15]);
      gredu     = new double[int(2 * param[15])];

      for (int i = 0; i < 2 * param[15]; i++)
        if (Tcl_GetDouble(interp, argv[i + 19], &gredu[i]) != TCL_OK) {
          opserr << "WARNING invalid " << arg[i - 3] << "\n";
          opserr << "nDMaterial PressureIndependMultiYield: " << tag << "\n";
          return TCL_ERROR;
        }
    }

    if (gredu != 0) {
      for (int i = 19 + int(2 * param[15]); i < argc; i++)
        if (Tcl_GetDouble(interp, argv[i], &param[i - 3 - int(2 * param[15])]) != TCL_OK) {
          opserr << "WARNING invalid " << arg[i - 3 - int(2 * param[15])] << "\n";
          opserr << "nDMaterial PressureDependMultiYield: " << tag << "\n";
          return TCL_ERROR;
        }
    } else {
      for (int i = 19; i < argc; i++)
        if (Tcl_GetDouble(interp, argv[i], &param[i - 3]) != TCL_OK) {
          opserr << "WARNING invalid " << arg[i - 3 - int(2 * param[15])] << "\n";
          opserr << "nDMaterial PressureDependMultiYield: " << tag << "\n";
          return TCL_ERROR;
        }
    }

    PressureDependMultiYield* temp = new PressureDependMultiYield(
        tag, param[0], param[1], param[2], param[3], param[4], param[5], param[6], param[7],
        param[8], param[9], param[10], param[11], param[12], param[13], param[14], param[15], gredu,
        param[16], param[17], param[18], param[19], param[20], param[21], param[22], param[23]);

    theMaterial = temp;
    if (gredu != 0) {
      delete[] gredu;
      gredu = 0;
    }
  }

  // Pressure Dependend Multi-yield, by ZHY
  else if (strcmp(argv[1], "PressureDependMultiYield02") == 0) {
    const int numParam = 13;
    const int totParam = 26;
    int tag;
    double param[totParam];
    param[numParam]      = 20;
    param[numParam + 1]  = 5.0;
    param[numParam + 2]  = 3.;
    param[numParam + 3]  = 1.;
    param[numParam + 4]  = 0.;
    param[numParam + 5]  = 0.6;
    param[numParam + 6]  = 0.9;
    param[numParam + 7]  = 0.02;
    param[numParam + 8]  = 0.7;
    param[numParam + 9]  = 101.;
    param[numParam + 10] = 0.1;
    param[numParam + 11] = 0.;
    param[numParam + 12] = 1.;

    char* arg[] = {"nd",
                   "rho",
                   "refShearModul",
                   "refBulkModul",
                   "frictionAng",
                   "peakShearStra",
                   "refPress",
                   "pressDependCoe",
                   "phaseTransformAngle",
                   "contractionParam1",
                   "contractionParam3",
                   "dilationParam1",
                   "dilationParam3",
                   "numberOfYieldSurf (=20)",
                   "contractionParam2=5.0",
                   "dilationParam2=3.0",
                   "liquefactionParam1=1.0",
                   "liquefactionParam2=0.0",
                   "e (=0.6)",
                   "volLimit1 (=0.9)",
                   "volLimit2 (=0.02)",
                   "volLimit3 (=0.7)",
                   "Atmospheric pressure (=101)",
                   "cohesi (=.1)",
                   "Hv (=0)",
                   "Pv (=1.)"};
    if (argc < (3 + numParam)) {
      opserr << "WARNING insufficient arguments\n";
      opserr << "Want: nDMaterial PressureDependMultiYield02 tag? " << arg[0];
      opserr << "? "
             << "\n";
      opserr << arg[1] << "? " << arg[2] << "? " << arg[3] << "? "
             << "\n";
      opserr << arg[4] << "? " << arg[5] << "? " << arg[6] << "? "
             << "\n";
      opserr << arg[7] << "? " << arg[8] << "? " << arg[9] << "? "
             << "\n";
      opserr << arg[10] << "? " << arg[11] << "? " << arg[12] << "? "
             << "\n";
      opserr << arg[13] << "? " << arg[14] << "? " << arg[15] << "? "
             << "\n";
      opserr << arg[16] << "? " << arg[17] << "? " << arg[18] << "? "
             << "\n";
      opserr << arg[19] << "? " << arg[20] << "? " << arg[21] << "? "
             << "\n";
      opserr << arg[22] << "? " << arg[23] << "? " << "\n";
      return TCL_ERROR;
    }

    if (Tcl_GetInt(interp, argv[2], &tag) != TCL_OK) {
      opserr << "WARNING invalid PressureDependMultiYield02 tag" << "\n";
      return TCL_ERROR;
    }

    int in = 17;
    for (int i = 3; (i < argc && i < in); i++)
      if (Tcl_GetDouble(interp, argv[i], &param[i - 3]) != TCL_OK) {
        opserr << "WARNING invalid " << arg[i - 3] << "\n";
        opserr << "nDMaterial PressureDependMultiYield02: " << tag << "\n";
        return TCL_ERROR;
      }

    static double* gredu = 0;

    // user defined yield surfaces
    if (param[numParam] < 0 && param[numParam] > -100) {
      param[numParam] = -int(param[numParam]);
      gredu           = new double[int(2 * param[numParam])];

      for (int i = 0; i < 2 * param[numParam]; i++)
        if (Tcl_GetDouble(interp, argv[i + in], &gredu[i]) != TCL_OK) {
          opserr << "WARNING invalid " << arg[i - 3] << "\n";
          opserr << "nDMaterial PressureIndependMultiYield: " << tag << "\n";
          return TCL_ERROR;
        }
    }

    if (gredu != 0) {
      for (int i = in + int(2 * param[numParam]); i < argc; i++)
        if (Tcl_GetDouble(interp, argv[i], &param[i - 3 - int(2 * param[numParam])]) != TCL_OK) {
          opserr << "WARNING invalid " << arg[i - 3 - int(2 * param[numParam])] << "\n";
          opserr << "nDMaterial PressureDependMultiYield02: " << tag << "\n";
          return TCL_ERROR;
        }
    } else {
      for (int i = in; i < argc; i++)
        if (Tcl_GetDouble(interp, argv[i], &param[i - 3]) != TCL_OK) {
          opserr << "WARNING invalid " << arg[i - 3 - int(2 * param[numParam])] << "\n";
          opserr << "nDMaterial PressureDependMultiYield02: " << tag << "\n";
          return TCL_ERROR;
        }
    }


    PressureDependMultiYield02* temp = new PressureDependMultiYield02(
        tag, param[0], param[1], param[2], param[3], param[4], param[5], param[6], param[7],
        param[8], param[9], param[10], param[11], param[12], param[13], gredu, param[14], param[15],
        param[16], param[17], param[18], param[19], param[20], param[21], param[22], param[23],
        param[24], param[25]);

    theMaterial = temp;
    if (gredu != 0) {
      delete[] gredu;
      gredu = 0;
    }
  }

  // Fluid Solid Porous, by ZHY
  else if (strcmp(argv[1], "FluidSolidPorous") == 0) {

    int tag;
    double param[4];
    char* arg[] = {"nd", "soilMatTag", "combinedBulkModul", "Atmospheric pressure"};
    if (argc < 6) {
      opserr << "WARNING insufficient arguments\n";
      opserr << "Want: nDMaterial FluidSolidPorous tag? " << arg[0];
      opserr << "? "
             << "\n";
      opserr << arg[1] << "? " << arg[2] << "? " << "\n";
      return TCL_ERROR;
    }

    if (Tcl_GetInt(interp, argv[2], &tag) != TCL_OK) {
      opserr << "WARNING invalid FluidSolidPorous tag" << "\n";
      return TCL_ERROR;
    }

    for (int i = 3; i < 6; i++)
      if (Tcl_GetDouble(interp, argv[i], &param[i - 3]) != TCL_OK) {
        opserr << "WARNING invalid " << arg[i - 3] << "\n";
        opserr << "nDMaterial FluidSolidPorous: " << tag << "\n";
        return TCL_ERROR;
      }

    NDMaterial* soil = builder->getTypedObject<NDMaterial>(param[1]);
    if (soil == 0) {
      opserr << "WARNING FluidSolidPorous: couldn't get soil material ";
      opserr << "tagged: " << param[1] << "\n";
      return TCL_ERROR;
    }

    param[3] = 101.;
    if (argc == 7) {
      if (Tcl_GetDouble(interp, argv[6], &param[3]) != TCL_OK) {
        opserr << "WARNING invalid " << arg[3] << "\n";
        opserr << "nDMaterial FluidSolidPorous: " << tag << "\n";
        return TCL_ERROR;
      }
    }

    theMaterial = new FluidSolidPorousMaterial(tag, param[0], *soil, param[2], param[3]);
  }

  else if (strcmp(argv[1], "Template3Dep") == 0) {
    theMaterial = TclModelBuilder_addTemplate3Dep(clientData, interp, argc, argv, theTclBuilder, 2);
  }

  else if (strcmp(argv[1], "NewTemplate3Dep") == 0) {
    theMaterial =
        TclModelBuilder_addNewTemplate3Dep(clientData, interp, argc, argv, theTclBuilder, 2);
  }

  else if (strcmp(argv[1], "FiniteDeformationElastic3D") == 0 ||
           strcmp(argv[1], "FDElastic3D") == 0) {
    theMaterial = TclModelBuilder_addFiniteDeformationElastic3D(clientData, interp, argc, argv,
                                                                theTclBuilder, 1);
  }

  else if (strcmp(argv[1], "FiniteDeformationEP3D") == 0 || strcmp(argv[1], "FDEP3D") == 0) {
    theMaterial =
        TclModelBuilder_addFiniteDeformationEP3D(clientData, interp, argc, argv, theTclBuilder, 2);
  }


  else if (strcmp(argv[1], "PlaneStressMaterial") == 0 || 
           strcmp(argv[1], "PlaneStress") == 0) {
    if (argc < 4) {
      opserr << "WARNING insufficient arguments\n";
      opserr << "Want: nDMaterial PlaneStress tag? matTag?" << "\n";
      return TCL_ERROR;
    }

    int tag, matTag;

    if (Tcl_GetInt(interp, argv[2], &tag) != TCL_OK) {
      opserr << "WARNING invalid nDMaterial PlaneStress tag" << "\n";
      return TCL_ERROR;
    }

    if (Tcl_GetInt(interp, argv[3], &matTag) != TCL_OK) {
      opserr << "WARNING invalid matTag" << "\n";
      opserr << "PlaneStress: " << matTag << "\n";
      return TCL_ERROR;
    }

    NDMaterial* threeDMaterial = builder->getTypedObject<NDMaterial>(matTag);
    if (threeDMaterial == 0) {
      opserr << "WARNING nD material does not exist\n";
      opserr << "nD material: " << matTag;
      opserr << "\nPlaneStress nDMaterial: " << tag << "\n";
      return TCL_ERROR;
    }

    theMaterial = new PlaneStressMaterial(tag, *threeDMaterial);
  }


  else if (strcmp(argv[1], "PlateFiberMaterial") == 0 || strcmp(argv[1], "PlateFiber") == 0) {

    if (argc < 4) {
      opserr << "WARNING insufficient arguments\n";
      opserr << "Want: nDMaterial PlateFiber tag? matTag?" << "\n";
      return TCL_ERROR;
    }

    int tag, matTag;

    if (Tcl_GetInt(interp, argv[2], &tag) != TCL_OK) {
      opserr << "WARNING invalid nDMaterial PlateFiber tag" << "\n";
      return TCL_ERROR;
    }

    if (Tcl_GetInt(interp, argv[3], &matTag) != TCL_OK) {
      opserr << "WARNING invalid matTag" << "\n";
      opserr << "PlateFiber: " << matTag << "\n";
      return TCL_ERROR;
    }

    NDMaterial* threeDMaterial = builder->getTypedObject<NDMaterial>(matTag);
    if (threeDMaterial == 0) {
      opserr << "WARNING nD material does not exist\n";
      opserr << "nD material: " << matTag;
      opserr << "\nPlateFiber nDMaterial: " << tag << "\n";
      return TCL_ERROR;
    }

    theMaterial = new PlateFiberMaterial(tag, *threeDMaterial);
  }

  else if (strcmp(argv[1], "BeamFiberMaterial") == 0 || 
           strcmp(argv[1], "BeamFiber") == 0) {
    if (argc < 4) {
      opserr << "WARNING insufficient arguments\n";
      opserr << "Want: nDMaterial BeamFiber tag? matTag?" << "\n";
      return TCL_ERROR;
    }

    int tag, matTag;

    if (Tcl_GetInt(interp, argv[2], &tag) != TCL_OK) {
      opserr << "WARNING invalid nDMaterial BeamFiber tag" << "\n";
      return TCL_ERROR;
    }

    if (Tcl_GetInt(interp, argv[3], &matTag) != TCL_OK) {
      opserr << "WARNING invalid matTag" << "\n";
      opserr << "BeamFiber: " << matTag << "\n";
      return TCL_ERROR;
    }

    NDMaterial* threeDMaterial = builder->getTypedObject<NDMaterial>(matTag);
    if (threeDMaterial == nullptr) {
      return TCL_ERROR;
    }

    theMaterial = new BeamFiberMaterial(tag, *threeDMaterial);
  }

  else {
    theMaterial = TclModelBuilder_addFeapMaterial(clientData, interp, argc, argv, theTclBuilder);
  }

  if (theMaterial == nullptr)
    return TCL_ERROR;

  // Now add the material to the modelBuilder
  if (theTclBuilder->addNDMaterial(*theMaterial) < 0) {
    opserr << "WARNING could not add material to the domain\n";
    opserr << *theMaterial << "\n";
    delete theMaterial;
    return TCL_ERROR;
  }

  return TCL_OK;
#endif

  return TCL_ERROR;

}




#if 0
Template3Dep* TclModelBuilder_addTemplate3Dep(ClientData clientData, Tcl_Interp* interp, int argc,
                                              TCL_Char** argv, TclModelBuilder* theTclBuilder,
                                              int eleArgStart);

NewTemplate3Dep* TclModelBuilder_addNewTemplate3Dep(ClientData clientData, Tcl_Interp* interp,
                                                    int argc, TCL_Char** argv,
                                                    TclModelBuilder* theTclBuilder,
                                                    int eleArgStart);

FiniteDeformationElastic3D*
TclModelBuilder_addFiniteDeformationElastic3D(ClientData clientData, Tcl_Interp* interp, int argc,
                                              TCL_Char** argv, TclModelBuilder* theTclBuilder,
                                              int eleArgStart);

FiniteDeformationEP3D* TclModelBuilder_addFiniteDeformationEP3D(ClientData clientData,
                                                                Tcl_Interp* interp, int argc,
                                                                TCL_Char** argv,
                                                                TclModelBuilder* theTclBuilder,
                                                                int eleArgStart);

NDMaterial* TclModelBuilder_addFeapMaterial(ClientData clientData, Tcl_Interp* interp, int argc,
                                            TCL_Char** argv, TclModelBuilder* theTclBuilder);




template <typename MatType>
int
TclCommand_newMinMaxND(ClientData clientData, Tcl_Interp* interp, int argc, const char**const argv)
{
    if (argc < 4) {
      opserr << "WARNING insufficient arguments\n";
      opserr << "Want: uniaxialMaterial MinMax tag? matTag?";
      opserr << " <-min min?> <-max max?>" << endln;
      return TCL_ERROR;
    }

    int tag, matTag;
    
    if (Tcl_GetInt(interp, argv[2], &tag) != TCL_OK) {
      opserr << "WARNING invalid uniaxialMaterial MinMax tag" << endln;
      return TCL_ERROR;               
    }

    if (Tcl_GetInt(interp, argv[3], &matTag) != TCL_OK) {
      opserr << "WARNING invalid component tag\n";
      opserr << "uniaxialMaterial MinMax: " << tag << endln;
      return TCL_ERROR;
    }

    // Search for min and max strains
    double epsmin = NEG_INF_STRAIN;
    double epsmax = POS_INF_STRAIN;
      
    for (int j = 4; j < argc; j++) {
      if (strcmp(argv[j],"-min") == 0) {
        if ((j+1) >= argc || Tcl_GetDouble (interp, argv[j+1], &epsmin) != TCL_OK) {
          opserr << "WARNING invalid min\n";
          opserr << "uniaxialMaterial MinMax: " << tag << endln;
          return TCL_ERROR;
        }
        j++;
      }
      if (strcmp(argv[j],"-max") == 0) {
        if ((j+1) >= argc || Tcl_GetDouble (interp, argv[j+1], &epsmax) != TCL_OK) {
          opserr << "WARNING invalid max\n";
          opserr << "uniaxialMaterial MinMax: " << tag << endln;
          return TCL_ERROR;
        }
        j++;
      }
    }
      
    UniaxialMaterial *theMat = theTclBuilder->getUniaxialMaterial(matTag);

    if (theMat == 0) {
      opserr << "WARNING component material does not exist\n";
      opserr << "Component material: " << matTag; 
      opserr << "\nuniaxialMaterial MinMax: " << tag << endln;
      return TCL_ERROR;
    }

    // Parsing was successful, allocate the material
    theMaterial = new MinMaxNDMaterial(tag, *theMat, epsmin, epsmax);    
}
#endif