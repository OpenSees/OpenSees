//===----------------------------------------------------------------------===//
//
//        OpenSees - Open System for Earthquake Engineering Simulation
//
//===----------------------------------------------------------------------===//
//
// Description: This file contains the function invoked when the user 
// invokes the nDMaterial command in the interpreter.
//
// Developed by:
//   Frank McKenna (fmckenna@ce.berkeley.edu)
//   Gregory L. Fenves (fenves@ce.berkeley.edu)
//   Filip C. Filippou (filippou@ce.berkeley.edu)
//
// With extensive contributions from
//   Boris Jeremic    (jeremic@ucdavis.edu)
//   Zaohui Yang      (zhyang@ucdavis.edu)
//   Zhao Cheng       (zcheng@ucdavis.edu)
//
#include "material.hpp"
#include <BasicModelBuilder.h>
#include <elementAPI.h>
#include <tcl.h>
#include <packages.h>
#include <string.h>

#include <PressureDependentElastic3D.h>
#include <J2Plasticity.h>
#include <MultiaxialCyclicPlasticity.h> // Gang Wang

#include <PlaneStressMaterial.h>
#include <PlaneStrainMaterial.h>        // Antonios Vytiniotis:
#include <PlateFiberMaterial.h>

// start Yuli Huang & Xinzheng Lu
#include <PlateRebarMaterial.h>
#include <PlateFromPlaneStressMaterial.h>
#include <ConcreteS.h>
#include <PlaneStressUserMaterial.h>
// end Yuli Huang & Xinzheng Lu

#include <CapPlasticity.h>           // Quan Gu & ZhiJian Qiu  2013
#include <SimplifiedJ2.h>            // Quan Gu & ZhiJian Qiu 2013-6-26
#include <PlaneStressSimplifiedJ2.h> // Quan Gu & ZhiJian Qiu 2013-6-26

#include <BeamFiberMaterial.h>
#include <ConcreteMcftNonLinear5.h>
#include <ConcreteMcftNonLinear7.h>

#include <PressureIndependMultiYield.h>
#include <PressureDependMultiYield.h>
#include <PressureDependMultiYield02.h>
#include <PressureDependMultiYield03.h>
#include <FluidSolidPorousMaterial.h>

#include <J2PlasticityThermal.h>                 // added by L.Jiang [SIF]
#include <PlateFiberMaterialThermal.h>           // L.Jiang [SIF]
#include <PlateFromPlaneStressMaterialThermal.h> // Liming Jiang [SIF]
#include <PlateRebarMaterialThermal.h>           // Liming Jiang [SIF]

#include <MultiYieldSurfaceClay.h>

#if 0
extern NDMaterial *Tcl_addWrapperNDMaterial(matObj *, ClientData, Tcl_Interp *,
                                            int, TCL_Char **);
#endif
extern OPS_Routine OPS_J2Plasticity;
extern OPS_Routine OPS_J2BeamFiber2dMaterial;
extern OPS_Routine OPS_J2BeamFiber3dMaterial;

#if defined(OPSDEF_Material_FEAP)
NDMaterial *TclBasicBuilder_addFeapMaterial(ClientData clientData,
                                            Tcl_Interp *interp, int argc,
                                            TCL_Char ** const argv);
#endif // _OPS_Material_FEAP


extern "C" int OPS_ResetInputNoBuilder(ClientData clientData, Tcl_Interp *interp, int cArg,
                          int mArg, TCL_Char ** const argv, Domain *domain);

typedef struct ndMaterialPackageCommand {
  char *funcName;
  void *(*funcPtr)();
  struct ndMaterialPackageCommand *next;
} NDMaterialPackageCommand;

static NDMaterialPackageCommand *theNDMaterialPackageCommands = nullptr;

int
TclCommand_addNDMaterial(ClientData clientData, Tcl_Interp *interp,
                                 int argc, TCL_Char ** const argv)
{
  G3_Runtime *rt = G3_getRuntime(interp);
  BasicModelBuilder *builder = static_cast<BasicModelBuilder*>(clientData);

  // Make sure there is a minimum number of arguments
  if (argc < 3) {
    opserr << "WARNING insufficient number of ND material arguments\n";
    opserr << "Want: nDMaterial type? tag? <specific material args>" << "\n";
    return TCL_ERROR;
  }

  OPS_ResetInputNoBuilder(clientData, interp, 2, argc, argv, 0);

  // Pointer to an ND material that will be added to the model builder
  NDMaterial *theMaterial = nullptr;

  auto tcl_cmd = material_dispatch.find(std::string(argv[1]));
  if (tcl_cmd != material_dispatch.end()) {
    void* theMat = (*tcl_cmd->second)(rt, argc, &argv[0]);
    if (theMat != 0)
      theMaterial = (NDMaterial *)theMat;
    else
      return TCL_ERROR;
  }

  // Check argv[1] for ND material type
  else if (strcmp(argv[1], "J2BeamFiber") == 0) {
    void *theMat = 0;
    if (builder->getNDM() == 2)
      theMat = OPS_J2BeamFiber2dMaterial(rt, argc, argv);
    else
      theMat = OPS_J2BeamFiber3dMaterial(rt, argc, argv);

    if (theMat != 0)
      theMaterial = (NDMaterial *)theMat;
    else
      return TCL_ERROR;
  }

  else if (strcmp(argv[1], "PressureDependentElastic3D") == 0) {
    if (argc < 6) {
      opserr << "WARNING insufficient arguments\n";
      opserr << "Want: nDMaterial PressureDependentElastic3D tag? E? v? rho?"
             << "\n";
      return TCL_ERROR;
    }

    int tag = 0;
    double E = 0.0;
    double v = 0.0;
    double rho = 0.0;
    double expp = 0.0;
    double prp = 0.0;
    double pop = 0.0;

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

    //////////////////////////////////////////////////////////////////////////////////
    if (argc == 6) {
      theMaterial = new PressureDependentElastic3D(tag, E, v, rho);
      // opserr << "nDMaterial PressureDependentElastic3D: expp =" << expp <<
      // "\n";
    }
    //////////////////////////////////////////////////////////////////////////////////
    else if (argc == 7) {
      // get the exponent of the pressure sensitive elastic material)
      if (Tcl_GetDouble(interp, argv[6], &expp) != TCL_OK) {
        opserr << "WARNING invalid v\n";
        opserr << "nDMaterial PressureDependentElastic3D: " << tag << "\n";
        return TCL_ERROR;
      }
      theMaterial = new PressureDependentElastic3D(tag, E, v, rho, expp);
      // opserr << "nDMaterial PressureDependentElastic3D: expp =" << expp <<
      // "\n";
    }
    //////////////////////////////////////////////////////////////////////////////////
    else if (argc == 8) {
      // get the exponent pressure of the pressure sensitive elastic material)
      if (Tcl_GetDouble(interp, argv[6], &expp) != TCL_OK) {
        opserr << "WARNING invalid v\n";
        opserr << "nDMaterial PressureDependentElastic3D: expp" << tag << "\n";
        return TCL_ERROR;
      }
      // get the reference pressure of the pressure sensitive elastic material)
      if (Tcl_GetDouble(interp, argv[7], &prp) != TCL_OK) {
        opserr << "WARNING invalid v\n";
        opserr << "nDMaterial PressureDependentElastic3D: prp " << tag << "\n";
        return TCL_ERROR;
      }
      // opserr << "nDMaterial ElasticIsotropic3D: prp =" << prp << "\n";
      theMaterial = new PressureDependentElastic3D(tag, E, v, rho, expp, prp);
    }
    //////////////////////////////////////////////////////////////////////////////////
    else if (argc >= 9) {
      // get the exponent of the pressure sensitive elastic material)
      if (Tcl_GetDouble(interp, argv[6], &expp) != TCL_OK) {
        opserr << "WARNING invalid v\n";
        opserr << "nDMaterial PressureDependentElastic3D: expp" << tag << "\n";
        return TCL_ERROR;
      }
      // get the reference pressure of the pressure sensitive elastic material)
      if (Tcl_GetDouble(interp, argv[7], &prp) != TCL_OK) {
        opserr << "WARNING invalid v\n";
        opserr << "nDMaterial PressureDependentElastic3D: prp" << tag << "\n";
        return TCL_ERROR;
      }
      // get the cutoff pressure po of the pressure sensitive elastic material)
      if (Tcl_GetDouble(interp, argv[8], &pop) != TCL_OK) {
        opserr << "WARNING invalid v\n";
        opserr << "nDMaterial PressureDependentElastic3D: pop" << tag << "\n";
        return TCL_ERROR;
      }
      // opserr << "nDMaterial PressureDependentElastic3D: pop =" << pop <<
      // "\n";
      theMaterial =
          new PressureDependentElastic3D(tag, E, v, rho, expp, prp, pop);
    }

  }

  // Check argv[1] for J2PlaneStrain material type
  else if ((strcmp(argv[1], "J2Plasticity") == 0) ||
           (strcmp(argv[1], "J2") == 0)) {

    void *theMat = OPS_J2Plasticity(rt, argc, argv);
    if (theMat != 0)
      theMaterial = (NDMaterial *)theMat;
    else
      return TCL_ERROR;

  }

  /////////////////////////////////////////////////////////////////
  /*
     nDmaterial PlaneStressJ2  $matTag  $G  $K  $sig0  $H_kin  $H_iso


       PlaneStress (int tag,
                                   int nd,
                                   NDMaterial &the3DMaterial);

  */

  else if ((strcmp(argv[1], "PlaneStressSimplifiedJ2") == 0)) {
    if (argc < 8) {
      opserr << "WARNING insufficient arguments\n";
      opserr << "Want: nDmaterial Simplified3DJ2  $matTag  $G  $K  $sig0  "
                "$H_kin  $H_iso"
             << "\n";
      return TCL_ERROR;
    }

    int tag;
    double K, G, sig0, H_kin, H_iso;

    if (Tcl_GetInt(interp, argv[2], &tag) != TCL_OK) {
      opserr << "WARNING invalid SimplifiedJ2 tag" << "\n";
      return TCL_ERROR;
    }

    if (Tcl_GetDouble(interp, argv[3], &G) != TCL_OK) {
      opserr << "WARNING invalid G\n";
      opserr << "nDMaterial SimplifiedJ2: " << tag << "\n";
      return TCL_ERROR;
    }

    if (Tcl_GetDouble(interp, argv[4], &K) != TCL_OK) {
      opserr << "WARNING invalid K\n";
      opserr << "nDMaterial SimplifiedJ2: " << tag << "\n";
      return TCL_ERROR;
    }

    if (Tcl_GetDouble(interp, argv[5], &sig0) != TCL_OK) {
      opserr << "WARNING invalid sig0\n";
      opserr << "nDMaterial SimplifiedJ2: " << tag << "\n";
      return TCL_ERROR;
    }

    if (Tcl_GetDouble(interp, argv[6], &H_kin) != TCL_OK) {
      opserr << "WARNING invalid H_kin\n";
      opserr << "nDMaterial SimplifiedJ2: " << tag << "\n";
      return TCL_ERROR;
    }

    if (Tcl_GetDouble(interp, argv[7], &H_iso) != TCL_OK) {
      opserr << "WARNING invalid H_iso\n";
      opserr << "nDMaterial SimplifiedJ2: " << tag << "\n";
      return TCL_ERROR;
    }

    NDMaterial *theMaterial2 =
        new SimplifiedJ2(tag, 3, G, K, sig0, H_kin, H_iso);

    theMaterial = new PlaneStressSimplifiedJ2(tag, 2, *theMaterial2);

    //	delete theMaterial2;

  }
  /////////////////////////////////////////////////////////////////

  //
  //  MultiAxialCyclicPlasticity Model   by Gang Wang
  //
  //  nDMaterial MultiaxialCyclicPlasticity $tag, $rho, $K, $G,
  //      $Su , $Ho , $h, $m, $beta, $KCoeff
  // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  else if ((strcmp(argv[1], "MultiaxialCyclicPlasticity") == 0) ||
           (strcmp(argv[1], "MCP") == 0)) {
    if (argc < 12) {
      opserr << "WARNING insufficient arguments\n";
      opserr << "Want: nDMaterial MultiaxialCyclicPlasticity tag? rho? K? G? "
                "Su? Ho? h? m? beta? KCoeff? <eta?>"
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

    theMaterial = new MultiaxialCyclicPlasticity(tag, 0, rho, K, G, Su, Ho, h,
                                                 m, beta, Kcoeff, eta);
  }

  // Pressure Independent Multi-yield, by ZHY
  else if (strcmp(argv[1], "PressureIndependMultiYield") == 0) {
    const int numParam = 6;
    const int totParam = 10;
    int tag;
    double param[totParam];
    param[6] = 0.0;
    param[7] = 100.;
    param[8] = 0.0;
    param[9] = 20;

    const char *arg[] = {"nd",
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

    for (int i = 3; (i < argc && i < 13); ++i)
      if (Tcl_GetDouble(interp, argv[i], &param[i - 3]) != TCL_OK) {
        opserr << "WARNING invalid " << arg[i - 3] << "\n";
        opserr << "nDMaterial PressureIndependMultiYield: " << tag << "\n";
        return TCL_ERROR;
      }

    static double *gredu = 0;
    // user defined yield surfaces
    if (param[9] < 0 && param[9] > -40) {
      param[9] = -int(param[9]);
      gredu = new double[int(2 * param[9])];
      for (int i = 0; i < 2 * param[9]; ++i)
        if (Tcl_GetDouble(interp, argv[i + 13], &gredu[i]) != TCL_OK) {
          opserr << "WARNING invalid " << arg[i - 3] << "\n";
          opserr << "nDMaterial PressureIndependMultiYield: " << tag << "\n";
          return TCL_ERROR;
        }
    }

    PressureIndependMultiYield *temp = new PressureIndependMultiYield(
        tag, param[0], param[1], param[2], param[3], param[4], param[5],
        param[6], param[7], param[8], param[9], gredu);
    theMaterial = temp;

    if (gredu != 0) {
      delete[] gredu;
      gredu = 0;
    }
  }

  // Pressure Independent Multi-yield, by Quan Gu
  else if (strcmp(argv[1], "MultiYieldSurfaceClay") == 0) {
    const int numParam = 6;
    const int totParam = 10;
    int tag;
    double param[totParam];
    param[6] = 0.0;
    param[7] = 100.;
    param[8] = 0.0;
    param[9] = 20;

    const char *arg[] = {"nd",
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
      opserr << "Want: nDMaterial MultiYieldSurfaceClay tag? " << arg[0];
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
      opserr << "WARNING invalid MultiYieldSurfaceClay tag" << "\n";
      return TCL_ERROR;
    }

    for (int i = 3; (i < argc && i < 13); ++i)
      if (Tcl_GetDouble(interp, argv[i], &param[i - 3]) != TCL_OK) {
        opserr << "WARNING invalid " << arg[i - 3] << "\n";
        opserr << "nDMaterial MultiYieldSurfaceClay: " << tag << "\n";
        return TCL_ERROR;
      }

    static double *gredu = 0;
    // user defined yield surfaces
    if (param[9] < 0 && param[9] > -40) {
      param[9] = -int(param[9]);
      gredu = new double[int(2 * param[9])];
      for (int i = 0; i < 2 * param[9]; ++i)
        if (Tcl_GetDouble(interp, argv[i + 13], &gredu[i]) != TCL_OK) {
          opserr << "WARNING invalid " << arg[i - 3] << "\n";
          opserr << "nDMaterial MultiYieldSurfaceClay: " << tag << "\n";
          return TCL_ERROR;
        }
    }

    MultiYieldSurfaceClay *temp = new MultiYieldSurfaceClay(
        tag, param[0], param[1], param[2], param[3], param[4], param[5],
        param[6], param[7], param[8], param[9], gredu);
    theMaterial = temp;

    if (gredu != 0) {
      delete[] gredu;
      gredu = 0;
    }
  }
  // ============

  // Pressure Dependent Multi-yield, by ZHY
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

    const char *arg[] = {"nd",
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

    for (int i = 3; (i < argc && i < 19); ++i)
      if (Tcl_GetDouble(interp, argv[i], &param[i - 3]) != TCL_OK) {
        opserr << "WARNING invalid " << arg[i - 3] << "\n";
        opserr << "nDMaterial PressureDependMultiYield: " << tag << "\n";
        return TCL_ERROR;
      }

    static double *gredu = 0;
    // user defined yield surfaces
    if (param[15] < 0 && param[15] > -40) {
      param[15] = -int(param[15]);
      gredu = new double[int(2 * param[15])];

      for (int i = 0; i < 2 * param[15]; ++i)
        if (Tcl_GetDouble(interp, argv[i + 19], &gredu[i]) != TCL_OK) {
          opserr << "WARNING invalid " << arg[i - 3] << "\n";
          opserr << "nDMaterial PressureIndependMultiYield: " << tag << "\n";
          return TCL_ERROR;
        }
    }

    if (gredu != 0) {
      for (int i = 19 + int(2 * param[15]); i < argc; ++i)
        if (Tcl_GetDouble(interp, argv[i],
                          &param[i - 3 - int(2 * param[15])]) != TCL_OK) {
          opserr << "WARNING invalid " << arg[i - 3 - int(2 * param[15])]
                 << "\n";
          opserr << "nDMaterial PressureDependMultiYield: " << tag << "\n";
          return TCL_ERROR;
        }
    } else {
      for (int i = 19; i < argc; ++i)
        if (Tcl_GetDouble(interp, argv[i], &param[i - 3]) != TCL_OK) {
          opserr << "WARNING invalid " << arg[i - 3 - int(2 * param[15])]
                 << "\n";
          opserr << "nDMaterial PressureDependMultiYield: " << tag << "\n";
          return TCL_ERROR;
        }
    }

    PressureDependMultiYield *temp = new PressureDependMultiYield(
        tag, param[0], param[1], param[2], param[3], param[4], param[5],
        param[6], param[7], param[8], param[9], param[10], param[11], param[12],
        param[13], param[14], param[15], gredu, param[16], param[17], param[18],
        param[19], param[20], param[21], param[22], param[23]);

    theMaterial = temp;
    if (gredu != 0) {
      delete[] gredu;
      gredu = 0;
    }
  }

  // Pressure Dependent Multi-yield, by ZHY
  else if (strcmp(argv[1], "PressureDependMultiYield02") == 0) {
    const int numParam = 13;
    const int totParam = 26;
    int tag;
    double param[totParam];
    param[numParam] = 20;
    param[numParam + 1] = 5.0;
    param[numParam + 2] = 3.;
    param[numParam + 3] = 1.;
    param[numParam + 4] = 0.;
    param[numParam + 5] = 0.6;
    param[numParam + 6] = 0.9;
    param[numParam + 7] = 0.02;
    param[numParam + 8] = 0.7;
    param[numParam + 9] = 101.;
    param[numParam +10] = 0.1;
    param[numParam +11] = 0.;
    param[numParam +12] = 1.;

    const char *arg[] = {"nd",
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
      return TCL_ERROR;
    }

    if (Tcl_GetInt(interp, argv[2], &tag) != TCL_OK) {
      opserr << "WARNING invalid PressureDependMultiYield02 tag" << "\n";
      return TCL_ERROR;
    }

    int in = 17;
    for (int i = 3; (i < argc && i < in); ++i)
      if (Tcl_GetDouble(interp, argv[i], &param[i - 3]) != TCL_OK) {
        opserr << "WARNING invalid " << arg[i - 3] << "\n";
        opserr << "nDMaterial PressureDependMultiYield02: " << tag << "\n";
        return TCL_ERROR;
      }

    static double *gredu = 0;

    // user defined yield surfaces
    if (param[numParam] < 0 && param[numParam] > -100) {
      param[numParam] = -int(param[numParam]);
      gredu = new double[int(2 * param[numParam])];

      for (int i = 0; i < 2 * param[numParam]; ++i)
        if (Tcl_GetDouble(interp, argv[i + in], &gredu[i]) != TCL_OK) {
          opserr << "WARNING invalid " << arg[i - 3] << "\n";
          opserr << "nDMaterial PressureIndependMultiYield: " << tag << "\n";
          return TCL_ERROR;
        }
    }

    if (gredu != 0) {
      for (int i = in + int(2 * param[numParam]); i < argc; ++i)
        if (Tcl_GetDouble(interp, argv[i],
                          &param[i - 3 - int(2 * param[numParam])]) != TCL_OK) {
          opserr << "WARNING invalid " << arg[i - 3 - int(2 * param[numParam])]
                 << "\n";
          opserr << "nDMaterial PressureDependMultiYield02: " << tag << "\n";
          return TCL_ERROR;
        }
    } else {
      for (int i = in; i < argc; ++i)
        if (Tcl_GetDouble(interp, argv[i], &param[i - 3]) != TCL_OK) {
          opserr << "WARNING invalid " << arg[i - 3 - int(2 * param[numParam])]
                 << "\n";
          opserr << "nDMaterial PressureDependMultiYield02: " << tag << "\n";
          return TCL_ERROR;
        }
    }

    PressureDependMultiYield02 *temp = new PressureDependMultiYield02(
        tag, param[0], param[1], param[2], param[3], param[4], param[5],
        param[6], param[7], param[8], param[9], param[10], param[11], param[12],
        param[13], gredu, param[14], param[15], param[16], param[17], param[18],
        param[19], param[20], param[21], param[22], param[23], param[24],
        param[25]);

    theMaterial = temp;
    if (gredu != 0) {
      delete[] gredu;
      gredu = 0;
    }
  }

  // nDMaterial PressureDependMultiYield03  $tag  $nd  $rho  $refShearModul
  // $refBulkModul $frictionAng  $peakShearStra  $refPress  $pressDependCoe
  // $PTAng $mType $ca  $cb $cc $cd $ce $da $db $dc <$noYieldSurf=20
  // <$r1 $Gs1 …>  $liquefac1=1. $liquefac2=0. $pa=101 <$c=1.73>>
  // PressureDependMultiYield03 (based on PressureDependMultiYield02).
  else if (strcmp(argv[1], "PressureDependMultiYield03") == 0) {
    const int numParam = 18;
    const int totParam = 23;
    int tag;
    double param[totParam];
    param[numParam] = 20;
    param[numParam + 1] = 1.;
    param[numParam + 2] = 0.;
    param[numParam + 3] = 101.;
    param[numParam + 4] = 1.73;

    const char *arg[] = {"nd",
                         "rho",
                         "refShearModul",
                         "refBulkModul",
                         "frictionAng",
                         "peakShearStra",
                         "refPress",
                         "pressDependCoe",
                         "phaseTransformAngle",
                         "mType",
                         "ca",
                         "cb",
                         "cc",
                         "cd",
                         "ce",
                         "da",
                         "db",
                         "dc",
                         "numberOfYieldSurf (=20)",
                         "liquefactionParam1=1.0",
                         "liquefactionParam2=0.0",
                         "Atmospheric pressure (=101)",
                         "cohesi (=1.73)"};

    if (argc < (3 + numParam)) { // 3 refers to "nDMaterial
                                 // PressureDependMultiYield03  $tag"
      opserr << "WARNING insufficient arguments\n";
      opserr << "Want: nDMaterial PressureDependMultiYield03 tag? " << arg[0];
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
      opserr << arg[19] << "? " << arg[20] << "? " << arg[21] << "? " << arg[22]
             << "? " << "\n";
      return TCL_ERROR;
    }

    if (Tcl_GetInt(interp, argv[2], &tag) != TCL_OK) {
      opserr << "WARNING invalid PressureDependMultiYield03 tag" << "\n";
      return TCL_ERROR;
    }

    int in = 22;
    for (int i = 3; (i < argc && i < in); ++i)
      if (Tcl_GetDouble(interp, argv[i], &param[i - 3]) != TCL_OK) {
        opserr << "WARNING invalid " << arg[i - 3] << "\n";
        opserr << "nDMaterial PressureDependMultiYield03: " << tag << "\n";
        return TCL_ERROR;
      }

    static double *gredu = 0;

    // user defined yield surfaces
    if (param[numParam] < 0 && param[numParam] > -100) {
      param[numParam] = -int(param[numParam]);
      gredu = new double[int(2 * param[numParam])];

      for (int i = 0; i < 2 * param[numParam]; ++i)
        if (Tcl_GetDouble(interp, argv[i + in], &gredu[i]) != TCL_OK) {
          opserr << "WARNING invalid " << arg[i - 3] << "\n";
          opserr << "nDMaterial PressureDependMultiYield03: " << tag << "\n";
          return TCL_ERROR;
        }
    }

    if (gredu != 0) {
      for (int i = in + int(2 * param[numParam]); i < argc; ++i)
        if (Tcl_GetDouble(interp, argv[i],
                          &param[i - 3 - int(2 * param[numParam])]) != TCL_OK) {
          opserr << "WARNING invalid " << arg[i - 3 - int(2 * param[numParam])]
                 << "\n";
          opserr << "nDMaterial PressureDependMultiYield03: " << tag << "\n";
          return TCL_ERROR;
        }
    } else {
      for (int i = in; i < argc; ++i)
        if (Tcl_GetDouble(interp, argv[i], &param[i - 3]) != TCL_OK) {
          opserr << "WARNING invalid " << arg[i - 3 - int(2 * param[numParam])]
                 << "\n";
          opserr << "nDMaterial PressureDependMultiYield03: " << tag << "\n";
          return TCL_ERROR;
        }
    }

    PressureDependMultiYield03 *temp = new PressureDependMultiYield03(
        tag, param[0], param[1], param[2], param[3], param[4], param[5],
        param[6], param[7], param[8], param[9], param[10], param[11], param[12],
        param[13], param[14], param[15], param[16], param[17], param[18], gredu,
        param[19], param[20], param[21], param[22]);

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
    const char *arg[] = {"nd", "soilMatTag", "combinedBulkModul",
                   "Atmospheric pressure"};
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

    for (int i = 3; i < 6; ++i)
      if (Tcl_GetDouble(interp, argv[i], &param[i - 3]) != TCL_OK) {
        opserr << "WARNING invalid " << arg[i - 3] << "\n";
        opserr << "nDMaterial FluidSolidPorous: " << tag << "\n";
        return TCL_ERROR;
      }

    NDMaterial *soil = builder->getTypedObject<NDMaterial>(param[1]);
    if (soil == nullptr)
      return TCL_ERROR;


    param[3] = 101.;
    if (argc == 7) {
      if (Tcl_GetDouble(interp, argv[6], &param[3]) != TCL_OK) {
        opserr << "WARNING invalid " << arg[3] << "\n";
        opserr << "nDMaterial FluidSolidPorous: " << tag << "\n";
        return TCL_ERROR;
      }
    }

    theMaterial =
        new FluidSolidPorousMaterial(tag, param[0], *soil, param[2], param[3]);
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

    NDMaterial *threeDMaterial = builder->getTypedObject<NDMaterial>(matTag);
    if (threeDMaterial == nullptr)
      return TCL_ERROR;

    theMaterial = new PlaneStressMaterial(tag, *threeDMaterial);
  }

  // PlaneStrainMaterial
  else if (strcmp(argv[1], "PlaneStrainMaterial") == 0 ||
           strcmp(argv[1], "PlaneStrain") == 0) {
    if (argc < 4) {
      opserr << "WARNING insufficient arguments\n";
      opserr << "Want: nDMaterial PlaneStrain tag? matTag?" << "\n";
      return TCL_ERROR;
    }

    int tag, matTag;

    if (Tcl_GetInt(interp, argv[2], &tag) != TCL_OK) {
      opserr << "WARNING invalid nDMaterial PlaneStrain tag" << "\n";
      return TCL_ERROR;
    }

    if (Tcl_GetInt(interp, argv[3], &matTag) != TCL_OK) {
      opserr << "WARNING invalid matTag" << "\n";
      opserr << "PlaneStrain: " << matTag << "\n";
      return TCL_ERROR;
    }

    NDMaterial *threeDMaterial = builder->getTypedObject<NDMaterial>(matTag);
    if (threeDMaterial == nullptr)
      return TCL_ERROR;

    theMaterial = new PlaneStrainMaterial(tag, *threeDMaterial);
  }

  // ----- Cap plasticity model ------    // Quan Gu & ZhiJian Qiu  2013

  // format nDmaterial CapPlasticity $tag $ndm $rho $G $K $X $D $W $R $lambda
  // $theta $beta $alpha $T $tol
  else if (strcmp(argv[1], "CapPlasticity") == 0) {

    int tag;
    int ndm = 3;
    double rho = 0.0;
    double G = 1.0e10;
    double K = 1.1e10;
    double X = 1.1032e8;
    double D = 4.6412e-10;
    double W = 0.42;
    double R = 4.43;
    double lambda = 7.9979e6;
    double theta = 0.11;
    double beta = 6.3816e-8;
    double alpha = 2.6614e7;
    double T = -2.0684e6;
    double tol = 1.0e-10;

    if (Tcl_GetInt(interp, argv[2], &tag) != TCL_OK) {
      opserr << "WARNING invalid CapPlasticity tag" << "\n";
      return TCL_ERROR;
    }

    if (Tcl_GetInt(interp, argv[3], &ndm) != TCL_OK) {
      opserr << "WARNING invalid CapPlasticity nd" << "\n";
      return TCL_ERROR;
    }

    if (Tcl_GetDouble(interp, argv[4], &rho) != TCL_OK) {
      opserr << "WARNING invalid CapPlasticity rho" << "\n";
      return TCL_ERROR;
    }

    if (Tcl_GetDouble(interp, argv[5], &G) != TCL_OK) {
      opserr << "WARNING invalid CapPlasticity G" << "\n";
      return TCL_ERROR;
    }

    if (Tcl_GetDouble(interp, argv[6], &K) != TCL_OK) {
      opserr << "WARNING invalid CapPlasticity K" << "\n";
      return TCL_ERROR;
    }

    if (argc > 7) {

      if (Tcl_GetDouble(interp, argv[7], &X) != TCL_OK) {
        opserr << "WARNING invalid CapPlasticity X" << "\n";
        return TCL_ERROR;
      }

      if (Tcl_GetDouble(interp, argv[8], &D) != TCL_OK) {
        opserr << "WARNING invalid CapPlasticity D" << "\n";
        return TCL_ERROR;
      }

      if (Tcl_GetDouble(interp, argv[9], &W) != TCL_OK) {
        opserr << "WARNING invalid CapPlasticity W" << "\n";
        return TCL_ERROR;
      }

      if (Tcl_GetDouble(interp, argv[10], &R) != TCL_OK) {
        opserr << "WARNING invalid CapPlasticity R" << "\n";
        return TCL_ERROR;
      }

      if (Tcl_GetDouble(interp, argv[11], &lambda) != TCL_OK) {
        opserr << "WARNING invalid CapPlasticity lambda" << "\n";
        return TCL_ERROR;
      }

      if (Tcl_GetDouble(interp, argv[12], &theta) != TCL_OK) {
        opserr << "WARNING invalid CapPlasticity theta" << "\n";
        return TCL_ERROR;
      }

      if (Tcl_GetDouble(interp, argv[13], &beta) != TCL_OK) {
        opserr << "WARNING invalid CapPlasticity beta" << "\n";
        return TCL_ERROR;
      }

      if (Tcl_GetDouble(interp, argv[14], &alpha) != TCL_OK) {
        opserr << "WARNING invalid CapPlasticity alpha" << "\n";
        return TCL_ERROR;
      }

      if (Tcl_GetDouble(interp, argv[15], &T) != TCL_OK) {
        opserr << "WARNING invalid CapPlasticity T" << "\n";
        return TCL_ERROR;
      }

      if (Tcl_GetDouble(interp, argv[16], &tol) != TCL_OK) {
        opserr << "WARNING invalid CapPlasticity tol" << "\n";
        return TCL_ERROR;
      }

    } // end if

    theMaterial = new CapPlasticity(tag, G, K, rho, X, D, W, R, lambda, theta,
                                    beta, alpha, T, ndm, tol);

  }

  /////////////////////////////////////////////////////////////////
  /*
     nDmaterial Simplified3DJ2  $matTag  $G  $K  $sig0  $H_kin  $H_iso


      SimplifiedJ2 (int tag,
                                   int nd,
                                   double G,
                                   double K,
                                   double sigmaY0,
                                   double H_kin,
                                   double H_iso);

  */

  // Check argv[1] for J2PlaneStrain material type
  else if ((strcmp(argv[1], "Simplified3DJ2") == 0) ||
           (strcmp(argv[1], "3DJ2") == 0)) {
    if (argc < 8) {
      opserr << "WARNING insufficient arguments\n";
      opserr << "Want: nDmaterial Simplified3DJ2  $matTag  $G  $K  $sig0  "
                "$H_kin  $H_iso"
             << "\n";
      return TCL_ERROR;
    }

    int tag;
    double K, G, sig0, H_kin, H_iso;

    if (Tcl_GetInt(interp, argv[2], &tag) != TCL_OK) {
      opserr << "WARNING invalid SimplifiedJ2 tag" << "\n";
      return TCL_ERROR;
    }

    if (Tcl_GetDouble(interp, argv[3], &G) != TCL_OK) {
      opserr << "WARNING invalid G\n";
      opserr << "nDMaterial SimplifiedJ2: " << tag << "\n";
      return TCL_ERROR;
    }

    if (Tcl_GetDouble(interp, argv[4], &K) != TCL_OK) {
      opserr << "WARNING invalid K\n";
      opserr << "nDMaterial SimplifiedJ2: " << tag << "\n";
      return TCL_ERROR;
    }

    if (Tcl_GetDouble(interp, argv[5], &sig0) != TCL_OK) {
      opserr << "WARNING invalid sig0\n";
      opserr << "nDMaterial SimplifiedJ2: " << tag << "\n";
      return TCL_ERROR;
    }

    if (Tcl_GetDouble(interp, argv[6], &H_kin) != TCL_OK) {
      opserr << "WARNING invalid H_kin\n";
      opserr << "nDMaterial SimplifiedJ2: " << tag << "\n";
      return TCL_ERROR;
    }

    if (Tcl_GetDouble(interp, argv[7], &H_iso) != TCL_OK) {
      opserr << "WARNING invalid H_iso\n";
      opserr << "nDMaterial SimplifiedJ2: " << tag << "\n";
      return TCL_ERROR;
    }

    theMaterial = new SimplifiedJ2(tag, 3, G, K, sig0, H_kin, H_iso);
  }

  else if (strcmp(argv[1], "PlateRebarMaterial") == 0 ||
           strcmp(argv[1], "PlateRebar") == 0) {
    if (argc < 5) {
      opserr << "WARNING insufficient arguments\n";
      opserr << "Want: nDMaterial PlateRebar tag? matTag? angle?" << "\n";
      return TCL_ERROR;
    }

    int tag, matTag;
    double angle;

    if (Tcl_GetInt(interp, argv[2], &tag) != TCL_OK) {
      opserr << "WARNING invalid nDMaterial PlateRebar tag" << "\n";
      return TCL_ERROR;
    }

    if (Tcl_GetInt(interp, argv[3], &matTag) != TCL_OK) {
      opserr << "WARNING invalid matTag" << "\n";
      opserr << "PlateRebar: " << tag << "\n";
      return TCL_ERROR;
    }

    UniaxialMaterial *theMat = builder->getTypedObject<UniaxialMaterial>(matTag);
    if (theMat == nullptr)
      return TCL_ERROR;

    if (Tcl_GetDouble(interp, argv[4], &angle) != TCL_OK) {
      opserr << "WARNING invalid angle" << "\n";
      opserr << "PlateRebar: " << tag << "\n";
      return TCL_ERROR;
    }

    theMaterial = new PlateRebarMaterial(tag, *theMat, angle);
  }

  // start Yuli Huang & Xinzheng Lu PlateFromPlaneStressMaterial
  else if (strcmp(argv[1], "PlateFromPlaneStressMaterial") == 0 ||
           strcmp(argv[1], "PlateFromPlaneStress") == 0) {
    if (argc < 5) {
      opserr << "WARNING insufficient arguments\n";
      opserr << "Want: nDMaterial PlateFromPlaneStress tag? matTag? gmod?"
             << "\n";
      return TCL_ERROR;
    }

    int tag, matTag;
    double gmod;

    if (Tcl_GetInt(interp, argv[2], &tag) != TCL_OK) {
      opserr << "WARNING invalid nDMaterial PlateFromPlaneStress tag" << "\n";
      return TCL_ERROR;
    }

    if (Tcl_GetInt(interp, argv[3], &matTag) != TCL_OK) {
      opserr << "WARNING invalid matTag" << "\n";
      opserr << "PlateFromPlaneStress: " << tag << "\n";
      return TCL_ERROR;
    }

    NDMaterial *theMat = builder->getTypedObject<NDMaterial>(matTag);
    if (theMat == nullptr)
      return TCL_ERROR;

    if (Tcl_GetDouble(interp, argv[4], &gmod) != TCL_OK) {
      opserr << "WARNING invalid gmod" << "\n";
      opserr << "PlateFromPlaneStress: " << tag << "\n";
      return TCL_ERROR;
    }

    theMaterial = new PlateFromPlaneStressMaterial(tag, *theMat, gmod);
  }

  // start Yuli Huang & Xinzheng Lu ConcreteS
  else if (strcmp(argv[1], "ConcreteS") == 0) {
    if (argc < 8) {
      opserr << "WARNING insufficient arguments\n";
      opserr << "Want: nDMaterial ConcreteS tag? E? nu? fc? ft? Es?" << "\n";
      return TCL_ERROR;
    }

    int tag;
    double E, nu, fc, ft, Es;

    if (Tcl_GetInt(interp, argv[2], &tag) != TCL_OK) {
      opserr << "WARNING invalid nDMaterial ConcreteS tag" << "\n";
      return TCL_ERROR;
    }

    if (Tcl_GetDouble(interp, argv[3], &E) != TCL_OK) {
      opserr << "WARNING invalid E" << "\n";
      opserr << "ConcreteS: " << tag << "\n";
      return TCL_ERROR;
    }

    if (Tcl_GetDouble(interp, argv[4], &nu) != TCL_OK) {
      opserr << "WARNING invalid nu" << "\n";
      opserr << "ConcreteS: " << tag << "\n";
      return TCL_ERROR;
    }

    if (Tcl_GetDouble(interp, argv[5], &fc) != TCL_OK) {
      opserr << "WARNING invalid fc" << "\n";
      opserr << "ConcreteS: " << tag << "\n";
      return TCL_ERROR;
    }

    if (Tcl_GetDouble(interp, argv[6], &ft) != TCL_OK) {
      opserr << "WARNING invalid ft" << "\n";
      opserr << "ConcreteS: " << tag << "\n";
      return TCL_ERROR;
    }

    if (Tcl_GetDouble(interp, argv[7], &Es) != TCL_OK) {
      opserr << "WARNING invalid Es" << "\n";
      opserr << "ConcreteS: " << tag << "\n";
      return TCL_ERROR;
    }

    theMaterial = new ConcreteS(tag, E, nu, fc, ft, Es);
  }
  // end Yuli Huang & Xinzheng Lu ConcreteS

  // start Yuli Huang & Xinzheng Lu PlaneStressUserMaterial
  else if (strcmp(argv[1], "PlaneStressUserMaterial") == 0) {
    if (argc < 6) {
      opserr << "WARNING insufficient arguments\n";
      opserr << "Want: nDMaterial PlaneStressUserMaterial tag? nstatevs? "
                "nprops? prop1? ... propn?"
             << "\n";
      return TCL_ERROR;
    }

    int tag, nstatevs, nprops;
    double *props, p;

    if (Tcl_GetInt(interp, argv[2], &tag) != TCL_OK) {
      opserr << "WARNING invalid nDMaterial PlaneStressUserMaterial tag"
             << "\n";
      return TCL_ERROR;
    }

    if (Tcl_GetInt(interp, argv[3], &nstatevs) != TCL_OK) {
      opserr << "WARNING invalid nDMaterial PlaneStressUserMaterial nstatevs"
             << "\n";
      return TCL_ERROR;
    }

    if (nstatevs < 1)
      nstatevs = 1;

    if (Tcl_GetInt(interp, argv[4], &nprops) != TCL_OK) {
      opserr << "WARNING invalid nDMaterial PlaneStressUserMaterial nprops"
             << "\n";
      return TCL_ERROR;
    }

    if (nprops < 1)
      nprops = 1;

    props = new double[nprops];
    for (int i = 0; i < nprops; ++i) {
      if (Tcl_GetDouble(interp, argv[5 + i], &p) != TCL_OK) {
        opserr << "WARNING invalid prop" << "\n";
        opserr << "PlaneStressUserMaterial: " << tag << "\n";
        return TCL_ERROR;
      }
      props[i] = p;
    }

    theMaterial = new PlaneStressUserMaterial(tag, nstatevs, nprops, props);

    if (props != nullptr)
      delete [] props;
  }
  // end Yuli Huang & Xinzheng Lu PlaneStressUserMaterial

  else if (strcmp(argv[1], "ConcreteMcftNonLinear7") == 0 ||
           strcmp(argv[1], "ConcreteMcftNonLinear5") == 0) {
    if (argc < 11) {
      opserr << "WARNING insufficient arguments\n";
      opserr << "Want: nDMaterial ConcreteMcftNonlinear7 tag? fcu? ecu? Ec? "
                "fcr? Esv? fyv? alphaV? RoV?"
             << "\n";
      return TCL_ERROR;
    }

    int tag = 0;
    double fcu = 0.0;
    double ecu = 0.0;
    double Ec = 0.0;
    double fcr = 0.0;
    double Esv = 0.0;
    double fyv = 0.0;
    double alphaV = 0.0;
    double RoV = 0.0;

    if (Tcl_GetInt(interp, argv[2], &tag) != TCL_OK) {
      opserr << "WARNING invalid ConcreteMcftNonlinear7: tag" << "\n";
      return TCL_ERROR;
    }

    if (Tcl_GetDouble(interp, argv[3], &fcu) != TCL_OK) {
      opserr << "WARNING invalid fcu\n";
      opserr << "nDMaterial ConcreteMcftNonLinearNonLinear5: fcu" << tag
             << "\n";
      return TCL_ERROR;
    }

    if (Tcl_GetDouble(interp, argv[4], &ecu) != TCL_OK) {
      opserr << "WARNING invalid ecu\n";
      opserr << "nDMaterial ConcreteMcftNonLinearNonLinear5: ecu" << tag
             << "\n";
      return TCL_ERROR;
    }

    if (Tcl_GetDouble(interp, argv[5], &Ec) != TCL_OK) {
      opserr << "WARNING invalid Ec\n";
      opserr << "nDMaterial ConcreteMcftNonlinear7: Ec" << tag << "\n";
      return TCL_ERROR;
    }

    if (Tcl_GetDouble(interp, argv[6], &fcr) != TCL_OK) {
      opserr << "WARNING invalid fcr\n";
      opserr << "nDMaterial ConcreteMcftNonlinear7: fcr" << tag << "\n";
      return TCL_ERROR;
    }

    if (Tcl_GetDouble(interp, argv[7], &Esv) != TCL_OK) {
      opserr << "WARNING invalid Esv\n";
      opserr << "nDMaterial ConcreteMcftNonlinear7: Esv" << tag << "\n";
      return TCL_ERROR;
    }

    if (Tcl_GetDouble(interp, argv[8], &fyv) != TCL_OK) {
      opserr << "WARNING invalid fyv\n";
      opserr << "nDMaterial ConcreteMcftNonlinear7: fyv" << tag << "\n";
      return TCL_ERROR;
    }

    if (Tcl_GetDouble(interp, argv[9], &alphaV) != TCL_OK) {
      opserr << "WARNING invalid alphaV\n";
      opserr << "nDMaterial ConcreteMcftNonlinear7: alphaV" << tag << "\n";
      return TCL_ERROR;
    }

    if (Tcl_GetDouble(interp, argv[10], &RoV) != TCL_OK) {
      opserr << "WARNING invalid RoV\n";
      opserr << "nDMaterial ConcreteMcftNonlinear7: RoV" << tag << "\n";
      return TCL_ERROR;
    }

    if (strcmp(argv[1], "ConcreteMcftNonLinear7") == 0)
      theMaterial = new ConcreteMcftNonLinear7(tag, fcu, ecu, Ec, fcr, Esv, fyv,
                                               alphaV, RoV);
    else
      theMaterial = new ConcreteMcftNonLinear5(tag, fcu, ecu, Ec, fcr, Esv, fyv,
                                               alphaV, RoV);
  }

  else if (strcmp(argv[1], "Bidirectional") == 0) {
    opserr << "nDMaterial Bidirectional is now a section model, please "
           << "change to \'section Bidirectional\'" << "\n";
    return TCL_ERROR;
  }

  //-------------------------------------------------------------
  else if (strcmp(argv[1], "PlateFromPlaneStressThermal") == 0) {
    if (argc < 5) {
      opserr << "WARNING insufficient arguments\n";
      opserr << "Want: nDMaterial PlateFromPlaneStress tag? matTag? gmod?"
             << "\n";
      return TCL_ERROR;
    }

    int tag, matTag;
    double gmod;

    if (Tcl_GetInt(interp, argv[2], &tag) != TCL_OK) {
      opserr << "WARNING invalid nDMaterial PlateFromPlaneStress tag" << "\n";
      return TCL_ERROR;
    }

    if (Tcl_GetInt(interp, argv[3], &matTag) != TCL_OK) {
      opserr << "WARNING invalid matTag" << "\n";
      opserr << "PlateFromPlaneStress: " << tag << "\n";
      return TCL_ERROR;
    }

    NDMaterial *theMat = builder->getTypedObject<NDMaterial>(matTag);
    if (theMat == nullptr)
      return TCL_ERROR;

    if (Tcl_GetDouble(interp, argv[4], &gmod) != TCL_OK) {
      opserr << "WARNING invalid gmod" << "\n";
      opserr << "PlateFromPlaneStress: " << tag << "\n";
      return TCL_ERROR;
    }

    theMaterial = new PlateFromPlaneStressMaterialThermal(tag, *theMat, gmod);
  } else if (strcmp(argv[1], "PlateRebarMaterialThermal") == 0 ||
             strcmp(argv[1], "PlateRebarThermal") == 0) {
    if (argc < 5) {
      opserr << "WARNING insufficient arguments\n";
      opserr << "Want: nDMaterial PlateRebar tag? matTag? angle?" << "\n";
      return TCL_ERROR;
    }

    int tag, matTag;
    double angle;

    if (Tcl_GetInt(interp, argv[2], &tag) != TCL_OK) {
      opserr << "WARNING invalid nDMaterial PlateRebar tag" << "\n";
      return TCL_ERROR;
    }

    if (Tcl_GetInt(interp, argv[3], &matTag) != TCL_OK) {
      opserr << "WARNING invalid matTag" << "\n";
      opserr << "PlateRebar: " << tag << "\n";
      return TCL_ERROR;
    }

    UniaxialMaterial *theMat = builder->getTypedObject<UniaxialMaterial>(matTag);
    if (theMat == nullptr)
      return TCL_ERROR;


    if (Tcl_GetDouble(interp, argv[4], &angle) != TCL_OK) {
      opserr << "WARNING invalid angle" << "\n";
      opserr << "PlateRebar: " << tag << "\n";
      return TCL_ERROR;
    }

    theMaterial = new PlateRebarMaterialThermal(tag, *theMat, angle);
//
  } else if ((strcmp(argv[1], "J2PlasticityThermal") == 0) ||
             (strcmp(argv[1], "J2Thermal") == 0)) {
    if (argc < 9) {
      opserr << "WARNING insufficient arguments\n";
      opserr << "Want: nDMaterial J2PlasticityThermal tag? K? G? sig0? sigInf? "
                "delta? H? <eta?>"
             << "\n";
      return TCL_ERROR;
    }

    int tag;
    double K, G, sig0, sigInf, delta, H;
    double eta = 0.0;

    if (Tcl_GetInt(interp, argv[2], &tag) != TCL_OK) {
      opserr << "WARNING invalid J2PlasticityThermal tag" << "\n";
      return TCL_ERROR;
    }

    if (Tcl_GetDouble(interp, argv[3], &K) != TCL_OK) {
      opserr << "WARNING invalid K\n";
      opserr << "nDMaterial J2PlasticityThermal: " << tag << "\n";
      return TCL_ERROR;
    }

    if (Tcl_GetDouble(interp, argv[4], &G) != TCL_OK) {
      opserr << "WARNING invalid G\n";
      opserr << "nDMaterial J2PlasticityThermal: " << tag << "\n";
      return TCL_ERROR;
    }

    if (Tcl_GetDouble(interp, argv[5], &sig0) != TCL_OK) {
      opserr << "WARNING invalid sig0\n";
      opserr << "nDMaterial J2PlasticityThermal: " << tag << "\n";
      return TCL_ERROR;
    }

    if (Tcl_GetDouble(interp, argv[6], &sigInf) != TCL_OK) {
      opserr << "WARNING invalid sigInf\n";
      opserr << "nDMaterial J2PlasticityThermal: " << tag << "\n";
      return TCL_ERROR;
    }

    if (Tcl_GetDouble(interp, argv[7], &delta) != TCL_OK) {
      opserr << "WARNING invalid delta\n";
      opserr << "nDMaterial J2PlasticityThermal: " << tag << "\n";
      return TCL_ERROR;
    }
    if (Tcl_GetDouble(interp, argv[8], &H) != TCL_OK) {
      opserr << "WARNING invalid H\n";
      opserr << "nDMaterial J2PlasticityThermal: " << tag << "\n";
      return TCL_ERROR;
    }
    if (argc > 9 && Tcl_GetDouble(interp, argv[9], &eta) != TCL_OK) {
      opserr << "WARNING invalid eta\n";
      opserr << "nDMaterial J2PlasticityThermal: " << tag << "\n";
      return TCL_ERROR;
    }

    theMaterial =
        new J2PlasticityThermal(tag, 0, K, G, sig0, sigInf, delta, H, eta);

  } else if (strcmp(argv[1], "PlateFiberMaterialThermal") == 0 ||
             strcmp(argv[1], "PlateFiberThermal") == 0) {
    if (argc < 4) {
      opserr << "WARNING insufficient arguments\n";
      opserr << "Want: nDMaterial PlateFiberThermal tag? matTag?" << "\n";
      return TCL_ERROR;
    }

    int tag, matTag;

    if (Tcl_GetInt(interp, argv[2], &tag) != TCL_OK) {
      opserr << "WARNING invalid nDMaterial PlateFiberThermal tag" << "\n";
      return TCL_ERROR;
    }

    if (Tcl_GetInt(interp, argv[3], &matTag) != TCL_OK) {
      opserr << "WARNING invalid matTag" << "\n";
      opserr << "PlateFiberThermal: " << matTag << "\n";
      return TCL_ERROR;
    }

    NDMaterial *threeDMaterial = builder->getTypedObject<NDMaterial>(matTag);
    if (threeDMaterial == nullptr)
      return TCL_ERROR;

    theMaterial = new PlateFiberMaterialThermal(tag, *threeDMaterial);
  }
  //--------End of adding PlateFiberMaterialThermal

  // end of adding thermo-mechanical nd materials-L.Jiang[SIF]

#if defined(OPSDEF_Material_FEAP)
  else {
    theMaterial = TclBasicBuilder_addFeapMaterial(clientData, interp, argc, argv);
  }
#endif // _OPS_Material_FEAP

  if (theMaterial == 0) {
    //
    // maybe element in a class package already loaded
    //  loop through linked list of loaded functions comparing names & if find
    //  call it
    //

    NDMaterialPackageCommand *matCommands = theNDMaterialPackageCommands;
    bool found = false;
    while (matCommands != NULL && found == false) {
      if (strcmp(argv[1], matCommands->funcName) == 0) {
        theMaterial = (NDMaterial *)(*(matCommands->funcPtr))();
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
    // maybe material in a routine
    //

    char *matType = new char[strlen(argv[1]) + 1];
    strcpy(matType, argv[1]);
    matObj *matObject = OPS_GetMaterialType(matType, strlen(matType));

    delete[] matType;

    if (matObject != 0) {

      theMaterial = Tcl_addWrapperNDMaterial(matObject, clientData, interp,
                                             argc, argv);

      if (theMaterial == 0)
        delete matObject;
    }
#endif
  }

  //
  // maybe material class exists in a package yet to be loaded
  //

  if (theMaterial == 0) {

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
      NDMaterialPackageCommand *theMatCommand = new NDMaterialPackageCommand;
      theMatCommand->funcPtr = funcPtr;
      theMatCommand->funcName = matName;
      theMatCommand->next = theNDMaterialPackageCommands;
      theNDMaterialPackageCommands = theMatCommand;

      theMaterial = (NDMaterial *)(*funcPtr)();
    }
  }

  if (theMaterial == nullptr) {
    opserr << "WARNING could not create nDMaterial: " << argv[1];
    return TCL_ERROR;
  }

  // Now add the material to the modelBuilder
  if (builder->addTaggedObject<NDMaterial>(*theMaterial) != TCL_OK ) {

    opserr << "WARNING could not add material to the domain\n";
    opserr << *theMaterial << "\n";
    delete theMaterial;
    return TCL_ERROR;
  }

  return TCL_OK;
}
