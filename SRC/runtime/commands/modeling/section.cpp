//===----------------------------------------------------------------------===//
//
//        OpenSees - Open System for Earthquake Engineering Simulation
//
//===----------------------------------------------------------------------===//
//
// Description: This file contains the function invoked when the user invokes
// the section command in the interpreter.
//
// Written: rms, mhs, cmp
// Created: 07/99
//
#include <assert.h>
#include <tcl.h>
#include <runtimeAPI.h>
#include <G3_Logging.h>
#include <elementAPI.h>
#include <BasicModelBuilder.h>

#include <string.h>
#include <fstream>
#include <iostream>

#include <packages.h>

extern "C" int OPS_ResetInputNoBuilder(ClientData clientData,
                                       Tcl_Interp *interp, int cArg, int mArg,
                                       TCL_Char ** const argv, Domain *domain);

#include <FrameSection.h>
#include <ElasticMaterial.h>
#include <ElasticSection2d.h>
#include <ElasticSection3d.h>
#include <ElasticShearSection2d.h>
#include <ElasticShearSection3d.h>
#include <ElasticWarpingShearSection2d.h>
// #include <ElasticTubeSection3d.h>
#include <ParallelSection.h>
#include <FiberSection2d.h>
#include <NDFiberSection2d.h>
#include <NDFiberSection3d.h>
#include <NDFiberSectionWarping2d.h>
#include <FiberSection2dInt.h>
#include <FiberSection3d.h>
#include <FrameFiberSection3d.h>
#include <FiberSectionAsym3d.h>
//#include <FiberSectionGJ.h>


#include <LayeredShellFiberSection.h> // Yuli Huang & Xinzheng Lu

#include <ElasticPlateSection.h>
#include <ElasticMembranePlateSection.h>
#include <MembranePlateFiberSection.h>

// SectionBuilder
#include <QuadPatch.h>
#include <CircPatch.h>
#include <QuadCell.h>
#include <StraightReinfLayer.h>
#include <CircReinfLayer.h>
#include <ReinfBar.h>
#include <SectionBuilder/FiberSectionBuilder.h>

//
#include <Bidirectional.h>
#include <Elliptical2.h>
#include <Isolator2spring.h>

//--- Adding Thermo-mechanical Sections:[BEGIN]   by UoE OpenSees Group ---//
#include <FiberSection2dThermal.h>
#include <FiberSection3dThermal.h> //Added by L.Jiang [SIF] 2017
#include <MembranePlateFiberSectionThermal.h> //Added by Liming, [SIF] 2017
#include <LayeredShellFiberSectionThermal.h>  //Added by Liming, [SIF] 2017
//--- Adding Thermo-mechanical Sections: [END]   by UoE OpenSees Group ---//

//#include <McftSection2dfiber.h>

extern OPS_Routine OPS_ElasticSection;
extern OPS_Routine OPS_ElasticWarpingShearSection2d;
// extern OPS_Routine OPS_ElasticTubeSection3d;
extern OPS_Routine OPS_UniaxialSection;
extern OPS_Routine OPS_ParallelSection;
extern OPS_Routine OPS_Bidirectional;
extern OPS_Routine OPS_Elliptical2;
extern OPS_Routine OPS_ReinforcedConcreteLayeredMembraneSection; // M. J. Nunez - UChile
extern OPS_Routine OPS_LayeredMembraneSection; // M. J. Nunez - UChile

// TODO: Make OPS_Routine
extern void *OPS_ElasticMembraneSection(); // M. J. Nunez - UChile

Tcl_CmdProc TclCommand_newElasticSection;
Tcl_CmdProc TclCommand_addFiberSection;
Tcl_CmdProc TclCommand_addFiberIntSection;
Tcl_CmdProc TclCommand_addUCFiberSection;
Tcl_CmdProc TclCommand_addSectionAggregator;

// extern OPS_Routine OPS_WFSection2d;
// extern OPS_Routine OPS_RCCircularSection;
// extern OPS_Routine OPS_RCSection2d;
// extern OPS_Routine OPS_RCTBeamSection2d;
// extern OPS_Routine OPS_RCTunnelSection;
// extern OPS_Routine OPS_TubeSection;

SectionForceDeformation *
TclBasicBuilderYS_SectionCommand(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char ** const argv);

int
TclCommand_addSection(ClientData clientData, Tcl_Interp *interp,
                              int argc, TCL_Char ** const argv)
{
  assert(clientData != nullptr);
  G3_Runtime *rt = G3_getRuntime(interp);
  BasicModelBuilder *builder = static_cast<BasicModelBuilder*>(clientData);
  Domain *theDomain = builder->getDomain();

  // Make sure there is a minimum number of arguments
  if (argc < 3) {
    opserr << G3_ERROR_PROMPT << "insufficient number of section arguments\n";
    opserr << "Want: section type? tag? <specific material args>" << endln;
    return TCL_ERROR;
  }

  OPS_ResetInputNoBuilder(clientData, interp, 2, argc, argv, theDomain);

  // Pointer to a section that will be added to the model builder
  SectionForceDeformation *theSection = nullptr;

  // Check argv[1] to dispatch section type

  if (strcmp(argv[1], "Fiber") == 0 || 
      strcmp(argv[1], "fiberSec") == 0 ||
      strcmp(argv[1], "FiberFrame") == 0 ||
      strcmp(argv[1], "FrameFiber") == 0 ||
      strcmp(argv[1], "FiberSection") == 0 ||
      // Shear
      strcmp(argv[1], "NDFiber") == 0 ||
      strcmp(argv[1], "NDFiberWarping") == 0 ||
      // Thermal
      strcmp(argv[1], "FiberThermal") == 0 ||
      strcmp(argv[1], "fiberSecThermal") == 0 || 
      // Asymmetric
      strcmp(argv[1], "FiberAsym") == 0 ||
      strcmp(argv[1], "fiberSecAsym") == 0)

    return TclCommand_addFiberSection(clientData, interp, argc, argv);

  else if (strcmp(argv[1], "FiberInt") == 0) {
    // TODO
    opserr << "FiberInt is currently broken\n";
    return TCL_ERROR;
    // return TclCommand_addFiberIntSection(clientData, interp, argc, argv);
  }

  else if (strcmp(argv[1], "UCFiber") == 0)
    return TclCommand_addUCFiberSection(clientData, interp, argc, argv);

  else if (strcmp(argv[1], "Parallel") == 0) {
    SectionForceDeformation *theSection = 
                 (SectionForceDeformation*)OPS_ParallelSection(rt, argc, argv);

    if (theSection == nullptr || builder->addTaggedObject<SectionForceDeformation>(*theSection) < 0) {
      if (theSection != nullptr)
        delete theSection;
      return TCL_ERROR;
    } else
      return TCL_OK;
  }

  else if ((strcmp(argv[1], "FrameElastic") == 0) ||
           (strcmp(argv[1], "ElasticFrame") == 0)) {
    return TclCommand_newElasticSection(clientData, interp, argc, argv);
  }

  else if (strcmp(argv[1], "Elastic") == 0) {
    if (getenv("SEC"))
      return TclCommand_newElasticSection(clientData, interp, argc, argv);

    else {
      FrameSection *theSection = (FrameSection *)OPS_ElasticSection(rt, argc, argv);
      // Now add the section to the modelBuilder
      if (theSection == nullptr || builder->addTaggedObject<FrameSection>(*theSection) < 0) {
        if (theSection != nullptr)
          delete theSection;
        return TCL_ERROR;
      } else
        return TCL_OK;
    }
  }

  else if (strcmp(argv[1], "ElasticWarpingShear") == 0) {
    FrameSection *theSection = (FrameSection *)OPS_ElasticWarpingShearSection2d(rt, argc, argv);
    // Now add the section to the modelBuilder
    if (theSection == nullptr || builder->addTaggedObject<FrameSection>(*theSection) < 0) {
      if (theSection != nullptr)
        delete theSection;
      return TCL_ERROR;
    } else
      return TCL_OK;
  }

  else if (strcmp(argv[1], "Generic1D") == 0 ||
           strcmp(argv[1], "Generic1d") == 0 ||
           strcmp(argv[1], "Uniaxial") == 0) {
    FrameSection *theSection = (FrameSection *)OPS_UniaxialSection(rt, argc, argv);
    // Now add the section to the modelBuilder
    if (theSection == nullptr || builder->addTaggedObject<FrameSection>(*theSection) < 0) {
      if (theSection != nullptr)
        delete theSection;
      return TCL_ERROR;
    } else
      return TCL_OK;
  }

  else if (strcmp(argv[1], "Bidirectional") == 0) {
    void *theMat = OPS_Bidirectional(rt, argc, argv);
    if (theMat != 0)
      theSection = (SectionForceDeformation *)theMat;
    else
      return TCL_ERROR;
  }

  else if (strcmp(argv[1], "Elliptical") == 0 ||
           strcmp(argv[1], "Elliptical2") == 0) {
    void *theMat = OPS_Elliptical2(rt, argc, argv);
    if (theMat != 0)
      theSection = (SectionForceDeformation *)theMat;
    else
      return TCL_ERROR;
  }

  else if (strcmp(argv[1], "WFSection2d") == 0 ||
           strcmp(argv[1], "WSection2d") == 0) {

      opserr << "WFSection2d has been removed; use the from_aisc utility to "
             << "generate AISC sections from Python.\n";
      return TCL_ERROR;
#if 0
    void *theMat = OPS_WFSection2d(rt, argc, argv);
    if (theMat != 0)
      theSection = (SectionForceDeformation *)theMat;
    else
      return TCL_ERROR;
#endif
  }


  //
  // Membrane
  //
  else if ((strcmp(argv[1], "ReinforcedConcreteLayeredMembraneSection") == 0) || 
           (strcmp(argv[1], "RCLayeredMembraneSection") == 0) || 
           (strcmp(argv[1], "RCLMS") == 0)) {
      void* theMat = OPS_ReinforcedConcreteLayeredMembraneSection(rt, argc, argv);
      if (theMat != 0)
          theSection = (SectionForceDeformation*)theMat;
      else
          return TCL_ERROR;
  }

  else if ((strcmp(argv[1], "LayeredMembraneSection") == 0) || (strcmp(argv[1], "LMS") == 0)) {
      void* theMat = OPS_LayeredMembraneSection(rt, argc, argv);
      if (theMat != 0)
          theSection = (SectionForceDeformation*)theMat;
      else
          return TCL_ERROR;
  }

  else if (strcmp(argv[1], "ElasticMembraneSection") == 0) {
      void* theMat = OPS_ElasticMembraneSection();
      if (theMat != 0)
          theSection = (SectionForceDeformation*)theMat;
      else
          return TCL_ERROR;
  }

  else if (strcmp(argv[1], "AddDeformation") == 0 ||
           strcmp(argv[1], "Aggregator") == 0  ||
           strcmp(argv[1], "Aggregate") == 0) 
    return TclCommand_addSectionAggregator(clientData, interp, argc, argv);


  else if (strcmp(argv[1], "ElasticPlateSection") == 0) {

    if (argc < 5) {
      opserr << G3_ERROR_PROMPT << "insufficient arguments\n";
      opserr << "Want: section ElasticPlateSection tag? E? nu? h? " << endln;
      return TCL_ERROR;
    }

    int tag;
    double E, nu, h;
    if (Tcl_GetInt(interp, argv[2], &tag) != TCL_OK) {
      opserr << G3_ERROR_PROMPT << "invalid section ElasticPlateSection tag" << endln;
      return TCL_ERROR;
    }

    if (Tcl_GetDouble(interp, argv[3], &E) != TCL_OK) {
      opserr << G3_ERROR_PROMPT << "invalid E" << endln;
      opserr << "ElasticPlateSection section: " << tag << endln;
      return TCL_ERROR;
    }

    if (Tcl_GetDouble(interp, argv[4], &nu) != TCL_OK) {
      opserr << G3_ERROR_PROMPT << "invalid nu" << endln;
      opserr << "ElasticPlateSection section: " << tag << endln;
      return TCL_ERROR;
    }

    if (Tcl_GetDouble(interp, argv[5], &h) != TCL_OK) {
      opserr << G3_ERROR_PROMPT << "invalid h" << endln;
      opserr << "ElasticPlateSection section: " << tag << endln;
      return TCL_ERROR;
    }

    theSection = new ElasticPlateSection(tag, E, nu, h);
    // Now add the material to the modelBuilder
    if (builder->addTaggedObject<SectionForceDeformation>(*theSection) < 0) {
      delete theSection; // invoke the material objects destructor, otherwise mem leak
      return TCL_ERROR;
    } else
      return TCL_OK;
  }

  else if (strcmp(argv[1], "ElasticMembranePlateSection") == 0) {
    if (argc < 5) {
      opserr << G3_ERROR_PROMPT << "insufficient arguments\n";
      opserr << "Want: section ElasticMembranePlateSection tag? E? nu? h? "
                "<rho?> <Ep_mod?>"
             << endln;
      return TCL_ERROR;
    }

    int tag;
    double E, nu, h;
    double rho = 0.0;
    double Ep_mod = 1.0;

    if (Tcl_GetInt(interp, argv[2], &tag) != TCL_OK) {
      opserr << G3_ERROR_PROMPT << "invalid section ElasticMembranePlateSection tag"
             << endln;
      return TCL_ERROR;
    }

    if (Tcl_GetDouble(interp, argv[3], &E) != TCL_OK) {
      opserr << G3_ERROR_PROMPT << "invalid E" << endln;
      opserr << "ElasticMembranePlateSection section: " << tag << endln;
      return TCL_ERROR;
    }

    if (Tcl_GetDouble(interp, argv[4], &nu) != TCL_OK) {
      opserr << G3_ERROR_PROMPT << "invalid nu" << endln;
      opserr << "ElasticMembranePlateSection section: " << tag << endln;
      return TCL_ERROR;
    }

    if (Tcl_GetDouble(interp, argv[5], &h) != TCL_OK) {
      opserr << G3_ERROR_PROMPT << "invalid h" << endln;
      opserr << "ElasticMembranePlateSection section: " << tag << endln;
      return TCL_ERROR;
    }

    if (argc > 6 && Tcl_GetDouble(interp, argv[6], &rho) != TCL_OK) {
      opserr << G3_ERROR_PROMPT << "invalid rho" << endln;
      opserr << "ElasticMembranePlateSection section: " << tag << endln;
      return TCL_ERROR;
    }

    if (argc > 7 && Tcl_GetDouble(interp, argv[7], &Ep_mod) != TCL_OK) {
      opserr << G3_ERROR_PROMPT << "invalid Ep_mod" << endln;
      opserr << "ElasticMembranePlateSection section: " << tag << endln;
      return TCL_ERROR;
    }

    theSection = new ElasticMembranePlateSection(tag, E, nu, h, rho, Ep_mod);
  }

  else if (strcmp(argv[1], "PlateFiber") == 0) {
    if (argc < 5) {
      opserr << G3_ERROR_PROMPT << "insufficient arguments\n";
      opserr << "Want: section PlateFiber tag? matTag? h? " << endln;
      return TCL_ERROR;
    }

    double h;
    int tag, matTag;
    if (Tcl_GetInt(interp, argv[2], &tag) != TCL_OK) {
      opserr << G3_ERROR_PROMPT << "invalid section PlateFiber tag" << endln;
      return TCL_ERROR;
    }

    if (Tcl_GetInt(interp, argv[3], &matTag) != TCL_OK) {
      opserr << G3_ERROR_PROMPT << "invalid matTag" << endln;
      opserr << "PlateFiber section: " << matTag << endln;
      return TCL_ERROR;
    }

    if (Tcl_GetDouble(interp, argv[4], &h) != TCL_OK) {
      opserr << G3_ERROR_PROMPT << "invalid h" << endln;
      opserr << "PlateFiber section: " << tag << endln;
      return TCL_ERROR;
    }

    NDMaterial *theMaterial = builder->getTypedObject<NDMaterial>(matTag);
    if (theMaterial == nullptr) {
      opserr << G3_ERROR_PROMPT << "nD material does not exist\n";
      opserr << "nD material: " << matTag;
      opserr << "\nPlateFiber section: " << tag << endln;
      return TCL_ERROR;
    }

    theSection = new MembranePlateFiberSection(tag, h, *theMaterial);
    // Now add the material to the modelBuilder
    if (builder->addTaggedObject<SectionForceDeformation>(*theSection) < 0) {
      delete theSection;
      return TCL_ERROR;
    } else
      return TCL_OK;
  }

  // start Yuli Huang & Xinzheng Lu LayeredShellFiberSection
  else if (strcmp(argv[1], "LayeredShell") == 0) {
    if (argc < 6) {
      opserr << G3_ERROR_PROMPT << "insufficient arguments" << endln;
      opserr << "Want: section LayeredShell tag? nLayers? matTag1? h1? ... "
                "matTagn? hn? "
             << endln;
      return TCL_ERROR;
    }

    int tag, nLayers, matTag;
    double h, *thickness;
    NDMaterial **theMats;

    if (Tcl_GetInt(interp, argv[2], &tag) != TCL_OK) {
      opserr << G3_ERROR_PROMPT << "invalid section LayeredShell tag" << endln;
      return TCL_ERROR;
    }

    if (Tcl_GetInt(interp, argv[3], &nLayers) != TCL_OK) {
      opserr << G3_ERROR_PROMPT << "invalid nLayers" << endln;
      opserr << "LayeredShell section: " << tag << endln;
      return TCL_ERROR;
    }

    if (nLayers < 3) {
      opserr << "ERROR number of layers must be larger than 2" << endln;
      opserr << "LayeredShell section: " << tag << endln;
      return TCL_ERROR;
    }

    theMats = new NDMaterial *[nLayers];
    thickness = new double[nLayers];

    if (argc < 3+2*nLayers) {
      opserr << G3_ERROR_PROMPT << "Must provide " << 2*nLayers << " layers\n";
      return TCL_ERROR;
    }

    for (int iLayer = 0; iLayer < nLayers; iLayer++) {

      if (Tcl_GetInt(interp, argv[4 + 2 * iLayer], &matTag) != TCL_OK) {
        opserr << G3_ERROR_PROMPT << "invalid matTag" << endln;
        opserr << "LayeredShell section: " << tag << endln;
        return TCL_ERROR;
      }

      theMats[iLayer] = builder->getTypedObject<NDMaterial>(matTag);
      if (theMats[iLayer] == 0) {
        opserr << G3_ERROR_PROMPT << "nD material does not exist" << endln;
        opserr << "nD material: " << matTag;
        return TCL_ERROR;
      }

      if (Tcl_GetDouble(interp, argv[5 + 2 * iLayer], &h) != TCL_OK) {
        opserr << G3_ERROR_PROMPT << "invalid h" << endln;
        opserr << "LayeredShell section: " << tag << endln;
        return TCL_ERROR;
      }

      if (h < 0) {
        opserr << G3_ERROR_PROMPT << "invalid h" << endln;
        opserr << "PlateFiber section: " << tag << endln;
        return TCL_ERROR;
      }

      thickness[iLayer] = h;
    }

    theSection = new LayeredShellFiberSection(tag, nLayers, thickness, theMats);
    if (thickness != nullptr)
      delete[] thickness;
    if (theMats != 0)
      delete[] theMats;

    // Now add the material to the modelBuilder
    if (builder->addTaggedObject<SectionForceDeformation>(*theSection) < 0) {
      delete theSection;
      return TCL_ERROR;
    } else
      return TCL_OK;
  }
  // end Yuli Huang & Xinzheng Lu LayeredShellFiberSection

  //-----Thermo-mechanical shell sections added by L.Jiang [SIF]
  else if (strcmp(argv[1], "PlateFiberThermal") == 0) {
    if (argc < 5) {
      opserr << G3_ERROR_PROMPT << "insufficient arguments\n";
      opserr << "Want: section PlateFiberThermal tag? matTag? h? " << endln;
      return TCL_ERROR;
    }

    int tag, matTag;
    double h;

    if (Tcl_GetInt(interp, argv[2], &tag) != TCL_OK) {
      opserr << G3_ERROR_PROMPT << "invalid section PlateFiberThermal tag" << endln;
      return TCL_ERROR;
    }

    if (Tcl_GetInt(interp, argv[3], &matTag) != TCL_OK) {
      opserr << G3_ERROR_PROMPT << "invalid matTag" << endln;
      opserr << "PlateFiberThermal section: " << matTag << endln;
      return TCL_ERROR;
    }

    if (Tcl_GetDouble(interp, argv[4], &h) != TCL_OK) {
      opserr << G3_ERROR_PROMPT << "invalid h" << endln;
      opserr << "PlateFiberThermal section: " << tag << endln;
      return TCL_ERROR;
    }

    NDMaterial *theMaterial = builder->getTypedObject<NDMaterial>(matTag);
    if (theMaterial == nullptr) {
      opserr << G3_ERROR_PROMPT << "nD material does not exist\n";
      opserr << "nD material: " << matTag;
      opserr << "\nPlateFiberThermal section: " << tag << endln;
      return TCL_ERROR;
    }

    theSection = new MembranePlateFiberSectionThermal(tag, h, *theMaterial);
    // Now add the material to the modelBuilder
    if (builder->addTaggedObject<SectionForceDeformation>(*theSection) < 0) {
      delete theSection; // invoke the material objects destructor, otherwise mem leak
      return TCL_ERROR;
    } else
      return TCL_OK;
  }

  // LayeredShellFiberSectionThermal based on the
  // LayeredShellFiberSectionThermal by Yuli Huang & Xinzheng Lu
  else if (strcmp(argv[1], "LayeredShellThermal") == 0) {
    if (argc < 6) {
      opserr << G3_ERROR_PROMPT << "insufficient arguments" << endln;
      opserr << "Want: section LayeredShellThermal tag? nLayers? matTag1? h1? "
                "... matTagn? hn? "
             << endln;
      return TCL_ERROR;
    }

    int tag, nLayers, matTag;
    double h, *thickness;
    NDMaterial **theMats;

    if (Tcl_GetInt(interp, argv[2], &tag) != TCL_OK) {
      opserr << G3_ERROR_PROMPT << "invalid section LayeredShellThermal tag" << endln;
      return TCL_ERROR;
    }

    if (Tcl_GetInt(interp, argv[3], &nLayers) != TCL_OK) {
      opserr << G3_ERROR_PROMPT << "invalid nLayers" << endln;
      opserr << "LayeredShellThermal section: " << tag << endln;
      return TCL_ERROR;
    }

    if (nLayers < 3) {
      opserr << "ERROR number of layers must be larger than 2" << endln;
      opserr << "LayeredShellThermal section: " << tag << endln;
      return TCL_ERROR;
    }

    theMats = new NDMaterial *[nLayers];
    thickness = new double[nLayers];

    for (int iLayer = 0; iLayer < nLayers; iLayer++) {
      if (Tcl_GetInt(interp, argv[4 + 2 * iLayer], &matTag) != TCL_OK) {
        opserr << G3_ERROR_PROMPT << "invalid matTag" << endln;
        opserr << "LayeredShellThermal section: " << tag << endln;
        return TCL_ERROR;
      }

      theMats[iLayer] = builder->getTypedObject<NDMaterial>(matTag);
      if (theMats[iLayer] == 0) {
        opserr << G3_ERROR_PROMPT << "nD material does not exist" << endln;
        ;
        opserr << "nD material: " << matTag;
        opserr << "LayeredShellThermal section: " << tag << endln;
        return TCL_ERROR;
      }

      if (Tcl_GetDouble(interp, argv[5 + 2 * iLayer], &h) != TCL_OK) {
        opserr << G3_ERROR_PROMPT << "invalid h" << endln;
        opserr << "LayeredShellThermal section: " << tag << endln;
        return TCL_ERROR;
      }

      if (h < 0) {
        opserr << G3_ERROR_PROMPT << "invalid h" << endln;
        opserr << "LayeredShellThermal section: " << tag << endln;
        return TCL_ERROR;
      }

      thickness[iLayer] = h;
    }

    theSection =
        new LayeredShellFiberSectionThermal(tag, nLayers, thickness, theMats);
    if (thickness != 0)
      delete[] thickness;
    if (theMats != 0)
      delete[] theMats;

    // Now add the material to the modelBuilder
    if (builder->addTaggedObject<SectionForceDeformation>(*theSection) < 0) {
      delete theSection; // invoke the material objects destructor, otherwise mem leak
      return TCL_ERROR;
    } else
      return TCL_OK;
  }
  // end L.Jiang [SIF] added based on LayeredShellFiberSectionThermal section
  //
  else if (strcmp(argv[1], "Iso2spring") == 0) {
    if (argc < 10) {
      opserr << G3_ERROR_PROMPT << "insufficient arguments\n";
      opserr
          << "Want: section Iso2spring tag? tol? k1? Fy? k2? kv? hb? Pe? <Po?>"
          << endln;
      return TCL_ERROR;
    }

    int tag;
    double tol, k1, Fy, kb, kvo, hb, Pe;
    double Po = 0.0;

    if (Tcl_GetInt(interp, argv[2], &tag) != TCL_OK) {
      opserr << G3_ERROR_PROMPT << "invalid Iso2spring tag" << endln;
      return TCL_ERROR;
    }

    if (Tcl_GetDouble(interp, argv[3], &tol) != TCL_OK) {
      opserr << G3_ERROR_PROMPT << "invalid tol\n";
      opserr << "section Iso2spring: " << tag << endln;
      return TCL_ERROR;
    }

    if (Tcl_GetDouble(interp, argv[4], &k1) != TCL_OK) {
      opserr << G3_ERROR_PROMPT << "invalid k1\n";
      opserr << "section Iso2spring: " << tag << endln;
      return TCL_ERROR;
    }

    if (Tcl_GetDouble(interp, argv[5], &Fy) != TCL_OK) {
      opserr << G3_ERROR_PROMPT << "invalid Fy\n";
      opserr << "section Iso2spring: " << tag << endln;
      return TCL_ERROR;
    }

    if (Tcl_GetDouble(interp, argv[6], &kb) != TCL_OK) {
      opserr << G3_ERROR_PROMPT << "invalid k2\n";
      opserr << "section Iso2spring: " << tag << endln;
      return TCL_ERROR;
    }

    if (Tcl_GetDouble(interp, argv[7], &kvo) != TCL_OK) {
      opserr << G3_ERROR_PROMPT << "invalid kv\n";
      opserr << "section Iso2spring: " << tag << endln;
      return TCL_ERROR;
    }
    if (Tcl_GetDouble(interp, argv[8], &hb) != TCL_OK) {
      opserr << G3_ERROR_PROMPT << "invalid hb\n";
      opserr << "section Iso2spring: " << tag << endln;
      return TCL_ERROR;
    }

    if (Tcl_GetDouble(interp, argv[9], &Pe) != TCL_OK) {
      opserr << G3_ERROR_PROMPT << "invalid Pe\n";
      opserr << "section Iso2spring: " << tag << endln;
      return TCL_ERROR;
    }
    if (argc > 10) {
      if (Tcl_GetDouble(interp, argv[10], &Po) != TCL_OK) {
        opserr << G3_ERROR_PROMPT << "invalid Po\n";
        opserr << "section Iso2spring: " << tag << endln;
        return TCL_ERROR;
      }
    }

    theSection = new Isolator2spring(tag, tol, k1, Fy, kb, kvo, hb, Pe, Po);

    // Now add the material to the modelBuilder
    if (builder->addTaggedObject<SectionForceDeformation>(*theSection) < 0) {
      delete theSection; // invoke the material objects destructor, otherwise mem leak
      return TCL_ERROR;
    } else
      return TCL_OK;
  }

  else {
    theSection = TclBasicBuilderYS_SectionCommand(clientData, interp, argc, argv);
  }

  // Ensure we have created the Material, out of memory if got here and no
  // section
  if (theSection == nullptr) {
    opserr << G3_ERROR_PROMPT << "could not create section " << argv[1] << endln;
    return TCL_ERROR;
  }

  // Now add the material to the modelBuilder
  if (builder->addTaggedObject<SectionForceDeformation>(*theSection) < 0) {
    opserr << *theSection << endln;
    delete theSection; // invoke the material objects destructor, otherwise mem leak
    return TCL_ERROR;
  } else
    return TCL_OK;
}


struct FiberSectionConfig {
   bool isND            = false;
   bool isAsym          = false;
   bool isWarping       = false;
   bool isThermal       = false;
   bool isNew           = false; // use new FrameFiberSection class
   bool computeCentroid = true;
   double xz[2];
   double alpha;
   double density;
   bool use_density = false;
};

static SectionBuilder* 
findSectionBuilder(BasicModelBuilder* builder, Tcl_Interp *interp, int argc, const char** const argv)
{
  int tag;
  bool section_passed = false;
  for (int i = 0; i<argc; ++i) {
    if (strcmp(argv[i], "-section") == 0) {
      if (Tcl_GetInt(interp, argv[i+1], &tag) != TCL_OK) {
        opserr << G3_ERROR_PROMPT << "failed to parse section tag \"" << argv[i+1] << "\"\n";
        return nullptr;
      } else {
        section_passed = true;
        break;
      }
    }
  }

  if (!section_passed)
   if (builder->getCurrentSectionBuilder(tag) != 0) {
     return nullptr;
   }

  if (tag == -1)
    return nullptr;

  return builder->getTypedObject<SectionBuilder>(tag);

}


// build the section
// This function assumes torsion is not NULL when num==3
static int
initSectionCommands(ClientData clientData, Tcl_Interp *interp,
                    int secTag, UniaxialMaterial *theTorsion, double Ys, double Zs, 
                    double alpha, const FiberSectionConfig& options)
{
  assert(clientData != nullptr);
  BasicModelBuilder *builder = static_cast<BasicModelBuilder*>(clientData);

  // dimension of the structure (1d, 2d, or 3d)
  int ndm = builder->getNDM();

  SectionBuilder  *sbuilder = nullptr;
  FrameSection    *section  = nullptr;
  // create 2d section
  if (ndm == 2) {
    if (options.isND) {
      if (options.isWarping) {
        auto sec = new NDFiberSectionWarping2d(secTag, 30, alpha);
        sbuilder = new FiberSectionBuilder<2, NDMaterial, NDFiberSectionWarping2d>(*builder, *sec);
        section = sec;
      } else {
        auto sec = new NDFiberSection2d(secTag, options.computeCentroid);
        sbuilder = new FiberSectionBuilder<2, NDMaterial, NDFiberSection2d>(*builder, *sec);
        section = sec;
      }
    } else {
      if (options.isThermal) {
        auto sec = new FiberSection2dThermal(secTag, options.computeCentroid);
        sbuilder = new FiberSectionBuilder<2, UniaxialMaterial, FiberSection2dThermal>(*builder, *sec);
        section = sec;
      } else {
        auto sec = new FiberSection2d(secTag, options.computeCentroid);
        sbuilder = new FiberSectionBuilder<2, UniaxialMaterial, FiberSection2d>(*builder, *sec);
        section = sec;
      }
    }
  } else if (ndm == 3) {
    // This function is not called when torsion is NULL and num==3
    assert(theTorsion != nullptr);

    if (options.isND) {
      auto sec = new NDFiberSection3d(secTag,
                                      options.computeCentroid);
      sbuilder = new FiberSectionBuilder<3, NDMaterial, NDFiberSection3d>(*builder, *sec);
      section = sec;

    } else {
      if (options.isThermal) {
        auto sec = new FiberSection3dThermal(secTag, /* TODO: torsion */
                                            options.computeCentroid);
        sbuilder = new FiberSectionBuilder<3, UniaxialMaterial, FiberSection3dThermal>(*builder, *sec);
        section = sec;
      } else if (options.isAsym) {
        auto sec = new FiberSectionAsym3d(secTag, 30, theTorsion, Ys, Zs);
        sbuilder = new FiberSectionBuilder<3, UniaxialMaterial, FiberSectionAsym3d>(*builder, *sec);
        section = sec;
      } else {
        if (options.isNew) {
          auto sec = new FrameFiberSection3d(secTag, 30, *theTorsion, options.computeCentroid, 
                                             options.density, options.use_density);
          sbuilder = new FiberSectionBuilder<3, UniaxialMaterial, FrameFiberSection3d>(*builder, *sec);
          section = sec;
        } else {
          auto sec = new FiberSection3d(secTag, 30, *theTorsion, options.computeCentroid);
          sbuilder = new FiberSectionBuilder<3, UniaxialMaterial, FiberSection3d>(*builder, *sec);
          section = sec;
        }
      }
    }

  } else {
    opserr << G3_ERROR_PROMPT << "Model dimension (ndm = " << ndm
           << ") is incompatible with available frame elements\n";
    return TCL_ERROR;
  }

  if (builder->addTaggedObject<FrameSection>(*section) < 0) {
    return TCL_ERROR;
  }
  if (builder->addTypedObject<SectionBuilder>(secTag, sbuilder) < 0) {
    opserr << G3_ERROR_PROMPT << "cannot add section\n";
    return TCL_ERROR;
  }

  return TCL_OK;
}

int
TclCommand_addFiberSection(ClientData clientData, Tcl_Interp *interp, int argc,
                           TCL_Char ** const argv)
{
  assert(clientData != nullptr);
  BasicModelBuilder* builder = static_cast<BasicModelBuilder*>(clientData);

  // Check if we are being invoked from Python or Tcl
  bool openseespy = false;
  if (Tcl_GetVar(interp, "opensees::pragma::openseespy", 0) != nullptr) { 
    openseespy = true;
  }

  int ndm = builder->getNDM();

  // cmp - Check argument counts; In Tcl we require the brace argument
  //       in Python it can be omitted, but only for legacy reasons
  if (argc < 4 && !openseespy) {
    opserr << "Insufficient arguments, expected at least 4\n";
    return TCL_ERROR;
  }
  else if (argc < 3 && openseespy) {
    opserr << "Insufficient arguments, expected at least 3\n";
    return TCL_ERROR;
  }


  int secTag;
  if (Tcl_GetInt(interp, argv[2], &secTag) != TCL_OK) {
    opserr << G3_ERROR_PROMPT << "bad command - want: \nsection fiberSec secTag { "
              "\n\tpatch <patch arguments> \n\tlayer <layer arguments> \n}\n";
    return TCL_ERROR;
  }

  builder->setCurrentSectionBuilder(secTag);

  FiberSectionConfig options;
  if (strcmp(argv[1], "NDFiber") == 0)
    options.isND = true;

  if (strcmp(argv[1], "NDFiberWarping") == 0) {
    options.isND = true;
    options.isWarping = true;
  }
  else if (strcmp(argv[1], "FrameFiber") == 0 ||
             strcmp(argv[1], "FiberFrame") == 0)
    options.isNew = true;

  else if (strcmp(argv[1], "FiberThermal") == 0 ||
            strcmp(argv[1], "fiberSecThermal") == 0)
    options.isThermal = true;

  else if (strstr(argv[1], "Asym") != nullptr)
    options.isAsym    = true;


  int iarg  = 3;
  double GJ;
  UniaxialMaterial *torsion = nullptr;
  bool deleteTorsion = false;
  bool shearParsed = false;
  double Ys=0.0, Zs=0.0; // coords of shear center relative to
                         // centroid

//// Interaction parameters
//int NStrip1, NStrip2, NStrip3;
//double t1, t2, t3;

  while (iarg < argc) {

    if (strcmp(argv[iarg], "-noCentroid") == 0) {
      options.computeCentroid = false;
      iarg += 1;
    }

    else if (strcmp(argv[iarg], "-mass") == 0 && iarg + 1 < argc) {
      if (argc < iarg + 2) {
        opserr << G3_ERROR_PROMPT << "not enough -mass args need -mass mass?\n";
        return TCL_ERROR;
      }
      if (Tcl_GetDouble(interp, argv[iarg + 1], &options.density) != TCL_OK) {
        opserr << G3_ERROR_PROMPT << "invalid density";
        return TCL_ERROR;
      }
      options.use_density = true;

      iarg  += 2;
    }

    else if (strcmp(argv[iarg], "-GJ") == 0 && iarg + 1 < argc) {
      if (Tcl_GetDouble(interp, argv[iarg + 1], &GJ) != TCL_OK) {
        opserr << G3_ERROR_PROMPT << "invalid GJ";
        return TCL_ERROR;
      }
      deleteTorsion = true;
      torsion = new ElasticMaterial(0, GJ);

      iarg  += 2;
    }

    else if (strcmp(argv[iarg], "-torsion") == 0 && iarg + 1 < argc) {
      int torsionTag = 0;
      if (Tcl_GetInt(interp, argv[iarg + 1], &torsionTag) != TCL_OK) {
        opserr << G3_ERROR_PROMPT << "invalid torsionTag";
        return TCL_ERROR;
      }

      torsion = builder->getTypedObject<UniaxialMaterial>(torsionTag);
      if (torsion == nullptr) {
        opserr << G3_ERROR_PROMPT << "uniaxial material does not exist\n";
        opserr << "uniaxial material: " << torsionTag;
        opserr << "\nFiberSection3d: " << secTag << endln;
        return TCL_ERROR;
      }

      iarg += 2;
    }

    else if (strstr(argv[1], "Asym") != nullptr && !shearParsed) {
      if (iarg + 1 >= argc) {
        opserr << G3_ERROR_PROMPT << "Asym sections require shear center before fiber block.\n";
        return TCL_ERROR;
      }
      if (Tcl_GetDouble(interp, argv[iarg], &Ys) != TCL_OK) {
        opserr << G3_ERROR_PROMPT << "invalid Ys";
        return TCL_ERROR;
      }
      if (Tcl_GetDouble(interp, argv[iarg+1], &Zs) != TCL_OK) {
        opserr << G3_ERROR_PROMPT << "invalid Zs";
        return TCL_ERROR;
      }
      shearParsed = true;
      iarg  += 2;
    }

    else {
      // braces; skip and handle later
//    iarg += 1;
      break;
    }
  }

  if (torsion == nullptr && ndm == 3) {
    opserr << G3_ERROR_PROMPT << "- no torsion specified for 3D fiber section, use -GJ or "
              "-torsion\n";
    opserr << "\nFiberSection3d: " << secTag << endln;
    return TCL_ERROR;
  }

  // initialize  the fiber section (for building)                 // TODO, alpha
  if (initSectionCommands(clientData, interp, secTag, torsion, Ys, Zs, 1.0, options) != TCL_OK) {
    opserr << G3_ERROR_PROMPT << "error constructing the section\n";
    return TCL_ERROR;
  }

  //
  // Execute the commands inside the braces (fibers, patches, and reinforcing layers)
  //
  if (iarg < argc && Tcl_Eval(interp, argv[iarg]) != TCL_OK) {
    // Assume the subcommands have printed a message regarding the error
    return TCL_ERROR;
  }

  if (deleteTorsion)
    delete torsion;

  return TCL_OK;
}

int
TclCommand_addFiberIntSection(ClientData clientData, Tcl_Interp *interp,
                              int argc, TCL_Char ** const argv)
{
  assert(clientData != nullptr);
  BasicModelBuilder *builder = static_cast<BasicModelBuilder*>(clientData);
  int NDM = builder->getNDM();

  if (argc < 4)
    return TCL_ERROR;

  int secTag;
  if (Tcl_GetInt(interp, argv[2], &secTag) != TCL_OK) {
    opserr << G3_ERROR_PROMPT << "bad command - want: \nsection fiberInt secTag -GJ <GJ> { "
              "\n\tpatch <patch arguments> \n\tlayer <layer arguments> \n}\n";
    return TCL_ERROR;
  }

  builder->setCurrentSectionBuilder(secTag);


  int brace = 3; // Start of recursive parse
  double GJ = 1.0;
  bool deleteTorsion = false;
  UniaxialMaterial *torsion = 0;
  if (strcmp(argv[3], "-GJ") == 0) {
    if (Tcl_GetDouble(interp, argv[4], &GJ) != TCL_OK) {
      opserr << G3_ERROR_PROMPT << "invalid GJ";
      return TCL_ERROR;
    }
    torsion = new ElasticMaterial(0, GJ); // Is this gonna be a memory leak? MHS

    brace = 5;
  }
  int torsionTag = 0;
  if (strcmp(argv[3], "-torsion") == 0) {
    if (Tcl_GetInt(interp, argv[4], &torsionTag) != TCL_OK) {
      opserr << G3_ERROR_PROMPT << "invalid torsionTag";
      return TCL_ERROR;
    }

    torsion = builder->getTypedObject<UniaxialMaterial>(torsionTag);
    if (torsion == 0) {
      opserr << G3_ERROR_PROMPT << "uniaxial material does not exist\n";
      opserr << "uniaxial material: " << torsionTag;
      opserr << "\nFiberSection3d: " << secTag << endln;
      return TCL_ERROR;
    }

    brace = 5;
  }

  int NStrip1, NStrip2, NStrip3;
  double t1, t2, t3;

  if (strcmp(argv[3], "-NStrip") == 0) {

    if (Tcl_GetInt(interp, argv[4], &NStrip1) != TCL_OK) {
      opserr << G3_ERROR_PROMPT << "invalid NStrip1";
      return TCL_ERROR;
    }

    if (Tcl_GetDouble(interp, argv[5], &t1) != TCL_OK) {
      opserr << G3_ERROR_PROMPT << "invalid t1";
      return TCL_ERROR;
    }

    if (Tcl_GetInt(interp, argv[6], &NStrip2) != TCL_OK) {
      opserr << G3_ERROR_PROMPT << "invalid NStrip2";
      return TCL_ERROR;
    }

    if (Tcl_GetDouble(interp, argv[7], &t2) != TCL_OK) {
      opserr << G3_ERROR_PROMPT << "invalid t2";
      return TCL_ERROR;
    }

    if (Tcl_GetInt(interp, argv[8], &NStrip3) != TCL_OK) {
      opserr << G3_ERROR_PROMPT << "invalid NStrip3";
      return TCL_ERROR;
    }

    if (Tcl_GetDouble(interp, argv[9], &t3) != TCL_OK) {
      opserr << G3_ERROR_PROMPT << "invalid t3";
      return TCL_ERROR;
    }

    brace = 10; // may be 5
  }

#if 0
  // init  the fiber section (for building)                           // TODO, alpha
  if (initSectionCommands(clientData, interp, secTag, *torsion, Ys, Zs, 1.0) != TCL_OK) {
    opserr << G3_ERROR_PROMPT << "- error constructing the section\n";
    return TCL_ERROR;
  }
#endif


  // parse the information inside the braces (patches and reinforcing layers)
  if (Tcl_Eval(interp, argv[brace]) != TCL_OK) {
    opserr << G3_ERROR_PROMPT << "- error reading information in { } \n";
    return TCL_ERROR;
  }

  if (NDM == 3 && torsion == 0) {
    opserr << G3_ERROR_PROMPT << "- no torsion specified for 3D fiber section, use -GJ or "
              "-torsion\n";
    opserr << "\nFiberSectionInt3d: " << secTag << endln;
    return TCL_ERROR;
  }

#if 0 // TODO !!!
  // build the fiber section (for analysis)
  if (buildSectionInt(clientData, interp, secTag, *torsion, NStrip1, t1,
                      NStrip2, t2, NStrip3, t3) != TCL_OK) {
    opserr << G3_ERROR_PROMPT << "- error constructing the section\n";
    return TCL_ERROR;
  }
#endif

  if (deleteTorsion)
    delete torsion;

  return TCL_OK;
}

//
// add patch to fiber section
//
int
TclCommand_addPatch(ClientData clientData, 
                    Tcl_Interp *interp, 
                    int argc,
                    TCL_Char ** const argv)
{
  assert(clientData != nullptr);
  BasicModelBuilder* builder = static_cast<BasicModelBuilder*>(clientData);

  SectionBuilder* fiberSectionRepr = findSectionBuilder(builder, interp, argc, argv);
  if (fiberSectionRepr == nullptr) {
    opserr << G3_ERROR_PROMPT << "cannot retrieve section\n";
    return TCL_ERROR;
  }


  // make sure at least one other argument to contain patch type
  if (argc < 2) {
    opserr << G3_ERROR_PROMPT << "need to specify a patch type \n";
    return TCL_ERROR;
  }

  // check argv[1] for type of patch  and create the object
  if (strcmp(argv[1], "quad") == 0 || strcmp(argv[1], "quadr") == 0) {
    int numSubdivIJ, numSubdivJK, matTag;
    double vertexCoordY, vertexCoordZ;
    Matrix vertexCoords(4, 2);

    if (argc < 13) {
      opserr << G3_ERROR_PROMPT << "invalid number of parameters: patch quad matTag "
                "numSubdivIJ numSubdivJK yVertI zVertI yVertJ zVertJ yVertK "
                "zVertK yVertL zVertL\n";
      return TCL_ERROR;
    }

    int argi = 2;

    if (Tcl_GetInt(interp, argv[argi++], &matTag) != TCL_OK) {
      opserr << G3_ERROR_PROMPT << "invalid matTag: patch quad matTag numSubdivIJ "
                "numSubdivJK yVertI zVertI yVertJ zVertJ yVertK zVertK yVertL "
                "zVertL\n";
      return TCL_ERROR;
    }

    if (Tcl_GetInt(interp, argv[argi++], &numSubdivIJ) != TCL_OK) {
      opserr << G3_ERROR_PROMPT << "invalid numSubdivIJ: patch quad matTag numSubdivIJ "
                "numSubdivJK yVertI zVertI yVertJ zVertJ yVertK zVertK yVertL "
                "zVertL\n";
      return TCL_ERROR;
    }

    if (Tcl_GetInt(interp, argv[argi++], &numSubdivJK) != TCL_OK) {
      opserr << G3_ERROR_PROMPT << "invalid numSubdivJK: patch quad matTag numSubdivIJ "
                "numSubdivJK yVertI zVertI yVertJ zVertJ yVertK zVertK yVertL "
                "zVertL\n";
      return TCL_ERROR;
    }

    for (int j = 0; j < 4; j++) {
      if (Tcl_GetDouble(interp, argv[argi++], &vertexCoordY) != TCL_OK) {
        opserr << G3_ERROR_PROMPT << "invalid Coordinate y: ...yVertI zVertI yVertJ "
                  "zVertJ yVertK zVertK yVertL zVertL\n";
        return TCL_ERROR;
      }

      if (Tcl_GetDouble(interp, argv[argi++], &vertexCoordZ) != TCL_OK) {
        opserr << G3_ERROR_PROMPT << "invalid Coordinate z: ...yVertI zVertI yVertJ "
                  "zVertJ yVertK zVertK yVertL zVertL\n";
        return TCL_ERROR;
      }

      vertexCoords(j, 0) = vertexCoordY;
      vertexCoords(j, 1) = vertexCoordZ;
    }

    // Done parsing
    QuadPatch patch(matTag, numSubdivIJ, numSubdivJK, vertexCoords);
    int error = fiberSectionRepr->addPatch(patch);
    if (error != 0)
      return TCL_ERROR;

  }

  // check argv[1] for type of patch  and create the object
  else if (strcmp(argv[1], "rect") == 0 ||
           strcmp(argv[1], "rectangular") == 0) {

    int numSubdivIJ, numSubdivJK, matTag;
    double vertexCoordY, vertexCoordZ;
    Matrix vertexCoords(4, 2);

    if (argc < 9) {
      opserr << G3_ERROR_PROMPT << "invalid number of parameters: patch quad matTag "
                "numSubdivIJ numSubdivJK yVertI zVertI yVertK zVertK\n";
      return TCL_ERROR;
    }

    int argi = 2;
    if (Tcl_GetInt(interp, argv[argi++], &matTag) != TCL_OK) {
      opserr << G3_ERROR_PROMPT << "invalid matTag: patch quad matTag numSubdivIJ "
                "numSubdivJK yVertI zVertI yVertJ zVertJ yVertK zVertK yVertL "
                "zVertL\n";
      return TCL_ERROR;
    }

    if (Tcl_GetInt(interp, argv[argi++], &numSubdivIJ) != TCL_OK) {
      opserr << G3_ERROR_PROMPT << "invalid numSubdivIJ: patch quad matTag numSubdivIJ "
                "numSubdivJK yVertI zVertI yVertJ zVertJ yVertK zVertK yVertL "
                "zVertL\n";
      return TCL_ERROR;
    }

    if (Tcl_GetInt(interp, argv[argi++], &numSubdivJK) != TCL_OK) {
      opserr << G3_ERROR_PROMPT << "invalid numSubdivJK: patch quad matTag numSubdivIJ "
                "numSubdivJK yVertI zVertI yVertJ zVertJ yVertK zVertK yVertL "
                "zVertL\n";
      return TCL_ERROR;
    }

    for (int j = 0; j < 2; j++) {
      if (Tcl_GetDouble(interp, argv[argi++], &vertexCoordY) != TCL_OK) {
        opserr << G3_ERROR_PROMPT << "invalid Coordinate y: ...yVertI zVertI yVertJ "
                  "zVertJ yVertK zVertK yVertL zVertL\n";
        return TCL_ERROR;
      }

      if (Tcl_GetDouble(interp, argv[argi++], &vertexCoordZ) != TCL_OK) {
        opserr << G3_ERROR_PROMPT << "invalid Coordinate z: ...yVertI zVertI yVertJ "
                  "zVertJ yVertK zVertK yVertL zVertL\n";
        return TCL_ERROR;
      }

      vertexCoords(j * 2, 0) = vertexCoordY;
      vertexCoords(j * 2, 1) = vertexCoordZ;
    }

    vertexCoords(1, 0) = vertexCoords(2, 0);
    vertexCoords(1, 1) = vertexCoords(0, 1);
    vertexCoords(3, 0) = vertexCoords(0, 0);
    vertexCoords(3, 1) = vertexCoords(2, 1);

    // create patch
    QuadPatch patch(matTag, numSubdivIJ, numSubdivJK, vertexCoords);

    // add patch to section representation
    int error = fiberSectionRepr->addPatch(patch);
    if (error) {
      return TCL_ERROR;
    }
  }

  else if (strcmp(argv[1], "circ") == 0) {
    int numSubdivRad, numSubdivCirc, matTag;
    double yCenter, zCenter;
    Vector centerPosition(2);
    double intRad, extRad;
    double startAng, endAng;

    int argi = 2;
    if (argc < 11) {
      opserr << G3_ERROR_PROMPT << "invalid number of parameters: patch circ matTag "
                "numSubdivCirc numSubdivRad yCenter zCenter intRad extRad "
                "startAng endAng\n";
      return TCL_ERROR;
    }

    if (Tcl_GetInt(interp, argv[argi++], &matTag) != TCL_OK) {
      opserr << G3_ERROR_PROMPT << "invalid matTag: patch circ matTag numSubdivCirc "
                "numSubdivRad yCenter zCenter intRad extRad startAng endAng\n";
      return TCL_ERROR;
    }

    if (Tcl_GetInt(interp, argv[argi++], &numSubdivCirc) != TCL_OK) {
      opserr
          << G3_ERROR_PROMPT << "invalid numSubdivCirc: patch circ matTag numSubdivCirc "
             "numSubdivRad yCenter zCenter intRad extRad startAng endAng\n";
      return TCL_ERROR;
    }

    if (Tcl_GetInt(interp, argv[argi++], &numSubdivRad) != TCL_OK) {
      opserr << G3_ERROR_PROMPT << "invalid numSubdivRad: patch circ matTag numSubdivCirc "
                "numSubdivRad yCenter zCenter intRad extRad startAng endAng\n";
      return TCL_ERROR;
    }

    if (Tcl_GetDouble(interp, argv[argi++], &yCenter) != TCL_OK) {
      opserr << G3_ERROR_PROMPT << "invalid yCenter: patch circ matTag numSubdivCirc "
                "numSubdivRad yCenter zCenter intRad extRad startAng endAng\n";
      return TCL_ERROR;
    }

    if (Tcl_GetDouble(interp, argv[argi++], &zCenter) != TCL_OK) {
      opserr << G3_ERROR_PROMPT << "invalid zCenter: patch circ matTag numSubdivCirc "
                "numSubdivRad yCenter zCenter intRad extRad startAng endAng\n";
      return TCL_ERROR;
    }

    if (Tcl_GetDouble(interp, argv[argi++], &intRad) != TCL_OK) {
      opserr << G3_ERROR_PROMPT << "invalid intRad: patch circ matTag numSubdivCirc "
                "numSubdivRad yCenter zCenter intRad extRad startAng endAng\n";
      return TCL_ERROR;
    }

    if (Tcl_GetDouble(interp, argv[argi++], &extRad) != TCL_OK) {
      opserr << G3_ERROR_PROMPT << "invalid extRad: patch circ matTag numSubdivCirc "
                "numSubdivRad yCenter zCenter intRad extRad startAng endAng\n";
      return TCL_ERROR;
    }

    if (Tcl_GetDouble(interp, argv[argi++], &startAng) != TCL_OK) {
      opserr << G3_ERROR_PROMPT << "invalid startAng: patch circ matTag numSubdivCirc "
                "numSubdivRad yCenter zCenter intRad extRad startAng endAng\n";
      return TCL_ERROR;
    }

    if (Tcl_GetDouble(interp, argv[argi++], &endAng) != TCL_OK) {
      opserr << G3_ERROR_PROMPT << "invalid endAng: patch circ matTag numSubdivCirc "
                "numSubdivRad yCenter zCenter intRad extRad startAng endAng\n";
      return TCL_ERROR;
    }

    centerPosition(0) = yCenter;
    centerPosition(1) = zCenter;

    // create patch
    CircPatch patch(matTag, numSubdivCirc, numSubdivRad, centerPosition,
                    intRad, extRad, startAng, endAng);

    // add patch to section
    int error = fiberSectionRepr->addPatch(patch);
    if (error) {
      return TCL_ERROR;
    }
  }

  else {
    opserr << G3_ERROR_PROMPT << "patch type is not available\n";
    return TCL_ERROR;
  }

  return TCL_OK;
}

// add fiber to fiber section
int
TclCommand_addFiber(ClientData clientData, Tcl_Interp *interp, int argc,
                    TCL_Char ** const argv)
{
  assert(clientData != nullptr);
  BasicModelBuilder* builder = static_cast<BasicModelBuilder*>(clientData);

  if (argc < 5) {
    opserr << G3_ERROR_PROMPT << "invalid num args: fiber yLoc zLoc area matTag\n";
    return TCL_ERROR;
  }

  SectionBuilder* fiberSectionRepr = findSectionBuilder(builder, interp, argc, argv);
  if (fiberSectionRepr == nullptr) {
    opserr << G3_ERROR_PROMPT << "cannot retrieve a section builder\n";
    return TCL_ERROR;
  }

//int numFibers = fiberSectionRepr->getNumFibers();


  double yLoc, zLoc, area;
  if (Tcl_GetDouble(interp, argv[1], &yLoc) != TCL_OK) {
    opserr << G3_ERROR_PROMPT << "invalid yLoc: fiber yLoc zLoc area matTag\n";
    return TCL_ERROR;
  }
  if (Tcl_GetDouble(interp, argv[2], &zLoc) != TCL_OK) {
    opserr << G3_ERROR_PROMPT << "invalid zLoc: fiber yLoc zLoc area matTag\n";
    return TCL_ERROR;
  }
  if (Tcl_GetDouble(interp, argv[3], &area) != TCL_OK) {
    opserr << G3_ERROR_PROMPT << "invalid area: fiber yLoc zLoc area matTag\n";
    return TCL_ERROR;
  }
  int matTag;
  if (Tcl_GetInt(interp, argv[4], &matTag) != TCL_OK) {
    opserr << G3_ERROR_PROMPT << "invalid matTag: fiber yLoc zLoc area matTag\n";
    return TCL_ERROR;
  }

  //
  // Add fiber to section builder
  //
  int ndm = builder->getNDM();
  int error = 0;
  if (ndm == 2) {
    Vector pos(2);
    pos(0) = yLoc;
    pos(1) = zLoc;
    error = fiberSectionRepr->addFiber(0, matTag, area, pos);
  } else if (ndm == 3) {
    Vector pos(2);
    pos(0) = yLoc;
    pos(1) = zLoc;
    error = fiberSectionRepr->addFiber(0, matTag, area, pos);
  }

  if (error) {
    opserr << G3_ERROR_PROMPT << "cannot add patch to section\n";
    return TCL_ERROR;
  }

  return TCL_OK;
}

// add Hfiber to fiber section
int
TclCommand_addHFiber(ClientData clientData, Tcl_Interp *interp, int argc,
                     TCL_Char ** const argv)
{
  assert(clientData != nullptr);
  BasicModelBuilder* builder = static_cast<BasicModelBuilder*>(clientData);


  SectionBuilder* fiberSectionRepr = findSectionBuilder(builder, interp, argc, argv);
  if (fiberSectionRepr == nullptr) {
    opserr << G3_ERROR_PROMPT << "cannot retrieve section\n";
    return TCL_ERROR;
  }

  // make sure at least one other argument to contain patch type
  if (argc < 5) {
    opserr << G3_ERROR_PROMPT << "invalid num args: Hfiber yLoc zLoc area matTag\n";
    return TCL_ERROR;
  }


  int matHTag;
  double yHLoc, zHLoc, Harea;

  if (Tcl_GetDouble(interp, argv[1], &yHLoc) != TCL_OK) {
    opserr << G3_ERROR_PROMPT << "invalid yLoc: Hfiber yLoc zLoc area matTag\n";
    return TCL_ERROR;
  }
  if (Tcl_GetDouble(interp, argv[2], &zHLoc) != TCL_OK) {
    opserr << G3_ERROR_PROMPT << "invalid zLoc: Hfiber yLoc zLoc area matTag\n";
    return TCL_ERROR;
  }
  if (Tcl_GetDouble(interp, argv[3], &Harea) != TCL_OK) {
    opserr << G3_ERROR_PROMPT << "invalid area: Hfiber yLoc zLoc area matTag\n";
    return TCL_ERROR;
  }

  if (Tcl_GetInt(interp, argv[4], &matHTag) != TCL_OK) {
    opserr << G3_ERROR_PROMPT << "invalid matTag: Hfiber yLoc zLoc area matTag\n";
    return TCL_ERROR;
  }

  Vector fiberHPosition(2);
  fiberHPosition(0) = yHLoc;
  fiberHPosition(1) = zHLoc;

  // add patch to section builder
  int error = fiberSectionRepr->addHFiber(0, matHTag, Harea, fiberHPosition);

  if (error)
    return TCL_ERROR;

  return TCL_OK;
}

// add layers of reinforcing bars to fiber section

int
TclCommand_addReinfLayer(ClientData clientData, Tcl_Interp *interp, int argc,
                         TCL_Char ** const argv)
{
  assert(clientData != nullptr);
  BasicModelBuilder *builder = static_cast<BasicModelBuilder*>(clientData);

  SectionBuilder* fiberSectionRepr = findSectionBuilder(builder, interp, argc, argv);
  if (fiberSectionRepr == nullptr) {
    opserr << G3_ERROR_PROMPT << "cannot retrieve section\n";
    return TCL_ERROR;
  }

  // make sure at least one other argument to contain layer type
  if (argc < 2) {
    opserr << G3_ERROR_PROMPT << "need to specify a layer type \n";
    return TCL_ERROR;
  }

  // check argv[1] for type of layer and create the object
  if (strcmp(argv[1], "straight") == 0 ||
      strcmp(argv[1], "line")     == 0) {
    if (argc < 9) {
      opserr << G3_ERROR_PROMPT << "invalid number of parameters: layer straight matTag "
                "numReinfBars reinfBarArea yStartPt zStartPt yEndPt zEndPt\n";
      return TCL_ERROR;
    }

    int matTag, numReinfBars;
    double reinfBarArea;
    double yStartPt, zStartPt, yEndPt, zEndPt;

    int argi = 2;

    if (Tcl_GetInt(interp, argv[argi++], &matTag) != TCL_OK) {
      opserr << G3_ERROR_PROMPT << "invalid matTag: layer straight matTag numReinfBars "
                "reinfBarArea  yStartPt zStartPt yEndPt zEndPt\n";
      return TCL_ERROR;
    }

    if (Tcl_GetInt(interp, argv[argi++], &numReinfBars) != TCL_OK) {
      opserr << G3_ERROR_PROMPT << "invalid numReinfBars: layer straight matTag "
                "numReinfBars reinfBarArea  yStartPt zStartPt yEndPt zEndPt\n";
      return TCL_ERROR;
    }

    if (Tcl_GetDouble(interp, argv[argi++], &reinfBarArea) != TCL_OK) {
      opserr << G3_ERROR_PROMPT << "invalid reinfBarArea: layer straight matTag "
                "numReinfBars reinfBarArea  yStartPt zStartPt yEndPt zEndPt\n";
      return TCL_ERROR;
    }

    if (Tcl_GetDouble(interp, argv[argi++], &yStartPt) != TCL_OK) {
      opserr << G3_ERROR_PROMPT << "invalid yStartPt: layer straight matTag numReinfBars "
                "reinfBarArea  yStartPt zStartPt yEndPt zEndPt\n";
      return TCL_ERROR;
    }

    if (Tcl_GetDouble(interp, argv[argi++], &zStartPt) != TCL_OK) {
      opserr << G3_ERROR_PROMPT << "invalid zStartPt: layer straight matTag numReinfBars "
                "reinfBarArea  yStartPt zStartPt yEndPt zEndPt\n";
      return TCL_ERROR;
    }

    if (Tcl_GetDouble(interp, argv[argi++], &yEndPt) != TCL_OK) {
      opserr << G3_ERROR_PROMPT << "invalid yEndPt: layer straight matTag numReinfBars "
                "reinfBarArea  yStartPt zStartPt yEndPt zEndPt\n";
      return TCL_ERROR;
    }

    if (Tcl_GetDouble(interp, argv[argi++], &zEndPt) != TCL_OK) {
      opserr << G3_ERROR_PROMPT << "invalid zEndPt: layer straight matTag numReinfBars "
                "reinfBarArea  yStartPt zStartPt yEndPt zEndPt\n";
      return TCL_ERROR;
    }

    // create the reinforcing layer

    static Vector startPt(2);
    static Vector endPt(2);

    startPt(0) = yStartPt;
    startPt(1) = zStartPt;
    endPt(0) = yEndPt;
    endPt(1) = zEndPt;

    StraightReinfLayer reinfLayer(matTag, numReinfBars, reinfBarArea, startPt, endPt);

    // add reinfLayer to section
    int error = fiberSectionRepr->addLayer(reinfLayer);
    if (error)
      return TCL_ERROR;

  } else if (strcmp(argv[1], "circ") == 0) {
    if (argc < 8) {
      opserr << G3_ERROR_PROMPT << "invalid number of parameters: layer circ matTag "
                "numReinfBars reinfBarArea yCenter zCenter arcRadius <startAng "
                "endAng>\n";
      return TCL_ERROR;
    }

    int  matTag, numReinfBars;
    double reinfBarArea;
    double yCenter, zCenter, radius, startAng, endAng;

    int argi = 2;

    if (Tcl_GetInt(interp, argv[argi++], &matTag) != TCL_OK) {
      opserr << G3_ERROR_PROMPT << "invalid matTag: layer circ matTag numReinfBars "
                "reinfBarArea yCenter zCenter radius startAng endAng\n";
      return TCL_ERROR;
    }

    if (Tcl_GetInt(interp, argv[argi++], &numReinfBars) != TCL_OK) {
      opserr << G3_ERROR_PROMPT << "invalid numReinfBars: layer circ matTag numReinfBars "
                "reinfBarArea yCenter zCenter radius startAng endAng\n";
      return TCL_ERROR;
    }

    if (Tcl_GetDouble(interp, argv[argi++], &reinfBarArea) != TCL_OK) {
      opserr << G3_ERROR_PROMPT << "invalid reinfBarArea: layer circ matTag numReinfBars "
                "reinfBarArea yCenter zCenter radius startAng endAng\n";
      return TCL_ERROR;
    }

    if (Tcl_GetDouble(interp, argv[argi++], &yCenter) != TCL_OK) {
      opserr << G3_ERROR_PROMPT << "invalid yCenter: layer circ matTag numReinfBars "
                "reinfBarArea yCenter zCenter radius startAng endAng\n";
      return TCL_ERROR;
    }

    if (Tcl_GetDouble(interp, argv[argi++], &zCenter) != TCL_OK) {
      opserr << G3_ERROR_PROMPT << "invalid zCenter: layer circ matTag numReinfBars "
                "reinfBarArea yCenter zCenter radius startAng endAng\n";
      return TCL_ERROR;
    }

    if (Tcl_GetDouble(interp, argv[argi++], &radius) != TCL_OK) {
      opserr << G3_ERROR_PROMPT << "invalid radius: layer circ matTag numReinfBars "
                "reinfBarArea yCenter zCenter radius startAng endAng\n";
      return TCL_ERROR;
    }

    bool anglesSpecified = false;

    if (argc > 9) {
      if (Tcl_GetDouble(interp, argv[argi++], &startAng) != TCL_OK) {
        opserr << G3_ERROR_PROMPT << "invalid startAng: layer circ matTag numReinfBars "
                  "reinfBarArea yCenter zCenter radius startAng endAng\n";
        return TCL_ERROR;
      }

      if (Tcl_GetDouble(interp, argv[argi++], &endAng) != TCL_OK) {
        opserr << G3_ERROR_PROMPT << "invalid endAng: layer circ matTag numReinfBars "
                  "reinfBarArea yCenter zCenter radius startAng endAng\n";
        return TCL_ERROR;
      }

      anglesSpecified = true;
    }

    // create the reinforcing layer

    static Vector center(2);

    center(0) = yCenter;
    center(1) = zCenter;

    CircReinfLayer *reinfLayer = nullptr;
    if (anglesSpecified)
      // Construct arc
      reinfLayer = new CircReinfLayer(matTag, numReinfBars, reinfBarArea,
                                      center, radius, startAng, endAng);
    else
      // Construct circle
      reinfLayer = new CircReinfLayer(matTag, numReinfBars, reinfBarArea,
                                      center, radius);

    // add reinfLayer to section
    int error = fiberSectionRepr->addLayer(*reinfLayer);
    delete reinfLayer;
    if (error)
      return TCL_ERROR;

  } else {
    opserr << G3_ERROR_PROMPT << "reinforcing layer type is not available\n";
    return TCL_ERROR;
  }

  return TCL_OK;
}


int
TclCommand_addUCFiberSection(ClientData clientData, Tcl_Interp *interp,
                             int argc, TCL_Char ** const argv)
{
  assert(clientData != nullptr);
  BasicModelBuilder* builder = static_cast<BasicModelBuilder*>(clientData);
  G3_Runtime *rt = G3_getRuntime(interp);
  int secTag;

  if (argc < 4)
    return TCL_ERROR;

  if (Tcl_GetInt(interp, argv[2], &secTag) != TCL_OK) {
    opserr << "could not read section tag\n";
    return TCL_ERROR;
  }

  builder->setCurrentSectionBuilder(secTag);

  // first create an empty FiberSection

  SectionForceDeformation *section = nullptr;
  FiberSection2d *section2d = nullptr;
  FiberSection3d *section3d = nullptr;
  bool computeCentroid = false;

  int NDM = builder->getNDM();
  if (NDM == 2) {
    section2d = new FiberSection2d(secTag, 30, computeCentroid);
    section = section2d;
    // SectionForceDeformation *section = new FiberSection(secTag, 0, 0);
  } else if (NDM == 3) {
    UniaxialMaterial *theGJ = new ElasticMaterial(0, 1e10);
    section3d =
        new FiberSection3d(secTag, 30, *theGJ, computeCentroid);
    section = section3d;
    delete theGJ;
  }

  if (section == nullptr) {
    return TCL_ERROR;
  }

  //
  // now parse the output file containing the fiber data,
  // create fibers and add them to the section
  //

  // open the file
  TCL_Char *fileName = argv[3];
  std::ifstream theFile;
  theFile.open(fileName, std::ios::in);
  if (!theFile) {
    opserr << "section UCFiber - could not open file named " << fileName;
    return TCL_ERROR;
  } else {
    int foundStart = 0;
    static char garbage[100];

    // parse through until find start of fiber data
    while (foundStart == 0 && theFile >> garbage)
      if (strcmp(garbage, "#FIBERS") == 0)
        foundStart = 1;

    if (foundStart == 0) {
      theFile.close();
      return TCL_ERROR;
    }

    // parse the fiber data until eof, creating a fiber and adding to section as
    // go
    double ycoord, zcoord, area, prestrain;
    int matTag;

    while (theFile >> ycoord >> zcoord >> area >> prestrain >> garbage >>
           matTag) {

      UniaxialMaterial *theMaterial = G3_getUniaxialMaterialInstance(rt,matTag);
      if (theMaterial == 0) {
        opserr << "section UCFiber - no material exists with tag << " << matTag
               << endln;
        return TCL_ERROR;
      }

      if (NDM == 2)
        section2d->addFiber(*theMaterial, area, zcoord);
      else
        section3d->addFiber(*theMaterial, area, ycoord, zcoord);

    }

    theFile.close();
  }

  // finally add the section to our modelbuilder
  if (builder->addTaggedObject<SectionForceDeformation>(*section) < 0) {
    return TCL_ERROR;
  }

  return TCL_OK;
}

#if 0
#include <UniaxialFiber2d.h>
#include <UniaxialFiber3d.h>
int
buildSectionInt(ClientData clientData, Tcl_Interp *interp, TclBasicBuilder *theTclBasicBuilder,
                int secTag, UniaxialMaterial &theTorsion, int NStrip1,
                double t1, int NStrip2, double t2, int NStrip3, double t3)
{
  assert(clientData != nullptr);
  BasicModelBuilder* builder = static_cast<BasicModelBuilder*>(clientData);
  SectionRepres *sectionRepres = theTclBasicBuilder->getSectionRepres(secTag);

  if (sectionRepres == nullptr) {
    opserr << G3_ERROR_PROMPT << "cannot retrieve section\n";
    return TCL_ERROR;
  }

  if (sectionRepres->getType() == SEC_TAG_FiberSection) {
    // build the section

    FiberSectionRepr *fiberSectionRepr = (FiberSectionRepr *)sectionRepres;

    int i, j, k;
    int numFibers;
    int numHFibers;

    int numPatches;
    Patch **patch;

    int numReinfLayers;
    ReinfLayer **reinfLayer;

    numPatches = fiberSectionRepr->getNumPatches();
    patch = fiberSectionRepr->getPatches();
    numReinfLayers = fiberSectionRepr->getNumReinfLayers();
    reinfLayer = fiberSectionRepr->getReinfLayers();

    int numSectionRepresFibers = fiberSectionRepr->getNumFibers();
    Fiber **sectionRepresFibers = fiberSectionRepr->getFibers();

    int numSectionRepresHFibers = fiberSectionRepr->getNumHFibers();
    Fiber **sectionRepresHFibers = fiberSectionRepr->getHFibers();

    numFibers = numSectionRepresFibers;
    for (int i = 0; i < numPatches; ++i)
      numFibers += patch[i]->getNumCells();

    for (int i = 0; i < numReinfLayers; ++i)
      numFibers += reinfLayer[i]->getNumReinfBars();

    numHFibers = numSectionRepresHFibers;

    static Vector fiberPosition(2);
    int matTag;

    ID fibersMaterial(numFibers - numSectionRepresFibers);
    Matrix fibersPosition(2, numFibers - numSectionRepresFibers);
    Vector fibersArea(numFibers - numSectionRepresFibers);

    int numCells;
    Cell **cell;

    k = 0;
    for (int i = 0; i < numPatches; ++i) {
      numCells = patch[i]->getNumCells();
      matTag = patch[i]->getMaterialID();

      cell = patch[i]->getCells();

      for (int j = 0; j < numCells; j++) {
        fibersMaterial(k) = matTag;
        fibersArea(k) = cell[j]->getArea();
        fiberPosition = cell[j]->getCentroidPosition();
        fibersPosition(0, k) = fiberPosition(0);
        fibersPosition(1, k) = fiberPosition(1);
        k++;
      }

      for (int j = 0; j < numCells; j++)
        delete cell[j];

      delete[] cell;
    }

    ReinfBar *reinfBar;
    int numReinfBars;

    for (int i = 0; i < numReinfLayers; ++i) {
      numReinfBars = reinfLayer[i]->getNumReinfBars();
      reinfBar = reinfLayer[i]->getReinfBars();
      matTag = reinfLayer[i]->getMaterialID();

      for (int j = 0; j < numReinfBars; j++) {
        fibersMaterial(k) = matTag;
        fibersArea(k) = reinfBar[j].getArea();
        fiberPosition = reinfBar[j].getPosition();

        fibersPosition(0, k) = fiberPosition(0);
        fibersPosition(1, k) = fiberPosition(1);

        k++;
      }
      delete[] reinfBar;
    }

    UniaxialMaterial *material = nullptr;

    Fiber **fiber = new Fiber *[numFibers];

    // copy the section repres fibers
    for (i = 0; i < numSectionRepresFibers; ++i)
      fiber[i] = sectionRepresFibers[i];

    Fiber **Hfiber = new Fiber *[numHFibers];

    // copy the section repres fibers
    for (int i = 0; i < numSectionRepresHFibers; ++i)
      Hfiber[i] = sectionRepresHFibers[i];

    // creates 2d section
    int NDM = builder->getNDM();
    if (NDM == 2) {
      k = 0;
      for (int i = numSectionRepresFibers; i < numFibers; ++i) {
        material = builder->getUniaxialMaterial(fibersMaterial(k));
        if (material == nullptr) {
          opserr << G3_ERROR_PROMPT << "invalid material ID for patch\n";
          return TCL_ERROR;
        }

        fiber[i] = new UniaxialFiber2d(k, *material, fibersArea(k),
                                       fibersPosition(0, k));

        k++;
      }

      SectionForceDeformation *section =
          new FiberSection2dInt(secTag, numFibers, fiber, numHFibers, Hfiber,
                                NStrip1, t1, NStrip2, t2, NStrip3, t3);

      // Delete fibers
      for (int i = 0; i < numFibers; ++i)
        delete fiber[i];

      for (int i = 0; i < numHFibers; ++i)
        delete Hfiber[i];

      if (section == nullptr) {
        opserr << G3_ERROR_PROMPT << "cannot construct section\n";
        return TCL_ERROR;
      }

      if (theTclBasicBuilder->addSection (*section) < 0) {
        opserr << G3_ERROR_PROMPT << "- cannot add section\n";
        return TCL_ERROR;
      }

    } else if (NDM == 3) {

      static Vector fiberPosition(2);
      k = 0;
      for (i = numSectionRepresFibers; i < numFibers; ++i) {
        material = builder->getUniaxialMaterial(fibersMaterial(k));
        if (material == nullptr) {
          opserr << G3_ERROR_PROMPT << "invalid material ID for patch\n";
          return TCL_ERROR;
        }

        fiberPosition(0) = fibersPosition(0, k);
        fiberPosition(1) = fibersPosition(1, k);

        fiber[i] =
            new UniaxialFiber3d(k, *material, fibersArea(k), fiberPosition);
        if (fibersArea(k) < 0)
          opserr << "ERROR: " << fiberPosition(0) << " " << fiberPosition(1)
                 << endln;
        k++;
      }

      SectionForceDeformation *section = 0;
      section = new FiberSection3d(secTag, numFibers, fiber, theTorsion,
                                   options.computeCentroid);

      // Delete fibers
      for (i = 0; i < numFibers; ++i)
        delete fiber[i];

      if (section == 0) {
        opserr << G3_ERROR_PROMPT << "- cannot construct section\n";
        return TCL_ERROR;
      }

      if (theTclBasicBuilder->addSection (*section) < 0) {
        opserr << G3_ERROR_PROMPT << "- cannot add section\n";
        return TCL_ERROR;
      }

    } else {
      opserr << G3_ERROR_PROMPT << "NDM = " << NDM
             << " is incompatible with available frame elements\n";
      return TCL_ERROR;
    }

    // Delete fiber array
    delete[] fiber;
    //   delete [] Hfiber;

  } else {
    opserr << G3_ERROR_PROMPT << "section invalid: can only build fiber sections\n";
    return TCL_ERROR;
  }

  return TCL_OK;
}
#endif


#if 0
SectionForceDeformation*
G3Parse_newTubeSection(G3_Runtime* rt, int argc, G3_Char ** const argv)
{
  SectionForceDeformation *theSection = nullptr;
  if (strcmp(argv[1], "Tube") == 0) {
    if (argc < 8) {
      opserr << G3_ERROR_PROMPT << "insufficient arguments\n";
      opserr << "Want: section Tube tag? matTag? D? t? nfw? nfr?" << endln;
      return nullptr;
    }

    int tag, matTag;
    double D, t;
    int nfw, nfr;

    if (Tcl_GetInt(interp, argv[2], &tag) != TCL_OK) {
      opserr << G3_ERROR_PROMPT << "invalid section Tube tag" << endln;
      return nullptr;
    }

    if (Tcl_GetInt(interp, argv[3], &matTag) != TCL_OK) {
      opserr << G3_ERROR_PROMPT << "invalid section Tube matTag" << endln;
      return nullptr;
    }

    if (Tcl_GetDouble(interp, argv[4], &D) != TCL_OK) {
      opserr << G3_ERROR_PROMPT << "invalid D" << endln;
      opserr << "Tube section: " << tag << endln;
      return nullptr;
    }

    if (Tcl_GetDouble(interp, argv[5], &t) != TCL_OK) {
      opserr << G3_ERROR_PROMPT << "invalid t" << endln;
      opserr << "Tube section: " << tag << endln;
      return nullptr;
    }

    if (Tcl_GetInt(interp, argv[6], &nfw) != TCL_OK) {
      opserr << G3_ERROR_PROMPT << "invalid nfw" << endln;
      opserr << "Tube section: " << tag << endln;
      return nullptr;
    }

    if (Tcl_GetInt(interp, argv[7], &nfr) != TCL_OK) {
      opserr << G3_ERROR_PROMPT << "invalid nfr" << endln;
      opserr << "Tube  section: " << tag << endln;
      return nullptr;
    }

    TubeSectionIntegration tubesect(D, t, nfw, nfr);

    int numFibers = tubesect.getNumFibers();

    if (argc > 8) {

      double shape = 1.0;
      if (argc > 9) {
        if (Tcl_GetDouble(interp, argv[9], &shape) != TCL_OK) {
          opserr << G3_ERROR_PROMPT << "invalid shape" << endln;
          opserr << "WFSection2d section: " << tag << endln;
          return nullptr;
        }
      }

      NDMaterial *theSteel = builder->getTypedObject<NDMaterial>(matTag);
      if (theSteel == 0)
        return nullptr;


      NDMaterial **theMats = new NDMaterial *[numFibers];

      tubesect.arrangeFibers(theMats, theSteel);

      // Parsing was successful, allocate the section
      theSection = 0;
      if (strcmp(argv[8], "-nd") == 0)
        theSection =
            new NDFiberSection3d(tag, numFibers, theMats, tubesect, shape);
      if (strcmp(argv[8], "-ndWarping") == 0)
        theSection = new NDFiberSectionWarping2d(tag, numFibers, theMats,
                                                 tubesect, shape);

      delete[] theMats;
    } else {
      UniaxialMaterial *theSteel = builder->getTypedObject<UniaxialMaterial>(matTag);
      if (theSteel == 0)
        return nullptr;

      UniaxialMaterial **theMats = new UniaxialMaterial *[numFibers];

      tubesect.arrangeFibers(theMats, theSteel);

      // Parsing was successful, allocate the section
      theSection = new FiberSection2d(tag, numFibers, theMats, tubesect);

      delete[] theMats;
    }
  }
}
#endif

