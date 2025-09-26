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
// Description: This file contains the function invoked when the user invokes
// the section command in the interpreter.
//===----------------------------------------------------------------------===//
//
// Membrane      nxx nyy nxy
// MembranePlate nxx nyy nxy mxx myy mxy vxz vyz
//
// Written: rms, mhs, cmp
// Created: 07/99
//
#include <set>
#include <assert.h>
#include <string.h>
#include <fstream>
#include <iostream>

#include <tcl.h>
#include <Parsing.h>
#include <Logging.h>
#include <ArgumentTracker.h>
#include <BasicModelBuilder.h>
#include <Parameter.h>

#include <packages.h>
#include <runtimeAPI.h>

using namespace OpenSees;

extern "C" int OPS_ResetInputNoBuilder(ClientData clientData,
                                       Tcl_Interp *interp, int cArg, int mArg,
                                       TCL_Char ** const argv, Domain *domain);

#include <PlaneStrainMaterial.h>
#include <PlaneStressMaterial.h>
#include <FrameSection.h>
#include <ElasticMaterial.h>
#include <ParallelSection.h>
#include <FiberSection2d.h>
#include <NDFiberSection2d.h>
#include <NDFiberSection3d.h>
#include <NDFiberSectionWarping2d.h>
#include <FiberSection2dInt.h>
#include <FiberSection3d.h>
#include <FrameFiberSection3d.h>
#include <FrameSolidSection3d.h>
#include <FiberSectionAsym3d.h>

// SectionBuilder
#include <QuadPatch.h>
#include <CircPatch.h>
#include <StraightReinfLayer.h>
#include <CircReinfLayer.h>
#include <FiberSectionBuilder.h>

//
#include <Bidirectional.h>
#include <Elliptical2.h>
#include <Isolator2spring.h>

#include <FiberSection2dThermal.h>
#include <FiberSection3dThermal.h> //Added by L.Jiang [SIF] 2017

//#include <McftSection2dfiber.h>

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
Tcl_CmdProc TclCommand_addPlaneSection;
Tcl_CmdProc TclCommand_addElasticShellSection;
Tcl_CmdProc TclCommand_addUniaxialSection;
Tcl_CmdProc TclCommand_ShellSection;

// extern OPS_Routine OPS_WFSection2d;
// extern OPS_Routine OPS_RCCircularSection;
// extern OPS_Routine OPS_RCSection2d;
// extern OPS_Routine OPS_RCTBeamSection2d;
// extern OPS_Routine OPS_RCTunnelSection;
// extern OPS_Routine OPS_TubeSection;

SectionForceDeformation *
TclBasicBuilderYS_SectionCommand(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char ** const argv);


#include <set>
#include <ArgumentTracker.h>
#include <UniaxialMaterial.h>
int
TclCommand_addTrussSection(ClientData clientData, Tcl_Interp *interp,
                              int argc, TCL_Char ** const argv)
{

  enum class Positions {
    Tag, MaterialTag, Area,
    EndRequired, End
  };
  
  ArgumentTracker<Positions> tracker;
  std::set<int> positional;
  int tag;
  int matTag;
  double area;

  for (int i=2; i<argc; i++) {
    if (strcmp(argv[i], "-material") == 0) {
      if (++i >= argc) {
        opserr << "Missing value for option " << argv[i-1] << "\n";
        return TCL_ERROR;
      }
      if (Tcl_GetInt(interp, argv[i], &matTag) != TCL_OK) {
        opserr << "Invalid value for option " << argv[i-1] << "\n";
        return TCL_ERROR;
      }
      tracker.consume(Positions::MaterialTag);
    }
    else if (strcmp(argv[i], "-area") == 0) {
      if (++i >= argc) {
        opserr << "Missing value for option " << argv[i-1] << "\n";
        return TCL_ERROR;
      }
      if (Tcl_GetDouble(interp, argv[i], &area) != TCL_OK) {
        opserr << "Invalid value for option " << argv[i-1] << "\n";
        return TCL_ERROR;
      }
      tracker.consume(Positions::Area);
    }
    else
      positional.insert(i);
  }

  //
  // Positional arguments
  //
  for (int i: positional) {
    if (tracker.current() == Positions::EndRequired)
      tracker.increment();

    switch (tracker.current()) {
      case Positions::Tag :
        if (Tcl_GetInt(interp, argv[i], &tag) != TCL_OK) {
          opserr << "Invalid tag " << argv[i] << "\n";
          return TCL_ERROR;
        }
        tracker.consume(Positions::Tag);
        break;

      case Positions::MaterialTag:
        if (Tcl_GetInt(interp, argv[i], &matTag) != TCL_OK) {
          opserr << "Invalid value for material tag " << argv[i] << "\n";
          return TCL_ERROR;
        }
        tracker.consume(Positions::MaterialTag);
        break;

      case Positions::Area:
        if (Tcl_GetDouble(interp, argv[i], &area) != TCL_OK) {
          opserr << "Invalid value for area " << argv[i] << "\n";
          return TCL_ERROR;
        }
        tracker.consume(Positions::Area);
        break;

      case Positions::EndRequired:
      case Positions::End:
        break;
    }
  }

  //
  if (tracker.current() < Positions::EndRequired) {
    opserr << "Missing required arguments: ";
    if (tracker.contains(Positions::Tag))
      opserr << "tag ";
    if (tracker.contains(Positions::MaterialTag))
      opserr << "material ";
    if (tracker.contains(Positions::Area))
      opserr << "area ";
    opserr << "\n";
    return TCL_ERROR;
  }

  BasicModelBuilder *builder = static_cast<BasicModelBuilder*>(clientData);

  UniaxialMaterial *material = builder->getTypedObject<UniaxialMaterial>(matTag);
  if (material == nullptr) {
    return TCL_ERROR;
  }
  auto fiber_section = new FrameFiberSection3d(tag, 1, nullptr, true, 0.0, 0);
  fiber_section->addFiber(*material, area, 0.0, 0.0);
  return builder->addTaggedObject<FrameSection>(*fiber_section);
}


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
    opserr << OpenSees::PromptValueError << "insufficient number of arguments\n";
    return TCL_ERROR;
  }

  OPS_ResetInputNoBuilder(clientData, interp, 2, argc, argv, theDomain);

  // Pointer to a section that will be added to the model builder
  SectionForceDeformation *theSection = nullptr;

  // Check argv[1] to dispatch section type

  if (strcmp(argv[1], "Fiber") == 0 || 
      strcmp(argv[1], "fiberSec") == 0 ||
      strcmp(argv[1], "FiberSec") == 0 ||
      strcmp(argv[1], "FiberFrame") == 0 ||
      strcmp(argv[1], "FrameFiber") == 0 ||
      strcmp(argv[1], "AxialFiber") == 0 ||
      strcmp(argv[1], "FiberSection") == 0 ||
      //
      strcmp(argv[1], "ShearFiber") == 0 ||
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


  else if (strcmp(argv[1], "Truss") == 0 ||
           strcmp(argv[1], "TrussSection") == 0 ||
           strcmp(argv[1], "TrussSection2d") == 0 ||
           strcmp(argv[1], "TrussSection3d") == 0) {
    return TclCommand_addTrussSection(clientData, interp, argc, argv);
  }



  else if (strcmp(argv[1], "FiberInt") == 0) {
    // TODO
    // return TclCommand_addFiberIntSection(clientData, interp, argc, argv);
    opserr << "FiberInt is currently broken\n";
    return TCL_ERROR;
  }

  else if ((strcmp(argv[1], "PlaneStrain") == 0) ||
           (strcmp(argv[1], "PlaneStress") == 0))
    return TclCommand_addPlaneSection(clientData, interp, argc, argv);

  else if (strcmp(argv[1], "UCFiber") == 0)
    return TclCommand_addUCFiberSection(clientData, interp, argc, argv);

  else if ((strcmp(argv[1], "Elastic") == 0) ||
    (strcmp(argv[1], "ElasticShear") == 0) ||
    (strcmp(argv[1], "ElasticWarpingShear") == 0) ||
    (strcmp(argv[1], "FrameElastic") == 0) ||
    (strcmp(argv[1], "ElasticFrame") == 0)) {
    return TclCommand_newElasticSection(clientData, interp, argc, argv);
  }

  else if (strcmp(argv[1], "Generic1D") == 0 ||
        strcmp(argv[1], "Generic1d") == 0 ||
        strcmp(argv[1], "Uniaxial") == 0) {
    return TclCommand_addUniaxialSection(clientData, interp, argc, argv);
  }

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
  }

  else if (strcmp(argv[1], "AddDeformation") == 0 ||
           strcmp(argv[1], "Aggregator") == 0  ||
           strcmp(argv[1], "Aggregate") == 0) 
    return TclCommand_addSectionAggregator(clientData, interp, argc, argv);


  //
  // Shell
  //
  else if ((strcmp(argv[1], "ElasticMembranePlateSection") == 0) ||
           (strcmp(argv[1], "ElasticShell") == 0) ||
           (strcmp(argv[1], "ElasticPlateSection") == 0) ||
           (strcmp(argv[1], "ElasticMembraneSection") == 0)) {
    return TclCommand_addElasticShellSection(clientData, interp, argc, argv);
  }

  else if ((strcmp(argv[1], "LayeredShell") == 0) || 
           (strcmp(argv[1], "LayeredShellThermal") == 0) ||
           (strcmp(argv[1], "PlateFiber") == 0) || 
           (strcmp(argv[1], "PlateFiberThermal") == 0)) {
    return TclCommand_ShellSection(clientData, interp, argc, argv);
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

  //
  //
  //
  else if (strcmp(argv[1], "Iso2spring") == 0) {
    if (argc < 10) {
      opserr << OpenSees::PromptValueError << "insufficient arguments\n";
      opserr
          << "Want: section Iso2spring tag? tol? k1? Fy? k2? kv? hb? Pe? <Po?>"
          << "\n";
      return TCL_ERROR;
    }

    int tag;
    double tol, k1, Fy, kb, kvo, hb, Pe;
    double Po = 0.0;

    if (Tcl_GetInt(interp, argv[2], &tag) != TCL_OK) {
      opserr << OpenSees::PromptValueError << "invalid Iso2spring tag" << "\n";
      return TCL_ERROR;
    }

    if (Tcl_GetDouble(interp, argv[3], &tol) != TCL_OK) {
      opserr << OpenSees::PromptValueError << "invalid tol\n";
      return TCL_ERROR;
    }

    if (Tcl_GetDouble(interp, argv[4], &k1) != TCL_OK) {
      opserr << OpenSees::PromptValueError << "invalid k1\n";
      return TCL_ERROR;
    }

    if (Tcl_GetDouble(interp, argv[5], &Fy) != TCL_OK) {
      opserr << OpenSees::PromptValueError << "invalid Fy\n";
      return TCL_ERROR;
    }

    if (Tcl_GetDouble(interp, argv[6], &kb) != TCL_OK) {
      opserr << OpenSees::PromptValueError << "invalid k2\n";
      return TCL_ERROR;
    }

    if (Tcl_GetDouble(interp, argv[7], &kvo) != TCL_OK) {
      opserr << OpenSees::PromptValueError << "invalid kv\n";
      return TCL_ERROR;
    }
    if (Tcl_GetDouble(interp, argv[8], &hb) != TCL_OK) {
      opserr << OpenSees::PromptValueError << "invalid hb\n";
      return TCL_ERROR;
    }

    if (Tcl_GetDouble(interp, argv[9], &Pe) != TCL_OK) {
      opserr << OpenSees::PromptValueError << "invalid Pe\n";
      return TCL_ERROR;
    }
    if (argc > 10) {
      if (Tcl_GetDouble(interp, argv[10], &Po) != TCL_OK) {
        opserr << OpenSees::PromptValueError << "invalid Po\n";
        opserr << "section Iso2spring: " << tag << "\n";
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
    opserr << OpenSees::PromptValueError << "could not create section " << argv[1] << "\n";
    return TCL_ERROR;
  }

  // Now add the material to the modelBuilder
  if (builder->addTaggedObject<SectionForceDeformation>(*theSection) < 0) {
    opserr << *theSection << "\n";
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
        opserr << OpenSees::PromptValueError << "failed to parse section tag \"" << argv[i+1] << "\"\n";
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
                    int secTag, 
                    UniaxialMaterial *theTorsion, 
                    double Ys, double Zs, 
                    double alpha, 
                    const FiberSectionConfig& options)
{
  assert(clientData != nullptr);
  BasicModelBuilder *builder = static_cast<BasicModelBuilder*>(clientData);

  // Dimension of the structure
  int ndm = builder->getNDM();

  SectionBuilder  *sbuilder = nullptr;
  FrameSection    *section  = nullptr;
  // Create 2d section
  if (ndm == 2) {
    if (options.isND) {
      if (options.isNew) {
        auto sec = new FrameSolidSection3d(secTag, 30);
        sbuilder = new FiberSectionBuilder<2, NDMaterial, FrameSolidSection3d>(*builder, *sec);
        section = sec;
      }
      else if (options.isWarping) {
        auto sec = new NDFiberSectionWarping2d(secTag, 30, alpha);
        sbuilder = new FiberSectionBuilder<2, NDMaterial, NDFiberSectionWarping2d>(*builder, *sec);
        section = sec;
      }
      else {
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
  }

  else if (ndm == 3) {

    if (options.isND) {
      if (options.isNew) {
        auto sec = new FrameSolidSection3d(secTag, 30);
        sbuilder = new FiberSectionBuilder<3, NDMaterial, FrameSolidSection3d>(*builder, *sec);
        section = sec;
      }
      else {
        auto sec = new NDFiberSection3d(secTag,
                                        options.computeCentroid);
        sbuilder = new FiberSectionBuilder<3, NDMaterial, NDFiberSection3d>(*builder, *sec);
        section = sec;
      }
    } else {


      if (options.isThermal) {
        if (theTorsion == nullptr) {
          opserr << OpenSees::PromptValueError 
                 << "FiberThermal section requires torsion\n";
          return TCL_ERROR;
        }
        auto sec = new FiberSection3dThermal(secTag, 30, *theTorsion,
                                             options.computeCentroid);
        sbuilder = new FiberSectionBuilder<3, UniaxialMaterial, FiberSection3dThermal>(*builder, *sec);
        section = sec;

      }
      else if (options.isAsym) {
        auto sec = new FiberSectionAsym3d(secTag, 30, theTorsion, Ys, Zs);
        sbuilder = new FiberSectionBuilder<3, UniaxialMaterial, FiberSectionAsym3d>(*builder, *sec);
        section = sec;

      }
      else {
        if (options.isNew) {
          auto sec = new FrameFiberSection3d(secTag, 30,  theTorsion, options.computeCentroid, 
                                             options.density, options.use_density);
          sbuilder = new FiberSectionBuilder<3, UniaxialMaterial, FrameFiberSection3d>(*builder, *sec);
          section = sec;
        } else {
          if (theTorsion == nullptr) {
            opserr << OpenSees::PromptValueError 
                  << "FiberThermal section requires torsion\n";
            return TCL_ERROR;
          }
          auto sec = new FiberSection3d(secTag, 30, *theTorsion, options.computeCentroid);
          sbuilder = new FiberSectionBuilder<3, UniaxialMaterial, FiberSection3d>(*builder, *sec);
          section = sec;
        }
      }
    }

  } else {
    opserr << OpenSees::PromptValueError << "Model dimension (ndm = " << ndm
           << ") is incompatible with available frame elements\n";
    return TCL_ERROR;
  }


  // In 2D truss elements still look for FrameSections
  if (builder->addTaggedObject<FrameSection>(*section) < 0) {
    return TCL_ERROR;
  }
  // if (ndm == 2) {
  //   if (builder->addTaggedObject<SectionForceDeformation>(*section->getCopy()) < 0) {
  //     return TCL_ERROR;
  //   }
  // }

  if (builder->addTypedObject<SectionBuilder>(secTag, sbuilder) < 0) {
    opserr << OpenSees::PromptValueError << "Faled to add section\n";
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
    opserr << OpenSees::PromptValueError 
           << "failed to parse section tag \"" << argv[2] << "\""
           << OpenSees::SignalMessageEnd;
    return TCL_ERROR;
  }

  builder->setCurrentSectionBuilder(secTag);

  FiberSectionConfig options;
  if (strcmp(argv[1], "NDFiber") == 0 ||
      strcmp(argv[1], "ShearFiber") == 0)
    options.isND = true;

  if (strcmp(argv[1], "NDFiberWarping") == 0) {
    options.isND = true;
    options.isWarping = true;
  }
  else if (strcmp(argv[1], "FrameFiber") == 0 ||
           strcmp(argv[1], "FiberFrame") == 0 ||
           strcmp(argv[1], "AxialFiber") == 0 ||
           strcmp(argv[1], "ShearFiber") == 0 
    )
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
        opserr << OpenSees::PromptValueError << "not enough -mass args need -mass mass?\n";
        return TCL_ERROR;
      }
      if (Tcl_GetDouble(interp, argv[iarg + 1], &options.density) != TCL_OK) {
        opserr << OpenSees::PromptValueError << "invalid density";
        return TCL_ERROR;
      }
      options.use_density = true;

      iarg  += 2;
    }

    else if (strcmp(argv[iarg], "-GJ") == 0 && iarg + 1 < argc) {
      if (Tcl_GetDouble(interp, argv[iarg + 1], &GJ) != TCL_OK) {
        opserr << OpenSees::PromptValueError << "invalid GJ";
        return TCL_ERROR;
      }
      deleteTorsion = true;
      torsion = new ElasticMaterial(0, GJ);

      iarg  += 2;
    }

    else if (strcmp(argv[iarg], "-torsion") == 0 && iarg + 1 < argc) {
      int torsionTag = 0;
      if (Tcl_GetInt(interp, argv[iarg + 1], &torsionTag) != TCL_OK) {
        opserr << OpenSees::PromptValueError << "invalid torsionTag";
        return TCL_ERROR;
      }

      torsion = builder->getTypedObject<UniaxialMaterial>(torsionTag);
      if (torsion == nullptr) {
        opserr << OpenSees::PromptValueError << "uniaxial material does not exist\n";
        opserr << "uniaxial material: " << torsionTag;
        opserr << "\nFiberSection3d: " << secTag << "\n";
        return TCL_ERROR;
      }

      iarg += 2;
    }

    else if (strstr(argv[1], "Asym") != nullptr && !shearParsed) {
      if (iarg + 1 >= argc) {
        opserr << OpenSees::PromptValueError << "Asym sections require shear center before fiber block.\n";
        return TCL_ERROR;
      }
      if (Tcl_GetDouble(interp, argv[iarg], &Ys) != TCL_OK) {
        opserr << OpenSees::PromptValueError << "invalid Ys";
        return TCL_ERROR;
      }
      if (Tcl_GetDouble(interp, argv[iarg+1], &Zs) != TCL_OK) {
        opserr << OpenSees::PromptValueError << "invalid Zs";
        return TCL_ERROR;
      }
      shearParsed = true;
      iarg  += 2;
    }

    else {
      // braces; skip and handle later
      break;
    }
  }

  if (torsion == nullptr && ndm == 3 && !options.isNew) {
    opserr << OpenSees::PromptValueError
           << "missing required torsion for 3D fiber section, use -GJ or "
              "-torsion\n";
    return TCL_ERROR;
  }

  // initialize  the fiber section (for building)            // TODO, pass alpha
  if (initSectionCommands(clientData, interp, secTag, torsion, Ys, Zs, 1.0, options) != TCL_OK) {
    opserr << OpenSees::PromptValueError << "error constructing the section\n";
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
    opserr << OpenSees::PromptValueError << "bad command - want: \nsection fiberInt secTag -GJ <GJ> { "
              "\n\tpatch <patch arguments> \n\tlayer <layer arguments> \n}\n";
    return TCL_ERROR;
  }

  builder->setCurrentSectionBuilder(secTag);

  //
  int brace = 3;
  double GJ = 1.0;
  bool deleteTorsion = false;
  UniaxialMaterial *torsion = nullptr;
  if (strcmp(argv[3], "-GJ") == 0) {
    if (Tcl_GetDouble(interp, argv[4], &GJ) != TCL_OK) {
      opserr << OpenSees::PromptValueError << "invalid GJ\n";
      return TCL_ERROR;
    }
    torsion = new ElasticMaterial(0, GJ); // Is this gonna be a memory leak? MHS

    brace = 5;
  }
  int torsionTag = 0;
  if (strcmp(argv[3], "-torsion") == 0) {
    if (Tcl_GetInt(interp, argv[4], &torsionTag) != TCL_OK) {
      opserr << OpenSees::PromptValueError << "invalid torsionTag\n";
      return TCL_ERROR;
    }

    torsion = builder->getTypedObject<UniaxialMaterial>(torsionTag);
    if (torsion == 0)
      return TCL_ERROR;

    brace = 5;
  }

  int NStrip1, NStrip2, NStrip3;
  double t1, t2, t3;

  if (strcmp(argv[3], "-NStrip") == 0) {

    if (Tcl_GetInt(interp, argv[4], &NStrip1) != TCL_OK) {
      opserr << OpenSees::PromptValueError << "invalid NStrip1\n";
      return TCL_ERROR;
    }

    if (Tcl_GetDouble(interp, argv[5], &t1) != TCL_OK) {
      opserr << OpenSees::PromptValueError << "invalid t1";
      return TCL_ERROR;
    }

    if (Tcl_GetInt(interp, argv[6], &NStrip2) != TCL_OK) {
      opserr << OpenSees::PromptValueError << "invalid NStrip2";
      return TCL_ERROR;
    }

    if (Tcl_GetDouble(interp, argv[7], &t2) != TCL_OK) {
      opserr << OpenSees::PromptValueError << "invalid t2";
      return TCL_ERROR;
    }

    if (Tcl_GetInt(interp, argv[8], &NStrip3) != TCL_OK) {
      opserr << OpenSees::PromptValueError << "invalid NStrip3";
      return TCL_ERROR;
    }

    if (Tcl_GetDouble(interp, argv[9], &t3) != TCL_OK) {
      opserr << OpenSees::PromptValueError << "invalid t3";
      return TCL_ERROR;
    }

    brace = 10; // may be 5
  }


  // parse the information inside the braces (patches and reinforcing layers)
  if (Tcl_Eval(interp, argv[brace]) != TCL_OK) {
    opserr << OpenSees::PromptValueError << "- error reading information in { } \n";
    return TCL_ERROR;
  }

  if (NDM == 3 && torsion == nullptr) {
    opserr << OpenSees::PromptValueError << "- no torsion specified for 3D fiber section, use -GJ or "
              "-torsion\n";
    opserr << "\nFiberSectionInt3d: " << secTag << "\n";
    return TCL_ERROR;
  }

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
    opserr << OpenSees::PromptValueError << "cannot retrieve section\n";
    return TCL_ERROR;
  }


  // make sure at least one other argument to contain patch type
  if (argc < 2) {
    opserr << OpenSees::PromptValueError << "need to specify a patch type \n";
    return TCL_ERROR;
  }

  // check argv[1] for type of patch  and create the object
  if (strcmp(argv[1], "quad") == 0 || strcmp(argv[1], "quadr") == 0) {
    int numSubdivIJ, numSubdivJK, matTag;
    double vertexCoordY, vertexCoordZ;
    MatrixND<4,2> vertexCoords{};

    if (argc < 13) {
      opserr << OpenSees::PromptValueError << "invalid number of parameters: patch quad matTag "
                "numSubdivIJ numSubdivJK yVertI zVertI yVertJ zVertJ yVertK "
                "zVertK yVertL zVertL\n";
      return TCL_ERROR;
    }

    int argi = 2;

    if (Tcl_GetInt(interp, argv[argi++], &matTag) != TCL_OK) {
      opserr << OpenSees::PromptValueError << "invalid matTag: patch quad matTag numSubdivIJ "
                "numSubdivJK yVertI zVertI yVertJ zVertJ yVertK zVertK yVertL "
                "zVertL\n";
      return TCL_ERROR;
    }

    if (Tcl_GetInt(interp, argv[argi++], &numSubdivIJ) != TCL_OK) {
      opserr << OpenSees::PromptValueError << "invalid numSubdivIJ: patch quad matTag numSubdivIJ "
                "numSubdivJK yVertI zVertI yVertJ zVertJ yVertK zVertK yVertL "
                "zVertL\n";
      return TCL_ERROR;
    }

    if (Tcl_GetInt(interp, argv[argi++], &numSubdivJK) != TCL_OK) {
      opserr << OpenSees::PromptValueError << "invalid numSubdivJK: patch quad matTag numSubdivIJ "
                "numSubdivJK yVertI zVertI yVertJ zVertJ yVertK zVertK yVertL "
                "zVertL\n";
      return TCL_ERROR;
    }

    for (int j = 0; j < 4; j++) {
      if (Tcl_GetDouble(interp, argv[argi++], &vertexCoordY) != TCL_OK) {
        opserr << OpenSees::PromptValueError << "invalid Coordinate y: ...yVertI zVertI yVertJ "
                  "zVertJ yVertK zVertK yVertL zVertL\n";
        return TCL_ERROR;
      }

      if (Tcl_GetDouble(interp, argv[argi++], &vertexCoordZ) != TCL_OK) {
        opserr << OpenSees::PromptValueError << "invalid Coordinate z: ...yVertI zVertI yVertJ "
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
    MatrixND<4,2> vertexCoords;

    if (argc < 9) {
      opserr << OpenSees::PromptValueError << "invalid number of parameters: patch quad matTag "
                "numSubdivIJ numSubdivJK yVertI zVertI yVertK zVertK\n";
      return TCL_ERROR;
    }

    int argi = 2;
    if (Tcl_GetInt(interp, argv[argi++], &matTag) != TCL_OK) {
      opserr << OpenSees::PromptValueError << "invalid matTag: patch quad matTag numSubdivIJ "
                "numSubdivJK yVertI zVertI yVertJ zVertJ yVertK zVertK yVertL "
                "zVertL\n";
      return TCL_ERROR;
    }

    if (Tcl_GetInt(interp, argv[argi++], &numSubdivIJ) != TCL_OK) {
      opserr << OpenSees::PromptValueError << "invalid numSubdivIJ: patch quad matTag numSubdivIJ "
                "numSubdivJK yVertI zVertI yVertJ zVertJ yVertK zVertK yVertL "
                "zVertL\n";
      return TCL_ERROR;
    }

    if (Tcl_GetInt(interp, argv[argi++], &numSubdivJK) != TCL_OK) {
      opserr << OpenSees::PromptValueError << "invalid numSubdivJK: patch quad matTag numSubdivIJ "
                "numSubdivJK yVertI zVertI yVertJ zVertJ yVertK zVertK yVertL "
                "zVertL\n";
      return TCL_ERROR;
    }

    for (int j = 0; j < 2; j++) {
      if (Tcl_GetDouble(interp, argv[argi++], &vertexCoordY) != TCL_OK) {
        opserr << OpenSees::PromptValueError << "invalid Coordinate y: ...yVertI zVertI yVertJ "
                  "zVertJ yVertK zVertK yVertL zVertL\n";
        return TCL_ERROR;
      }

      if (Tcl_GetDouble(interp, argv[argi++], &vertexCoordZ) != TCL_OK) {
        opserr << OpenSees::PromptValueError << "invalid Coordinate z: ...yVertI zVertI yVertJ "
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
      opserr << OpenSees::PromptValueError << "invalid number of parameters: patch circ matTag "
                "numSubdivCirc numSubdivRad yCenter zCenter intRad extRad "
                "startAng endAng\n";
      return TCL_ERROR;
    }

    if (Tcl_GetInt(interp, argv[argi++], &matTag) != TCL_OK) {
      opserr << OpenSees::PromptValueError << "invalid matTag: patch circ matTag numSubdivCirc "
                "numSubdivRad yCenter zCenter intRad extRad startAng endAng\n";
      return TCL_ERROR;
    }

    if (Tcl_GetInt(interp, argv[argi++], &numSubdivCirc) != TCL_OK) {
      opserr
          << OpenSees::PromptValueError << "invalid numSubdivCirc: patch circ matTag numSubdivCirc "
             "numSubdivRad yCenter zCenter intRad extRad startAng endAng\n";
      return TCL_ERROR;
    }

    if (Tcl_GetInt(interp, argv[argi++], &numSubdivRad) != TCL_OK) {
      opserr << OpenSees::PromptValueError << "invalid numSubdivRad: patch circ matTag numSubdivCirc "
                "numSubdivRad yCenter zCenter intRad extRad startAng endAng\n";
      return TCL_ERROR;
    }

    if (Tcl_GetDouble(interp, argv[argi++], &yCenter) != TCL_OK) {
      opserr << OpenSees::PromptValueError << "invalid yCenter: patch circ matTag numSubdivCirc "
                "numSubdivRad yCenter zCenter intRad extRad startAng endAng\n";
      return TCL_ERROR;
    }

    if (Tcl_GetDouble(interp, argv[argi++], &zCenter) != TCL_OK) {
      opserr << OpenSees::PromptValueError << "invalid zCenter: patch circ matTag numSubdivCirc "
                "numSubdivRad yCenter zCenter intRad extRad startAng endAng\n";
      return TCL_ERROR;
    }

    if (Tcl_GetDouble(interp, argv[argi++], &intRad) != TCL_OK) {
      opserr << OpenSees::PromptValueError << "invalid intRad: patch circ matTag numSubdivCirc "
                "numSubdivRad yCenter zCenter intRad extRad startAng endAng\n";
      return TCL_ERROR;
    }

    if (Tcl_GetDouble(interp, argv[argi++], &extRad) != TCL_OK) {
      opserr << OpenSees::PromptValueError << "invalid extRad: patch circ matTag numSubdivCirc "
                "numSubdivRad yCenter zCenter intRad extRad startAng endAng\n";
      return TCL_ERROR;
    }

    if (Tcl_GetDouble(interp, argv[argi++], &startAng) != TCL_OK) {
      opserr << OpenSees::PromptValueError << "invalid startAng: patch circ matTag numSubdivCirc "
                "numSubdivRad yCenter zCenter intRad extRad startAng endAng\n";
      return TCL_ERROR;
    }

    if (Tcl_GetDouble(interp, argv[argi++], &endAng) != TCL_OK) {
      opserr << OpenSees::PromptValueError << "invalid endAng: patch circ matTag numSubdivCirc "
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
    opserr << OpenSees::PromptValueError << "patch type is not available\n";
    return TCL_ERROR;
  }

  return TCL_OK;
}

// Add a fiber to a fiber section
int
TclCommand_addFiber(ClientData clientData, Tcl_Interp *interp, int argc,
                    TCL_Char ** const argv)
{
  enum class Position : int {
    Y, Z, Area, Material, End
  };
  ArgumentTracker<Position> tracker;
  std::set<int> positional;

  assert(clientData != nullptr);
  BasicModelBuilder* builder = static_cast<BasicModelBuilder*>(clientData);

  SectionBuilder* fiberSectionRepr = findSectionBuilder(builder, interp, argc, argv);
  if (fiberSectionRepr == nullptr) {
    opserr << OpenSees::PromptValueError << "cannot retrieve a section builder\n";
    return TCL_ERROR;
  }

  double yLoc, zLoc=0, area;
  int matTag;
  bool warn_2d_z = false;
  static constexpr int WarpModeCount = 3;
  double warp[WarpModeCount][3]{};
  int warp_arg = -1;
  for (int i=1; i<argc; i++) {
    if (strcmp(argv[i], "-section") == 0) {
      ++i;
    }
    else if (strcmp(argv[i], "-warp") == 0) {
      if (i + 1 >= argc) {
        opserr << OpenSees::PromptValueError << "missing warp argument\n";
        return TCL_ERROR;
      }
      warp_arg = i+1;
      i++;
    }
    else if (strcmp(argv[i], "-warn-2d-z") == 0) {
        warn_2d_z = true;
        i++;
    }
    else if (strcmp(argv[i], "-material") == 0) {
      if (argc == ++i || Tcl_GetInt(interp, argv[i], &matTag) != TCL_OK) {
        opserr << OpenSees::PromptValueError << "invalid material tag\n";
        return TCL_ERROR;
      }
      tracker.consume(Position::Material);
    }
    else if (strcmp(argv[i], "-area") == 0) {
      if (argc == ++i || Tcl_GetDouble(interp, argv[i], &area) != TCL_OK) {
        opserr << OpenSees::PromptValueError << "invalid area\n";
        return TCL_ERROR;
      }
      tracker.consume(Position::Area);
    }
    else if (strcmp(argv[i], "-y") == 0) {
      if (argc == ++i || Tcl_GetDouble(interp, argv[i], &yLoc) != TCL_OK) {
        opserr << OpenSees::PromptValueError << "invalid y coordinate\n";
        return TCL_ERROR;
      }
      tracker.consume(Position::Y);
    }
    else if (strcmp(argv[i], "-z") == 0) {
      if (argc == ++i || Tcl_GetDouble(interp, argv[i], &zLoc) != TCL_OK) {
        opserr << OpenSees::PromptValueError << "invalid z coordinate\n";
        return TCL_ERROR;
      }
      tracker.consume(Position::Z);
    }
    else {
      positional.insert(i);
    }
  }

  for (int i: positional) {
    switch (tracker.current()) {
      case Position::Y:
        if (Tcl_GetDouble(interp, argv[i], &yLoc) != TCL_OK) {
          opserr << OpenSees::PromptValueError << "invalid y coordinate\n";
          return TCL_ERROR;
        }
        tracker.consume(Position::Y);
        break;
      case Position::Z:
        if (Tcl_GetDouble(interp, argv[i], &zLoc) != TCL_OK) {
          opserr << OpenSees::PromptValueError << "invalid z coordinate\n";
          return TCL_ERROR;
        }
        tracker.consume(Position::Z);
        break;
      case Position::Area:
        if (Tcl_GetDouble(interp, argv[i], &area) != TCL_OK) {
          opserr << OpenSees::PromptValueError << "invalid area\n";
          return TCL_ERROR;
        }
        tracker.consume(Position::Area);
        break;
      case Position::Material:
        if (Tcl_GetInt(interp, argv[i], &matTag) != TCL_OK) {
          opserr << OpenSees::PromptValueError << "invalid material tag\n";
          return TCL_ERROR;
        }
        tracker.consume(Position::Material);
        break;
      default:
        opserr << OpenSees::PromptValueError << "unexpected argument at position " << i << "\n";
        return TCL_ERROR;
    }
  }

  if (builder->getNDM() == 2) {
    if (warn_2d_z && zLoc != 0.0) {
      opswrn << OpenSees::SignalWarning << "z coordinate ignored in 2D\n";
    }
  }
  if (tracker.current() != Position::End) {
    opserr << OpenSees::PromptValueError << "missing required arguments: ";
    while (tracker.current() != Position::End) {
      switch (tracker.current()) {
        case Position::Y:
          opserr << "y ";
          break;
        case Position::Z:
          opserr << "z ";
          break;
        case Position::Area:
          opserr << "area ";
          break;
        case Position::Material:
          opserr << "material ";
          break;
        case Position::End:
          break;
      }
      if (tracker.current() == Position::End)
        break;
      tracker.consume(tracker.current());
    }
    opserr << "\n";
    return TCL_ERROR;
  }

  //
  // process warping
  //
  int i_warp = 0;
  int          split_1_argc;
  const char **split_1_argv;
  if (warp_arg >= 0 && Tcl_SplitList(interp, argv[warp_arg], &split_1_argc, &split_1_argv) == TCL_OK) {
    int argi = 0;
    for (; i_warp<WarpModeCount; i_warp++) {
      if (argi >= split_1_argc)
        break;

      //
      int          split_argc;
      const char **split_argv;
      if (Tcl_SplitList(interp, split_1_argv[argi], &split_argc, &split_argv) != TCL_OK) {
        opserr << OpenSees::PromptValueError << "invalid warp\n";
        return TCL_ERROR;
      }

      if (split_argc != 3) {
        opserr << "WARNING warp parameter expected list of 3 floats\n";
          Tcl_Free((char *) split_argv);
          return TCL_ERROR;
      }

      for (int j = 0; j < 3; j++) {
        if (Tcl_GetDouble(interp, split_argv[j], &warp[i_warp][j]) != TCL_OK) {
          opserr << OpenSees::PromptValueError << "invalid warp\n";
          Tcl_Free((char *) split_argv);
          return TCL_ERROR;
        }
      }

      // Free memory allocated by Tcl_SplitList.
      Tcl_Free((char *) split_argv);
      argi++;
    }
}
  //
  // Add fiber to section builder
  //
  int ndm = builder->getNDM();
  int id = 0;
  if (ndm == 2) {
    Vector pos(2);
    pos(0) = yLoc;
    pos(1) = zLoc;
    id = fiberSectionRepr->addFiber(0, matTag, area, pos);
  } else if (ndm == 3) {
    Vector pos(2);
    pos(0) = yLoc;
    pos(1) = zLoc;
    id = fiberSectionRepr->addFiber(0, matTag, area, pos);
  }
  // set warping
  while (i_warp > 0) {
    if (0 > fiberSectionRepr->setWarping(id, i_warp-1, warp[i_warp-1])) {
      opserr << OpenSees::PromptValueError << "failed to set warping for fiber\n";
      return TCL_ERROR;
    }
    i_warp--;
  }

  if (id < 0) {
    opserr << OpenSees::PromptValueError << "Failed to add fiber to section\n";
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
    opserr << OpenSees::PromptValueError << "cannot retrieve section\n";
    return TCL_ERROR;
  }

  // make sure at least one other argument to contain patch type
  if (argc < 5) {
    opserr << OpenSees::PromptValueError << "invalid num args: Hfiber yLoc zLoc area matTag\n";
    return TCL_ERROR;
  }


  int matHTag;
  double yHLoc, zHLoc, Harea;

  if (Tcl_GetDouble(interp, argv[1], &yHLoc) != TCL_OK) {
    opserr << OpenSees::PromptValueError << "invalid yLoc: Hfiber yLoc zLoc area matTag\n";
    return TCL_ERROR;
  }
  if (Tcl_GetDouble(interp, argv[2], &zHLoc) != TCL_OK) {
    opserr << OpenSees::PromptValueError << "invalid zLoc: Hfiber yLoc zLoc area matTag\n";
    return TCL_ERROR;
  }
  if (Tcl_GetDouble(interp, argv[3], &Harea) != TCL_OK) {
    opserr << OpenSees::PromptValueError << "invalid area: Hfiber yLoc zLoc area matTag\n";
    return TCL_ERROR;
  }

  if (Tcl_GetInt(interp, argv[4], &matHTag) != TCL_OK) {
    opserr << OpenSees::PromptValueError << "invalid matTag: Hfiber yLoc zLoc area matTag\n";
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
    opserr << OpenSees::PromptValueError << "cannot retrieve section\n";
    return TCL_ERROR;
  }

  // make sure at least one other argument to contain layer type
  if (argc < 2) {
    opserr << OpenSees::PromptValueError << "need to specify a layer type \n";
    return TCL_ERROR;
  }

  // check argv[1] for type of layer and create the object
  if (strcmp(argv[1], "straight") == 0 ||
      strcmp(argv[1], "line")     == 0) {
    if (argc < 9) {
      opserr << OpenSees::PromptValueError << "invalid number of parameters: layer straight matTag "
                "numReinfBars reinfBarArea yStartPt zStartPt yEndPt zEndPt\n";
      return TCL_ERROR;
    }

    int matTag, numReinfBars;
    double reinfBarArea;
    double yStartPt, zStartPt, yEndPt, zEndPt;

    int argi = 2;

    if (Tcl_GetInt(interp, argv[argi++], &matTag) != TCL_OK) {
      opserr << OpenSees::PromptValueError << "invalid matTag: layer straight matTag numReinfBars "
                "reinfBarArea  yStartPt zStartPt yEndPt zEndPt\n";
      return TCL_ERROR;
    }

    if (Tcl_GetInt(interp, argv[argi++], &numReinfBars) != TCL_OK) {
      opserr << OpenSees::PromptValueError << "invalid numReinfBars: layer straight matTag "
                "numReinfBars reinfBarArea  yStartPt zStartPt yEndPt zEndPt\n";
      return TCL_ERROR;
    }

    if (Tcl_GetDouble(interp, argv[argi++], &reinfBarArea) != TCL_OK) {
      opserr << OpenSees::PromptValueError << "invalid reinfBarArea: layer straight matTag "
                "numReinfBars reinfBarArea  yStartPt zStartPt yEndPt zEndPt\n";
      return TCL_ERROR;
    }

    if (Tcl_GetDouble(interp, argv[argi++], &yStartPt) != TCL_OK) {
      opserr << OpenSees::PromptValueError << "invalid yStartPt: layer straight matTag numReinfBars "
                "reinfBarArea  yStartPt zStartPt yEndPt zEndPt\n";
      return TCL_ERROR;
    }

    if (Tcl_GetDouble(interp, argv[argi++], &zStartPt) != TCL_OK) {
      opserr << OpenSees::PromptValueError << "invalid zStartPt: layer straight matTag numReinfBars "
                "reinfBarArea  yStartPt zStartPt yEndPt zEndPt\n";
      return TCL_ERROR;
    }

    if (Tcl_GetDouble(interp, argv[argi++], &yEndPt) != TCL_OK) {
      opserr << OpenSees::PromptValueError << "invalid yEndPt: layer straight matTag numReinfBars "
                "reinfBarArea  yStartPt zStartPt yEndPt zEndPt\n";
      return TCL_ERROR;
    }

    if (Tcl_GetDouble(interp, argv[argi++], &zEndPt) != TCL_OK) {
      opserr << OpenSees::PromptValueError << "invalid zEndPt: layer straight matTag numReinfBars "
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
      opserr << OpenSees::PromptValueError << "invalid number of parameters: layer circ matTag "
                "numReinfBars reinfBarArea yCenter zCenter arcRadius <startAng "
                "endAng>\n";
      return TCL_ERROR;
    }

    int  matTag, numReinfBars;
    double reinfBarArea;
    double yCenter, zCenter, radius, startAng, endAng;

    int argi = 2;

    if (Tcl_GetInt(interp, argv[argi++], &matTag) != TCL_OK) {
      opserr << OpenSees::PromptValueError << "invalid matTag\n";
      return TCL_ERROR;
    }

    if (Tcl_GetInt(interp, argv[argi++], &numReinfBars) != TCL_OK) {
      opserr << OpenSees::PromptValueError << "invalid numReinfBars\n";
      return TCL_ERROR;
    }

    if (Tcl_GetDouble(interp, argv[argi++], &reinfBarArea) != TCL_OK) {
      opserr << OpenSees::PromptValueError << "invalid reinfBarArea\n";
      return TCL_ERROR;
    }

    if (Tcl_GetDouble(interp, argv[argi++], &yCenter) != TCL_OK) {
      opserr << OpenSees::PromptValueError << "invalid yCenter\n";
      return TCL_ERROR;
    }

    if (Tcl_GetDouble(interp, argv[argi++], &zCenter) != TCL_OK) {
      opserr << OpenSees::PromptValueError << "invalid zCenter\n";
      return TCL_ERROR;
    }

    if (Tcl_GetDouble(interp, argv[argi++], &radius) != TCL_OK) {
      opserr << OpenSees::PromptValueError << "invalid radius\n";
      return TCL_ERROR;
    }

    bool anglesSpecified = false;

    if (argc > 9) {
      if (Tcl_GetDouble(interp, argv[argi++], &startAng) != TCL_OK) {
        opserr << OpenSees::PromptValueError << "invalid startAng\n";
        return TCL_ERROR;
      }

      if (Tcl_GetDouble(interp, argv[argi++], &endAng) != TCL_OK) {
        opserr << OpenSees::PromptValueError << "invalid endAng\n";
        return TCL_ERROR;
      }

      anglesSpecified = true;
    }

    // create the reinforcing layer

    static Vector center(2);

    center(0) = yCenter;
    center(1) = zCenter;

    int error = -1;
    // construct and add to section
    if (anglesSpecified) {
      // Construct arc
      CircReinfLayer reinfLayer(matTag, numReinfBars, reinfBarArea,
                                center, radius, startAng, endAng);
      error = fiberSectionRepr->addLayer(reinfLayer);
    } else {
      // Construct circle
      CircReinfLayer reinfLayer(matTag, numReinfBars, reinfBarArea,
                                center, radius);
      error = fiberSectionRepr->addLayer(reinfLayer);
    }

    if (error)
      return TCL_ERROR;

  } else {
    opserr << OpenSees::PromptValueError << "reinforcing layer type is not available\n";
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
  // Now parse the output file containing the fiber data,
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
               << "\n";
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
