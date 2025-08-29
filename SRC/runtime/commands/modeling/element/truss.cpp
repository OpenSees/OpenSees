//===----------------------------------------------------------------------===//
//
//                                   xara
//                              https://xara.so
//
//===----------------------------------------------------------------------===//
//
// Copyright (c) 2025, Claudio M. Perez
// All rights reserved.  No warranty, explicit or implicit, is provided.
//
// This source code is licensed under the BSD 2-Clause License.
// See LICENSE file or https://opensource.org/licenses/BSD-2-Clause
//
//===----------------------------------------------------------------------===//
//
// Written: cmp 
// April 2025
//
#include <set>
#include <assert.h>
#include <string.h>
// OpenSees
#include <tcl.h>
#include <Domain.h>
#include <Logging.h>
#include <Parsing.h>
#include <ArgumentTracker.h>
#include <FrameSection.h>
#include <UniaxialMaterial.h>
#include <BasicModelBuilder.h>
// Elements
#include <TrussSection.h>
#include <CorotTruss.h>
#include <CorotTrussSection.h>
#include <FrameFiberSection3d.h>
#include <FrameSolidSection3d.h>

#ifdef _MSC_VER
#  include <string.h>
#  define strcasecmp _stricmp
#else
#  include <strings.h>
#endif
#define strcmp strcasecmp

//
//  element Truss        $tag $iNode $jNode $A $matTag <-rho $rho> <-cMass $flag> <-doRayleigh $flag> <-useInitialDisp $flag>
//  element Truss        $tag $iNode $jNode $sectTag <-rho $rho> <-cMass $flag> <-doRayleigh $flag>
//  element TrussSection $tag $iNode $jNode $sectTag <-rho $rho> <-cMass $flag> <-doRayleigh $flag>
//


template <typename Positions>
static int 
CreateTruss(ClientData clientData, Tcl_Interp *interp, int argc, 
                TCL_Char **argv)
{
  assert(clientData != nullptr);

  // Parsing is performed in three steps:
  // 1. Collect nodes
  // 2. Parse out keywords
  // 3. Parse remaining positional arguments.

  ArgumentTracker<Positions> tracker;
  std::set<int> positional;


  BasicModelBuilder* builder = static_cast<BasicModelBuilder*>(clientData);

  int tag;
  if (argc < 3 || (Tcl_GetInt(interp, argv[2], &tag) != TCL_OK)) {
    opserr << OpenSees::PromptValueError
           << "failed to read integer tag\n";
    return TCL_ERROR;
  }
  tracker.consume(Positions::Tag);

  double rho     = 0.0;
  int doRayleigh = 0; // by default rayleigh not done
  int cMass      = 0; // by default use lumped mass matrix
  int ndm        = builder->getNDM();
  double area    = 0;
  UniaxialMaterial *material = nullptr;
  FrameSection *section = nullptr;

  //
  // 1. Parse nodes
  //
  int nodes[2] = {0, 0};
  int node_end = 5;
  if (argc < 4 || (Tcl_GetInt(interp, argv[3], &nodes[0]) != TCL_OK)) {
    opserr << OpenSees::PromptValueError
           << "failed to read integer node tag\n";
    return TCL_ERROR;
  }
  tracker.consume(Positions::iNode);
  if (argc < 5 || (Tcl_GetInt(interp, argv[4], &nodes[1]) != TCL_OK)) {
    opserr << OpenSees::PromptValueError
           << "failed to read integer node tag\n";
    return TCL_ERROR;
  }
  tracker.consume(Positions::jNode);

  //
  // 2. Keywords 
  //
  for (int i=node_end; i<argc; i++) {
    if (strcmp(argv[i], "-rho") == 0) {
      if (argc == ++i || Tcl_GetDouble(interp, argv[i], &rho) != TCL_OK) {
        opserr << OpenSees::PromptValueError
               << "failed to read density\n";
        return TCL_ERROR;
      }
      tracker.consume(Positions::Density);
    }

    else if (strcmp(argv[i], "-section") == 0) {
      int sec;
      if (argc == ++i || Tcl_GetInt(interp, argv[i], &sec) != TCL_OK) {
        opserr << OpenSees::PromptValueError
               << "failed to read section tag\n";
        return TCL_ERROR;
      }
      section = builder->getTypedObject<FrameSection>(sec);
      if (section == nullptr)
        return TCL_ERROR;

      tracker.consume(Positions::Section);
      tracker.consume(Positions::Area);
      tracker.consume(Positions::Material);
    }

    else if (strcmp(argv[i], "-material") == 0) {
      int mat;
      if (argc == ++i || Tcl_GetInt(interp, argv[i], &mat) != TCL_OK) {
        opserr << OpenSees::PromptValueError
               << "failed to read material tag\n";
        return TCL_ERROR;
      }
      
      material = builder->getTypedObject<UniaxialMaterial>(mat);
      if (material == nullptr)
        return TCL_ERROR;
      tracker.consume(Positions::Material);
    }
    else if ((strcmp(argv[i], "-A") == 0) ||
             (strcmp(argv[i], "-area") == 0)) {
      if (argc == ++i || Tcl_GetDouble(interp, argv[i], &area) != TCL_OK) {
        opserr << OpenSees::PromptValueError
               << "failed to read area\n";
        return TCL_ERROR;
      }
      tracker.consume(Positions::Area);
    }
    else if (strcmp(argv[i], "-cMass") == 0) {
      if (argc == ++i || Tcl_GetInt(interp, argv[i], &cMass) != TCL_OK) {
        opserr << OpenSees::PromptValueError
               << "failed to read mass flag\n";
        return TCL_ERROR;
      }
      tracker.consume(Positions::MassFlag);
    }
    else if (strcmp(argv[i], "-doRayleigh") == 0) {
      if (argc == ++i || Tcl_GetInt(interp, argv[i], &doRayleigh) != TCL_OK) {
        opserr << OpenSees::PromptValueError
               << "failed to read rayleigh flag\n";
        return TCL_ERROR;
      }
      tracker.consume(Positions::RayleighFlag);
    }
    else if (strcmp(argv[i], "-useInitialDisp") == 0) {
      if (argc == ++i || Tcl_GetInt(interp, argv[i], &doRayleigh) != TCL_OK) {
        opserr << OpenSees::PromptValueError
               << "failed to read initial disp flag\n";
        return TCL_ERROR;
      }
      tracker.consume(Positions::UseInitialDisp);
    }
    else {
      positional.insert(i);
    }
  }

  //
  // 3. Positional arguments
  //
  for (int i : positional) {

    if (tracker.current() == Positions::EndRequired)
      tracker.increment();

    switch (tracker.current()) {
      case Positions::Tag:
        if (Tcl_GetInt(interp, argv[i], &tag) != TCL_OK) {
          opserr << OpenSees::PromptValueError
                 << "failed to read integer tag\n";
          return TCL_ERROR;
        }
        tracker.increment();
        break;
      
      case Positions::iNode:
        if (Tcl_GetInt(interp, argv[i], &nodes[0]) != TCL_OK) {
          opserr << OpenSees::PromptValueError
                 << "failed to read integer node tag\n";
          return TCL_ERROR;
        }
        tracker.increment();
        break;
      case Positions::jNode:
        if (Tcl_GetInt(interp, argv[i], &nodes[1]) != TCL_OK) {
          opserr << OpenSees::PromptValueError
                 << "failed to read integer node tag\n";
          return TCL_ERROR;
        }
        tracker.increment();
        break;

      case Positions::Material: {
        int mat;
        if (Tcl_GetInt(interp, argv[i], &mat) != TCL_OK) {
          opserr << OpenSees::PromptValueError
                 << "failed to read integer material tag\n";
          return TCL_ERROR;
        }
        material = builder->getTypedObject<UniaxialMaterial>(mat);
        if (material == nullptr) {
          return TCL_ERROR;
        }
        tracker.increment();
        break;
      }

      case Positions::Density:
        if (Tcl_GetDouble(interp, argv[i], &rho) != TCL_OK) {
          opserr << OpenSees::PromptValueError
                 << "failed to read density\n";
          return TCL_ERROR;
        }
        tracker.increment();
        break;

      case Positions::MassFlag:
        if (Tcl_GetInt(interp, argv[i], &cMass) != TCL_OK) {
          opserr << OpenSees::PromptValueError
                 << "failed to read mass flag\n";
          return TCL_ERROR;
        }
        tracker.increment();
        break;
      case Positions::RayleighFlag:
        if (Tcl_GetInt(interp, argv[i], &doRayleigh) != TCL_OK) {
          opserr << OpenSees::PromptValueError
                 << "failed to read rayleigh flag\n";
          return TCL_ERROR;
        }
        tracker.increment();
        break;
      case Positions::UseInitialDisp:
        if (Tcl_GetInt(interp, argv[i], &doRayleigh) != TCL_OK) {
          opserr << OpenSees::PromptValueError
                 << "failed to read initial disp flag\n";
          return TCL_ERROR;
        }
        tracker.increment();
        break;
      
      case Positions::Area:
        if (Tcl_GetDouble(interp, argv[i], &area) != TCL_OK) {
          opserr << OpenSees::PromptValueError
                 << "failed to read area\n";
          return TCL_ERROR;
        }
        tracker.increment();
        break; 
    
      case Positions::Section:
        int sec;
        if (Tcl_GetInt(interp, argv[i], &sec) != TCL_OK) {
          opserr << OpenSees::PromptValueError
                 << "failed to read section tag\n";
          return TCL_ERROR;
        }
        section = builder->getTypedObject<FrameSection>(sec);
        if (section == nullptr)
          return TCL_ERROR;
        tracker.increment();
        break;

      case Positions::EndRequired:
        // This will not be reached
        break;

      case Positions::End:
        opserr << OpenSees::PromptValueError
               << "unexpected argument " << argv[i] << "\n";
        return TCL_ERROR;
    }
  }

  //
  // 4. Check required positional arguments
  //
  if (tracker.current() < Positions::EndRequired) {
    opserr << OpenSees::PromptValueError
           << "missing required arguments: ";
    while (tracker.current() != Positions::EndRequired) {
      switch (tracker.current()) {
        case Positions::Tag:
          opserr << "tag ";
          break;
        case Positions::Material:
          opserr << "material ";
          break;
        default:
          break;
      }
      tracker.increment();
    }
    return TCL_ERROR;
  }

  //
  //
  //
  if (section == nullptr) {
    if (tracker.contains(Positions::Area)) {
        opserr << OpenSees::PromptValueError
               << "missing required argument area\n";
        return TCL_ERROR;
    }
    if (material == nullptr) {
      opserr << OpenSees::PromptValueError
             << "missing required argument material\n";
      return TCL_ERROR;
    }

    auto fiber_section = new FrameFiberSection3d(0, 1, nullptr, true, 0.0, 0);
    fiber_section->addFiber(*material, area, 0, 0);
    section = fiber_section;
  }

  //
  //
  //
  if (strstr(argv[1], "Corot") != nullptr && strstr(argv[1], "2") == nullptr) {
    if (ndm != 3) {
      opserr << OpenSees::PromptValueError
             << "CorotTruss only valid in 3D\n";
      return TCL_ERROR;
    }

    builder->getDomain()->addElement(new CorotTrussSection(tag, ndm, nodes[0], nodes[1], *section, rho, cMass, doRayleigh));
  }

  else if (strstr(argv[1], "Corot") == nullptr && strstr(argv[1], "2") == nullptr) {
    builder->getDomain()->addElement(new TrussSection(tag, ndm, nodes[0], nodes[1], *section, rho, cMass, doRayleigh));
  }
  return TCL_OK;
}


int
TclCommand_addTruss(ClientData clientData, 
                    Tcl_Interp *interp,  
                    Tcl_Size argc, 
                    TCL_Char ** const argv)
{

  if (strstr(argv[1], "ection") != nullptr || 
     ((strcasecmp(argv[1], "Truss") == 0) && (argc == 6))) {
    enum class Arguments : int {
        Tag,
        iNode,
        jNode,
        Section,
      EndRequired,
      End,
        Density,
        MassFlag,
        RayleighFlag,
        UseInitialDisp,
        Area,
        Material,
    };
    return CreateTruss<Arguments>(clientData, interp, argc, argv);
  }
  else {
    enum class Arguments : int {
        Tag,
        iNode,
        jNode,
        Area,
        Material,
      EndRequired,
      End,
        Density,
        MassFlag,
        RayleighFlag,
        UseInitialDisp,
        Section,
    };
    return CreateTruss<Arguments>(clientData, interp, argc, argv);
  }

}