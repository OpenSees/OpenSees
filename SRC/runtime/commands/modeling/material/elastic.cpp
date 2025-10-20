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
// cmp
//
#include <tcl.h>
#include <string.h>
#include <BasicModelBuilder.h>

#include <set>
#include <cstdlib>
#include "BasicModelBuilder.h"
#include "Logging.h"
#include "Parsing.h"
#include "ArgumentTracker.h"
#include "isotropy.h"

#include <ElasticMaterial.h>
#include <ElasticIsotropic.h>
#include <ElasticIsotropicMaterial.h>
#include <ElasticOrthotropicMaterial.h>

#include <ElasticIsotropicMaterial.h>
#include "ElasticIsotropicPlaneStress2D.h"
#include "ElasticIsotropicPlaneStrain2D.h"
#include <material/elastic/ElasticIsotropic3DThermal.h>
#include <material/elastic/ElasticIsotropicAxiSymm.h>
#include <material/elastic/ElasticIsotropicBeamFiber.h>
#include <material/elastic/ElasticIsotropicBeamFiber2d.h>
#include <material/elastic/ElasticIsotropicMaterial.h>
#include <material/elastic/ElasticIsotropicMaterialThermal.h>
#include <material/elastic/ElasticIsotropicPlateFiber.h>
#include <material/elastic/ElasticIsotropicThreeDimensional.h>

#include <PlaneStressMaterial.h>
#include <BeamFiberMaterial.h>



template <typename Position>
int
TclCommand_newElasticParser(ClientData clientData, Tcl_Interp *interp,
                            int argc, TCL_Char ** const argv)
{
  assert(clientData != nullptr);

  ArgumentTracker<Position> tracker;
  std::set<int> positional;

  int niso = (
    (Position::E      < Position::EndRequired) +
    (Position::G      < Position::EndRequired) +
    (Position::Nu     < Position::EndRequired) +
    (Position::K      < Position::EndRequired) +
    (Position::Lambda < Position::EndRequired)
  );

  int tag;
  double density = 0.0;
  // Isotropy
  IsotropicConstants consts {};
  // Viscosity
  double eta=0;


  // Isotropy
  IsotropicParse iso {consts, niso};
  if (TclCommand_setIsotropicParameters((ClientData)&iso, interp, argc, argv) == TCL_OK) {
    tracker.consume(Position::E);
    tracker.consume(Position::G);
    tracker.consume(Position::Nu);
    tracker.consume(Position::K);
    tracker.consume(Position::Lambda);
  }

  //
  // Keywords
  //
  for (int i=2; i<argc; i++) {
    if (iso.positions.find(i) != iso.positions.end()) {
      continue;
    }

    else if (strcmp(argv[i], "-rho") == 0 || strcmp(argv[i], "-density") == 0) {
      if (++i >= argc) {
          opserr << "Missing value for option " << argv[i-1] << "\n";
          return TCL_ERROR;
      }
      if (Tcl_GetDouble(interp, argv[i], &density) != TCL_OK) {
          opserr << "Invalid density value for option " << argv[i-1] << "\n";
          return TCL_ERROR;
      }
    }
    else
      positional.insert(i);
  }

  //
  // Positional arguments
  //
  for (int i : positional) {
  
    if (tracker.current() == Position::EndRequired)
      tracker.increment();

    switch (tracker.current()) {
      case Position::Tag :
        if (Tcl_GetInt(interp, argv[i], &tag) != TCL_OK) {
            opserr << OpenSees::PromptParseError << "invalid section Elastic tag.\n";
            return TCL_ERROR;           
        } else {
          tracker.increment();
          break;
        }
      case Position::E:
        if (Tcl_GetDouble (interp, argv[i], &consts.E) != TCL_OK) {
            opserr << OpenSees::PromptParseError << "invalid E.\n";
            return TCL_ERROR;
        } else {
          tracker.increment();
          break;
        }

      case Position::G:
        if (Tcl_GetDouble (interp, argv[i], &consts.G) != TCL_OK) {
            opserr << OpenSees::PromptParseError << "invalid G.\n";
            return TCL_ERROR;
        } else {
          tracker.increment();
          break;
        }

      case Position::K:
        if (Tcl_GetDouble (interp, argv[i], &consts.K) != TCL_OK) {
            opserr << OpenSees::PromptParseError << "invalid K.\n";
            return TCL_ERROR;
        } 
        tracker.increment();
        break;

      case Position::Nu:
        if (Tcl_GetDouble (interp, argv[i], &consts.nu) != TCL_OK) {
            opserr << OpenSees::PromptParseError << "invalid nu.\n";
            return TCL_ERROR;
        } 
        tracker.increment();
        break;
      
      case Position::Lambda:
        if (Tcl_GetDouble (interp, argv[i], &consts.lambda) != TCL_OK) {
            opserr << OpenSees::PromptParseError << "invalid Lame lambda.\n";
            return TCL_ERROR;
        } 
        tracker.increment();
        break;
      
      case Position::Eta:
        if (Tcl_GetDouble (interp, argv[i], &eta) != TCL_OK) {
            opserr << OpenSees::PromptParseError << "invalid eta.\n";
            return TCL_ERROR;
        } 
        tracker.increment();
        break;

      case Position::Density:
        if (Tcl_GetDouble (interp, argv[i], &density) != TCL_OK) {
            opserr << OpenSees::PromptParseError << "invalid density.\n";
            return TCL_ERROR;
        } 
        tracker.increment();
        break;

      case Position::EndRequired:
        // This will not be reached
        break;

      case Position::End:
        opserr << OpenSees::PromptParseError << "unexpected argument " << argv[i] << ".\n";
        return TCL_ERROR;
    }
  }

  if (tracker.current() < Position::EndRequired) {
    opserr << OpenSees::PromptParseError
            << "missing required arguments: ";

    while (tracker.current() != Position::EndRequired) {
      switch (tracker.current()) {
        case Position::Tag :
          opserr << "tag ";
          break;
        case Position::E:
          opserr << "E ";
          break;
        case Position::G:
          opserr << "G ";
          break;
        case Position::K:
          opserr << "K ";
          break;
        case Position::Nu:
          opserr << "nu ";
          break;
        case Position::Lambda:
          opserr << "lambda ";
          break;

        case Position::EndRequired:
        case Position::End:
        default:
          break;
      }

      if (tracker.current() == Position::EndRequired)
        break;

      tracker.consume(tracker.current());
    }

    opserr << "\n";

    return TCL_ERROR;
  }

  //
  // Create the material
  //
  BasicModelBuilder *builder = static_cast<BasicModelBuilder*>(clientData);

  if ((strcmp(argv[1], "ElasticIsotropic") == 0) ||
      (strcmp(argv[1], "Elastic") == 0)) {
    double E = consts.E;
    double nu = consts.nu;
    if (builder->addTaggedObject<NDMaterial>(*new ElasticIsotropicMaterial(tag, E, nu, density)) != TCL_OK ) {
      return TCL_ERROR;
    }
    if (strcmp(argv[0], "material") == 0) {
      if (builder->addTaggedObject<UniaxialMaterial>(*new ElasticMaterial(tag, E, eta, E, density)) != TCL_OK ) {
        return TCL_ERROR;
      }
      if (builder->addTaggedObject<Mate<3>>(*new ElasticIsotropic<3>(tag, E, nu, density)) != TCL_OK ) {
        return TCL_ERROR;
      }
    }
    return TCL_OK;
  }
  else if (strcmp(argv[1], "ElasticIsotropic3D") == 0) {
    double E = consts.E;
    double nu = consts.nu;
    if (builder->addTaggedObject<NDMaterial>(*new ElasticIsotropicThreeDimensional(tag, E, nu, density)) != TCL_OK ) {
      return TCL_ERROR;
    }
    return TCL_OK;
  }

  return TCL_ERROR;
}

int
TclCommand_newElasticMaterial(ClientData clientData, Tcl_Interp *interp,
                              Tcl_Size argc, TCL_Char ** const argv)
{
  //
  if (strcmp(argv[1], "ElasticIsotropic") == 0 ||
      strcmp(argv[1], "Elastic") == 0 ||
      strcmp(argv[1], "PlaneStressSimplifiedJ2") == 0) {

    // "ElasticIsotropic" tag?  E?  nu?  rho?
    enum class Position : int {
      Tag, E, Nu,  EndRequired, 
      Density, End,
      G, K, Lambda, Eta
    };
    return TclCommand_newElasticParser<Position>(clientData, interp, argc, argv);
  }

  return TCL_ERROR;
}

int
TclCommand_newElasticUniaxialMaterial(ClientData clientData, Tcl_Interp *interp,
                                      Tcl_Size argc, TCL_Char ** const argv)
{
  enum class Position : int {
    Tag, Epos, EndRequired, 
    Eta, Eneg, Density, End
  };

  ArgumentTracker<Position> tracker;
  std::set<int> positional;

  int tag;
  double Epos = 0.0;
  double Eneg = 0.0;
  double density = 0.0;
  double eta = 0.0;

  if (Tcl_GetInt(interp, argv[2], &tag) != TCL_OK) {
    opserr << OpenSees::PromptParseError 
           << "invalid tag."
           << OpenSees::SignalMessageEnd;
    return TCL_ERROR; 
  }
  tracker.consume(Position::Tag);

  for (int i=3; i<argc; i++) {
    if (strcmp(argv[i], "-Epos") == 0 || strcmp(argv[i], "-E") == 0) {
      if (++i >= argc) {
          opserr << "Missing value for option " << argv[i-1] << "\n";
          return TCL_ERROR;
      }
      if (Tcl_GetDouble(interp, argv[i], &Epos) != TCL_OK) {
          opserr << "Invalid value for option " << argv[i-1] << "\n";
          return TCL_ERROR;
      }
      tracker.consume(Position::Epos);
      if (!tracker.contains(Position::Eneg)) {
        Eneg = Epos; // Default to Epos if Eneg is not specified
      }
    }
    else if (strcmp(argv[i], "-Eneg") == 0) {
      if (++i >= argc) {
          opserr << "Missing value for option " << argv[i-1] << "\n";
          return TCL_ERROR;
      }
      if (Tcl_GetDouble(interp, argv[i], &Eneg) != TCL_OK) {
          opserr << "Invalid value for option " << argv[i-1] << "\n";
          return TCL_ERROR;
      }
      tracker.consume(Position::Eneg);
    }
    else if (strcmp(argv[i], "-rho") == 0 || strcmp(argv[i], "-density") == 0) {
      if (++i >= argc) {
          opserr << "Missing value for option " << argv[i-1] << "\n";
          return TCL_ERROR;
      }
      if (Tcl_GetDouble(interp, argv[i], &density) != TCL_OK) {
          opserr << "Invalid density value for option " << argv[i-1] << "\n";
          return TCL_ERROR;
      }
      tracker.consume(Position::Density);
    }
    else if (strcmp(argv[i], "-eta") == 0 || strcmp(argv[i], "-viscosity") == 0) {
      if (++i >= argc) {
          opserr << "Missing value for option " << argv[i-1] << "\n";
          return TCL_ERROR;
      }
      if (Tcl_GetDouble(interp, argv[i], &eta) != TCL_OK) {
          opserr << "Invalid value for option " << argv[i-1] << "\n";
          return TCL_ERROR;
      }
      tracker.consume(Position::Eta);
    }
    else
      positional.insert(i);
  }

  for (int i : positional) {
    if (tracker.current() == Position::EndRequired)
      tracker.increment();

    switch (tracker.current()) {
      case Position::Tag :
        if (Tcl_GetInt(interp, argv[i], &tag) != TCL_OK) {
            opserr << OpenSees::PromptParseError << "invalid tag.\n";
            return TCL_ERROR;           
        } else {
          tracker.increment();
          break;
        }
      case Position::Epos:
        if (Tcl_GetDouble (interp, argv[i], &Epos) != TCL_OK) {
            opserr << OpenSees::PromptParseError << "invalid Epos.\n";
            return TCL_ERROR;
        } else {
          Eneg = Epos;
          tracker.increment();
          break;
        }
      case Position::Eneg:
        if (Tcl_GetDouble (interp, argv[i], &Eneg) != TCL_OK) {
            opserr << OpenSees::PromptParseError << "invalid Eneg.\n";
            return TCL_ERROR;
        } else {
          tracker.increment();
          break;
        }
      case Position::Density:
        if (Tcl_GetDouble (interp, argv[i], &density) != TCL_OK) {
            opserr << OpenSees::PromptParseError << "invalid density.\n";
            return TCL_ERROR;
        } else {
          tracker.increment();
          break;
        }
      case Position::Eta:
        if (Tcl_GetDouble (interp, argv[i], &eta) != TCL_OK) {
            opserr << OpenSees::PromptParseError << "invalid eta.\n";
            return TCL_ERROR;
        } else {
          tracker.increment();
          break;
        }
      case Position::EndRequired:
      case Position::End:
        opserr << OpenSees::PromptParseError 
               << "unexpected argument " << argv[i] << ".\n";
        return TCL_ERROR;
    }
  }
  

  if (tracker.current() < Position::EndRequired) {
    opserr << OpenSees::PromptParseError
            << "missing required arguments: ";

    while (tracker.current() != Position::EndRequired) {
      switch (tracker.current()) {
        case Position::Tag :
          opserr << "tag ";
          break;
        case Position::Epos:
          opserr << "Epos ";
          break;
        case Position::Eneg:
          opserr << "Eneg ";
          break;
        case Position::Density:
          opserr << "density ";
          break;
        case Position::Eta:
          opserr << "eta ";
          break;

        case Position::EndRequired:
        case Position::End:
        default:
          break;
      }

      if (tracker.current() == Position::EndRequired)
        break;

      tracker.consume(tracker.current());
    }

    opserr << "\n";

    return TCL_ERROR;
  }

  BasicModelBuilder *builder = static_cast<BasicModelBuilder*>(clientData);
  if (strcmp(argv[1], "Elastic") == 0) {
    if (builder->addTaggedObject<UniaxialMaterial>(*new ElasticMaterial(tag, Epos, eta, Eneg, density)) != TCL_OK ) {
      return TCL_ERROR;
    }
    return TCL_OK;
  }

  return TCL_ERROR;
}


enum class MaterialSymmetry {
  Triclinic     , // 21
  Monoclinic    , // 13
  Orthorhombic  , //  9
  Tetragonal    , //  7
//Tetragonal    , //  6
  Rhombohedral  , //  7
//Rhombohedral  , //  6
  Hexagonal     , //  5
  Cubic         , //  3
  Isotropic       //  2
};

