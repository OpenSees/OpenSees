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
#include <tcl.h>
#include <set>
#include <string.h>
#ifdef _MSC_VER 
#  include <string.h>
#  define strcasecmp _stricmp
#else
#  include <strings.h>
#endif
#include <Logging.h>
#include <Parsing.h>
#include <ArgumentTracker.h>
#include <BasicModelBuilder.h>
#include "isotropy.h"

#include <damage/FariaPlasticDamage3d.h>

int
TclCommand_newConcreteMaterial(ClientData clientData, Tcl_Interp *interp,
                                int argc, TCL_Char ** const argv)
{

  assert(clientData != nullptr);
  enum class Position {
    Tag, E, Nu, PeakTension, PeakCompression, EndRequired,
    Beta, Ap, An, Bn,
    Density,
    G, K, Lambda,
    End,
  };
  ArgumentTracker<Position> tracker;
  std::set<int> positional;


  // Count the number of required isotropic parameters.
  // This is needed to accommodate uniaxial materials, which generally
  // only require one isotropic parameter (ie, E).
  int niso = (
    (Position::E      < Position::EndRequired) +
    (Position::G      < Position::EndRequired) +
    (Position::Nu     < Position::EndRequired) +
    (Position::K      < Position::EndRequired) +
    (Position::Lambda < Position::EndRequired)
  );

  // 
  // Values we're parsing for
  //
  int tag;
  double density = 0.0;
  // Isotropy
  IsotropicConstants consts {};
  // Plasticity
  double Fc, Ft=0;
  // Hardening
  double beta = 0.6;
  double Ap = 0.5,
         An = 2.0,
         Bn = 0.75;
  

  //
  // 1. Keyword arguments
  //

  // Isotropy
  IsotropicParse iso {consts, niso};
  if (TclCommand_setIsotropicParameters((ClientData)&iso, interp, argc, argv) == TCL_OK) {
    tracker.consume(Position::E);
    tracker.consume(Position::G);
    tracker.consume(Position::Nu);
    tracker.consume(Position::K);
    tracker.consume(Position::Lambda);
  }

  // Other arguments
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

    // Compression
    else if (strcmp(argv[i], "-Fc") == 0 || 
             strcmp(argv[i], "-fc") == 0) {
      if (++i >= argc) {
          opserr << "Missing value for option " << argv[i-1] << "\n";
          return TCL_ERROR;
      }
      if (Tcl_GetDouble(interp, argv[i], &Fc) != TCL_OK) {
          opserr << "Invalid " << &argv[i-1][1] << " value " << argv[i] << "\n";
          return TCL_ERROR;
      }
      tracker.consume(Position::PeakCompression);
    }
    // Tension
    else if ((strcmp(argv[i], "-Ft") == 0)) {
      if (++i >= argc) {
          opserr << "Missing value for option " << argv[i-1] << "\n";
          return TCL_ERROR;
      }
      if (Tcl_GetDouble(interp, argv[i], &Ft) != TCL_OK) {
          opserr << "Invalid " << &argv[i-1][1] << " value " << argv[i] << "\n";
          return TCL_ERROR;
      }
      tracker.consume(Position::PeakTension);
    }
    //
    // Hardening
    //
    else if (strcmp(argv[i], "-beta") == 0) {
      if (++i >= argc) {
          opserr << "Missing value for option " << argv[i-1] << "\n";
          return TCL_ERROR;
      }
      if (Tcl_GetDouble(interp, argv[i], &beta) != TCL_OK) {
          opserr << "Invalid value for option " << argv[i-1] << "\n";
          return TCL_ERROR;
      }
      tracker.consume(Position::Beta);
    }
    //
    //
    else if (strcmp(argv[i], "-Ap") == 0) {
      if (++i >= argc) {
          opserr << "Missing value for option " << argv[i-1] << "\n";
          return TCL_ERROR;
      }
      if (Tcl_GetDouble(interp, argv[i], &Ap) != TCL_OK) {
          opserr << "Invalid value for option " << argv[i-1] << "\n";
          return TCL_ERROR;
      }
      tracker.consume(Position::Ap);
    }
    else if (strcmp(argv[i], "-An") == 0) {
      if (++i >= argc) {
          opserr << "Missing value for option " << argv[i-1] << "\n";
          return TCL_ERROR;
      }
      if (Tcl_GetDouble(interp, argv[i], &An) != TCL_OK) {
          opserr << "Invalid value for option " << argv[i-1] << "\n";
          return TCL_ERROR;
      }
      tracker.consume(Position::An);
    }
    else if (strcmp(argv[i], "-Bn") == 0) {
      if (++i >= argc) {
          opserr << "Missing value for option " << argv[i-1] << "\n";
          return TCL_ERROR;
      }
      if (Tcl_GetDouble(interp, argv[i], &Bn) != TCL_OK) {
          opserr << "Invalid value for option " << argv[i-1] << "\n";
          return TCL_ERROR;
      }
      tracker.consume(Position::Bn);
    }
    else
      positional.insert(i);
  }

  //
  // 2) Positional arguments
  //
  for (int i : positional) {
  
    if (tracker.current() == Position::EndRequired)
      tracker.increment();

    switch (tracker.current()) {
      // General
      case Position::Tag :
        if (Tcl_GetInt(interp, argv[i], &tag) != TCL_OK) {
            opserr << OpenSees::PromptParseError << "invalid section Elastic tag.\n";
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

      // Isotropy
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
        } else {
          tracker.increment();
          break;
        }

      case Position::Nu:
        if (Tcl_GetDouble (interp, argv[i], &consts.nu) != TCL_OK) {
            opserr << OpenSees::PromptParseError << "invalid nu.\n";
            return TCL_ERROR;
        } else {
          tracker.increment();
          break;
        }
      
      case Position::Lambda:
        if (Tcl_GetDouble (interp, argv[i], &consts.lambda) != TCL_OK) {
            opserr << OpenSees::PromptParseError << "invalid Lame lambda.\n";
            return TCL_ERROR;
        } else {
          tracker.increment();
          break;
        }
    
      // Yielding
      case Position::PeakCompression:
        if (Tcl_GetDouble(interp, argv[i], &Fc) != TCL_OK) {
            opserr << OpenSees::PromptParseError << "invalid Fc.\n";
            return TCL_ERROR;           
        } else {
          tracker.increment();
          break;
        }

      case Position::PeakTension:
        if (Tcl_GetDouble (interp, argv[i], &Ft) != TCL_OK) {
            opserr << OpenSees::PromptParseError << "invalid Ft.\n";
            return TCL_ERROR;
        } else {
          tracker.increment();
          break;
        }
      case Position::Beta:
        if (Tcl_GetDouble (interp, argv[i], &beta) != TCL_OK) {
            opserr << OpenSees::PromptParseError << "invalid beta.\n";
            return TCL_ERROR;
        } else {
          tracker.increment();
          break;
        }
      case Position::Ap:
        if (Tcl_GetDouble (interp, argv[i], &Ap) != TCL_OK) {
            opserr << OpenSees::PromptParseError << "invalid Ap.\n";
            return TCL_ERROR;
        } else {
          tracker.increment();
          break;
        }
      case Position::An:
        if (Tcl_GetDouble (interp, argv[i], &An) != TCL_OK) {
            opserr << OpenSees::PromptParseError << "invalid An.\n";
            return TCL_ERROR;
        } else {
          tracker.increment();
          break;
        }
      case Position::Bn:
        if (Tcl_GetDouble (interp, argv[i], &Bn) != TCL_OK) {
            opserr << OpenSees::PromptParseError << "invalid Bn.\n";
            return TCL_ERROR;
        } else {
          tracker.increment();
          break;
        }
    

      case Position::EndRequired:
        // This will not be reached
        break;

      case Position::End:
        opserr << OpenSees::PromptParseError << "unexpected argument " << argv[i] << ".\n";
        return TCL_ERROR;
    }
  }


  //
  // 3. Check for required arguments
  //
  if (tracker.current() < Position::EndRequired) {
    opserr << OpenSees::PromptParseError
            << "missing required arguments: ";
    while (tracker.current() != Position::EndRequired) {
      switch (tracker.current()) {
        case Position::Tag :
          opserr << "tag ";
          break;
        // Isotropy
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
        // Yielding
        case Position::PeakCompression:
          opserr << "Fc ";
          break;
        case Position::PeakTension:
          opserr << "Ft ";
          break;
        case Position::Ap:
          opserr << "Ap ";
          break;
        // Hardening
        case Position::An:
          opserr << "An ";
          break;
        case Position::Bn:
          opserr << "Bn ";
          break;
        case Position::Beta:
          opserr << "beta ";
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
  // Create the material (TODO)
  //
  BasicModelBuilder *builder = static_cast<BasicModelBuilder*>(clientData);
  if ((strcmp(argv[1], "PlasticDamageConcrete3d") == 0) ||
      (strcasecmp(argv[1], "PlasticDamageConcrete") == 0) ||
      (strcmp(argv[1], "FariaPlasticDamage") == 0)) {

    NDMaterial* theMaterial = new FariaPlasticDamage3d(tag, consts.E, consts.nu, Ft, Fc, 
                                                        beta, Ap, An, Bn, density);
    if (builder->addTaggedObject<NDMaterial>(*theMaterial) != TCL_OK ) {
      delete theMaterial;
      return TCL_ERROR;
    }
    return TCL_OK;
  }

  return TCL_ERROR;
}
