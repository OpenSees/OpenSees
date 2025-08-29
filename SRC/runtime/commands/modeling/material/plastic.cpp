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
// Description: This file implements a unified parser for plasticity materials.
//
// Written: cmp
// April 2025
//
#include <tcl.h>
#include <set>
#include <string.h>
#include <Logging.h>
#include <Parsing.h>
#include <ArgumentTracker.h>
#include <BasicModelBuilder.h>

#include "isotropy.h"

#include <HardeningMaterial.h>
#include <SimplifiedJ2.h>
#include <PlaneStressSimplifiedJ2.h>
#include <J2Plasticity.h>
#include <J2PlasticityThermal.h>
#include <MultiaxialCyclicPlasticity.h>
#include <DruckerPrager.h>
#include <J2BeamFiber2d.h>
#include <J2BeamFiber3d.h>


template <typename Position>
static inline int
TclCommand_newPlasticParser(ClientData clientData, Tcl_Interp *interp,
                                  int argc, TCL_Char ** const argv)
{

  assert(clientData != nullptr);

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
  double Fy, Fsat = 0, Fo = 0;
  // Hardening
  double Hiso=0,
         Hkin=0;
  struct {
    double theta = 1.0;
    double Hsat  = 0;
    double Hmix  = 0;
  } hard{};
  bool mix = Position::Theta < Position::End;
  // Viscosity
  double eta=0;
  // Drucker-Prager
  double rho = 0, rho_bar = 0;
  double atm = 101.0;
  double delta2 = 0.0;

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
    // Yielding
    else if (strcmp(argv[i], "-Fy") == 0 || 
             strcmp(argv[i], "-fy") == 0 ||
             strcmp(argv[i], "-yield-stress") == 0) {
      if (++i >= argc) {
          opserr << "Missing value for option " << argv[i-1] << "\n";
          return TCL_ERROR;
      }
      if (Tcl_GetDouble(interp, argv[i], &Fy) != TCL_OK) {
          opserr << "Invalid yield stress value for option " << argv[i-1] << "\n";
          return TCL_ERROR;
      }
      tracker.consume(Position::YieldStress);
    }
    else if ((strcmp(argv[i], "-Fo") == 0) || 
             (strcmp(argv[i], "-Ko") == 0)) {
      if (++i >= argc) {
          opserr << "Missing value for option " << argv[i-1] << "\n";
          return TCL_ERROR;
      }
      if (Tcl_GetDouble(interp, argv[i], &Fo) != TCL_OK) {
          opserr << "Invalid initial saturation stress value for option " << argv[i-1] << "\n";
          return TCL_ERROR;
      }
      tracker.consume(Position::SatStress0);
    }
    else if ((strcmp(argv[i], "-Fs") == 0) || 
             (strcmp(argv[i], "-Fsat") == 0) ||
             (strcmp(argv[i], "-fsat") == 0)) {
      if (++i >= argc) {
          opserr << "Missing value for option " << argv[i-1] << "\n";
          return TCL_ERROR;
      }
      if (Tcl_GetDouble(interp, argv[i], &Fsat) != TCL_OK) {
          opserr << "Invalid saturation stress value for option " << argv[i-1] << "\n";
          return TCL_ERROR;
      }
      tracker.consume(Position::SatStress);
    }
    //
    // Hardening
    //
    else if (strcmp(argv[i], "-Hiso") == 0 || strcmp(argv[i], "-isotropic-hardening") == 0) {
      if (++i >= argc) {
          opserr << "Missing value for option " << argv[i-1] << "\n";
          return TCL_ERROR;
      }
      if (Tcl_GetDouble(interp, argv[i], &Hiso) != TCL_OK) {
          opserr << "Invalid value for option " << argv[i-1] << "\n";
          return TCL_ERROR;
      }
      mix = false;
      tracker.consume(Position::Hiso);
    }
    else if (strcmp(argv[i], "-Hkin") == 0 || strcmp(argv[i], "-kinematic-hardening") == 0) {
      if (++i >= argc) {
          opserr << "Missing value for option " << argv[i-1] << "\n";
          return TCL_ERROR;
      }
      if (Tcl_GetDouble(interp, argv[i], &Hkin) != TCL_OK) {
          opserr << "Invalid value for option " << argv[i-1] << "\n";
          return TCL_ERROR;
      }
      mix = false;
      tracker.consume(Position::Hkin);
    }
    else if (strcmp(argv[i], "-H") == 0 || strcmp(argv[i], "-Hmix") == 0) {
      if (++i >= argc) {
          opserr << "Missing value for option " << argv[i-1] << "\n";
          return TCL_ERROR;
      }
      if (Tcl_GetDouble(interp, argv[i], &hard.Hmix) != TCL_OK) {
          opserr << "Invalid value for option " << argv[i-1] << "\n";
          return TCL_ERROR;
      }
      mix = true;
      tracker.consume(Position::Hmix);
      tracker.consume(Position::Hiso);
      tracker.consume(Position::Hkin);
    }

    else if (strcmp(argv[i], "-theta") == 0 || strcmp(argv[i], "-mix") == 0) {
      if (++i >= argc) {
          opserr << "Missing value for option " << argv[i-1] << "\n";
          return TCL_ERROR;
      }
      if (Tcl_GetDouble(interp, argv[i], &hard.theta) != TCL_OK) {
          opserr << "Invalid value for option " << argv[i-1] << "\n";
          return TCL_ERROR;
      }
      tracker.consume(Position::Theta);
    }
    else if (strcmp(argv[i], "-Hsat") == 0 || 
             strcmp(argv[i], "-delta") == 0 || 
             strcmp(argv[i], "-delta1") == 0) {
      if (++i >= argc) {
          opserr << "Missing value for option " << argv[i-1] << "\n";
          return TCL_ERROR;
      }
      if (Tcl_GetDouble(interp, argv[i], &hard.Hsat) != TCL_OK) {
          opserr << "Invalid value for option " << argv[i-1] << "\n";
          return TCL_ERROR;
      }
      tracker.consume(Position::Hsat);
    }
    // Drucker-Prager
    else if (strcmp(argv[i], "-delta2") == 0 ||
             strcmp(argv[i], "-Hten") == 0) {
      if (++i >= argc) {
          opserr << "Missing value for option " << argv[i-1] << "\n";
          return TCL_ERROR;
      }
      if (Tcl_GetDouble(interp, argv[i], &delta2) != TCL_OK) {
          opserr << "Invalid value for option " << argv[i-1] << "\n";
          return TCL_ERROR;
      }
      tracker.consume(Position::Delta2);
    }
    else if (strcmp(argv[i], "-atm") == 0) {
      if (++i >= argc) {
          opserr << "Missing value for option " << argv[i-1] << "\n";
          return TCL_ERROR;
      }
      if (Tcl_GetDouble(interp, argv[i], &atm) != TCL_OK) {
          opserr << "Invalid value for option " << argv[i-1] << "\n";
          return TCL_ERROR;
      }
      tracker.consume(Position::Atm);
    }
    else if (strcmp(argv[i], "-Rvol") == 0) {
      if (++i >= argc) {
          opserr << "Missing value for option " << argv[i-1] << "\n";
          return TCL_ERROR;
      }
      if (Tcl_GetDouble(interp, argv[i], &rho) != TCL_OK) {
          opserr << "Invalid value for option " << argv[i-1] << "\n";
          return TCL_ERROR;
      }
      tracker.consume(Position::Rho);
    }
    else if (strcmp(argv[i], "-Rbar") == 0 ||
             strcmp(argv[i], "-rhoBar") == 0) {
      if (++i >= argc) {
          opserr << "Missing value for option " << argv[i-1] << "\n";
          return TCL_ERROR;
      }
      if (Tcl_GetDouble(interp, argv[i], &rho_bar) != TCL_OK) {
          opserr << "Invalid value for option " << argv[i-1] << "\n";
          return TCL_ERROR;
      }
      tracker.consume(Position::RhoBar);
    }
    // Viscosity
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
      case Position::YieldStress:
        if (Tcl_GetDouble(interp, argv[i], &Fy) != TCL_OK) {
            opserr << OpenSees::PromptParseError << "invalid yield stress.\n";
            return TCL_ERROR;           
        } else {
          tracker.increment();
          break;
        }
      case Position::SatStress:
        if (Tcl_GetDouble (interp, argv[i], &Fsat) != TCL_OK) {
            opserr << OpenSees::PromptParseError << "invalid saturation stress.\n";
            return TCL_ERROR;
        } else {
          tracker.increment();
          break;
        }
      
      case Position::SatStress0:
        if (Tcl_GetDouble (interp, argv[i], &Fo) != TCL_OK) {
            opserr << OpenSees::PromptParseError << "invalid initial saturation stress.\n";
            return TCL_ERROR;
        } else {
          tracker.increment();
          break;
        }
      // Hardening
      case Position::Hiso:
        if (Tcl_GetDouble (interp, argv[i], &Hiso) != TCL_OK) {
            opserr << OpenSees::PromptParseError << "invalid Hiso.\n";
            return TCL_ERROR;
        } else {
          tracker.consume(Position::Hiso);
          tracker.consume(Position::Theta);
          break;
        }
      case Position::Hkin:
        if (Tcl_GetDouble (interp, argv[i], &Hkin) != TCL_OK) {
            opserr << OpenSees::PromptParseError << "invalid Hkin.\n";
            return TCL_ERROR;
        } else {
          tracker.consume(Position::Hkin);
          tracker.consume(Position::Theta);
          break;
        }
      case Position::Hsat:
        if (Tcl_GetDouble (interp, argv[i], &hard.Hsat) != TCL_OK) {
            opserr << OpenSees::PromptParseError << "invalid Hsat.\n";
            return TCL_ERROR;
        } else {
          tracker.increment();
          break;
        }
      case Position::Hmix:
        if (Tcl_GetDouble (interp, argv[i], &hard.Hmix) != TCL_OK) {
            opserr << OpenSees::PromptParseError << "invalid Hmix.\n";
            return TCL_ERROR;
        } else {
          mix = true;
          tracker.consume(Position::Hmix);
          tracker.consume(Position::Hiso);
          tracker.consume(Position::Hkin);
          break;
        }
      case Position::Theta:
        if (Tcl_GetDouble (interp, argv[i], &hard.theta) != TCL_OK) {
            opserr << OpenSees::PromptParseError << "invalid hardening theta.\n";
            return TCL_ERROR;
        } else {
          tracker.consume(Position::Theta);
          break;
        }

      // Drucker 
      case Position::Delta2:
        if (Tcl_GetDouble (interp, argv[i], &delta2) != TCL_OK) {
            opserr << OpenSees::PromptParseError << "invalid delta2.\n";
            return TCL_ERROR;
        } else {
          tracker.increment();
          break;
        }
      case Position::Rho:
        if (Tcl_GetDouble (interp, argv[i], &rho) != TCL_OK) {
            opserr << OpenSees::PromptParseError << "invalid Rvol.\n";
            return TCL_ERROR;
        } else {
          tracker.increment();
          break;
        }
      case Position::Atm:
        if (Tcl_GetDouble (interp, argv[i], &atm) != TCL_OK) {
            opserr << OpenSees::PromptParseError << "invalid atm.\n";
            return TCL_ERROR;
        } else {
          tracker.increment();
          break;
        }
      case Position::RhoBar:
        if (Tcl_GetDouble (interp, argv[i], &rho_bar) != TCL_OK) {
            opserr << OpenSees::PromptParseError << "invalid Rbar.\n";
            return TCL_ERROR;
        } else {
          tracker.increment();
          break;
        }
      // Viscosity
      case Position::Eta:
        if (Tcl_GetDouble (interp, argv[i], &eta) != TCL_OK) {
            opserr << OpenSees::PromptParseError << "invalid eta.\n";
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

  if (mix) {
    Hiso =        hard.theta  * hard.Hmix;
    Hkin = (1.0 - hard.theta) * hard.Hmix;
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
        case Position::YieldStress:
          opserr << "Fy ";
          break;
        case Position::SatStress:
          opserr << "Fsat ";
          break;
        case Position::SatStress0:
          opserr << "Fo ";
          break;
        // Hardening
        case Position::Hiso:
          opserr << "Hiso ";
          break;
        case Position::Hkin:
          opserr << "Hkin ";
          break;
        case Position::Hmix:
          opserr << "Hmix ";
          break;
        case Position::Theta:
          opserr << "theta ";
          break;

        // Drucker-Prager
        case Position::Delta2:
          opserr << "Hten ";
          break;
        case Position::Rho:
          opserr << "Rvol ";
          break;
        case Position::RhoBar:
          opserr << "Rbar ";
          break;
        
        // Viscosity
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

  //
  // Create the material (TODO)
  //
  BasicModelBuilder *builder = static_cast<BasicModelBuilder*>(clientData);
  if ((strcmp(argv[1], "Hardening") == 0) ||
      (strcmp(argv[1], "Steel") == 0)) {

    UniaxialMaterial* theMaterial = new HardeningMaterial(tag, consts.E, Fy, Hiso, Hkin, eta);
    if (builder->addTaggedObject<UniaxialMaterial>(*theMaterial) != TCL_OK ) {
      delete theMaterial;
      return TCL_ERROR;
    }
    return TCL_OK;
  }

  else if (strcmp(argv[1], "J2BeamFiber") == 0) {
    NDMaterial* theMaterial = nullptr;
    if (builder->getNDM() == 2)
      theMaterial = new J2BeamFiber2d(tag, consts.E, consts.G, Fy, Hkin, Hiso);
    else 
      theMaterial = new J2BeamFiber3d(tag, consts.E, consts.G, Fy, Hkin, Hiso);

    if (builder->addTaggedObject<NDMaterial>(*theMaterial) != TCL_OK ) {
      delete theMaterial;
      return TCL_ERROR;
    }
    return TCL_OK;
  }

  else if ((strcmp(argv[1], "Simplified3DJ2") == 0) ||
           (strcmp(argv[1], "SimplifiedJ2") == 0) ||
           (strcmp(argv[1], "J2Simplified") == 0) ||
           (strcmp(argv[1], "J2L") == 0) ||
           (strcmp(argv[1], "PlaneStressSimplifiedJ2") == 0) ||
           (strcmp(argv[1], "3DJ2") == 0)) {

    NDMaterial* theMaterial = new SimplifiedJ2(tag, 3, consts.G, consts.K, Fy, Hkin, Hiso, density);
    if (strcmp(argv[1], "PlaneStressSimplifiedJ2") == 0) {
      theMaterial = new PlaneStressSimplifiedJ2(tag, 2, *theMaterial);
    }
    if (builder->addTaggedObject<NDMaterial>(*theMaterial) != TCL_OK ) {
      delete theMaterial;
      return TCL_ERROR;
    }
    return TCL_OK;
  }

  else if ((strcmp(argv[1], "J2") == 0) ||
           (strcmp(argv[1], "J2N") == 0) ||
           (strcmp(argv[1], "J2Plasticity")  == 0)) {

    NDMaterial* theMaterial = new J2Plasticity(tag, 0, consts.K, consts.G, 
                                               Fy, Fsat, hard.Hsat, Hiso, eta, density);
    if (builder->addTaggedObject<NDMaterial>(*theMaterial) != TCL_OK ) {
      delete theMaterial;
      return TCL_ERROR;
    }
    return TCL_OK;
  }

  else if ((strcmp(argv[1], "J2PlasticityThermal") == 0) ||
           (strcmp(argv[1], "J2Thermal") == 0)) {
    NDMaterial* theMaterial = new J2PlasticityThermal(tag, 0, consts.K, consts.G, 
              Fy, Fsat, hard.Hsat, Hiso, eta, density);
    if (builder->addTaggedObject<NDMaterial>(*theMaterial) != TCL_OK ) {
    delete theMaterial;
    return TCL_ERROR;
    }
    return TCL_OK;
  }
  else if (strcmp(argv[1], "DruckerPrager") == 0 ||
           strcmp(argv[1], "DP") == 0) {

    NDMaterial* theMaterial = new DruckerPrager(tag, 0, consts.K, consts.G,
                  Fy, rho, rho_bar, Fsat, Fo,
                  hard.Hsat, delta2, hard.Hmix, hard.theta, density, atm);
    if (builder->addTaggedObject<NDMaterial>(*theMaterial) != TCL_OK ) {
      delete theMaterial;
      return TCL_ERROR;
    }
    return TCL_OK;
  }
  return TCL_ERROR;
}


int
TclCommand_newPlasticMaterial(ClientData clientData, Tcl_Interp *interp,
                              int argc, TCL_Char ** const argv)
{
  // 
  if (strcmp(argv[1], "Simplified3DJ2") == 0 ||
      strcmp(argv[1], "SimplifiedJ2") == 0 ||
      strcmp(argv[1], "J2Simplified") == 0 ||
      strcmp(argv[1], "J2L") == 0 ||
      strcmp(argv[1], "3DJ2") == 0 ||
      strcmp(argv[1], "PlaneStressSimplifiedJ2") == 0) {

    // "SimplifiedJ2"  tag?  G?  K?  Fy? Hkin?  Hiso?
    enum class Position : int {
      Tag, G, K, YieldStress, EndRequired, 
      Hkin, Hiso,
      End,
      E, Nu, Lambda, Eta, Theta, Hmix, Hsat,
      SatStress, SatStress0,
      Delta2, Rho, RhoBar, Atm,
      Density
    };
    return TclCommand_newPlasticParser<Position>(clientData, interp, argc, argv);
  }
  else if (strcmp(argv[1], "J2BeamFiber") == 0) {
    // J2BeamFiber $tag $E $v $sigmaY $Hiso $Hkin <$rho>
    enum class Position : int {
      Tag, E, G, YieldStress, EndRequired,
      Hkin, Hiso,
      Density,
      End,
      Nu, K, Eta, Lambda, Theta, Hmix, Hsat,
      SatStress, SatStress0,
      Delta2, Rho, RhoBar, Atm
    };
    return TclCommand_newPlasticParser<Position>(clientData, interp, argc, argv);
  }

  // "UniaxialJ2Plasticity" tag? E? sigmaY? Hkin? <Hiso?>
  else if (strcmp(argv[1], "UniaxialJ2Plasticity") == 0) {
  }

  else if (strcmp(argv[1], "HardeningMaterial") == 0 ||
           strcmp(argv[1], "Hardening")  == 0 ||
           strcmp(argv[1], "Hardening2") == 0 ||
           strcmp(argv[1], "Steel") == 0) {

    // "Hardening"  tag?  E?  Y?  Hiso?  Hkin?
    enum class Position : int {
      Tag, E, YieldStress, Hiso, EndRequired, 
      Hkin,
      End,
      // Keyword-only arguments
      Density,
      // Unused
      Eta, G, K, Nu, Lambda, Theta, Hmix, Hsat,
      SatStress, SatStress0, 
      Delta2, Rho, RhoBar, Atm,
    };
    return TclCommand_newPlasticParser<Position>(clientData, interp, argc, argv);
  }

  else if (strcmp(argv[1], "J2") == 0 ||
           strcmp(argv[1], "J2N") == 0 ||
           strcmp(argv[1], "J2Plasticity")  == 0) {

    // "J2Plasticity" tag? K? G? sig0? sigInf? delta? Hiso? <eta?>
    enum class Position : int {
      Tag, K, G, YieldStress, SatStress, Hsat, Hiso, EndRequired, 
      Eta,                                           End,
      E, Nu, Lambda, Hkin, Theta, Hmix, SatStress0, 
      Delta2, Rho, RhoBar, Atm,
      Density
    };
    return TclCommand_newPlasticParser<Position>(clientData, interp, argc, argv);
  }
  else if (strcmp(argv[1], "DP") == 0 ||
           strcmp(argv[1], "DruckerPrager")  == 0) {

    // DruckerPrager tag? K? G? sigma_y? rho? rho_bar? Kinf? Ko? delta1? delta2? H? theta? <massDensity? atm?>
    enum class Position : int {
      Tag, K, G, YieldStress, Rho, RhoBar, 
        SatStress, SatStress0, Hsat, Delta2, Hmix, Theta, EndRequired, 
      Density, Atm, End,
      Eta, E, Nu, Lambda, Hiso, Hkin
    };
    return TclCommand_newPlasticParser<Position>(clientData, interp, argc, argv);
  }
  return TCL_ERROR;
}


#include <UniaxialJ2Plasticity.h>
int
TclCommand_newUniaxialJ2Plasticity(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char** const argv)
{
    // ----- 1D J2 Plasticity ----
    if (argc < 7) {
      opserr << "WARNING invalid number of arguments\n";
      opserr << "Want: uniaxialMaterial UniaxialJ2Plasticity tag? E? sigmaY? Hkin? <Hiso?>"
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
      return TCL_ERROR;
    }

    if (Tcl_GetDouble(interp, argv[4], &sigmaY) != TCL_OK) {
      opserr << "WARNING invalid sigmaY\n";
      return TCL_ERROR;
    }

    if (Tcl_GetDouble(interp, argv[5], &Hkin) != TCL_OK) {
      opserr << "WARNING invalid Hkin\n";
      return TCL_ERROR;
    }

    if (argc >= 7)
      if (Tcl_GetDouble(interp, argv[6], &Hiso) != TCL_OK) {
        opserr << "WARNING invalid Hiso\n";
        return TCL_ERROR;
      }

    // Parsing was successful, allocate the material
    UniaxialMaterial* theMaterial = new UniaxialJ2Plasticity(tag, E, sigmaY, Hkin, Hiso);

   assert(clientData != nullptr);
   BasicModelBuilder *builder = static_cast<BasicModelBuilder*>(clientData);
   builder->addTaggedObject<UniaxialMaterial>(*theMaterial);
   return TCL_OK;
}


int
TclCommand_newJ2Simplified(ClientData clientData, Tcl_Interp* interp, int argc, const char** const argv)
{

  BasicModelBuilder* builder = static_cast<BasicModelBuilder*>(clientData);

  if (argc < 8) {
    opserr << "WARNING insufficient arguments\n";
    opserr << "Want: nDmaterial Simplified3DJ2  $matTag  $G  $K  $sig0 $Hkin  $Hiso"
            << "\n";
    return TCL_ERROR;
  }

  int tag;
  double K, G, sig0, H_kin, H_iso;

  if (Tcl_GetInt(interp, argv[2], &tag) != TCL_OK) {
    opserr << "WARNING invalid tag" << "\n";
    return TCL_ERROR;
  }

  if (Tcl_GetDouble(interp, argv[3], &G) != TCL_OK) {
    opserr << "WARNING invalid G\n";
    return TCL_ERROR;
  }

  if (Tcl_GetDouble(interp, argv[4], &K) != TCL_OK) {
    opserr << "WARNING invalid K\n";
    return TCL_ERROR;
  }

  if (Tcl_GetDouble(interp, argv[5], &sig0) != TCL_OK) {
    opserr << "WARNING invalid sig0\n";
    return TCL_ERROR;
  }

  if (Tcl_GetDouble(interp, argv[6], &H_kin) != TCL_OK) {
    opserr << "WARNING invalid Hkin\n";
    return TCL_ERROR;
  }

  if (Tcl_GetDouble(interp, argv[7], &H_iso) != TCL_OK) {
    opserr << "WARNING invalid Hiso\n";
    return TCL_ERROR;
  }

  NDMaterial* theMaterial = nullptr;
  
  if ((strcmp(argv[1], "Simplified3DJ2") == 0) ||
      (strcmp(argv[1], "3DJ2") == 0) ||
      (strcmp(argv[1], "SimplifiedJ2") == 0) ||
      (strcmp(argv[1], "PlaneStressSimplifiedJ2") == 0)) {
    double density = 0.0;
    theMaterial = new SimplifiedJ2(tag, 3, G, K, sig0, H_kin, H_iso, density);

    if (strcmp(argv[1], "PlaneStressSimplifiedJ2") == 0) {
      theMaterial = new PlaneStressSimplifiedJ2(tag, 2, *theMaterial);
    }
  }

  //
  if (builder->addTaggedObject<NDMaterial>(*theMaterial) != TCL_OK ) {
    delete theMaterial;
    return TCL_ERROR;
  }
  return TCL_OK;
}


#if 0

int
TclCommand_newJ2Material(ClientData clientData,
                         Tcl_Interp* interp,
                         int argc,
                         const char** const argv)
{
  BasicModelBuilder* builder = static_cast<BasicModelBuilder*>(clientData);

  if (argc < 9) {
    opserr << "WARNING insufficient arguments\n";
    opserr << "Want: nDMaterial J2Plasticity tag? K? G? sig0? sigInf? delta? H? <eta?>" << "\n";
    return TCL_ERROR;
  }

  int tag;
  double K, G, sig0, sigInf, delta, H;
  double eta = 0.0;

  if (Tcl_GetInt(interp, argv[2], &tag) != TCL_OK) {
    opserr << "WARNING invalid J2Plasticity tag" << "\n";
    return TCL_ERROR;
  }

  if (Tcl_GetDouble(interp, argv[3], &K) != TCL_OK) {
    opserr << "WARNING invalid K\n";
    return TCL_ERROR;
  }

  if (Tcl_GetDouble(interp, argv[4], &G) != TCL_OK) {
    opserr << "WARNING invalid G\n";
    return TCL_ERROR;
  }

  if (Tcl_GetDouble(interp, argv[5], &sig0) != TCL_OK) {
    opserr << "WARNING invalid sig0\n";
    return TCL_ERROR;
  }

  if (Tcl_GetDouble(interp, argv[6], &sigInf) != TCL_OK) {
    opserr << "WARNING invalid sigInf\n";
    return TCL_ERROR;
  }

  if (Tcl_GetDouble(interp, argv[7], &delta) != TCL_OK) {
    opserr << "WARNING invalid delta\n";
    return TCL_ERROR;
  }
  if (Tcl_GetDouble(interp, argv[8], &H) != TCL_OK) {
    opserr << "WARNING invalid H\n";
    return TCL_ERROR;
  }
  if (argc > 9 && Tcl_GetDouble(interp, argv[9], &eta) != TCL_OK) {
    opserr << "WARNING invalid eta\n";
    return TCL_ERROR;
  }


  //
  NDMaterial* theMaterial = nullptr;

  if ((strcmp(argv[1], "J2Plasticity") == 0) || 
      (strcmp(argv[1], "J2") == 0)) {
    theMaterial = new J2Plasticity(tag, 0, K, G, sig0, sigInf, delta, H, eta);
  }
  
  if (theMaterial == nullptr)
    return TCL_ERROR;

  if (builder->addTaggedObject<NDMaterial>(*theMaterial) != TCL_OK ) {
    delete theMaterial;
    return TCL_ERROR;
  }
  return TCL_OK;
}

int
TclCommand_newPlasticMaterial(ClientData clientData, Tcl_Interp* interp, int argc, const char**const argv)
{
  if ((strcmp(argv[1], "J2Plasticity") == 0) || (strcmp(argv[1], "J2") == 0))
  {
    if (argc < 9) {
      opserr << "WARNING insufficient arguments\n";
      opserr << "Want: nDMaterial J2Plasticity tag? K? G? sig0? sigInf? delta? H? <eta?>" << "\n";
      return TCL_ERROR;
    }

    int tag;
    double K, G, sig0, sigInf, delta, H;
    double eta = 0.0;

    if (Tcl_GetInt(interp, argv[2], &tag) != TCL_OK) {
      opserr << "WARNING invalid J2Plasticity tag" << "\n";
      return TCL_ERROR;
    }

    if (Tcl_GetDouble(interp, argv[3], &K) != TCL_OK) {
      opserr << "WARNING invalid K\n";
      return TCL_ERROR;
    }

    if (Tcl_GetDouble(interp, argv[4], &G) != TCL_OK) {
      opserr << "WARNING invalid G\n";
      return TCL_ERROR;
    }

    if (Tcl_GetDouble(interp, argv[5], &sig0) != TCL_OK) {
      opserr << "WARNING invalid sig0\n";
      return TCL_ERROR;
    }

    if (Tcl_GetDouble(interp, argv[6], &sigInf) != TCL_OK) {
      opserr << "WARNING invalid sigInf\n";
      return TCL_ERROR;
    }

    if (Tcl_GetDouble(interp, argv[7], &delta) != TCL_OK) {
      opserr << "WARNING invalid delta\n";
      return TCL_ERROR;
    }
    if (Tcl_GetDouble(interp, argv[8], &H) != TCL_OK) {
      opserr << "WARNING invalid H\n";
      return TCL_ERROR;
    }
    if (argc > 9 && Tcl_GetDouble(interp, argv[9], &eta) != TCL_OK) {
      opserr << "WARNING invalid eta\n";
      return TCL_ERROR;
    }

    theMaterial = new J2Plasticity(tag, 0, K, G, sig0, sigInf, delta, H, eta);
  }
}
#endif
