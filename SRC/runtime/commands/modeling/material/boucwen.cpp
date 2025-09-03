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
// Written: cmp
//
#include <tcl.h>
#include <Logging.h>
#include <Parsing.h>
#include <ArgumentTracker.h>
#include <BasicModelBuilder.h>
#include <string.h>

#ifdef _MSC_VER
#  define strcasecmp _stricmp
#else
#  include <strings.h>
#endif
#define strcmp strcasecmp

#include <BoucWen/BWBF.h>
#include <BoucWen/BWBN.h>
#include <BoucWen/BoucWenMaterial.h>
#include <BoucWen/BoucWenOriginal.h>

template <typename Positions>
static int
ParseBoucWen(ClientData clientData, Tcl_Interp *interp,
                  int argc, TCL_Char ** const argv)
{

  BasicModelBuilder *builder = static_cast<BasicModelBuilder *>(clientData);

  ArgumentTracker<Positions> tracker;
  std::set<int> positional;

  int tag;
  double E, Fy, alpha;
  double beta = 0.5, gamma=0.5;
  double Ao=1, 
         delta_a = 0,
         delta_n = 0, 
         delta_v = 0;
  
  double pinch_slope = 0.0;
  double pinch_slip  = 0.0;
  double pinch_start = 0.0;
  double pinch_rate  = 0.0;
  double pinch_size  = 0.0;
  double pinch_lamda = 0.5;

  int iterations = 25;
  double tolerance = 1e-8;
  double n = 1.0;

  // BoucWenOriginal-only
  double mu = 2.0, alphaNL = 0.0;

  //
  // Begin parse
  //

  if (Tcl_GetInt(interp, argv[2], &tag) != TCL_OK) {
    opserr << "WARNING invalid uniaxialMaterial tag\n";
    return TCL_ERROR;
  }

  for (int i=2; i<argc; i++) {
    if ((strcmp(argv[i], "-fy") == 0) || (strcmp(argv[i], "-Fy") == 0)) {
      if (++i >= argc) {
        opserr << "Missing value for option " << argv[i-1] << "\n";
        return TCL_ERROR;
      }
      if (Tcl_GetDouble(interp, argv[i], &Fy) != TCL_OK) {
        opserr << "Invalid value for option " << argv[i-1] << "\n";
        return TCL_ERROR;
      }
      tracker.consume(Positions::Fy);
    }
    else if ((strcasecmp(argv[i], "-E") == 0) || (strcmp(argv[i], "-ko") == 0)) {
      if (++i >= argc) {
        opserr << "Missing value for option " << argv[i-1] << "\n";
        return TCL_ERROR;
      }
      if (Tcl_GetDouble(interp, argv[i], &E) != TCL_OK) {
        opserr << "Invalid value for option " << argv[i-1] << "\n";
        return TCL_ERROR;
      }
      tracker.consume(Positions::E);
    }
    else if (strcasecmp(argv[i], "-alpha") == 0 || strcasecmp(argv[i], "-b") == 0) {
      if (++i >= argc) {
        opserr << "Missing value for option " << argv[i-1] << "\n";
        return TCL_ERROR;
      }
      if (Tcl_GetDouble(interp, argv[i], &alpha) != TCL_OK) {
        opserr << "Invalid value for option " << argv[i-1] << "\n";
        return TCL_ERROR;
      }
      tracker.consume(Positions::Alpha);
    }
    else if (strcasecmp(argv[i], "-n") == 0 || strcasecmp(argv[i], "-N") == 0) {
      if (++i >= argc) {
        opserr << "Missing value for option " << argv[i-1] << "\n";
        return TCL_ERROR;
      }
      if (Tcl_GetDouble(interp, argv[i], &n) != TCL_OK) {
        opserr << "Invalid value for option " << argv[i-1] << "\n";
        return TCL_ERROR;
      }
      tracker.consume(Positions::N);
    }

    else if (strcmp(argv[i], "-beta") == 0) {
      if (++i >= argc) {
        opserr << "Missing value for option " << argv[i-1] << "\n";
        return TCL_ERROR;
      }
      if (Tcl_GetDouble(interp, argv[i], &beta) != TCL_OK) {
        opserr << "Invalid value for option " << argv[i-1] << "\n";
        return TCL_ERROR;
      }
      tracker.consume(Positions::Beta);
    }
    else if (strcmp(argv[i], "-gamma") == 0) {
      if (++i >= argc) {
        opserr << "Missing value for option " << argv[i-1] << "\n";
        return TCL_ERROR;
      }
      if (Tcl_GetDouble(interp, argv[i], &gamma) != TCL_OK) {
        opserr << "Invalid value for option " << argv[i-1] << "\n";
        return TCL_ERROR;
      }
      tracker.consume(Positions::Gamma);
    }
    //
    // Damage
    //
    else if (strcmp(argv[i], "-Ao") == 0) {
      if (++i >= argc) {
        opserr << "Missing value for option " << argv[i-1] << "\n";
        return TCL_ERROR;
      }
      if (Tcl_GetDouble(interp, argv[i], &Ao) != TCL_OK) {
        opserr << "Invalid value for option " << argv[i-1] << "\n";
        return TCL_ERROR;
      }
      tracker.consume(Positions::Ao);
    }
    else if (strcmp(argv[i], "-delta_a") == 0) {
      if (++i >= argc) {
        opserr << "Missing value for option " << argv[i-1] << "\n";
        return TCL_ERROR;
      }
      if (Tcl_GetDouble(interp, argv[i], &delta_a) != TCL_OK) {
        opserr << "Invalid value for option " << argv[i-1] << "\n";
        return TCL_ERROR;
      }
      tracker.consume(Positions::DeltaA);
    }
    else if (strcmp(argv[i], "-delta_n") == 0) {
      if (++i >= argc) {
        opserr << "Missing value for option " << argv[i-1] << "\n";
        return TCL_ERROR;
      }
      if (Tcl_GetDouble(interp, argv[i], &delta_n) != TCL_OK) {
        opserr << "Invalid value for option " << argv[i-1] << "\n";
        return TCL_ERROR;
      }
      tracker.consume(Positions::DeltaN);
    }
    else if (strcmp(argv[i], "-delta_v") == 0) {
      if (++i >= argc) {
        opserr << "Missing value for option " << argv[i-1] << "\n";
        return TCL_ERROR;
      }
      if (Tcl_GetDouble(interp, argv[i], &delta_v) != TCL_OK) {
        opserr << "Invalid value for option " << argv[i-1] << "\n";
        return TCL_ERROR;
      }
      tracker.consume(Positions::DeltaV);
    }
    //
    // Pinching
    //
    else if (strcmp(argv[i], "-pinch_slope") == 0 || strcmp(argv[i], "-p") == 0) {
      if (++i >= argc) {
        opserr << "Missing value for option " << argv[i-1] << "\n";
        return TCL_ERROR;
      }
      if (Tcl_GetDouble(interp, argv[i], &pinch_slope) != TCL_OK) {
        opserr << "Invalid value for option " << argv[i-1] << "\n";
        return TCL_ERROR;
      }
      tracker.consume(Positions::PinchSlope);
    }
    else if (strcmp(argv[i], "-pinch_start") == 0) {
      if (++i >= argc) {
        opserr << "Missing value for option " << argv[i-1] << "\n";
        return TCL_ERROR;
      }
      if (Tcl_GetDouble(interp, argv[i], &pinch_start) != TCL_OK) {
        opserr << "Invalid value for option " << argv[i-1] << "\n";
        return TCL_ERROR;
      }
      tracker.consume(Positions::PinchStart);
    }
    else if (strcmp(argv[i], "-pinch_rate") == 0) {
      if (++i >= argc) {
        opserr << "Missing value for option " << argv[i-1] << "\n";
        return TCL_ERROR;
      }
      if (Tcl_GetDouble(interp, argv[i], &pinch_rate) != TCL_OK) {
        opserr << "Invalid value for option " << argv[i-1] << "\n";
        return TCL_ERROR;
      }
      tracker.consume(Positions::PinchRate);
    }
    else if (strcmp(argv[i], "-pinch_slip") == 0) {
      if (++i >= argc) {
        opserr << "Missing value for option " << argv[i-1] << "\n";
        return TCL_ERROR;
      }
      if (Tcl_GetDouble(interp, argv[i], &pinch_slip) != TCL_OK) {
        opserr << "Invalid value for option " << argv[i-1] << "\n";
        return TCL_ERROR;
      }
      tracker.consume(Positions::PinchSlip);
    }
    else if (strcmp(argv[i], "-pinch_size") == 0) {
      if (++i >= argc) {
        opserr << "Missing value for option " << argv[i-1] << "\n";
        return TCL_ERROR;
      }
      if (Tcl_GetDouble(interp, argv[i], &pinch_size) != TCL_OK) {
        opserr << "Invalid value for option " << argv[i-1] << "\n";
        return TCL_ERROR;
      }
      tracker.consume(Positions::PinchSize);
    }
    else if (strcmp(argv[i], "-pinch_lamda") == 0) {
      if (++i >= argc) {
        opserr << "Missing value for option " << argv[i-1] << "\n";
        return TCL_ERROR;
      }
      if (Tcl_GetDouble(interp, argv[i], &pinch_lamda) != TCL_OK) {
        opserr << "Invalid value for option " << argv[i-1] << "\n";
        return TCL_ERROR;
      }
      tracker.consume(Positions::PinchLamda);
    }
    //
    // Numerics
    //
    else if (strcmp(argv[i], "-tolerance") == 0) {
      if (++i >= argc) {
        opserr << "Missing value for option " << argv[i-1] << "\n";
        return TCL_ERROR;
      }
      if (Tcl_GetDouble(interp, argv[i], &tolerance) != TCL_OK) {
        opserr << "Invalid value for option " << argv[i-1] << "\n";
        return TCL_ERROR;
      }
      tracker.consume(Positions::Tolerance);
    }
    else if (strcmp(argv[i], "-iterations") == 0) {
      if (++i >= argc) {
        opserr << "Missing value for option " << argv[i-1] << "\n";
        return TCL_ERROR;
      }
      if (Tcl_GetInt(interp, argv[i], &iterations) != TCL_OK) {
        opserr << "Invalid value for option " << argv[i-1] << "\n";
        return TCL_ERROR;
      }
      tracker.consume(Positions::Iterations);
    }
    else 
      positional.insert(i);
  }


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
      
      case Positions::Fy:
        if (Tcl_GetDouble(interp, argv[i], &Fy) != TCL_OK) {
          opserr << "Invalid value for Fy " << argv[i] << "\n";
          return TCL_ERROR;
        }
        tracker.consume(Positions::Fy);
        break;
      case Positions::E:
        if (Tcl_GetDouble(interp, argv[i], &E) != TCL_OK) {
          opserr << "Invalid value for E " << argv[i] << "\n";
          return TCL_ERROR;
        }
        tracker.consume(Positions::E);
        break;
      case Positions::Alpha:
        if (Tcl_GetDouble(interp, argv[i], &alpha) != TCL_OK) {
          opserr << "Invalid value for Alpha " << argv[i] << "\n";
          return TCL_ERROR;
        }
        tracker.consume(Positions::Alpha);
        break;
      case Positions::N:
        if (Tcl_GetDouble(interp, argv[i], &n) != TCL_OK) {
          opserr << "Invalid value for N " << argv[i] << "\n";
          return TCL_ERROR;
        }
        tracker.consume(Positions::N);
        break;
      case Positions::Beta:
        if (Tcl_GetDouble(interp, argv[i], &beta) != TCL_OK) {
          opserr << "Invalid value for beta " << argv[i] << "\n";
          return TCL_ERROR;
        }
        tracker.consume(Positions::Beta);
        break;
      case Positions::Gamma:
        if (Tcl_GetDouble(interp, argv[i], &gamma) != TCL_OK) {
          opserr << "Invalid value for gamma " << argv[i] << "\n";
          return TCL_ERROR;
        }
        tracker.consume(Positions::Gamma);
        break;
      case Positions::Tolerance:
        if (Tcl_GetDouble(interp, argv[i], &tolerance) != TCL_OK) {
          opserr << "Invalid value for tolerance " << argv[i] << "\n";
          return TCL_ERROR;
        }
        tracker.consume(Positions::Tolerance);
        break;
      case Positions::Iterations:
        if (Tcl_GetInt(interp, argv[i], &iterations) != TCL_OK) {
          opserr << "Invalid value for iterations " << argv[i] << "\n";
          return TCL_ERROR;
        }
        tracker.consume(Positions::Iterations);
        break;
      //
      // Degradation
      //
      case Positions::Ao:
        if (Tcl_GetDouble(interp, argv[i], &Ao) != TCL_OK) {
          opserr << "Invalid value for Ao " << argv[i] << "\n";
          return TCL_ERROR;
        }
        tracker.consume(Positions::Ao);
        break;
      case Positions::DeltaA:
        if (Tcl_GetDouble(interp, argv[i], &delta_a) != TCL_OK) {
          opserr << "Invalid value for delta_a " << argv[i] << "\n";
          return TCL_ERROR;
        }
        tracker.consume(Positions::DeltaA);
        break;
      case Positions::DeltaN:
        if (Tcl_GetDouble(interp, argv[i], &delta_n) != TCL_OK) {
          opserr << "Invalid value for delta_n " << argv[i] << "\n";
          return TCL_ERROR;
        }
        tracker.consume(Positions::DeltaN);
        break;
      case Positions::DeltaV:
        if (Tcl_GetDouble(interp, argv[i], &delta_v) != TCL_OK) {
          opserr << "Invalid value for delta_v " << argv[i] << "\n";
          return TCL_ERROR;
        }
        tracker.consume(Positions::DeltaV);
        break;
      //
      // Schellenberg extras
      //
      case Positions::AlphaNL:
        if (Tcl_GetDouble(interp, argv[i], &alphaNL) != TCL_OK) {
          opserr << "Invalid value for alphaNL " << argv[i] << "\n";
          return TCL_ERROR;
        }
        tracker.consume(Positions::AlphaNL);
        break;
      case Positions::Mu:
        if (Tcl_GetDouble(interp, argv[i], &mu) != TCL_OK) {
          opserr << "Invalid value for mu " << argv[i] << "\n";
          return TCL_ERROR;
        }
        tracker.consume(Positions::Mu);
        break;
      
      case Positions::PinchRate:
      case Positions::PinchSlip:
      case Positions::PinchSize:
      case Positions::PinchSlope:
      case Positions::PinchStart:
      case Positions::PinchLamda:
        break;
      //
      //
      case Positions::EndRequired:
        // This will not be reached
        break;

      case Positions::End:
        opserr << "Invalid value for option " << argv[i] << "\n";
        return TCL_ERROR;
    }
  }

  if (tracker.current() < Positions::EndRequired) {
    opserr << "Missing required arguments: ";
    while (tracker.current() != Positions::End) {
      switch (tracker.current()) {
        case Positions::Fy:
          opserr << "Fy ";
          break;
        case Positions::E:
            opserr << "E ";
            break;
        case Positions::Alpha:
            opserr << "Alpha ";
            break;
        case Positions::N:
            opserr << "N ";
            break;
        case Positions::Beta:
            opserr << "Beta ";
            break;
        case Positions::Gamma:
            opserr << "Gamma ";
            break;
        case Positions::Ao:
            opserr << "Ao ";
            break;
        case Positions::DeltaA:
            opserr << "delta_a ";
            break;
        case Positions::DeltaV:
            opserr << "delta_v ";
            break;
        case Positions::DeltaN:
            opserr << "delta_n ";
            break;
        
        case Positions::Tolerance:
            opserr << "tolerance ";
            break;
        case Positions::Iterations:
            opserr << "iterations ";
            break;

        case Positions::PinchSlope:
            opserr << "PinchSlope ";
            break;
        case Positions::PinchStart:
            opserr << "PinchStart ";
            break;
        case Positions::PinchRate:
            opserr << "PinchRate ";
            break;
        case Positions::PinchSize:
            opserr << "PinchSize ";
            break;
        
        case Positions::Mu:
        case Positions::AlphaNL:
        case Positions::EndRequired:
        case Positions::End:
        default:
          break;
      }
    
      if (tracker.current() == Positions::EndRequired)
        break;

      tracker.consume(tracker.current());
    }
    opserr << "\n";
    return TCL_ERROR;
  }

  //
  //
  //

  //
  //
  //
  UniaxialMaterial *theMaterial = nullptr;
  if ((strcmp(argv[1], "BWBF") == 0) ||
      (strcmp(argv[1], "Bouc") == 0)) {

    theMaterial =
        new BWBF(tag, 
                 E, 
                 Fy, 
                 alpha, 
                 n, 
                 beta,
                 delta_a, delta_v, delta_n,
                 //
                 pinch_slope, 
                 pinch_slip, 
                 pinch_start, 
                 pinch_rate, 
                 pinch_size,
                 pinch_lamda,
                 //
                 tolerance, 
                 iterations);
  }
  else if ((strcmp(argv[1], "BoucWen") == 0)) {

    theMaterial =
        new BoucWenMaterial(tag, alpha, E, n, gamma, beta, Ao, 
                            delta_a, delta_v, delta_n, 
                            tolerance, iterations);
  }
  else if (strcmp(argv[1], "BoucWenOriginal") == 0) {
    theMaterial =
        new BoucWenOriginal(tag, E, Fy, alpha, alphaNL, mu, 
                            n, beta, gamma,
                            tolerance, iterations);
  }
  else if (strcmp(argv[1], "BWBN") == 0) {
    theMaterial =
        new BWBN(tag, alpha, E, n, gamma, beta, Ao, 
                 pinch_start, pinch_slip, pinch_slope, 
                 pinch_size, pinch_rate, pinch_lamda,
                 tolerance, iterations);
  }
  else {
    opserr << "WARNING BoucWen: invalid material type\n";
    return TCL_ERROR;
  }

  if (theMaterial == nullptr)
    return TCL_ERROR;

  return builder->addTaggedObject<UniaxialMaterial>(*theMaterial);
}


int
TclCommand_newBoucWen(ClientData clientData, Tcl_Interp *interp,
                          int argc, TCL_Char ** const argv)
{

  if (strcmp(argv[1], "BWBF") == 0 ||
      strcmp(argv[1], "Bouc") == 0) {

    enum class Positions: int {
      Tag,
        E, Fy, Alpha, N, 
      EndRequired,
        Beta,
        Tolerance, Iterations,
        Ao, DeltaA, DeltaV, DeltaN,
        PinchSlope, PinchStart, PinchRate, PinchSize, PinchSlip, PinchLamda,
      EndAllowed,
        Mu, AlphaNL,
        Gamma,
      End
    };

    return ParseBoucWen<Positions>(clientData, interp, argc, argv);
  }

  else if ((strcmp(argv[1], "BoucWen") == 0)) {
    
    enum class Positions: int {
      Tag,
        Alpha, E, N, Gamma, Beta, Ao, DeltaA, DeltaV, DeltaN,
      EndRequired,
      EndAllowed,
        Fy,
        Tolerance, Iterations,
        PinchSlope, PinchStart, PinchRate, PinchSize, PinchSlip, PinchLamda,
        Mu, AlphaNL,
      End
    };
    return ParseBoucWen<Positions>(clientData, interp, argc, argv);
  }

  else if ((strcmp(argv[1], "BWBN") == 0)) {
    
    enum class Positions: int {
      Tag,
        Alpha, E, N, Gamma, Beta, Ao,
        PinchStart, PinchSlip, PinchSlope, PinchSize, PinchRate, PinchLamda,
      EndRequired,
      EndAllowed,
        Fy,
        Tolerance, Iterations,
        DeltaA, DeltaV, DeltaN,
        Mu, AlphaNL,
      End
    };
    return ParseBoucWen<Positions>(clientData, interp, argc, argv);
  }

  else if (strcmp(argv[1], "BoucWenOriginal") == 0) {
    // A. Schellenberg
    enum class Positions: int {
      Tag,
        E, Fy, Alpha,
      EndRequired,
        AlphaNL,
        Mu,
        N, Beta, Gamma,
        Tolerance, Iterations,
      End,
        Ao, DeltaA, DeltaV, DeltaN,
        PinchSlope, PinchStart, PinchRate, PinchSize, PinchSlip, PinchLamda,
    };
    return ParseBoucWen<Positions>(clientData, interp, argc, argv);
  }

  return TCL_ERROR;
}
