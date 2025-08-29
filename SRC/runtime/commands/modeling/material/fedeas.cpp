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
//
#include <tcl.h>
#include <Logging.h>
#include <Parsing.h>
#include <BasicModelBuilder.h>
#include <ArgumentTracker.h>
#include <string.h>

#ifdef _MSC_VER
#  define strcasecmp _stricmp
#else
#  include <strings.h>
#endif
#define strcmp strcasecmp

#include <Steel01.h>
#include <Steel01Thermal.h>
#include <Steel2.h>
#include <Steel02.h>
#include <Steel02Thermal.h>
#include <Concrete01.h>
#include <Concrete02.h>


template <typename Positions>
static int
FedeasConcrParse(ClientData clientData, Tcl_Interp *interp,
                  int argc, TCL_Char ** const argv)
{

  BasicModelBuilder *builder = static_cast<BasicModelBuilder *>(clientData);

  ArgumentTracker<Positions> tracker;
  std::set<int> positional;


  int tag;
  double fpc, epsc0, fpcu, epscu;
  double rat=0.1, ft=0, Ets=0;

  if (Tcl_GetInt(interp, argv[2], &tag) != TCL_OK) {
    opserr << "WARNING invalid uniaxialMaterial tag\n";
    return TCL_ERROR;
  }

  for (int i=2; i<argc; i++) {
    if ((strcmp(argv[i], "-fpc") == 0) || (strcmp(argv[i], "-Fc") == 0)) {
      if (++i >= argc) {
        opserr << "Missing value for option " << argv[i-1] << "\n";
        return TCL_ERROR;
      }
      if (Tcl_GetDouble(interp, argv[i], &fpc) != TCL_OK) {
        opserr << "Invalid value for option " << argv[i-1] << "\n";
        return TCL_ERROR;
      }
      if (tracker.contains(Positions::ft))
        ft = 0.1*fpc;
      tracker.consume(Positions::fpc);
    }
    else if ((strcasecmp(argv[i], "-epsc0") == 0) || (strcmp(argv[i], "-ec0") == 0)) {
      if (++i >= argc) {
        opserr << "Missing value for option " << argv[i-1] << "\n";
        return TCL_ERROR;
      }
      if (Tcl_GetDouble(interp, argv[i], &epsc0) != TCL_OK) {
        opserr << "Invalid value for option " << argv[i-1] << "\n";
        return TCL_ERROR;
      }
      tracker.consume(Positions::epsc0);
    }
    else if (strcasecmp(argv[i], "-fpcu") == 0 || strcasecmp(argv[i], "-Fcu") == 0) {
      if (++i >= argc) {
        opserr << "Missing value for option " << argv[i-1] << "\n";
        return TCL_ERROR;
      }
      if (Tcl_GetDouble(interp, argv[i], &fpcu) != TCL_OK) {
        opserr << "Invalid value for option " << argv[i-1] << "\n";
        return TCL_ERROR;
      }
      tracker.consume(Positions::fpcu);
    }
    else if (strcasecmp(argv[i], "-epscu") == 0 || strcmp(argv[i], "-ecu") == 0) {
      if (++i >= argc) {
        opserr << "Missing value for option " << argv[i-1] << "\n";
        return TCL_ERROR;
      }
      if (Tcl_GetDouble(interp, argv[i], &epscu) != TCL_OK) {
        opserr << "Invalid value for option " << argv[i-1] << "\n";
        return TCL_ERROR;
      }
      tracker.consume(Positions::epscu);
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
      
      case Positions::fpc:
        if (Tcl_GetDouble(interp, argv[i], &fpc) != TCL_OK) {
          opserr << "Invalid value for Fc " << argv[i] << "\n";
          return TCL_ERROR;
        }
        tracker.consume(Positions::fpc);
        break;
      case Positions::epsc0:
        if (Tcl_GetDouble(interp, argv[i], &epsc0) != TCL_OK) {
          opserr << "Invalid value for ec0 " << argv[i] << "\n";
          return TCL_ERROR;
        }
        tracker.consume(Positions::epsc0);
        break;
      case Positions::fpcu:
        if (Tcl_GetDouble(interp, argv[i], &fpcu) != TCL_OK) {
          opserr << "Invalid value for Fcu " << argv[i] << "\n";
          return TCL_ERROR;
        }
        tracker.consume(Positions::fpcu);
        break;
      case Positions::epscu:
        if (Tcl_GetDouble(interp, argv[i], &epscu) != TCL_OK) {
          opserr << "Invalid value for ecu " << argv[i] << "\n";
          return TCL_ERROR;
        }
        tracker.consume(Positions::epscu);
        break;
      case Positions::rat:
        if (Tcl_GetDouble(interp, argv[i], &rat) != TCL_OK) {
          opserr << "Invalid value for option " << argv[i] << "\n";
          return TCL_ERROR;
        }
        tracker.consume(Positions::rat);
        break;
      case Positions::ft:
        if (Tcl_GetDouble(interp, argv[i], &ft) != TCL_OK) {
          opserr << "Invalid value for option " << argv[i] << "\n";
          return TCL_ERROR;
        }
        tracker.consume(Positions::ft);
        break;
      case Positions::Ets:
        if (Tcl_GetDouble(interp, argv[i], &Ets) != TCL_OK) {
          opserr << "Invalid value for option " << argv[i] << "\n";
          return TCL_ERROR;
        }
        tracker.consume(Positions::Ets);
        break;

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
        case Positions::fpc:
          opserr << "fpc ";
          break;
        case Positions::epsc0:
          opserr << "epsc0 ";
          break;
        case Positions::fpcu:
          opserr << "Fcu ";
          break;
        case Positions::epscu:
          opserr << "epscu ";
          break;
        case Positions::rat:
          opserr << "rat ";
          break;
        case Positions::ft:
          opserr << "ft ";
          break;
        case Positions::Ets:
          opserr << "Ets ";
          break;
        case Positions::EndRequired:
        case Positions::End:
        default:
          break;
      }
    
      if (tracker.current() == Positions::End)
        break;

      tracker.consume(tracker.current());
    }
    opserr << "\n";
    return TCL_ERROR;
  }

  //
  //
  //
  if (fpcu > 0.0) {
    fpcu *= -1;
    // opswrn << OpenSees::SignalWarning <<  "Fcu should be negative\n";
  }

  //
  //
  //
  UniaxialMaterial *theMaterial = nullptr;
  if (strcmp(argv[1], "Concrete1") == 0 ||
      strcasecmp(argv[1], "Concrete01") == 0) {

    theMaterial = new Concrete01(tag, fpc, epsc0, fpcu, epscu);
  }

  else if ((strcmp(argv[1], "concr2") == 0) ||
           (strcmp(argv[1], "Concrete02") == 0)) {

    theMaterial =
        new Concrete02(tag, fpc, epsc0, fpcu, epscu, rat, ft, Ets);
  }

  if (theMaterial == nullptr)
    return TCL_ERROR;

  return builder->addTaggedObject<UniaxialMaterial>(*theMaterial);
}


template <typename Positions>
static int
FedeasSteelParse(ClientData clientData, Tcl_Interp *interp,
                                  int argc, TCL_Char ** const argv)
{
  BasicModelBuilder *builder = static_cast<BasicModelBuilder *>(clientData);

  ArgumentTracker<Positions> tracker;
  std::set<int> positional;

  UniaxialMaterial *theMaterial = nullptr;

  int tag;
  double fy, E, b;

  double   a1 = 0.0,
           a2 = 1.0,
           a3 = 0.0,
           a4 = 1.0;
  double   R0 = 15.0,
          cR1 = 0.925,
          cR2 = 0.15;

  for (int i=2; i<argc; i++) {
    if ((strcmp(argv[i], "-fy") == 0) ||
        strcmp(argv[i], "-Fy") == 0) {
      if (++i >= argc) {
        opserr << "Missing value for option " << argv[i-1] << "\n";
        return TCL_ERROR;
      }
      if (Tcl_GetDouble(interp, argv[i], &fy) != TCL_OK) {
        opserr << "Invalid value for option " << argv[i-1] << "\n";
        return TCL_ERROR;
      }
      tracker.consume(Positions::fy);
    }
    else if (strcmp(argv[i], "-E") == 0) {
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
    else if (strcmp(argv[i], "-b") == 0) {
      if (++i >= argc) {
        opserr << "Missing value for option " << argv[i-1] << "\n";
        return TCL_ERROR;
      }
      if (Tcl_GetDouble(interp, argv[i], &b) != TCL_OK) {
        opserr << "Invalid value for option " << argv[i-1] << "\n";
        return TCL_ERROR;
      }
      tracker.consume(Positions::b);
    }
    else if (strcmp(argv[i], "-R0") == 0) {
      if (++i >= argc) {
        opserr << "Missing value for option " << argv[i-1] << "\n";
        return TCL_ERROR;
      }
      if (Tcl_GetDouble(interp, argv[i], &R0) != TCL_OK) {
        opserr << "Invalid value for option " << argv[i-1] << "\n";
        return TCL_ERROR;
      }
      tracker.consume(Positions::R0);
    }
    else if (strcmp(argv[i], "-cR1") == 0) {
      if (++i >= argc) {
        opserr << "Missing value for option " << argv[i-1] << "\n";
        return TCL_ERROR;
      }
      if (Tcl_GetDouble(interp, argv[i], &cR1) != TCL_OK) {
        opserr << "Invalid value for option " << argv[i-1] << "\n";
        return TCL_ERROR;
      }
      tracker.consume(Positions::cR1);
    }
    else if (strcmp(argv[i], "-cR2") == 0) {
      if (++i >= argc) {
        opserr << "Missing value for option " << argv[i-1] << "\n";
        return TCL_ERROR;
      }
      if (Tcl_GetDouble(interp, argv[i], &cR2) != TCL_OK) {
        opserr << "Invalid value for option " << argv[i-1] << "\n";
        return TCL_ERROR;
      }
      tracker.consume(Positions::cR2);
    }
    else if (strcmp(argv[i], "-a1") == 0) {
      if (++i >= argc) {
        opserr << "Missing value for option " << argv[i-1] << "\n";
        return TCL_ERROR;
      }
      if (Tcl_GetDouble(interp, argv[i], &a1) != TCL_OK) {
        opserr << "Invalid value for option " << argv[i-1] << "\n";
        return TCL_ERROR;
      }
      tracker.consume(Positions::a1);
    }
    else if (strcmp(argv[i], "-a2") == 0) {
      if (++i >= argc) {
        opserr << "Missing value for option " << argv[i-1] << "\n";
        return TCL_ERROR;
      }
      if (Tcl_GetDouble(interp, argv[i], &a2) != TCL_OK) {
        opserr << "Invalid value for option " << argv[i-1] << "\n";
        return TCL_ERROR;
      }
      tracker.consume(Positions::a2);
    }
    else if (strcmp(argv[i], "-a3") == 0) {
      if (++i >= argc) {
        opserr << "Missing value for option " << argv[i-1] << "\n";
        return TCL_ERROR;
      }
      if (Tcl_GetDouble(interp, argv[i], &a3) != TCL_OK) {
        opserr << "Invalid value for option " << argv[i-1] << "\n";
        return TCL_ERROR;
      }
      tracker.consume(Positions::a3);
    }
    else if (strcmp(argv[i], "-a4") == 0) {
      if (++i >= argc) {
        opserr << "Missing value for option " << argv[i-1] << "\n";
        return TCL_ERROR;
      }
      if (Tcl_GetDouble(interp, argv[i], &a4) != TCL_OK) {
        opserr << "Invalid value for option " << argv[i-1] << "\n";
        return TCL_ERROR;
      }
      tracker.consume(Positions::a4);
    }
    else
      positional.insert(i);
  }

  //
  // Positional arguments
  //
  for (int i : positional) {

    if (tracker.current() == Positions::EndRequired)
      tracker.increment();

    switch (tracker.current()) {
      case Positions::Tag:
        if (Tcl_GetInt(interp, argv[i], &tag) != TCL_OK) {
          opserr << "invalid tag.\n";
          return TCL_ERROR;
        } else {
          tracker.increment();
          break;
        }
      case Positions::fy :
        if (Tcl_GetDouble(interp, argv[i], &fy) != TCL_OK) {
          opserr << "invalid Fy.\n";
          return TCL_ERROR;
        } else {
          tracker.increment();
          break;
        }
      case Positions::E:
        if (Tcl_GetDouble(interp, argv[i], &E) != TCL_OK) {
          opserr << "invalid E.\n";
          return TCL_ERROR;
        } else {
          tracker.increment();
          break;
        }
      case Positions::b:
        if (Tcl_GetDouble(interp, argv[i], &b) != TCL_OK) {
          opserr << "invalid b.\n";
          return TCL_ERROR;
        } else {
          tracker.increment();
          break;
        }
      case Positions::R0:
        if (Tcl_GetDouble(interp, argv[i], &R0) != TCL_OK) {
          opserr << "invalid R0.\n";
          return TCL_ERROR;
        } else {
          tracker.increment();
          break;
        }
      case Positions::cR1:
        if (Tcl_GetDouble(interp, argv[i], &cR1) != TCL_OK) {
          opserr << "invalid cR1.\n";
          return TCL_ERROR;
        } else {
          tracker.increment();
          break;
        }
      case Positions::cR2:
        if (Tcl_GetDouble(interp, argv[i], &cR2) != TCL_OK) {
          opserr << "invalid cR2.\n";
          return TCL_ERROR;
        } else {
          tracker.increment();
          break;
        }
      case Positions::a1:
        if (Tcl_GetDouble(interp, argv[i], &a1) != TCL_OK) {
          opserr << "invalid a1.\n";
          return TCL_ERROR;
        } else {
          tracker.increment();
          break;
        }
      case Positions::a2:
        if (Tcl_GetDouble(interp, argv[i], &a2) != TCL_OK) {
          opserr << "invalid a2.\n";
          return TCL_ERROR;
        } else {
          tracker.increment();
          break;
        }
      case Positions::a3:
        if (Tcl_GetDouble(interp, argv[i], &a3) != TCL_OK) {
          opserr << "invalid a3.\n";
          return TCL_ERROR;
        } else {
          tracker.increment();
          break;
        }
      case Positions::a4:
        if (Tcl_GetDouble(interp, argv[i], &a4) != TCL_OK) {
          opserr << "invalid a4.\n";
          return TCL_ERROR;
        } else {
          tracker.increment();
          break;
        }
      case Positions::sig0:
        if (Tcl_GetDouble(interp, argv[i], &a4) != TCL_OK) {
          opserr << "invalid sig0.\n";
          return TCL_ERROR;
        } else {
          tracker.increment();
          break;
        }

      case Positions::EndRequired:
        // This will not be reached
        break;

      case Positions::End:
        opserr << "unexpected argument " << argv[i] << ".\n";
        return TCL_ERROR;
    }
  }

  // Check all required arguments are present
  if (tracker.current() < Positions::EndRequired) {
    opserr << "missing required arguments: ";
    while (tracker.current() != Positions::EndRequired) {
      switch (tracker.current()) {
        case Positions::Tag :
          opserr << "tag ";
          break;
        case Positions::fy :
          opserr << "Fy ";
          break;
        case Positions::E:
          opserr << "E ";
          break;
        case Positions::b:
          opserr << "b ";
          break;
        case Positions::R0:
          opserr << "R0 ";
          break;
        case Positions::cR1:
          opserr << "cR1 ";
          break;
        case Positions::cR2:
          opserr << "cR2 ";
          break;
        case Positions::a1:
          opserr << "a1 ";
          break;
        case Positions::a2:
          opserr << "a2 ";
          break;
        case Positions::a3:
          opserr << "a3 ";
          break;
        case Positions::a4:
          opserr << "a4 ";
          break;
        case Positions::sig0:
          opserr << "sig0 ";
          break;

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


  if (strcmp(argv[1], "Steel1") == 0 || 
      strcmp(argv[1], "Steel01") == 0) {
    theMaterial = new Steel01(tag, fy, E, b, a1, a2, a3, a4);
  }

  else if (strcmp(argv[1], "Steel01Thermal") == 0) {
    theMaterial = new Steel01Thermal(tag, fy, E, b, a1, a2, a3, a4);
  }

  else if ((strcmp(argv[1], "Steel2") == 0)) {
    theMaterial = new Steel2(tag, fy, E, b, R0, cR1, cR2, a1, a2, a3, a4);
  }

  else if ((strcmp(argv[1], "Steel02") == 0)) {
    theMaterial = new Steel02(tag, fy, E, b, R0, cR1, cR2, a1, a2, a3, a4);
  }

  else if ((strcmp(argv[1], "Steel02Thermal") == 0)) {
    theMaterial = new Steel02Thermal(tag, fy, E, b, R0, cR1, cR2);
  }

  if (theMaterial == nullptr)
    return TCL_ERROR;

  return builder->addTaggedObject<UniaxialMaterial>(*theMaterial);
}


int
TclCommand_newFedeasSteel(ClientData clientData, Tcl_Interp *interp,
                          int argc, TCL_Char ** const argv)
{

  if (strcmp(argv[1], "Steel01") == 0 ||
      strcmp(argv[1], "Steel01Thermal") == 0 ||
      strcmp(argv[1], "Steel1") == 0) {

    // uniaxialMaterial Steel01 tag? fy? E? b? <a1? a2? a3? a4?>
    enum class Positions: int {
      Tag,
      fy, E, b,       EndRequired, 
      a1, a2, a3, a4, End,
      R0, cR1, cR2, sig0
    };

    return FedeasSteelParse<Positions>(clientData, interp, argc, argv);
  }

  else if ((strcmp(argv[1], "Steel02") == 0) || 
           (strcmp(argv[1], "Steel2") == 0) || 
           (strcmp(argv[1], "Steel02Thermal") == 0) || 
           (strcmp(argv[1], "SteelMP") == 0)
  ) {
    
    // uniaxialMaterial Steel02 $tag $Fy $E $b $R0 $cR1 $cR2 <$a1 $a2 $a3 $a4 $sigInit>
    enum class Positions: int {
      Tag,
      fy, E, b,                           EndRequired, 
      R0, cR1, cR2, a1, a2, a3, a4, sig0, End
    };
    return FedeasSteelParse<Positions>(clientData, interp, argc, argv);
  }

  return TCL_ERROR;
}


int
TclCommand_newFedeasConcrete(ClientData clientData, Tcl_Interp *interp,
                          int argc, TCL_Char ** const argv)
{

  if (strcmp(argv[1], "Concrete01") == 0 ||
      strcmp(argv[1], "Concrete1") == 0) {

    // uniaxialMaterial Concrete01 tag? fpc? epsc0? fpcu? epscu?
    enum class Positions: int {
      Tag,
      fpc, epsc0, fpcu, epscu,  EndRequired, 
      End,
      rat, ft, Ets
    };

    return FedeasConcrParse<Positions>(clientData, interp, argc, argv);
  }

  else if ((strcmp(argv[1], "Concrete02") == 0) || 
           (strcmp(argv[1], "Concrete2") == 0) || 
           (strcmp(argv[1], "Concrete02Thermal") == 0)
  ) {
    
    // uniaxialMaterial Concrete02 tag? fpc? epsc0? fpcu? epscu? rat? ft? Ets?
    enum class Positions: int {
      Tag,
      fpc, epsc0, fpcu, epscu, EndRequired,
      rat, ft, Ets, 
      End
    };
    return FedeasConcrParse<Positions>(clientData, interp, argc, argv);
  }

  return TCL_ERROR;
}




#if 0


int
TclBasicBuilder_addUniaxialConcrete(ClientData clientData, Tcl_Interp *interp,
                                  int argc, TCL_Char ** const argv)
{

  BasicModelBuilder *builder = static_cast<BasicModelBuilder *>(clientData);

  if (argc < 3) {
    opserr << "WARNING insufficient number of arguments\n";
    return TCL_ERROR;
  }

  // enum Positions {
  //   E, End
  // };


  UniaxialMaterial *theMaterial = nullptr;

  double fpc, epsc0, fpcu, epscu;
  double rat, ft, Ets;
  int tag;

  if (Tcl_GetInt(interp, argv[2], &tag) != TCL_OK) {
    opserr << "WARNING invalid uniaxialMaterial tag\n";
    return TCL_ERROR;
  }

  if (strcmp(argv[1], "Concrete1") == 0 ||
      strcmp(argv[1], "concrete01") == 0) {

    if (argc < 7) {
      opserr << "WARNING invalid number of arguments\n";
      opserr
          << "Want: uniaxialMaterial Concrete01 tag? fpc? epsc0? fpcu? epscu?"
          << endln;
      return TCL_ERROR;
    }

    if (Tcl_GetDouble(interp, argv[3], &fpc) != TCL_OK) {
      opserr << "WARNING invalid fpc\n";
      return TCL_ERROR;
    }
    if (Tcl_GetDouble(interp, argv[4], &epsc0) != TCL_OK) {
      opserr << "WARNING invalid epsc0\n";
      return TCL_ERROR;
    }
    if (Tcl_GetDouble(interp, argv[5], &fpcu) != TCL_OK) {
      opserr << "WARNING invalid fpcu\n";
      return TCL_ERROR;
    }
    if (Tcl_GetDouble(interp, argv[6], &epscu) != TCL_OK) {
      opserr << "WARNING invalid epscu\n";
      return TCL_ERROR;
    }

    theMaterial = new Concrete01(tag, fpc, epsc0, fpcu, epscu);
  }

  else if (strcmp(argv[1], "concr2") == 0) {
    if (argc < 10) {
      opserr << "WARNING invalid number of arguments\n";
      opserr << "Want: uniaxialMaterial Concrete02 tag? fpc? epsc0? fpcu? epscu? rat? ft? Ets?"
             << endln;
      return TCL_ERROR;
    }

    if (Tcl_GetDouble(interp, argv[3], &fpc) != TCL_OK) {
      opserr << "WARNING invalid fpc\n";
      return TCL_ERROR;
    }
    if (Tcl_GetDouble(interp, argv[4], &epsc0) != TCL_OK) {
      opserr << "WARNING invalid epsc0\n";
      return TCL_ERROR;
    }
    if (Tcl_GetDouble(interp, argv[5], &fpcu) != TCL_OK) {
      opserr << "WARNING invalid fpcu\n";
      return TCL_ERROR;
    }
    if (Tcl_GetDouble(interp, argv[6], &epscu) != TCL_OK) {
      opserr << "WARNING invalid epscu\n";
      return TCL_ERROR;
    }
    if (Tcl_GetDouble(interp, argv[7], &rat) != TCL_OK) {
      opserr << "WARNING invalid rat\n";
      return TCL_ERROR;
    }
    if (Tcl_GetDouble(interp, argv[8], &ft) != TCL_OK) {
      opserr << "WARNING invalid Ft\n";
      return TCL_ERROR;
    }
    if (Tcl_GetDouble(interp, argv[9], &Ets) != TCL_OK) {
      opserr << "WARNING invalid Ets\n";
      return TCL_ERROR;
    }

    theMaterial =
        new Concrete02(tag, fpc, epsc0, fpcu, epscu, rat, ft, Ets);
  }

  return builder->addTaggedObject<UniaxialMaterial>(*theMaterial);
}



#include <FedeasHardeningMaterial.h>
#include <FedeasBond1Material.h>
#include <FedeasBond2Material.h>
#include <FedeasConcr3Material.h>
#include <FedeasHyster1Material.h>
#include <FedeasHyster2Material.h>
#include <PlasticDamageMaterial.h>
int
Cmd(ClientData clientData, Tcl_Interp *interp,
    int argc, TCL_Char ** const argv)
  {
  if (strcmp(argv[1], "Hardening1") == 0 ||
  strcmp(argv[1], "Hardening01") == 0) {
  if (argc < 7) {
  opserr << "WARNING invalid number of arguments\n";
  opserr << "Want: uniaxialMaterial Hardening01 tag? E? sigY? Hiso? Hkin?"
        << endln;
  return TCL_ERROR;
  }

  double E, sigY, Hiso, Hkin;

  if (Tcl_GetDouble(interp, argv[3], &E) != TCL_OK) {
  opserr << "WARNING invalid E\n";
  return TCL_ERROR;
  }
  if (Tcl_GetDouble(interp, argv[4], &sigY) != TCL_OK) {
  opserr << "WARNING invalid sigY\n";
  return TCL_ERROR;
  }
  if (Tcl_GetDouble(interp, argv[5], &Hiso) != TCL_OK) {
  opserr << "WARNING invalid Hiso\n";
  return TCL_ERROR;
  }
  if (Tcl_GetDouble(interp, argv[6], &Hkin) != TCL_OK) {
  opserr << "WARNING invalid Hkin\n";
  return TCL_ERROR;
  }

  theMaterial = new FedeasHardeningMaterial(tag, E, sigY, Hiso, Hkin);
  }

  else if (strcmp(argv[1], "Bond1") == 0 || strcmp(argv[1], "Bond01") == 0) {
  if (argc < 15) {
  opserr << "WARNING invalid number of arguments\n";
  opserr << "Want: uniaxialMaterial Bond01 tag? u1p? q1p? u2p? u3p? q3p? "
            "u1n? q1n? u2n? u3n? q3n? s0? bb?"
        << endln;
  return TCL_ERROR;
  }

  double u1p, q1p, u2p, u3p, q3p;
  double u1n, q1n, u2n, u3n, q3n;
  double s0, bb;

  if (Tcl_GetDouble(interp, argv[3], &u1p) != TCL_OK) {
  opserr << "WARNING invalid u1p\n";
  return TCL_ERROR;
  }
  if (Tcl_GetDouble(interp, argv[4], &q1p) != TCL_OK) {
  opserr << "WARNING invalid q1p\n";
  return TCL_ERROR;
  }
  if (Tcl_GetDouble(interp, argv[5], &u2p) != TCL_OK) {
  opserr << "WARNING invalid u2p\n";
  return TCL_ERROR;
  }
  if (Tcl_GetDouble(interp, argv[6], &u3p) != TCL_OK) {
  opserr << "WARNING invalid u3p\n";
  return TCL_ERROR;
  }
  if (Tcl_GetDouble(interp, argv[7], &q3p) != TCL_OK) {
  opserr << "WARNING invalid q3p\n";
  return TCL_ERROR;
  }

  if (Tcl_GetDouble(interp, argv[8], &u1n) != TCL_OK) {
  opserr << "WARNING invalid u1n\n";
  return TCL_ERROR;
  }
  if (Tcl_GetDouble(interp, argv[9], &q1n) != TCL_OK) {
  opserr << "WARNING invalid q1n\n";
  return TCL_ERROR;
  }
  if (Tcl_GetDouble(interp, argv[10], &u2n) != TCL_OK) {
  opserr << "WARNING invalid u2n\n";
  return TCL_ERROR;
  }
  if (Tcl_GetDouble(interp, argv[11], &u3n) != TCL_OK) {
  opserr << "WARNING invalid u3n\n";
  return TCL_ERROR;
  }
  if (Tcl_GetDouble(interp, argv[12], &q3n) != TCL_OK) {
  opserr << "WARNING invalid q3n\n";
  return TCL_ERROR;
  }

  if (Tcl_GetDouble(interp, argv[13], &s0) != TCL_OK) {
  opserr << "WARNING invalid s0\n";
  return TCL_ERROR;
  }
  if (Tcl_GetDouble(interp, argv[14], &bb) != TCL_OK) {
  opserr << "WARNING invalid bb\n";
  return TCL_ERROR;
  }

  theMaterial = new FedeasBond1Material(tag, u1p, q1p, u2p, u3p, q3p, u1n,
                                      q1n, u2n, u3n, q3n, s0, bb);
  }

  else if (strcmp(argv[1], "Bond2") == 0 || strcmp(argv[1], "Bond02") == 0) {
  if (argc < 17) {

  opserr << "WARNING invalid number of arguments\n";
  opserr << "Want: uniaxialMaterial Bond02 tag? u1p? q1p? u2p? u3p? q3p? "
            "u1n? q1n? u2n? u3n? q3n? s0? bb? alp? aln?"
        << endln;
  return TCL_ERROR;
  }

  double u1p, q1p, u2p, u3p, q3p;
  double u1n, q1n, u2n, u3n, q3n;
  double s0, bb, alp, aln;

  if (Tcl_GetDouble(interp, argv[3], &u1p) != TCL_OK) {
  opserr << "WARNING invalid u1p\n";
  return TCL_ERROR;
  }
  if (Tcl_GetDouble(interp, argv[4], &q1p) != TCL_OK) {
  opserr << "WARNING invalid q1p\n";
  return TCL_ERROR;
  }
  if (Tcl_GetDouble(interp, argv[5], &u2p) != TCL_OK) {
  opserr << "WARNING invalid u2p\n";
  return TCL_ERROR;
  }
  if (Tcl_GetDouble(interp, argv[6], &u3p) != TCL_OK) {
  opserr << "WARNING invalid u3p\n";
  return TCL_ERROR;
  }
  if (Tcl_GetDouble(interp, argv[7], &q3p) != TCL_OK) {
  opserr << "WARNING invalid q3p\n";
  return TCL_ERROR;
  }

  if (Tcl_GetDouble(interp, argv[8], &u1n) != TCL_OK) {
  opserr << "WARNING invalid u1n\n";
  return TCL_ERROR;
  }
  if (Tcl_GetDouble(interp, argv[9], &q1n) != TCL_OK) {
  opserr << "WARNING invalid q1n\n";
  return TCL_ERROR;
  }
  if (Tcl_GetDouble(interp, argv[10], &u2n) != TCL_OK) {
  opserr << "WARNING invalid u2n\n";
  return TCL_ERROR;
  }
  if (Tcl_GetDouble(interp, argv[11], &u3n) != TCL_OK) {
  opserr << "WARNING invalid u3n\n";
  return TCL_ERROR;
  }
  if (Tcl_GetDouble(interp, argv[12], &q3n) != TCL_OK) {
  opserr << "WARNING invalid q3n\n";
  return TCL_ERROR;
  }

  if (Tcl_GetDouble(interp, argv[13], &s0) != TCL_OK) {
  opserr << "WARNING invalid s0\n";
  return TCL_ERROR;
  }
  if (Tcl_GetDouble(interp, argv[14], &bb) != TCL_OK) {
  opserr << "WARNING invalid bb\n";
  return TCL_ERROR;
  }
  if (Tcl_GetDouble(interp, argv[15], &alp) != TCL_OK) {
  opserr << "WARNING invalid alp\n";
  return TCL_ERROR;
  }
  if (Tcl_GetDouble(interp, argv[16], &aln) != TCL_OK) {
  opserr << "WARNING invalid aln\n";
  return TCL_ERROR;
  }

  theMaterial = new FedeasBond2Material(tag, u1p, q1p, u2p, u3p, q3p, u1n,
                                      q1n, u2n, u3n, q3n, s0, bb, alp, aln);
  }

  else if (strcmp(argv[1], "Concrete3") == 0 ||
  strcmp(argv[1], "Concrete03") == 0) {
  if (argc < 13) {
  opserr << "WARNING invalid number of arguments\n";
  opserr << "Want: uniaxialMaterial Concrete03 tag? fpc? epsc0? fpcu? "
      "epscu? rat? ft? epst0? ft0? beta? epstu?"
    << endln;
  return TCL_ERROR;
  }

  double fpc, epsc0, fpcu, epscu;
  double rat, ft, epst0, ft0, beta, epstu;

  if (Tcl_GetDouble(interp, argv[3], &fpc) != TCL_OK) {
  opserr << "WARNING invalid fpc\n";
  return TCL_ERROR;
  }
  if (Tcl_GetDouble(interp, argv[4], &epsc0) != TCL_OK) {
  opserr << "WARNING invalid epsc0\n";
  return TCL_ERROR;
  }
  if (Tcl_GetDouble(interp, argv[5], &fpcu) != TCL_OK) {
  opserr << "WARNING invalid fpcu\n";
  return TCL_ERROR;
  }
  if (Tcl_GetDouble(interp, argv[6], &epscu) != TCL_OK) {
  opserr << "WARNING invalid epscu\n";
  return TCL_ERROR;
  }
  if (Tcl_GetDouble(interp, argv[7], &rat) != TCL_OK) {
  opserr << "WARNING invalid rat\n";
  return TCL_ERROR;
  }
  if (Tcl_GetDouble(interp, argv[8], &ft) != TCL_OK) {
  opserr << "WARNING invalid ft\n";
  return TCL_ERROR;
  }
  if (Tcl_GetDouble(interp, argv[9], &epst0) != TCL_OK) {
  opserr << "WARNING invalid epst0\n";
  return TCL_ERROR;
  }
  if (Tcl_GetDouble(interp, argv[10], &ft0) != TCL_OK) {
  opserr << "WARNING invalid ft0\n";
  return TCL_ERROR;
  }
  if (Tcl_GetDouble(interp, argv[11], &beta) != TCL_OK) {
  opserr << "WARNING invalid beta\n";
  return TCL_ERROR;
  }
  if (Tcl_GetDouble(interp, argv[12], &epstu) != TCL_OK) {
  opserr << "WARNING invalid epstu\n";
  return TCL_ERROR;
  }

  theMaterial = new FedeasConcr3Material(tag, fpc, epsc0, fpcu, epscu, rat,
                                  ft, epst0, ft0, beta, epstu);
  }

  else if (strcmp(argv[1], "Hysteretic1") == 0 ||
  strcmp(argv[1], "Hysteretic01") == 0) {
  if (argc < 15) {
  opserr << "WARNING invalid number of arguments\n";
  opserr << "Want: uniaxialMaterial Hysteretic01 tag? s1p? e1p? s2p? e2p? "
      "s1n? e1n? s2n? e1n? px? py? d1? d2?"
    << endln;
  return TCL_ERROR;
  }

  double s1p, e1p, s2p, e2p;
  double s1n, e1n, s2n, e2n;
  double px, py, d1, d2;

  if (Tcl_GetDouble(interp, argv[3], &s1p) != TCL_OK) {
  opserr << "WARNING invalid s1p\n";
  return TCL_ERROR;
  }
  if (Tcl_GetDouble(interp, argv[4], &e1p) != TCL_OK) {
  opserr << "WARNING invalid e1p\n";
  return TCL_ERROR;
  }
  if (Tcl_GetDouble(interp, argv[5], &s2p) != TCL_OK) {
  opserr << "WARNING invalid s2p\n";
  return TCL_ERROR;
  }
  if (Tcl_GetDouble(interp, argv[6], &e2p) != TCL_OK) {
  opserr << "WARNING invalid e2p\n";
  return TCL_ERROR;
  }
  if (Tcl_GetDouble(interp, argv[7], &s1n) != TCL_OK) {
  opserr << "WARNING invalid s1n\n";
  return TCL_ERROR;
  }
  if (Tcl_GetDouble(interp, argv[8], &e1n) != TCL_OK) {
  opserr << "WARNING invalid e1n\n";
  return TCL_ERROR;
  }
  if (Tcl_GetDouble(interp, argv[9], &s2n) != TCL_OK) {
  opserr << "WARNING invalid s2n\n";
  return TCL_ERROR;
  }
  if (Tcl_GetDouble(interp, argv[10], &e2n) != TCL_OK) {
  opserr << "WARNING invalid e2n\n";
  return TCL_ERROR;
  }
  if (Tcl_GetDouble(interp, argv[11], &px) != TCL_OK) {
  opserr << "WARNING invalid px\n";
  return TCL_ERROR;
  }
  if (Tcl_GetDouble(interp, argv[12], &py) != TCL_OK) {
  opserr << "WARNING invalid py\n";
  return TCL_ERROR;
  }
  if (Tcl_GetDouble(interp, argv[13], &d1) != TCL_OK) {
  opserr << "WARNING invalid d1\n";
  return TCL_ERROR;
  }
  if (Tcl_GetDouble(interp, argv[14], &d2) != TCL_OK) {
  opserr << "WARNING invalid d2\n";
  return TCL_ERROR;
  }

  theMaterial = new FedeasHyster1Material(tag, s1p, e1p, s2p, e2p, s1n, e1n,
                                  s2n, e2n, px, py, d1, d2);
  }

  else if (strcmp(argv[1], "Hysteretic2") == 0 ||
  strcmp(argv[1], "Hysteretic02") == 0) {
  if (argc < 19) {
  opserr << "WARNING invalid number of arguments\n";
  opserr << "Want: uniaxialMaterial Hysteretic02 tag? s1p? e1p? s2p? e2p? "
      "s3p? e3p? s1n? e1n? s2n? e1n? s3n? e3n? px? py? d1? d2?"
    << endln;
  return TCL_ERROR;
  }

  double s1p, e1p, s2p, e2p, s3p, e3p;
  double s1n, e1n, s2n, e2n, s3n, e3n;
  double px, py, d1, d2;

  if (Tcl_GetDouble(interp, argv[3], &s1p) != TCL_OK) {
  opserr << "WARNING invalid s1p\n";
  return TCL_ERROR;
  }
  if (Tcl_GetDouble(interp, argv[4], &e1p) != TCL_OK) {
  opserr << "WARNING invalid e1p\n";
  return TCL_ERROR;
  }
  if (Tcl_GetDouble(interp, argv[5], &s2p) != TCL_OK) {
  opserr << "WARNING invalid s2p\n";
  return TCL_ERROR;
  }
  if (Tcl_GetDouble(interp, argv[6], &e2p) != TCL_OK) {
  opserr << "WARNING invalid e2p\n";
  return TCL_ERROR;
  }
  if (Tcl_GetDouble(interp, argv[7], &s3p) != TCL_OK) {
  opserr << "WARNING invalid s2p\n";
  return TCL_ERROR;
  }
  if (Tcl_GetDouble(interp, argv[8], &e3p) != TCL_OK) {
  opserr << "WARNING invalid e2p\n";
  return TCL_ERROR;
  }
  if (Tcl_GetDouble(interp, argv[9], &s1n) != TCL_OK) {
  opserr << "WARNING invalid s1n\n";
  return TCL_ERROR;
  }
  if (Tcl_GetDouble(interp, argv[10], &e1n) != TCL_OK) {
  opserr << "WARNING invalid e1n\n";
  return TCL_ERROR;
  }
  if (Tcl_GetDouble(interp, argv[11], &s2n) != TCL_OK) {
  opserr << "WARNING invalid s2n\n";
  return TCL_ERROR;
  }
  if (Tcl_GetDouble(interp, argv[12], &e2n) != TCL_OK) {
  opserr << "WARNING invalid e2n\n";
  return TCL_ERROR;
  }
  if (Tcl_GetDouble(interp, argv[13], &s3n) != TCL_OK) {
  opserr << "WARNING invalid s2n\n";
  return TCL_ERROR;
  }
  if (Tcl_GetDouble(interp, argv[14], &e3n) != TCL_OK) {
  opserr << "WARNING invalid e2n\n";
  return TCL_ERROR;
  }
  if (Tcl_GetDouble(interp, argv[15], &px) != TCL_OK) {
  opserr << "WARNING invalid px\n";
  return TCL_ERROR;
  }
  if (Tcl_GetDouble(interp, argv[16], &py) != TCL_OK) {
  opserr << "WARNING invalid py\n";
  return TCL_ERROR;
  }
  if (Tcl_GetDouble(interp, argv[17], &d1) != TCL_OK) {
  opserr << "WARNING invalid d1\n";
  return TCL_ERROR;
  }
  if (Tcl_GetDouble(interp, argv[18], &d2) != TCL_OK) {
  opserr << "WARNING invalid d2\n";
  return TCL_ERROR;
  }

  theMaterial =
  new FedeasHyster2Material(tag, s1p, e1p, s2p, e2p, s3p, e3p, s1n, e1n,
                        s2n, e2n, s3n, e3n, px, py, d1, d2);
  }

  else if (strcmp(argv[1], "ConcretePlasticDamage") == 0 ||
           strcmp(argv[1], "PlasticDamage") == 0) {
    if (argc < 11) {
      opserr << "WARNING invalid number of arguments\n";
      opserr << "Want: uniaxialMaterial ConcretePlasticDamage tag? $Ec $Gf $Gc "
                "$ft $fcy $fc $ktcr $relax"
             << endln;
      return TCL_ERROR;
    }

    double Ec, Ft, Fc, ft_max, fcy, fc, ktcr, relax;

    if (Tcl_GetDouble(interp, argv[3], &Ec) != TCL_OK) {
      opserr << OpenSees::PromptValueError
             << "invalid Ec\n";
      return TCL_ERROR;
    }
    if (Tcl_GetDouble(interp, argv[4], &Ft) != TCL_OK) {
      opserr << OpenSees::PromptValueError
             << "invalid Ft\n";
      return TCL_ERROR;
    }
    if (Tcl_GetDouble(interp, argv[5], &Fc) != TCL_OK) {
      opserr << OpenSees::PromptValueError
             << "invalid Fc\n";
      return TCL_ERROR;
    }
    if (Tcl_GetDouble(interp, argv[6], &ft_max) != TCL_OK) {
      opserr << OpenSees::PromptValueError
             << "invalid ft_max\n";
      return TCL_ERROR;
    }
    if (Tcl_GetDouble(interp, argv[7], &fcy) != TCL_OK) {
      opserr << OpenSees::PromptValueError
             << "invalid fcy\n";
      return TCL_ERROR;
    }
    if (Tcl_GetDouble(interp, argv[8], &fc) != TCL_OK) {
      opserr << OpenSees::PromptValueError
             << "invalid fc\n";
      return TCL_ERROR;
    }
    if (Tcl_GetDouble(interp, argv[9], &ktcr) != TCL_OK) {
      opserr << OpenSees::PromptValueError
             << "invalid Ktcr\n";
      return TCL_ERROR;
    }
    if (Tcl_GetDouble(interp, argv[10], &relax) != TCL_OK) {
      opserr << OpenSees::PromptValueError
             << "invalid relax\n";
      return TCL_ERROR;
    }

    theMaterial = new PlasticDamageMaterial(tag, Ec, Ft, Fc, ft_max, fcy, fc,
                                            ktcr, relax);
  }
  else {
    opserr << "WARNING invalid uniaxialMaterial type\n";
    return TCL_ERROR;
  }
}
#endif
