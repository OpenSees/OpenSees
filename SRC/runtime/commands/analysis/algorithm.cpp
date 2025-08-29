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
//
// Description: This file implements commands that allow for construction
// and interaction with Algorithm objects. Any command which requires
// access to specific Algorithm types (from the standard library) should
// be implemented here.
//
#include <set>
#include <stdio.h>
#include <assert.h>
#include <unordered_map>

#include <tcl.h>
#include <Logging.h>
#include <Parsing.h>
#include <ArgumentTracker.h>
#include "BasicAnalysisBuilder.h"

// Algorithms
#include <Linear.h>
#include <NewtonRaphson.h>
#include <ModifiedNewton.h>
#include <NewtonHallM.h>
#include <Broyden.h>
#include <BFGS.h>
#include <AcceleratedNewton.h>
#include <ExpressNewton.h>

// LineSearch
#include <NewtonLineSearch.h>
#include <BisectionLineSearch.h>
#include <InitialInterpolatedLineSearch.h>
#include <RegulaFalsiLineSearch.h>
#include <SecantLineSearch.h>

// Accelerators
#include <RaphsonAccelerator.h>
#include <PeriodicAccelerator.h>
#include <KrylovAccelerator.h>
#include <SecantAccelerator1.h>
#include <SecantAccelerator2.h>
#include <SecantAccelerator3.h>
#include <MillerAccelerator.h>


#include <runtimeAPI.h>
class G3_Runtime;

extern "C" int OPS_ResetInputNoBuilder(ClientData clientData,
                                       Tcl_Interp *interp, int cArg, int mArg,
                                       TCL_Char ** const argv, Domain *domain);

typedef EquiSolnAlgo *(TclEquiSolnAlgo)(ClientData, Tcl_Interp *, int, TCL_Char **);


OPS_Routine OPS_ExpressNewton;
OPS_Routine OPS_NewtonHallM;

TclEquiSolnAlgo G3Parse_newEquiSolnAlgo;
static TclEquiSolnAlgo G3_newBroyden;
static TclEquiSolnAlgo G3_newBFGS;

Tcl_CmdProc TclCommand_newLinearAlgorithm;
Tcl_CmdProc TclCommand_newNewtonRaphson;
Tcl_CmdProc TclCommand_newNewtonHallM;
Tcl_CmdProc TclCommand_newAcceleratedNewton;
static Tcl_CmdProc TclCommand_newNewtonLineSearch;

namespace  OpenSees {
std::unordered_map<std::string, Tcl_CmdProc*> Algorithms {
  {"Linear",            TclCommand_newLinearAlgorithm},

  {"Newton",            TclCommand_newNewtonRaphson},
  {"ModifiedNewton",    TclCommand_newNewtonRaphson},
  {"NewtonHall",        TclCommand_newNewtonHallM},
  {"NewtonLineSearch",  TclCommand_newNewtonLineSearch},

  {"SecantNewton",      TclCommand_newAcceleratedNewton},
  {"MillerAccelerator", TclCommand_newAcceleratedNewton},
  {"KrylovNewton",      TclCommand_newAcceleratedNewton},
  {"PeriodicNewton",    TclCommand_newAcceleratedNewton},
  {"RaphsonNewton",     TclCommand_newAcceleratedNewton},
};
}

//
// command invoked to allow the SolnAlgorithm object to be built
//
int
TclCommand_specifyAlgorithm(ClientData clientData, Tcl_Interp *interp, Tcl_Size argc,
                 TCL_Char ** const argv)
{

  BasicAnalysisBuilder *builder = (BasicAnalysisBuilder *)clientData;
  assert(builder != nullptr);

  // Make sure at least one other argument to contain numberer
  if (argc < 2) {
    opserr << OpenSees::PromptValueError << "Need to specify an Algorithm type.\n";
    return TCL_ERROR;
  }


  auto command = OpenSees::Algorithms.find(std::string{argv[1]});
  if (command != OpenSees::Algorithms.end())
    return (*command->second)(clientData, interp, argc, argv);


  OPS_ResetInputNoBuilder(nullptr, interp, 2, argc, argv, nullptr);

  EquiSolnAlgo *theNewAlgo = G3Parse_newEquiSolnAlgo(clientData, interp, argc, argv);

  if (theNewAlgo == nullptr) {
    // Leave it to parsing routine to print error info, this way
    // we get more detail.
    return TCL_ERROR;

  } else {
    builder->set(theNewAlgo);
  }
  return TCL_OK;
}


EquiSolnAlgo *
G3Parse_newEquiSolnAlgo(ClientData clientData, Tcl_Interp *interp, Tcl_Size argc,
                        TCL_Char ** const argv)
{

  // check for type of Algorithm and create the object
  if (strcmp(argv[1], "Broyden") == 0)
    return G3_newBroyden(clientData, interp, argc, argv);

  else if (strcmp(argv[1], "BFGS") == 0)
    return G3_newBFGS(clientData, interp, argc, argv);

  EquiSolnAlgo *theNewAlgo = nullptr;
  G3_Runtime *rt = G3_getRuntime(interp);

  if ((strcmp(argv[1], "NewtonHallM") == 0) ||
           (strcmp(argv[1], "NewtonHall") == 0)) {
    void *theNewtonAlgo = OPS_NewtonHallM(rt, argc, argv);
    theNewAlgo = (EquiSolnAlgo *)theNewtonAlgo;
  }

  else if (strcmp(argv[1], "ExpressNewton") == 0) {
    void *theNewtonAlgo = OPS_ExpressNewton(rt, argc, argv);
    theNewAlgo = (EquiSolnAlgo *)theNewtonAlgo;
  }

  else {
    opserr << OpenSees::PromptValueError
           << "Unknown algorithm type '" << argv[1] << "'\n";
    return nullptr;
  }

  return theNewAlgo;
}


int
TclCommand_newLinearAlgorithm(ClientData clientData, Tcl_Interp *interp, Tcl_Size argc,
                           TCL_Char ** const argv)
{
  BasicAnalysisBuilder *builder = (BasicAnalysisBuilder *)clientData;
  assert(builder != nullptr);

  IncrementalIntegrator::TangentFlagType correction_tangent = CURRENT_TANGENT;
  int factorOnce = 0;
  int count = 2;
  while (count < argc) {
    if ((strcmp(argv[count], "-secant") == 0) ||
        (strcmp(argv[count], "-Secant") == 0)) {
      correction_tangent = CURRENT_SECANT;

    } else if ((strcmp(argv[count], "-initial") == 0) ||
               (strcmp(argv[count], "-Initial") == 0)) {
      correction_tangent = INITIAL_TANGENT;

    } else if ((strcmp(argv[count], "-factorOnce") == 0) ||
               (strcmp(argv[count], "-FactorOnce") == 0)) {
      factorOnce = 1;
    }
    count++;
  }

  builder->set(new Linear(correction_tangent, factorOnce));
  return TCL_OK;
}

int
TclCommand_newNewtonRaphson(ClientData clientData, 
                            Tcl_Interp* interp, 
                            Tcl_Size argc, TCL_Char**const argv)
{
  BasicAnalysisBuilder *builder = (BasicAnalysisBuilder *)clientData;
  assert(builder != nullptr);

  IncrementalIntegrator::TangentFlagType 
     correction_tangent = CURRENT_TANGENT,
     prediction_tangent = CURRENT_TANGENT;
  double iFactor = 0;
  double cFactor = 1;

  for (int i=2; i<argc; i++) {
    if (strcmp(argv[i],"-current")==0 ||
        strcmp(argv[i],"-currentTangent")==0) {
      correction_tangent = CURRENT_TANGENT;
      iFactor = 0;
      cFactor = 1.0;
    } else if (strcmp(argv[i],"-secant")==0 ||
        strcmp(argv[i],"-Secant")==0) {
      correction_tangent = CURRENT_SECANT;
      iFactor = 0;
      cFactor = 1.0;

    } else if (strcmp(argv[i],"-initial")==0 || 
               strcmp(argv[i],"-initialTangent")==0 ||
               strcmp(argv[i],"-Initial")==0) {
      correction_tangent = INITIAL_TANGENT;
      iFactor = 1.0;
      cFactor = 0.0;
    } else if (strcmp(argv[i],"-intialThenCurrent")==0 || 
               strcmp(argv[i],"-intialCurrent")==0) {
      prediction_tangent = INITIAL_TANGENT;
      iFactor = 0.0;
      cFactor = 1.0;

    } else if (strcmp(argv[i],"-hall")==0 || 
               strcmp(argv[i],"-Hall")==0) {
      correction_tangent = HALL_TANGENT;
      iFactor = 0.1;
      cFactor = 0.9;
      if (argc-i >= 2) {
        if (Tcl_GetDouble(interp, argv[i+1], &iFactor) != TCL_OK) {
          opserr << "WARNING invalid data reading ifactor\n";
          return TCL_ERROR;
        }
        if (Tcl_GetDouble(interp, argv[i+2], &cFactor) != TCL_OK) {
          opserr << "WARNING invalid data reading cfactor\n";
          return TCL_ERROR;
        }
        i += 2;
      }
    }
  }

  //
  //
  if (strcmp(argv[1], "Newton") == 0 ||
      strcmp(argv[1], "NewtonRaphson") == 0) {
    builder->set(new NewtonRaphson(prediction_tangent, 
                                   correction_tangent, iFactor, cFactor));
    return TCL_OK;
  }
  else if (strcmp(argv[1], "ModifiedNewton") == 0) {
    builder->set(new ModifiedNewton(correction_tangent, iFactor, cFactor));
    return TCL_OK;
  }
  return TCL_ERROR;
}


int
TclCommand_newNewtonHallM(ClientData clientData, Tcl_Interp* interp, Tcl_Size argc, TCL_Char**const argv)
{
  BasicAnalysisBuilder *builder = (BasicAnalysisBuilder *)clientData;
  assert(builder != nullptr);

  int method = 0;
  double iFactor = .1;
  double alpha = .01;
  double c = 100;
  double data[2];
  
  int numData = 1;
  if (OPS_GetDoubleInput(&numData,&data[0]) < 0) {
    opserr << "WARNING invalid data reading 2 hall factors\n";
    return 0;
  }
  iFactor = data[0];

  for (int i=2; i<argc; i++) {

      if (strcmp(argv[i], "-exp")==0 || 
          strcmp(argv[i], "-Exp")==0) {
        numData = 1;
        if (OPS_GetDoubleInput(&numData,&data[0]) < 0) {
          opserr << "WARNING invalid data reading 2 hall factors\n";
          return 0;
        } else 
          alpha = data[0];

      } else if (strcmp(argv[i], "-sigmoid")==0 || 
                 strcmp(argv[i], "-Sigmoid")==0) {
        method = 1;
        int numData = 2;
        if (OPS_GetDoubleInput(&numData,&data[0]) < 0) {
          opserr << "WARNING invalid data reading 2 hall factors\n";
          return 0;
        } else {
          alpha = data[0];
          c = data[1];
        }

      } else if (strcmp(argv[i], "-constant")==0 || 
                 strcmp(argv[i], "-Constant")==0) {
        method = 2;
        int numData = 1;
        if (OPS_GetDoubleInput(&numData,&data[0]) < 0) {
          opserr << "WARNING invalid data reading 2 hall factors\n";
          return 0;
        } else {
          c = data[0];
        }
      }
    }

  auto algorithm = new NewtonHallM(iFactor, method, alpha, c);
  builder->set(algorithm);
  return TCL_OK;
}



static EquiSolnAlgo *
G3_newBFGS(ClientData clientData, Tcl_Interp *interp, Tcl_Size argc, TCL_Char ** const argv)
{
  assert(clientData != nullptr);

  int correction_tangent = CURRENT_TANGENT;
  int count = -1;
  for (int i = 2; i < argc; ++i) {
    if (strcmp(argv[i], "-secant") == 0) {
      correction_tangent = CURRENT_SECANT;
    }
    else if (strcmp(argv[i], "-initial") == 0) {
      correction_tangent = INITIAL_TANGENT;
    }

    else if (strcmp(argv[i++], "-count") == 0 && i < argc) {
      if (i >= argc) {
        opserr << OpenSees::PromptValueError 
               << "Flag -count requires follow up argument\n";
        return nullptr;
      }
      if (Tcl_GetInt(interp, argv[i], &count) != TCL_OK) {
        opserr << OpenSees::PromptValueError 
               << "Invalid value for -count: " << argv[i] << "\n";
        return nullptr;
      }
    }
  }

  EquiSolnAlgo *theNewAlgo = nullptr;
  if (count == -1)
    theNewAlgo = new BFGS(correction_tangent, 10);
  else
    theNewAlgo = new BFGS(correction_tangent, count);

  return theNewAlgo;
}

static int
TclCommand_newNewtonLineSearch(ClientData clientData, Tcl_Interp *interp, Tcl_Size argc,
                       TCL_Char ** const argv)
{
  assert(clientData != nullptr);
  BasicAnalysisBuilder *builder = (BasicAnalysisBuilder *)clientData;

  enum class Positions: int {
    EndRequired,
    Tolerance,
    End,
    MaxIter,
    MaxEta,
    MinEta,
    PFlag,
    TypeSearch
  };

  enum class LineSearchType: int {
    InitialInterpolated,
    Bisection,
    Secant,
    RegulaFalsi
  };

  ArgumentTracker<Positions> tracker;
  std::set<int> positions;

  // set some default variable
  double tol = 0.8;
  int maxIter = 10;
  double maxEta = 10.0;
  double minEta = 0.1;
  int pFlag = 1;
  int typeSearch = 0;

  for (int i=2; i<argc; i++) {
    if (strcmp(argv[i], "-tol") == 0) {
      if (++i >= argc) {
        opserr << OpenSees::PromptValueError 
               << "Flag -tol requires follow up argument\n";
        return TCL_ERROR;
      }
      if (Tcl_GetDouble(interp, argv[i], &tol) != TCL_OK) {
        opserr << OpenSees::PromptValueError 
               << "Invalid value for -tol: " << argv[i] << "\n";
        return TCL_ERROR;
      }
      tracker.consume(Positions::Tolerance);
    }
    else if (strcmp(argv[i], "-maxIter") == 0) {
      if (++i >= argc) {
        opserr << OpenSees::PromptValueError 
               << "Flag -maxIter requires follow up argument\n";
        return TCL_ERROR;
      }
      if (Tcl_GetInt(interp, argv[i], &maxIter) != TCL_OK) {
        opserr << OpenSees::PromptValueError 
               << "Invalid value for -maxIter: " << argv[i] << "\n";
        return TCL_ERROR;
      }
      tracker.consume(Positions::MaxIter);

    } else if (strcmp(argv[i], "-pFlag") == 0) {
      if (++i >= argc) {
        opserr << OpenSees::PromptValueError 
               << "Flag -pFlag requires follow up argument"
               << OpenSees::SignalMessageEnd;
        return TCL_ERROR;
      }
      if (Tcl_GetInt(interp, argv[i], &pFlag) != TCL_OK) {
        opserr << OpenSees::PromptValueError 
               << "Invalid value for -pFlag: " << argv[i] << "\n";
        return TCL_ERROR;
      }
      tracker.consume(Positions::PFlag);

    } else if (strcmp(argv[i], "-minEta") == 0) {
      if (++i >= argc) {
        opserr << OpenSees::PromptValueError 
               << "Flag -minEta requires follow up argument\n";
        return TCL_ERROR;
      }
      if (Tcl_GetDouble(interp, argv[i], &minEta) != TCL_OK) {
        opserr << OpenSees::PromptValueError 
               << "Invalid value for -minEta: " << argv[i] << "\n";
        return TCL_ERROR;
      }
      tracker.consume(Positions::MinEta);

    } else if (strcmp(argv[i], "-maxEta") == 0) {
      if (++i >= argc) {
        opserr << OpenSees::PromptValueError 
               << "Flag -maxEta requires follow up argument\n";
        return TCL_ERROR;
      }
      if (Tcl_GetDouble(interp, argv[i], &maxEta) != TCL_OK) {
        opserr << OpenSees::PromptValueError 
               << "Invalid value for -maxEta: " << argv[i] << "\n";
        return TCL_ERROR;
      }
      tracker.consume(Positions::MaxEta);
    } 

    else if (strcmp(argv[i], "-type") == 0) {
      if (++i >= argc) {
        opserr << OpenSees::PromptValueError 
               << "Flag -type requires follow up argument\n";
        return TCL_ERROR;
      }
      if (strcmp(argv[i], "Bisection") == 0) {
        typeSearch = 1;
      } else if (strcmp(argv[i], "Secant") == 0) {
        typeSearch = 2;
      } else if (strcmp(argv[i], "RegulaFalsi") == 0 ||
                 strcmp(argv[i], "LinearInterpolated") == 0) {
        typeSearch = 3;
      } else if (strcmp(argv[i], "InitialInterpolated") == 0) {
        typeSearch = 0;
      } else {
        opserr << OpenSees::PromptValueError 
               << "Unknown line search type: " << argv[i] << "\n";
        return TCL_ERROR;
      }
      tracker.consume(Positions::TypeSearch);
    }
    else {
      positions.insert(i);
    }
  }

  for (int i : positions) {
    if (tracker.current() == Positions::EndRequired) {
      tracker.increment();
    }
    switch (tracker.current()) {
      case Positions::Tolerance:
        if (Tcl_GetDouble(interp, argv[i], &tol) != TCL_OK) {
          opserr << OpenSees::PromptValueError 
                 << "Invalid value for -tol: " << argv[i] << "\n";
          return TCL_ERROR;
        }
        tracker.consume(Positions::Tolerance);
        break;
      case Positions::MaxIter:
        if (Tcl_GetInt(interp, argv[i], &maxIter) != TCL_OK) {
          opserr << OpenSees::PromptValueError 
                 << "Invalid value for -maxIter: " << argv[i] << "\n";
          return TCL_ERROR;
        }
        tracker.consume(Positions::MaxIter);
        break;
      case Positions::PFlag:
        if (Tcl_GetInt(interp, argv[i], &pFlag) != TCL_OK) {
          opserr << OpenSees::PromptValueError 
                 << "Invalid value for -pFlag: " << argv[i] << "\n";
          return TCL_ERROR;
        }
        tracker.consume(Positions::PFlag);
        break;
      case Positions::MinEta:
        if (Tcl_GetDouble(interp, argv[i], &minEta) != TCL_OK) {
          opserr << OpenSees::PromptValueError 
                 << "Invalid value for -minEta: " << argv[i] << "\n";
          return TCL_ERROR;
        }
        tracker.consume(Positions::MinEta);
        break;
      case Positions::MaxEta:
        if (Tcl_GetDouble(interp, argv[i], &maxEta) != TCL_OK) {
          opserr << OpenSees::PromptValueError 
                 << "Invalid value for -maxEta: " << argv[i] << "\n";
          return TCL_ERROR;
        }
        tracker.consume(Positions::MaxEta);
        break;
      default:
        opserr << OpenSees::PromptValueError 
               << "Unexpected argument: " << argv[i] << "\n";
        return TCL_ERROR;
    }
  }

  LineSearch *theLineSearch = nullptr;
  if (typeSearch == 0)
    theLineSearch = new InitialInterpolatedLineSearch(tol, maxIter, minEta, maxEta, pFlag);
  else if (typeSearch == 1)
    theLineSearch = new BisectionLineSearch(tol, maxIter, minEta, maxEta, pFlag);
  else if (typeSearch == 2)
    theLineSearch = new SecantLineSearch(tol, maxIter, minEta, maxEta, pFlag);
  else if (typeSearch == 3)
    theLineSearch = new RegulaFalsiLineSearch(tol, maxIter, minEta, maxEta, pFlag);


  builder->set(new NewtonLineSearch(theLineSearch));
  return TCL_OK;
}



int
TclCommand_newAcceleratedNewton(ClientData clientData, Tcl_Interp *interp, Tcl_Size argc,
                             TCL_Char ** const argv)
{
  assert(clientData != nullptr);
  BasicAnalysisBuilder *builder = (BasicAnalysisBuilder *)clientData;
  
  IncrementalIntegrator::TangentFlagType
    incrementTangent = CURRENT_TANGENT,
    correction_tangent = CURRENT_TANGENT;

  int maxDim   = 3;
  int numTerms = 2; // Default for SecantAccelerator
  double iFactor = 0.0; // Initial factor for RaphsonAccelerator
  double cFactor = 1.0; // Current factor for RaphsonAccelerator

  enum class AcceleratorType {
    Miller,
    Secant,
    Raphson,
    Krylov,
    Periodic
  } type = AcceleratorType::Krylov;

  for (int i = 2; i < argc; ++i) {
    if ((strcmp(argv[i], "-iterate") == 0) ||
        (strcmp(argv[i], "-correction") == 0)) {
      if (++i >= argc) {
        opserr << OpenSees::PromptValueError 
               << "Flag -iterate requires follow up argument"
               << OpenSees::SignalMessageEnd;
        return TCL_ERROR;
      }
      if ((strcmp(argv[i], "current") == 0) || (strcmp(argv[i], "tangent") == 0))
        correction_tangent = CURRENT_TANGENT;
      else if (strcmp(argv[i], "initial") == 0)
        correction_tangent = INITIAL_TANGENT;
      else if ((strcmp(argv[i], "noTangent") == 0) || strcmp(argv[i], "None") == 0)
        correction_tangent = NO_TANGENT;
    }

    else if ((strcmp(argv[i], "-increment") == 0) ||
             (strcmp(argv[i], "-prediction") == 0)) {
      if (++i >= argc) {
        opserr << OpenSees::PromptValueError 
               << "Flag -increment requires follow up argument"
               << OpenSees::SignalMessageEnd;
        return TCL_ERROR;
      }
      if ((strcmp(argv[i], "current") == 0) || (strcmp(argv[i], "tangent") == 0))
        incrementTangent = CURRENT_TANGENT;
      else if (strcmp(argv[i], "initial") == 0)
        incrementTangent = INITIAL_TANGENT;
      else if ((strcmp(argv[i], "noTangent") == 0) || strcmp(argv[i], "None") == 0)
        incrementTangent = NO_TANGENT;
    }
    
    else if (strcmp(argv[i], "-maxDim") == 0 && i + 1 < argc) {
      i++;
      if (i >= argc || Tcl_GetInt(interp, argv[i], &maxDim) != TCL_OK) {
        opserr << OpenSees::PromptValueError 
               << "Flag -maxDim requires follow up argument"
               << OpenSees::SignalMessageEnd;
        return TCL_ERROR;
      }
    }

    else if (strcmp(argv[i], "-accelerator") == 0 && i + 1 < argc) {
      i++;
      if (i >= argc) {
        opserr << OpenSees::PromptValueError 
               << "Flag -accelerator requires follow up argument"
               << OpenSees::SignalMessageEnd;
        return TCL_ERROR;
      }
      if (strcmp(argv[i], "Miller") == 0) {
        type = AcceleratorType::Miller;
      } else if (strcmp(argv[i], "Secant") == 0 || strcmp(argv[i], "SecantNewton") == 0) {
        type = AcceleratorType::Secant;
      } else if (strcmp(argv[i], "Raphson") == 0 || strcmp(argv[i], "RaphsonNewton") == 0) {
        type = AcceleratorType::Raphson;
      } else if (strcmp(argv[i], "Krylov") == 0 || strcmp(argv[i],  "KrylovNewton") == 0) {
        type = AcceleratorType::Krylov;
      } else if (strcmp(argv[i], "Periodic") == 0 || strcmp(argv[i], "PeriodicNewton") == 0) {
        type = AcceleratorType::Periodic;
      } else {
        opserr << OpenSees::PromptValueError 
               << "Unknown accelerator '" << argv[i] << "'"
               << OpenSees::SignalMessageEnd;
        return TCL_ERROR;
      }
    }

    // SecantAccelerator
    else if (strcmp(argv[i], "-numTerms") == 0) {
      if (i+1 < argc)
        numTerms = atoi(argv[++i]);
      else {
        opserr << OpenSees::PromptValueError 
               << "Flag -numTerms requires follow up argument"
               << OpenSees::SignalMessageEnd;
        return TCL_ERROR;
      }
    }
  }

  Accelerator *accel = nullptr;
  if (strcmp(argv[1], "MillerNewton") == 0) {
    accel = new MillerAccelerator(maxDim, 0.01, correction_tangent);
  }

  else if (type == AcceleratorType::Secant ||
           strcmp(argv[1], "SecantNewton") == 0) {
    if (numTerms <= 1)
      accel = new SecantAccelerator1(maxDim, correction_tangent);
    else if (numTerms >= 3)
      accel = new SecantAccelerator3(maxDim, correction_tangent);
    else
      accel = new SecantAccelerator2(maxDim, correction_tangent);
  }

  else if (type == AcceleratorType::Raphson ||
           strcmp(argv[1], "RaphsonNewton") == 0) {
    accel = new RaphsonAccelerator(correction_tangent, iFactor, cFactor);
  }

  else if (type == AcceleratorType::Krylov ||
           strcmp(argv[1], "KrylovNewton") == 0) {
    accel = new KrylovAccelerator(maxDim, correction_tangent);
  }

  else if (type ==  AcceleratorType::Periodic ||
           strcmp(argv[1], "PeriodicNewton") == 0) {
    accel = new PeriodicAccelerator(maxDim, correction_tangent);
  }

  else {
    opserr << OpenSees::PromptValueError 
           << "Unknown accelerator type '" << argv[1] << "'\n";
    return TCL_ERROR;
  }

  builder->set(new AcceleratedNewton(accel, incrementTangent));

  return TCL_OK;
}


static EquiSolnAlgo *
G3_newBroyden(ClientData clientData, Tcl_Interp *interp, Tcl_Size argc,
              TCL_Char ** const argv)
{
  assert(clientData != nullptr);
  BasicAnalysisBuilder *builder = (BasicAnalysisBuilder *)clientData;

  int correction_tangent = CURRENT_TANGENT;
  int count = -1;

  for (int i = 2; i < argc; ++i) {
    if (strcmp(argv[i], "-secant") == 0) {
      correction_tangent = CURRENT_SECANT;

    } else if (strcmp(argv[i], "-initial") == 0) {
      correction_tangent = INITIAL_TANGENT;

    } else if (strcmp(argv[i++], "-count") == 0 && i < argc) {
      count = atoi(argv[i]);
    }
  }

  EquiSolnAlgo *theNewAlgo = nullptr;
  if (count == -1)
    theNewAlgo = new Broyden(correction_tangent);
  else
    theNewAlgo = new Broyden(correction_tangent, count);

  return theNewAlgo;
}

//
// Other commands
//
int
printAlgorithm(ClientData clientData, Tcl_Interp *interp, Tcl_Size argc,
               TCL_Char ** const argv, OPS_Stream &output)
{
  assert(clientData != nullptr);
  BasicAnalysisBuilder *builder = (BasicAnalysisBuilder *)clientData;
  const EquiSolnAlgo* theAlgorithm = builder->getAlgorithm();

  int eleArg = 0;
  if (theAlgorithm == nullptr) {
    opserr << OpenSees::PromptValueError << "No algorithm has been set.\n";
    return TCL_ERROR;
  }

  // if just 'print <filename> algorithm'- no flag
  if (argc == 0) {
    theAlgorithm->Print(output, 0);
    return TCL_OK;
  }

  // if 'print <filename> Algorithm flag' get the flag
  int flag;
  if (Tcl_GetInt(interp, argv[eleArg], &flag) != TCL_OK) {
    opserr << "WARNING print algorithm failed to get integer flag: \n";
    opserr << argv[eleArg] << endln;
    return TCL_ERROR;
  }
  theAlgorithm->Print(output, flag);
  return TCL_OK;
}

int
TclCommand_accelCPU(ClientData clientData, Tcl_Interp *interp, Tcl_Size argc, TCL_Char ** const argv)
{
  assert(clientData != nullptr);
  BasicAnalysisBuilder *builder = (BasicAnalysisBuilder *)clientData;
  const EquiSolnAlgo* algo = builder->getAlgorithm();

  if (algo == nullptr)
    return TCL_ERROR;

  Tcl_SetObjResult(interp, Tcl_NewDoubleObj(algo->getAccelTimeCPU()));

  return TCL_OK;
}


int
TclCommand_numFact(ClientData clientData, Tcl_Interp *interp, Tcl_Size argc, TCL_Char ** const argv)
{
  assert(clientData != nullptr);
  BasicAnalysisBuilder *builder = (BasicAnalysisBuilder *)clientData;
  const EquiSolnAlgo* algo = builder->getAlgorithm();

  if (algo == nullptr)
    return TCL_ERROR;

  Tcl_SetObjResult(interp, Tcl_NewIntObj(algo->getNumFactorizations()));

  return TCL_OK;
}

int
TclCommand_algorithmRecorder(ClientData clientData, Tcl_Interp *interp, Tcl_Size argc,
                        TCL_Char ** const argv)
{
#if 1
  return TCL_ERROR;
#else
  Recorder *theRecorder = nullptr;
  BasicAnalysisBuilder* builder = (BasicAnalysisBuilder*)clientData;
  Domain* domain = builder->getDomain();
  EquiSolnAlgo *theAlgo;
  TclCreateRecorder(clientData, interp, argc, argv, *domain, &theRecorder);

  if (theRecorder == nullptr) {
    char buffer[] = "-1";
    Tcl_SetResult(interp, buffer, TCL_VOLATILE);
    return TCL_ERROR;
  }

  // add the recorder to the domain,
  // NOTE: will not be called with theALgo == 0
  if (theAlgo != nullptr) {
    if ((theAlgo->addRecorder(*theRecorder)) < 0) {
      opserr << "WARNING could not add to domain - recorder " << argv[1]
             << endln;
      delete theRecorder;
      return TCL_ERROR;
    }
  }

  int recorderTag = theRecorder->getTag();
  Tcl_SetObjResult(interp, Tcl_NewIntObj(recorderTag));
  return TCL_OK;
#endif
}

int
TclCommand_totalCPU(ClientData clientData, Tcl_Interp *interp, Tcl_Size argc, TCL_Char ** const argv)
{
  assert(clientData != nullptr);
  const EquiSolnAlgo *algo = ((BasicAnalysisBuilder *)clientData)->getAlgorithm();

  if (algo == nullptr)
    return TCL_ERROR;

  Tcl_SetObjResult(interp, Tcl_NewDoubleObj(algo->getTotalTimeCPU()));

  return TCL_OK;
}

int
TclCommand_solveCPU(ClientData clientData, Tcl_Interp *interp, Tcl_Size argc, TCL_Char ** const argv)
{
  assert(clientData != nullptr);
  const EquiSolnAlgo *algo = ((BasicAnalysisBuilder *)clientData)->getAlgorithm();


  if (algo == nullptr)
    return TCL_ERROR;

  Tcl_SetObjResult(interp, Tcl_NewDoubleObj(algo->getSolveTimeCPU()));

  return TCL_OK;
}


int
TclCommand_numIter(ClientData clientData, Tcl_Interp *interp, Tcl_Size argc, TCL_Char ** const argv)
{
  assert(clientData != nullptr);
  const EquiSolnAlgo *algo = ((BasicAnalysisBuilder *)clientData)->getAlgorithm();

  if (algo == nullptr)
    return TCL_ERROR;

  Tcl_SetObjResult(interp, Tcl_NewIntObj(algo->getNumIterations()));

  return TCL_OK;
}

