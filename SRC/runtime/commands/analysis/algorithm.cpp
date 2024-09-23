//===----------------------------------------------------------------------===//
//
//        OpenSees - Open System for Earthquake Engineering Simulation
//
//===----------------------------------------------------------------------===//
//
// Description: This file implements commands that allow for construction
// and interaction with Algorithm objects. Any command which requires
// access to specific Algorithm types (from the standard library) should
// be implemented here.
//
#include <stdio.h>
#include <assert.h>
#include <unordered_map>

#include <tcl.h>
#include <Logging.h>
#include <Parsing.h>
#include "BasicAnalysisBuilder.h"

// Algorithms
#include <Linear.h>
#include <NewtonRaphson.h>
#include <ModifiedNewton.h>
#include <NewtonHallM.h>
#include <Broyden.h>
#include <BFGS.h>
#include <KrylovNewton.h>
#include <PeriodicNewton.h>
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


OPS_Routine OPS_NewtonRaphsonAlgorithm;
OPS_Routine OPS_ExpressNewton;
OPS_Routine OPS_ModifiedNewton;
OPS_Routine OPS_NewtonHallM;

TclEquiSolnAlgo G3Parse_newEquiSolnAlgo;
TclEquiSolnAlgo G3Parse_newSecantNewtonAlgorithm;
TclEquiSolnAlgo G3_newNewtonLineSearch;
static TclEquiSolnAlgo G3_newKrylovNewton;
static TclEquiSolnAlgo G3_newBroyden;
static TclEquiSolnAlgo G3_newBFGS;

Tcl_CmdProc TclCommand_newLinearAlgorithm;
Tcl_CmdProc TclCommand_newNewtonRaphson;
Tcl_CmdProc TclCommand_newModifiedNewton;
Tcl_CmdProc TclCommand_newNewtonHallM;

namespace  OpenSees {
std::unordered_map<std::string, Tcl_CmdProc*> Algorithms {
  {"Linear",         TclCommand_newLinearAlgorithm},
  {"Newton",         TclCommand_newNewtonRaphson},
  {"NewtonHall",     TclCommand_newNewtonHallM},
  {"ModifiedNewton", TclCommand_newModifiedNewton}
};
}

//
// command invoked to allow the SolnAlgorithm object to be built
//
int
TclCommand_specifyAlgorithm(ClientData clientData, Tcl_Interp *interp, int argc,
                 TCL_Char ** const argv)
{

  BasicAnalysisBuilder *builder = (BasicAnalysisBuilder *)clientData;
  assert(builder != nullptr);

  // Make sure at least one other argument to contain numberer
  if (argc < 2) {
    opserr << G3_ERROR_PROMPT << "Need to specify an Algorithm type.\n";
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
G3Parse_newEquiSolnAlgo(ClientData clientData, Tcl_Interp *interp, int argc,
                        TCL_Char ** const argv)
{

  // check for type of Algorithm and create the object
  if (strcmp(argv[1], "Broyden") == 0)
    return G3_newBroyden(clientData, interp, argc, argv);

  else if (strcmp(argv[1], "BFGS") == 0)
    return G3_newBFGS(clientData, interp, argc, argv);

  else if (strcmp(argv[1], "SecantNewton") == 0)
    return G3Parse_newSecantNewtonAlgorithm(clientData, interp, argc, argv);

  else if (strcmp(argv[1], "NewtonLineSearch") == 0)
    return G3_newNewtonLineSearch(clientData, interp, argc, argv);

  else if (strcmp(argv[1], "KrylovNewton") == 0)
    return G3_newKrylovNewton(clientData, interp, argc, argv);


  EquiSolnAlgo *theNewAlgo = nullptr;
  G3_Runtime *rt = G3_getRuntime(interp);

  if (strcmp(argv[1], "Newton") == 0) {
    void *theNewtonAlgo = OPS_NewtonRaphsonAlgorithm(rt, argc, argv);
    theNewAlgo = (EquiSolnAlgo *)theNewtonAlgo;
  }

  else if ((strcmp(argv[1], "NewtonHallM") == 0) ||
           (strcmp(argv[1], "NewtonHall") == 0)) {
    void *theNewtonAlgo = OPS_NewtonHallM(rt, argc, argv);
    theNewAlgo = (EquiSolnAlgo *)theNewtonAlgo;
  }

  else if (strcmp(argv[1], "ModifiedNewton") == 0) {
    void *theNewtonAlgo = OPS_ModifiedNewton(rt, argc, argv);
    theNewAlgo = (EquiSolnAlgo *)theNewtonAlgo;
  }

  else if (strcmp(argv[1], "ExpressNewton") == 0) {
    void *theNewtonAlgo = OPS_ExpressNewton(rt, argc, argv);
    theNewAlgo = (EquiSolnAlgo *)theNewtonAlgo;
  }

  else {
    opserr << G3_ERROR_PROMPT << "No EquiSolnAlgo of type '" << argv[1] << "' exists\n";
    return nullptr;
  }

  return theNewAlgo;
}

int
TclCommand_newLinearAlgorithm(ClientData clientData, Tcl_Interp *interp, int argc,
                           TCL_Char ** const argv)
{
  BasicAnalysisBuilder *builder = (BasicAnalysisBuilder *)clientData;
  assert(builder != nullptr);

  int formTangent = CURRENT_TANGENT;
  int factorOnce = 0;
  int count = 2;
  while (count < argc) {
    if ((strcmp(argv[count], "-secant") == 0) ||
        (strcmp(argv[count], "-Secant") == 0)) {
      formTangent = CURRENT_SECANT;

    } else if ((strcmp(argv[count], "-initial") == 0) ||
               (strcmp(argv[count], "-Initial") == 0)) {
      formTangent = INITIAL_TANGENT;

    } else if ((strcmp(argv[count], "-factorOnce") == 0) ||
               (strcmp(argv[count], "-FactorOnce") == 0)) {
      factorOnce = 1;
    }
    count++;
  }

  builder->set(new Linear(formTangent, factorOnce));
  return TCL_OK;
}

int
TclCommand_newNewtonRaphson(ClientData clientData, Tcl_Interp* interp, int argc, TCL_Char**const argv)
{
  BasicAnalysisBuilder *builder = (BasicAnalysisBuilder *)clientData;
  assert(builder != nullptr);

  int formTangent = CURRENT_TANGENT;
  double iFactor = 0;
  double cFactor = 1;

  for (int i=2; i<argc; i++) {
    if (strcmp(argv[i],"-secant")==0 || 
        strcmp(argv[i],"-Secant")==0) {
      formTangent = CURRENT_SECANT;
      iFactor = 0;
      cFactor = 1.0;

    } else if (strcmp(argv[i],"-initial")==0 || 
               strcmp(argv[i],"-Initial")==0) {
      formTangent = INITIAL_TANGENT;
      iFactor = 1.;
      cFactor = 0;

    } else if (strcmp(argv[i],"-intialThenCurrent")==0 || 
               strcmp(argv[i],"-intialCurrent")==0) {
      formTangent = INITIAL_THEN_CURRENT_TANGENT;
      iFactor = 0;
      cFactor = 1.0;

    } else if (strcmp(argv[i],"-hall")==0 || 
               strcmp(argv[i],"-Hall")==0) {

      formTangent = HALL_TANGENT;
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

  builder->set(new NewtonRaphson(formTangent, iFactor, cFactor));

  return TCL_OK;
}


int
TclCommand_newModifiedNewton(ClientData clientData, Tcl_Interp* interp, int argc, TCL_Char**const argv)
{
  BasicAnalysisBuilder *builder = (BasicAnalysisBuilder *)clientData;
  assert(builder != nullptr);

  int formTangent = CURRENT_TANGENT;
  double iFactor = 0;
  double cFactor = 1;

  for (int i=2; i<argc; i++) {
    if (strcmp(argv[i],"-secant") == 0) {
      formTangent = CURRENT_SECANT;

    } else if (strcmp(argv[i],"-initial") == 0) {
      formTangent = INITIAL_TANGENT;

    } else if (strcmp(argv[i],"-hall")==0 || 
              strcmp(argv[i],"-Hall")==0) {
      formTangent = HALL_TANGENT;
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

  auto algorithm = new ModifiedNewton(formTangent, iFactor, cFactor);
  builder->set(algorithm);
  return TCL_OK;
}

int
TclCommand_newNewtonHallM(ClientData clientData, Tcl_Interp* interp, int argc, TCL_Char**const argv)
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



EquiSolnAlgo *
G3Parse_newSecantNewtonAlgorithm(ClientData clientData, Tcl_Interp *interp,
                                 int argc, TCL_Char ** const argv)
{
  assert(clientData != nullptr);
  BasicAnalysisBuilder *builder = (BasicAnalysisBuilder *)clientData;

  ConvergenceTest *theTest = builder->getConvergenceTest();

  if (theTest == nullptr) {
    opserr << G3_ERROR_PROMPT << "No ConvergenceTest yet specified\n";
    return nullptr;
  }

  int incrementTangent = CURRENT_TANGENT;
  int iterateTangent = CURRENT_TANGENT;
  int maxDim = 3;
  int numTerms = 2;
  for (int i = 2; i < argc; ++i) {
    if (strcmp(argv[i], "-iterate") == 0 && i + 1 < argc) {
      i++;
      if (strcmp(argv[i], "current") == 0)
        iterateTangent = CURRENT_TANGENT;
      if (strcmp(argv[i], "initial") == 0)
        iterateTangent = INITIAL_TANGENT;
      if (strcmp(argv[i], "noTangent") == 0)
        iterateTangent = NO_TANGENT;

    } else if (strcmp(argv[i], "-increment") == 0 && i + 1 < argc) {
      i++;
      if (strcmp(argv[i], "current") == 0)
        incrementTangent = CURRENT_TANGENT;
      if (strcmp(argv[i], "initial") == 0)
        incrementTangent = INITIAL_TANGENT;
      if (strcmp(argv[i], "noTangent") == 0)
        incrementTangent = NO_TANGENT;

    } else if (strcmp(argv[i], "-maxDim") == 0 && i + 1 < argc) {
      i++;
      maxDim = atoi(argv[i]);

    } else if (strcmp(argv[i], "-numTerms") == 0) {
      if (i+1 < argc)
        numTerms = atoi(argv[++i]);
      else {
        opserr << G3_ERROR_PROMPT << "Flag -numTerms requires follow up argument\n";
        return nullptr;
      }
    }
  }

  Accelerator *theAccel = nullptr;
  if (numTerms <= 1)
    theAccel = new SecantAccelerator1(maxDim, iterateTangent);
  if (numTerms >= 3)
    theAccel = new SecantAccelerator3(maxDim, iterateTangent);
  if (numTerms == 2)
    theAccel = new SecantAccelerator2(maxDim, iterateTangent);
  return new AcceleratedNewton(*theTest, theAccel, incrementTangent);
}

static EquiSolnAlgo *
G3_newBFGS(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char ** const argv)
{
  assert(clientData != nullptr);
  BasicAnalysisBuilder *builder = (BasicAnalysisBuilder *)clientData;

  ConvergenceTest *theTest = builder->getConvergenceTest();

  if (theTest == nullptr) {
    opserr << G3_ERROR_PROMPT << "No ConvergenceTest yet specified\n";
    return nullptr;
  }


  int formTangent = CURRENT_TANGENT;
  int count = -1;
  for (int i = 2; i < argc; ++i) {
    if (strcmp(argv[i], "-secant") == 0) {
      formTangent = CURRENT_SECANT;
    } else if (strcmp(argv[i], "-initial") == 0) {
      formTangent = INITIAL_TANGENT;
    } else if (strcmp(argv[i++], "-count") == 0 && i < argc) {
      count = atoi(argv[i]);
    }
  }

  EquiSolnAlgo *theNewAlgo = nullptr;
  if (count == -1)
    theNewAlgo = new BFGS(*theTest, formTangent);
  else
    theNewAlgo = new BFGS(*theTest, formTangent, count);

  return theNewAlgo;
}

EquiSolnAlgo *
G3_newNewtonLineSearch(ClientData clientData, Tcl_Interp *interp, int argc,
                       TCL_Char ** const argv)
{

  assert(clientData != nullptr);
  BasicAnalysisBuilder *builder = (BasicAnalysisBuilder *)clientData;

  ConvergenceTest *theTest = builder->getConvergenceTest();

  if (theTest == nullptr) {
    opserr << G3_ERROR_PROMPT << " No ConvergenceTest yet specified\n";
    return nullptr;
  }

  int count = 2;

  // set some default variable
  double tol = 0.8;
  int maxIter = 10;
  double maxEta = 10.0;
  double minEta = 0.1;
  int pFlag = 1;
  int typeSearch = 0;

  while (count < argc) {
    if (strcmp(argv[count], "-tol") == 0) {
      count++;
      if (Tcl_GetDouble(interp, argv[count], &tol) != TCL_OK)
        return nullptr;
      count++;

    } else if (strcmp(argv[count], "-maxIter") == 0) {
      count++;
      if (Tcl_GetInt(interp, argv[count], &maxIter) != TCL_OK)
        return nullptr;

      count++;
    } else if (strcmp(argv[count], "-pFlag") == 0) {
      count++;
      if (Tcl_GetInt(interp, argv[count], &pFlag) != TCL_OK)
        return nullptr;
      count++;

    } else if (strcmp(argv[count], "-minEta") == 0) {
      count++;
      if (Tcl_GetDouble(interp, argv[count], &minEta) != TCL_OK)
        return nullptr;
      count++;

    } else if (strcmp(argv[count], "-maxEta") == 0) {
      count++;
      if (Tcl_GetDouble(interp, argv[count], &maxEta) != TCL_OK)
        return nullptr;
      count++;

    } else if (strcmp(argv[count], "-type") == 0) {
      count++;
      if (strcmp(argv[count], "Bisection") == 0)
        typeSearch = 1;
      else if (strcmp(argv[count], "Secant") == 0)
        typeSearch = 2;
      else if (strcmp(argv[count], "RegulaFalsi") == 0)
        typeSearch = 3;
      else if (strcmp(argv[count], "LinearInterpolated") == 0)
        typeSearch = 3;
      else if (strcmp(argv[count], "InitialInterpolated") == 0)
        typeSearch = 0;
      count++;

    } else
      count++;
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


  EquiSolnAlgo *theNewAlgo = nullptr;
  theNewAlgo = new NewtonLineSearch(*theTest, theLineSearch);
  return theNewAlgo;
}

static EquiSolnAlgo *
G3_newKrylovNewton(ClientData clientData, Tcl_Interp *interp, int argc,
                   TCL_Char ** const argv)
{
  assert(clientData != nullptr);
  BasicAnalysisBuilder *builder = (BasicAnalysisBuilder *)clientData;
  ConvergenceTest *theTest = builder->getConvergenceTest();

  if (theTest == nullptr) {
    opserr << G3_ERROR_PROMPT << "No ConvergenceTest yet specified\n";
    return nullptr;
  }

  int incrementTangent = CURRENT_TANGENT;
  int iterateTangent = CURRENT_TANGENT;
  int maxDim = 3;

  for (int i = 2; i < argc; ++i) {
    if (strcmp(argv[i], "-iterate") == 0 && i + 1 < argc) {
      i++;
      if (strcmp(argv[i], "current") == 0)
        iterateTangent = CURRENT_TANGENT;
      if (strcmp(argv[i], "initial") == 0)
        iterateTangent = INITIAL_TANGENT;
      if (strcmp(argv[i], "noTangent") == 0)
        iterateTangent = NO_TANGENT;

    } else if (strcmp(argv[i], "-increment") == 0 && i + 1 < argc) {
      i++;
      if (strcmp(argv[i], "current") == 0)
        incrementTangent = CURRENT_TANGENT;
      if (strcmp(argv[i], "initial") == 0)
        incrementTangent = INITIAL_TANGENT;
      if (strcmp(argv[i], "noTangent") == 0)
        incrementTangent = NO_TANGENT;

    } else if (strcmp(argv[i], "-maxDim") == 0 && i + 1 < argc) {
      i++;
      maxDim = atoi(argv[i]);
    }
  }

  Accelerator *theAccel = new KrylovAccelerator(maxDim, iterateTangent);

  EquiSolnAlgo *theNewAlgo = new AcceleratedNewton(*theTest, theAccel, incrementTangent);
  return theNewAlgo;
}

EquiSolnAlgo *
G3_newRaphsonNewton(ClientData clientData, Tcl_Interp *interp, int argc,
                    TCL_Char ** const argv)
{
  assert(clientData != nullptr);
  BasicAnalysisBuilder *builder = (BasicAnalysisBuilder *)clientData;
  ConvergenceTest *theTest = builder->getConvergenceTest();

  if (theTest == nullptr) {
    opserr << G3_ERROR_PROMPT << "No ConvergenceTest yet specified\n";
    return nullptr;
  }

  int incrementTangent = CURRENT_TANGENT;
  int iterateTangent = CURRENT_TANGENT;
  for (int i = 2; i < argc; ++i) {
    if (strcmp(argv[i], "-iterate") == 0 && i + 1 < argc) {
      i++;
      if (strcmp(argv[i], "current") == 0)
        iterateTangent = CURRENT_TANGENT;
      if (strcmp(argv[i], "initial") == 0)
        iterateTangent = INITIAL_TANGENT;
      if (strcmp(argv[i], "noTangent") == 0)
        iterateTangent = NO_TANGENT;

    } else if (strcmp(argv[i], "-increment") == 0 && i + 1 < argc) {
      i++;
      if (strcmp(argv[i], "current") == 0)
        incrementTangent = CURRENT_TANGENT;
      if (strcmp(argv[i], "initial") == 0)
        incrementTangent = INITIAL_TANGENT;
      if (strcmp(argv[i], "noTangent") == 0)
        incrementTangent = NO_TANGENT;
    }
  }

  Accelerator *theAccel;
  theAccel = new RaphsonAccelerator(iterateTangent);

  EquiSolnAlgo *theNewAlgo = new AcceleratedNewton(*theTest, theAccel, incrementTangent);
  return theNewAlgo;
}

EquiSolnAlgo *
G3_newMillerNewton(ClientData clientData, Tcl_Interp *interp, int argc,
                   TCL_Char ** const argv)
{
  assert(clientData != nullptr);
  BasicAnalysisBuilder *builder = (BasicAnalysisBuilder *)clientData;
  ConvergenceTest *theTest = builder->getConvergenceTest();

  if (theTest == nullptr) {
    opserr << G3_ERROR_PROMPT << "No ConvergenceTest yet specified\n";
    return nullptr;
  }

  int incrementTangent = CURRENT_TANGENT;
  int iterateTangent = CURRENT_TANGENT;
  int maxDim = 3;

  for (int i = 2; i < argc; ++i) {
    if (strcmp(argv[i], "-iterate") == 0 && i + 1 < argc) {
      i++;
      if (strcmp(argv[i], "current") == 0)
        iterateTangent = CURRENT_TANGENT;
      if (strcmp(argv[i], "initial") == 0)
        iterateTangent = INITIAL_TANGENT;
      if (strcmp(argv[i], "noTangent") == 0)
        iterateTangent = NO_TANGENT;
    } else if (strcmp(argv[i], "-increment") == 0 && i + 1 < argc) {
      i++;
      if (strcmp(argv[i], "current") == 0)
        incrementTangent = CURRENT_TANGENT;
      if (strcmp(argv[i], "initial") == 0)
        incrementTangent = INITIAL_TANGENT;
      if (strcmp(argv[i], "noTangent") == 0)
        incrementTangent = NO_TANGENT;
    } else if (strcmp(argv[i], "-maxDim") == 0 && i + 1 < argc) {
      i++;
      maxDim = atoi(argv[i]);
    }
  }

  Accelerator *theAccel = new MillerAccelerator(maxDim, 0.01, iterateTangent);

  EquiSolnAlgo *theNewAlgo = new AcceleratedNewton(*theTest, theAccel, incrementTangent);
  return theNewAlgo;
}

EquiSolnAlgo *
G3_newPeriodicNewton(ClientData clientData, Tcl_Interp *interp, int argc,
                     TCL_Char ** const argv)
{
  assert(clientData != nullptr);
  BasicAnalysisBuilder *builder = (BasicAnalysisBuilder *)clientData;
  ConvergenceTest *theTest = builder->getConvergenceTest();

  if (theTest == nullptr) {
    opserr << G3_ERROR_PROMPT << "No ConvergenceTest yet specified\n";
    return nullptr;
  }

  int incrementTangent = CURRENT_TANGENT;
  int iterateTangent = CURRENT_TANGENT;
  int maxDim = 3;
  for (int i = 2; i < argc; ++i) {
    if (strcmp(argv[i], "-iterate") == 0 && i + 1 < argc) {
      i++;
      if (strcmp(argv[i], "current") == 0)
        iterateTangent = CURRENT_TANGENT;
      if (strcmp(argv[i], "initial") == 0)
        iterateTangent = INITIAL_TANGENT;
      if (strcmp(argv[i], "noTangent") == 0)
        iterateTangent = NO_TANGENT;

    } else if (strcmp(argv[i], "-increment") == 0 && i + 1 < argc) {
      i++;
      if (strcmp(argv[i], "current") == 0)
        incrementTangent = CURRENT_TANGENT;
      if (strcmp(argv[i], "initial") == 0)
        incrementTangent = INITIAL_TANGENT;
      if (strcmp(argv[i], "noTangent") == 0)
        incrementTangent = NO_TANGENT;

    } else if (strcmp(argv[i], "-maxDim") == 0 && i + 1 < argc) {
      i++;
      maxDim = atoi(argv[i]);
    }
  }

  Accelerator *theAccel;
  theAccel = new PeriodicAccelerator(maxDim, iterateTangent);

  EquiSolnAlgo *theNewAlgo = nullptr;
  theNewAlgo = new AcceleratedNewton(*theTest, theAccel, incrementTangent);
  return theNewAlgo;
}

static EquiSolnAlgo *
G3_newBroyden(ClientData clientData, Tcl_Interp *interp, int argc,
              TCL_Char ** const argv)
{
  assert(clientData != nullptr);
  BasicAnalysisBuilder *builder = (BasicAnalysisBuilder *)clientData;
  ConvergenceTest *theTest = builder->getConvergenceTest();

  int formTangent = CURRENT_TANGENT;
  int count = -1;

  if (theTest == nullptr) {
    opserr << G3_ERROR_PROMPT << "No ConvergenceTest yet specified\n";
    return nullptr;
  }
  for (int i = 2; i < argc; ++i) {
    if (strcmp(argv[i], "-secant") == 0) {
      formTangent = CURRENT_SECANT;

    } else if (strcmp(argv[i], "-initial") == 0) {
      formTangent = INITIAL_TANGENT;

    } else if (strcmp(argv[i++], "-count") == 0 && i < argc) {
      count = atoi(argv[i]);
    }
  }

  EquiSolnAlgo *theNewAlgo = nullptr;
  if (count == -1)
    theNewAlgo = new Broyden(*theTest, formTangent);
  else
    theNewAlgo = new Broyden(*theTest, formTangent, count);

  return theNewAlgo;
}

//
// Other commands
//
int
printAlgorithm(ClientData clientData, Tcl_Interp *interp, int argc,
               TCL_Char ** const argv, OPS_Stream &output)
{
  assert(clientData != nullptr);
  BasicAnalysisBuilder *builder = (BasicAnalysisBuilder *)clientData;
  EquiSolnAlgo* theAlgorithm = builder->getAlgorithm();

  int eleArg = 0;
  if (theAlgorithm == nullptr) {
    opserr << G3_ERROR_PROMPT << "No algorithm has been set.\n";
    return TCL_ERROR;
  }

  // if just 'print <filename> algorithm'- no flag
  if (argc == 0) {
    theAlgorithm->Print(output);
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
TclCommand_accelCPU(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char ** const argv)
{
  assert(clientData != nullptr);
  BasicAnalysisBuilder *builder = (BasicAnalysisBuilder *)clientData;
  EquiSolnAlgo* algo = builder->getAlgorithm();

  if (algo == nullptr)
    return TCL_ERROR;

  Tcl_SetObjResult(interp, Tcl_NewDoubleObj(algo->getAccelTimeCPU()));

  return TCL_OK;
}

int
TclCommand_numFact(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char ** const argv)
{
  assert(clientData != nullptr);
  BasicAnalysisBuilder *builder = (BasicAnalysisBuilder *)clientData;
  EquiSolnAlgo* algo = builder->getAlgorithm();

  if (algo == nullptr)
    return TCL_ERROR;

  Tcl_SetObjResult(interp, Tcl_NewIntObj(algo->getNumFactorizations()));

  return TCL_OK;
}

int
TclCommand_algorithmRecorder(ClientData clientData, Tcl_Interp *interp, int argc,
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
TclCommand_totalCPU(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char ** const argv)
{
  assert(clientData != nullptr);
  EquiSolnAlgo *algo = ((BasicAnalysisBuilder *)clientData)->getAlgorithm();

  if (algo == nullptr)
    return TCL_ERROR;

  Tcl_SetObjResult(interp, Tcl_NewDoubleObj(algo->getTotalTimeCPU()));

  return TCL_OK;
}

int
TclCommand_solveCPU(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char ** const argv)
{
  assert(clientData != nullptr);
  EquiSolnAlgo *algo = ((BasicAnalysisBuilder *)clientData)->getAlgorithm();


  if (algo == nullptr)
    return TCL_ERROR;

  Tcl_SetObjResult(interp, Tcl_NewDoubleObj(algo->getSolveTimeCPU()));

  return TCL_OK;
}


int
TclCommand_numIter(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char ** const argv)
{
  assert(clientData != nullptr);
  EquiSolnAlgo *algo = ((BasicAnalysisBuilder *)clientData)->getAlgorithm();

  if (algo == nullptr)
    return TCL_ERROR;

  Tcl_SetObjResult(interp, Tcl_NewIntObj(algo->getNumIterations()));

  return TCL_OK;
}

