//===----------------------------------------------------------------------===//
//
//        OpenSees - Open System for Earthquake Engineering Simulation    
//
//===----------------------------------------------------------------------===//
//
// Description: This file contains functions that are responsible
// for orchestrating an analysis.
//
#include <tcl.h>
#include <assert.h>
#include <runtimeAPI.h>
#include <G3_Logging.h>
#include <StandardStream.h>
#include <FileStream.h>

#include <Matrix.h>
#include <Domain.h> // for modal damping
#include <AnalysisModel.h>

#include "BasicAnalysisBuilder.h"

#include <StaticAnalysis.h>
#include <DirectIntegrationAnalysis.h>
#include <VariableTimeStepDirectIntegrationAnalysis.h>

#include <EigenSOE.h>
#include <LinearSOE.h>
// for printA
#include <FullGenLinLapackSolver.h>
#include <FullGenLinSOE.h>
#include <GimmeMCK.h>

#include <LoadControl.h>
#include <EquiSolnAlgo.h>

#include <TransientIntegrator.h>
#include <StaticIntegrator.h>

// constraint handlers
#include <PlainHandler.h>
#include <PenaltyConstraintHandler.h>
#include <LagrangeConstraintHandler.h>
#include <TransformationConstraintHandler.h>

// numberers
#include <PlainNumberer.h>
#include <DOF_Numberer.h>
#include "analysis.h"


// for response spectrum analysis
extern void OPS_DomainModalProperties(G3_Runtime*);
extern int OPS_ResponseSpectrumAnalysis(G3_Runtime*);
extern "C" int OPS_ResetInputNoBuilder(ClientData clientData,
                                       Tcl_Interp *interp, int cArg, int mArg,
                                       TCL_Char ** const argv, Domain *domain);

Tcl_CmdProc TclCommand_clearAnalysis;
Tcl_CmdProc TclCommand_setNumberer;


//
// Add commands to the interpreter that take the AnalysisBuilder as clientData.
//
int
G3_AddTclAnalysisAPI(Tcl_Interp *interp, Domain* domain)
{

  BasicAnalysisBuilder *builder = new BasicAnalysisBuilder(domain);
  Tcl_CreateCommand(interp, "wipeAnalysis", &wipeAnalysis, builder, nullptr);
  Tcl_CreateCommand(interp, "_clearAnalysis", &TclCommand_clearAnalysis, builder, nullptr);

  Tcl_CreateCommand(interp, "numberer",   TclCommand_setNumberer, builder, nullptr);


  static int ncmd = sizeof(tcl_analysis_cmds)/sizeof(char_cmd);
  for (int i = 0; i < ncmd; ++i)
    Tcl_CreateCommand(interp, 
        tcl_analysis_cmds[i].name, 
        tcl_analysis_cmds[i].func, 
        (ClientData) builder, nullptr);

  return TCL_OK;
}

//
// command invoked to build an Analysis object
//
static int
specifyAnalysis(ClientData clientData, Tcl_Interp *interp, int argc,
                TCL_Char ** const argv)
{
  assert(clientData != nullptr);
  BasicAnalysisBuilder *builder = (BasicAnalysisBuilder*)clientData;

  if (argc < 2) {
    opserr << G3_ERROR_PROMPT << "need to specify an analysis type (Static, Transient)\n";
    return TCL_ERROR;
  }

  int argi = 1;

  if (strcmp(argv[argi], "-linear") == 0) {
    if (argc < 3) {
      opserr << G3_ERROR_PROMPT << "need to specify an analysis type (Static, Transient)\n";
      return TCL_ERROR;
    }
    Tcl_Eval(interp, "algorithm Linear\n"
                     "test FixedNumIter 1\n"
    );
    argi++;
  }
  if (strcmp(argv[argi], "Static") == 0) {
    builder->setStaticAnalysis();
    return TCL_OK;

  } else if (strcmp(argv[argi], "Transient") == 0) {
    builder->setTransientAnalysis();
    return TCL_OK;
  }

  else if (((strcmp(argv[1], "VariableTimeStepTransient") == 0) ||
          (strcmp(argv[1], "TransientWithVariableTimeStep") == 0) ||
          (strcmp(argv[1], "VariableTransient") == 0))) {
    opserr << "Unimplemented\n";
    return TCL_ERROR;

  } else {
    opserr << G3_ERROR_PROMPT << "Analysis type '" << argv[1]
      << "' does not exists (Static or Transient only). \n";
    return TCL_ERROR;
  }

  return TCL_OK;
}

//
// Command invoked to build the model, i.e. to invoke analyze()
// on the Analysis object
//
static int
analyzeModel(ClientData clientData, Tcl_Interp *interp, int argc,
             TCL_Char ** const argv)
{
  assert(clientData != nullptr);
  BasicAnalysisBuilder *builder = (BasicAnalysisBuilder*)clientData;

  VariableTimeStepDirectIntegrationAnalysis* theVariableTimeStepTransientAnalysis =
      builder->getVariableTimeStepDirectIntegrationAnalysis();

  int result = 0;
  switch (builder->CurrentAnalysisFlag) {
    case BasicAnalysisBuilder::STATIC_ANALYSIS: {
      int numIncr;
      if (argc < 2) {
        opserr << G3_ERROR_PROMPT << "static analysis: analysis numIncr?\n";
        return TCL_ERROR;
      }

      if (Tcl_GetInt(interp, argv[1], &numIncr) != TCL_OK)
        return TCL_ERROR;

      result = builder->analyze(numIncr, 0.0);
      break;
    }
    case BasicAnalysisBuilder::TRANSIENT_ANALYSIS: {
      double dT;
      int numIncr;
      if (argc < 3) {
        opserr << G3_ERROR_PROMPT << "transient analysis: analysis numIncr? deltaT?\n";
        return TCL_ERROR;
      }
      if (Tcl_GetInt(interp, argv[1], &numIncr) != TCL_OK)
        return TCL_ERROR;
      if (Tcl_GetDouble(interp, argv[2], &dT) != TCL_OK)
        return TCL_ERROR;

      if (argc == 6) {
        int Jd;
        double dtMin, dtMax;
        if (Tcl_GetDouble(interp, argv[3], &dtMin) != TCL_OK)
          return TCL_ERROR;
        if (Tcl_GetDouble(interp, argv[4], &dtMax) != TCL_OK)
          return TCL_ERROR;
        if (Tcl_GetInt(interp, argv[5], &Jd) != TCL_OK)
          return TCL_ERROR;

        if (theVariableTimeStepTransientAnalysis != nullptr)
          result = theVariableTimeStepTransientAnalysis->analyze(
              numIncr, dT, dtMin, dtMax, Jd);
        else {
          opserr << G3_ERROR_PROMPT << "analyze - no variable time step transient analysis "
                    "object constructed\n";
          return TCL_ERROR;
        }

      } else {
  //    result = theTransientAnalysis->analyze(numIncr, dT);
        result = builder->analyze(numIncr, dT);
      }
      break;
    }
    default:
      opserr << G3_ERROR_PROMPT << "No Analysis type has been specified \n";
      return TCL_ERROR;
  }

  char buffer[10];
  sprintf(buffer, "%d", result);
  Tcl_SetResult(interp, buffer, TCL_VOLATILE);

  return TCL_OK;
}


static int
initializeAnalysis(ClientData clientData, Tcl_Interp *interp, int argc,
                   TCL_Char ** const argv)
{
  assert(clientData != nullptr);
  BasicAnalysisBuilder *builder = (BasicAnalysisBuilder*)clientData;
  if (builder->initialize() < 0)
    return TCL_ERROR;
  return TCL_OK;
}



static int
eigenAnalysis(ClientData clientData, Tcl_Interp *interp, int argc,
              TCL_Char ** const argv)
{
  static  char *resDataPtr = 0;
  static  int resDataSize = 0;

  BasicAnalysisBuilder *builder = (BasicAnalysisBuilder*)clientData;

  Domain *domain = builder->getDomain();


  // make sure at least one other argument to contain type of system
  if (argc < 2) {
    opserr << G3_ERROR_PROMPT << "eigen <type> numModes?\n";
    return TCL_ERROR;
  }

  bool generalizedAlgo = true; 
      // 0 - frequency/generalized (default),
      // 1 - standard, 
      // 2 - buckling
  int typeSolver = EigenSOE_TAGS_ArpackSOE;
  int loc = 1;
  double shift = 0.0;
  bool findSmallest = true;
  int numEigen = 0;

  // Check type of eigenvalue analysis
  while (loc < (argc - 1)) {
    if ((strcmp(argv[loc], "frequency") == 0) ||
        (strcmp(argv[loc], "-frequency") == 0) ||
        (strcmp(argv[loc], "generalized") == 0) ||
        (strcmp(argv[loc], "-generalized") == 0))
      generalizedAlgo = true;

    else if ((strcmp(argv[loc], "standard") == 0) ||
             (strcmp(argv[loc], "-standard") == 0)) {
      generalizedAlgo = false;
      typeSolver = EigenSOE_TAGS_SymBandEigenSOE;
    }

    else if ((strcmp(argv[loc], "-findLargest") == 0))
      findSmallest = false;

    else if ((strcmp(argv[loc], "genBandArpack") == 0) ||
             (strcmp(argv[loc], "-genBandArpack") == 0) ||
             (strcmp(argv[loc], "genBandArpackEigen") == 0) ||
             (strcmp(argv[loc], "-genBandArpackEigen") == 0))
      typeSolver = EigenSOE_TAGS_ArpackSOE;

    else if ((strcmp(argv[loc], "symmBandLapack") == 0) ||
             (strcmp(argv[loc], "-symmBandLapack") == 0) ||
             (strcmp(argv[loc], "symmBandLapackEigen") == 0) ||
             (strcmp(argv[loc], "-symmBandLapackEigen") == 0))
      typeSolver = EigenSOE_TAGS_SymBandEigenSOE;

    else if ((strcmp(argv[loc], "fullGenLapack") == 0) ||
             (strcmp(argv[loc], "-fullGenLapack") == 0) ||
             (strcmp(argv[loc], "fullGenLapackEigen") == 0) ||
             (strcmp(argv[loc], "-fullGenLapackEigen") == 0))
      typeSolver = EigenSOE_TAGS_FullGenEigenSOE;

    else {
      opserr << "eigen - unknown option: " << argv[loc] << endln;
    }

    loc++;
  }

  // check argv[loc] for number of modes
  if ((Tcl_GetInt(interp, argv[loc], &numEigen) != TCL_OK) || numEigen < 0) {
    opserr << G3_ERROR_PROMPT << "eigen numModes?  - invalid numModes\n";
    return TCL_ERROR;
  }

  int requiredDataSize = 40 * numEigen;
  if (requiredDataSize > resDataSize) {
    if (resDataPtr != nullptr)
      delete[] resDataPtr;

    resDataPtr = new char[requiredDataSize];
    resDataSize = requiredDataSize;
  }

  for (int i = 0; i < requiredDataSize; ++i)
    resDataPtr[i] = '\n';

  //
  // create a transient analysis if no analysis exists
  // 
  builder->newEigenAnalysis(typeSolver, shift);

  int result = builder->eigen(numEigen,generalizedAlgo,findSmallest);

  if (result == 0) {
    const Vector &eigenvalues = domain->getEigenvalues();
    int cnt = 0;
    for (int i = 0; i < numEigen; ++i) {
      cnt += sprintf(&resDataPtr[cnt], "%35.20f  ", eigenvalues[i]);
    }

    Tcl_SetResult(interp, resDataPtr, TCL_STATIC);
  }

  return TCL_OK;
}

static int
modalProperties(ClientData clientData, Tcl_Interp *interp, int argc,
                TCL_Char ** const argv)
{
  G3_Runtime *rt = G3_getRuntime(interp);
  OPS_ResetInputNoBuilder(clientData, interp, 1, argc, argv, nullptr);
  OPS_DomainModalProperties(rt);
  return TCL_OK;
}

static int
responseSpectrum(ClientData clientData, Tcl_Interp *interp, int argc,
                 TCL_Char ** const argv)
{
  OPS_ResetInputNoBuilder(clientData, interp, 1, argc, argv, nullptr);
  G3_Runtime *rt = G3_getRuntime(interp);
  OPS_ResponseSpectrumAnalysis(rt);
  return TCL_OK;
}

// TODO: Move this to commands/modeling/damping.cpp? ...but it uses and
// AnalysisBuilder
static int
modalDamping(ClientData clientData, Tcl_Interp *interp, int argc,
             TCL_Char ** const argv)
{
  assert(clientData != nullptr);
  BasicAnalysisBuilder *builder = (BasicAnalysisBuilder*)clientData;
  int numEigen = builder->getNumEigen();

  if (argc < 2) {
    opserr
        << G3_ERROR_PROMPT << argv[0] << " ?factor - not enough arguments to command\n";
    return TCL_ERROR;
  }

  int numModes = argc - 1;

  if (numEigen == 0) {
    opserr << G3_WARN_PROMPT 
           << "- " << argv[0] << " - eigen command needs to be called first\n";

    numEigen = numModes;
    builder->newEigenAnalysis(EigenSOE_TAGS_ArpackSOE, 0.0);
    builder->eigen(numModes, true, true);
    // return TCL_ERROR;
  }

  /* 
   * "quick" modal damping adds modal damping forces to the right-hand side,
   * but does not add modal damping terms to the dynamic tangent.
   *
   * see https://portwooddigital.com/2022/11/08/quick-and-dirty-modal-damping/
   */
  bool do_tangent = true;
  if (strcmp(argv[0], "modalDampingQ") == 0)
    do_tangent = false;

  double factor = 0;
  Vector modalDampingValues(numEigen);

  if (numModes != 1 && numModes != numEigen) {
    // TODO: Just call eigen again?
    opserr << G3_ERROR_PROMPT << "modalDampingQ - same number of damping factors as modes must be "
              "specified\n";
//  opserr << "                    - same damping ratio will be applied to all\n";
    return TCL_ERROR;
  }

  //
  // read in values and set factors
  //
  if (numModes == numEigen) {

    // read in all factors one at a time
    for (int i = 0; i < numEigen; ++i) {
      if (Tcl_GetDouble(interp, argv[1 + i], &factor) != TCL_OK) {
        opserr << G3_ERROR_PROMPT << argv[0] << " - could not read factor at position "
               << i << "\n";
        return TCL_ERROR;
      }
      modalDampingValues[i] = factor;
    }

  } else {
    //  read in one & set all factors to that value
    if (Tcl_GetDouble(interp, argv[1], &factor) != TCL_OK) {
      opserr << G3_ERROR_PROMPT 
             << "rayleigh alphaM? betaK? betaK0? betaKc? - could not "
                "read betaK? \n";
      return TCL_ERROR;
    }

    for (int i = 0; i < numEigen; ++i)
      modalDampingValues[i] = factor;
  }

  // set factors in domain
  Domain *theDomain = builder->getDomain();
  assert(theDomain != nullptr);

  theDomain->setModalDampingFactors(&modalDampingValues, do_tangent);
  return TCL_OK;
}

static int
resetModel(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char ** const argv)
{
  assert(clientData != nullptr);
  BasicAnalysisBuilder *builder = (BasicAnalysisBuilder*)clientData;

  Domain* domain = builder->getDomain();
  assert(domain != nullptr);
  domain->revertToStart();

  TransientIntegrator *theTransientIntegrator =  builder->getTransientIntegrator();
  if (theTransientIntegrator != nullptr) {
    theTransientIntegrator->revertToStart();
  }
  return TCL_OK;
}


int
printIntegrator(ClientData clientData, Tcl_Interp *interp, int argc,
                TCL_Char ** const argv, OPS_Stream &output)
{
  assert(clientData != nullptr);
  BasicAnalysisBuilder *builder = (BasicAnalysisBuilder*)clientData;

  TransientIntegrator *theTransientIntegrator =  builder->getTransientIntegrator();
  StaticIntegrator *the_static_integrator = builder->getStaticIntegrator();

  int eleArg = 0;
  if (the_static_integrator == nullptr && theTransientIntegrator == nullptr)
    return TCL_OK;

  Integrator *theIntegrator;
  if (the_static_integrator != 0)
    theIntegrator = the_static_integrator;
  else
    theIntegrator = theTransientIntegrator;

  // if just 'print <filename> integrator'- no flag
  if (argc == 0) {
    theIntegrator->Print(output);
    return TCL_OK;
  }

  // if 'print <filename> Algorithm flag' get the flag
  int flag;
  if (Tcl_GetInt(interp, argv[eleArg], &flag) != TCL_OK) {
    opserr << G3_ERROR_PROMPT << "print algorithm failed to get integer flag: \n";
    opserr << argv[eleArg] << endln;
    return TCL_ERROR;
  }
  theIntegrator->Print(output, flag);
  return TCL_OK;
}

static int
printA(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char ** const argv)
{
  assert(clientData != nullptr);
  BasicAnalysisBuilder *builder = (BasicAnalysisBuilder*)clientData;

  int res = 0;

  FileStream outputFile;
  OPS_Stream *output = &opserr;
  LinearSOE  *oldSOE = builder->getLinearSOE();


  // Cant allocate theSolver on stack because theSOE is going to 
  // delete it
  FullGenLinLapackSolver *theSolver = new FullGenLinLapackSolver();
  FullGenLinSOE theSOE(*theSolver);

  builder->set(&theSOE, false);
  // invoke domainChange which constructs a graph and passes
  // it to the SOE. Otherwise, getA() returns null
  builder->domainChanged();


  bool ret = false;
  int currentArg = 1;
  double m = 0.0, c = 0.0, k = 0.0;
  bool do_mck = false;
  while (currentArg < argc) {
    if ((strcmp(argv[currentArg], "file") == 0) ||
        (strcmp(argv[currentArg], "-file") == 0)) {
      currentArg++;
      if (currentArg == argc) {
        opserr << G3_WARN_PROMPT << "-file missing argument\n";
        return TCL_ERROR;
      }

      if (outputFile.setFile(argv[currentArg]) != 0) {
        opserr << "printA <filename> .. - failed to open file: "
               << argv[currentArg] << endln;
        return TCL_ERROR;
      }
      output = &outputFile;
    } else if ((strcmp(argv[currentArg], "ret") == 0) ||
               (strcmp(argv[currentArg], "-ret") == 0)) {
      ret = true;

    } else if ((strcmp(argv[currentArg], "-m") == 0)) {
      currentArg++;
      if (Tcl_GetDouble(interp, argv[currentArg], &m) != TCL_OK) {
        opserr << G3_ERROR_PROMPT << "failed to read float following flag -m\n";
        return TCL_ERROR;
      }
      do_mck = true;

    } else if ((strcmp(argv[currentArg], "-c") == 0)) {
      currentArg++;
      if (Tcl_GetDouble(interp, argv[currentArg], &c) != TCL_OK) {
        opserr << G3_ERROR_PROMPT << "failed to read float following flag -c\n";
        return TCL_ERROR;
      }
      do_mck = true;

    } else if ((strcmp(argv[currentArg], "-k") == 0)) {
      currentArg++;
      if (Tcl_GetDouble(interp, argv[currentArg], &k) != TCL_OK) {
        opserr << G3_ERROR_PROMPT << "failed to read float following flag -k\n";
        return TCL_ERROR;
      }
      do_mck = true;
    }
    currentArg++;
  }
  
  //
  // Form the tangent
  //
  TransientIntegrator *oldint = nullptr;

  // construct integrator here so that it is not
  // destructed when the `if` scope ends
  GimmeMCK integrator(m, c, k, 0.0);
  if (do_mck) {
    oldint = builder->getTransientIntegrator();
    builder->set(integrator, false);
    integrator.formTangent(0);
    integrator.revertToLastStep();
  }

  else if (builder->getStaticIntegrator() != nullptr) {
    builder->getStaticIntegrator()->formTangent();
    builder->getStaticIntegrator()->revertToLastStep();
  }

  else if (builder->getTransientIntegrator() != nullptr) {
    builder->getTransientIntegrator()->formTangent(0);
    builder->getTransientIntegrator()->revertToLastStep();
  }
  builder->getDomain()->revertToLastCommit();

  const Matrix *A = theSOE.getA();
  if (A == nullptr) {
    opserr << OpenSees::PromptValueError 
           << "Could not get matrix from linear system\n";
    return TCL_ERROR;
  }

  if (ret) {
    int n = A->noRows();
    int m = A->noCols();
    if (n*m == 0) {
      opserr << OpenSees::PromptValueError 
             << "linear system is empty\n";
      return TCL_ERROR;
    }

    // Create an empty list with space preallocated for
    // n*m elements. This is not formally documented, but
    // it is mentioned here 
    //   https://wiki.tcl-lang.org/page/Tcl_NewListObj
    //
    // and evident from the source code here:
    //   https://github.com/enthought/tcl/blob/master/generic/tclListObj.c
    //
    Tcl_Obj* list = Tcl_NewListObj(n*m, nullptr);


    for (int i = 0; i < n; ++i) {
      for (int j = 0; j < m; j++)
        Tcl_ListObjAppendElement(interp, list, Tcl_NewDoubleObj((*A)(i, j)));
    }
    Tcl_SetObjResult(interp, list);

  } else {
    *output << *A;
    outputFile.close();
  }

  // put the original SOE back.
  if (oldSOE != nullptr)
    builder->set(oldSOE, true);

  if (oldint != nullptr)
    builder->set(*oldint, true);

  return res;
}

static int
printB(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char ** const argv)
{
  assert(clientData != nullptr);
  BasicAnalysisBuilder *builder = (BasicAnalysisBuilder*)clientData;

  int res = 0;

  FileStream outputFile;
  OPS_Stream *output = &opserr;

  bool ret = false;
  int currentArg = 1;
  while (currentArg < argc) {
    if ((strcmp(argv[currentArg], "file") == 0) ||
        (strcmp(argv[currentArg], "-file") == 0)) {
      currentArg++;
      if (currentArg == argc) {
        opserr << G3_WARN_PROMPT << "-file missing argument\n";
        return TCL_ERROR;
      }

      if (outputFile.setFile(argv[currentArg]) != 0) {
        opserr << "print <filename> .. - failed to open file: "
               << argv[currentArg] << endln;
        return TCL_ERROR;
      }
      output = &outputFile;
    } else if ((strcmp(argv[currentArg], "ret") == 0) ||
               (strcmp(argv[currentArg], "-ret") == 0)) {
      ret = true;
    }
    currentArg++;
  }

  LinearSOE  *theSOE = builder->getLinearSOE();
  if (theSOE != nullptr) {

    // TODO
    builder->formUnbalance();

    if (theSOE->getNumEqn() == 0) {
      opserr << OpenSees::PromptValueError << "System of equations is empty\n";
      return TCL_ERROR;
    }

    const Vector &b = theSOE->getB();

    if (ret) {
      const int size = b.Size();
      Tcl_Obj* list = Tcl_NewListObj(size, nullptr);
      for (int i = 0; i < size; ++i)
        Tcl_ListObjAppendElement(interp, list, Tcl_NewDoubleObj(b[i]));

      Tcl_SetObjResult(interp, list);
    } else {
      *output << b;
      outputFile.close();
    }
  }

  return res;
}

// This is removed in model.cpp
extern int
TclCommand_clearAnalysis(ClientData cd, Tcl_Interp *interp, int argc, TCL_Char ** const argv)
{

  if (cd != nullptr) {
    BasicAnalysisBuilder *builder = (BasicAnalysisBuilder *)cd;
    builder->wipe();
    delete builder;

    static int ncmd = sizeof(tcl_analysis_cmds)/sizeof(char_cmd);
    for (int i = 0; i < ncmd; ++i)
      Tcl_DeleteCommand(interp, tcl_analysis_cmds[i].name);

    Tcl_CreateCommand(interp, "wipeAnalysis",  &wipeAnalysis, nullptr, nullptr);
    Tcl_CreateCommand(interp, "_clearAnalysis", &TclCommand_clearAnalysis, nullptr, nullptr);
  }

  return TCL_OK;
}

static int
wipeAnalysis(ClientData cd, Tcl_Interp *interp, int argc, TCL_Char ** const argv)
{

  if (cd != nullptr) {
    BasicAnalysisBuilder *builder = (BasicAnalysisBuilder *)cd;
    builder->wipe();
  }
  return TCL_OK;
}


//
// command invoked to allow the ConstraintHandler object to be built
//
static int
specifyConstraintHandler(ClientData clientData, Tcl_Interp *interp, int argc,
                         TCL_Char ** const argv)
{
  
  BasicAnalysisBuilder *builder = (BasicAnalysisBuilder*)clientData;

  // make sure at least one other argument to contain type name
  if (argc < 2) {
    opserr << G3_ERROR_PROMPT << "need to specify a constraint type \n";
    return TCL_ERROR;
  }

  ConstraintHandler *theHandler = nullptr;
  // check argv[1] for type of handler and create the object
  if (strcmp(argv[1], "Plain") == 0)
    theHandler = new PlainHandler();

  else if (strcmp(argv[1], "Transformation") == 0) {
    theHandler = new TransformationConstraintHandler();
  }

  else if (strcmp(argv[1], "Penalty") == 0) {
    if (argc < 4) {
      opserr << "WARNING: need to specify alpha: handler Penalty alpha \n";
      return TCL_ERROR;
    }
    double alpha1, alpha2;
    if (Tcl_GetDouble(interp, argv[2], &alpha1) != TCL_OK)
      return TCL_ERROR;
    if (Tcl_GetDouble(interp, argv[3], &alpha2) != TCL_OK)
      return TCL_ERROR;
    theHandler = new PenaltyConstraintHandler(alpha1, alpha2);
  }

  else if (strcmp(argv[1], "Lagrange") == 0) {
    double alpha1 = 1.0;
    double alpha2 = 1.0;
    if (argc == 4) {
      if (Tcl_GetDouble(interp, argv[2], &alpha1) != TCL_OK)
        return TCL_ERROR;
      if (Tcl_GetDouble(interp, argv[3], &alpha2) != TCL_OK)
        return TCL_ERROR;
    }
    theHandler = new LagrangeConstraintHandler(alpha1, alpha2);
  }

  else {
    opserr << G3_ERROR_PROMPT << "ConstraintHandler type '" << argv[1]
      << "' does not exists \n\t(Plain, Penalty, Lagrange, Transformation) only\n";
    return TCL_ERROR;
  }

  builder->set(theHandler);
  return TCL_OK;
}
