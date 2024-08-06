//===----------------------------------------------------------------------===//
//
//        OpenSees - Open System for Earthquake Engineering Simulation
//
//===----------------------------------------------------------------------===//
//
// Description: This file implements commands that configure the linear
// solver.
//
#include <string>
#include <algorithm>
#ifdef _MSC_VER 
#  include <string.h>
#  define strcasecmp _stricmp
#else
#  include <strings.h>
#endif
// #define strcmp strcasecmp

#include <tcl.h>
#include <G3_Logging.h>
#include <runtimeAPI.h>
// #include "analysis.h"
#include <OPS_Globals.h>
#include "solver.hpp"
#include "BasicAnalysisBuilder.h"

// analysis
#include <StaticAnalysis.h>
#include <DirectIntegrationAnalysis.h>
#include <VariableTimeStepDirectIntegrationAnalysis.h>

// system of eqn and solvers
#include <SProfileSPDLinSolver.h>
#include <SProfileSPDLinSOE.h>
#include <ProfileSPDLinDirectThreadSolver.h>
#include <SparseGenColLinSOE.h>
#include <SparseGenRowLinSOE.h>
#include <SymSparseLinSOE.h>
#include <SymSparseLinSolver.h>

#ifdef _CUDA
#  include <BandGenLinSOE_Single.h>
#  include <BandGenLinLapackSolver_Single.h>
#endif

#ifdef _CULAS4
#  include <CulaSparseSolverS4.h>
#endif

#ifdef _CULAS5
#  include <CulaSparseSolverS5.h>
#endif

#if defined(_PARALLEL_PROCESSING)
//  parallel soe & solvers
#  include <DistributedBandSPDLinSOE.h>
#  include <DistributedSparseGenColLinSOE.h>
#  include <DistributedSparseGenRowLinSOE.h>
#  include <DistributedBandGenLinSOE.h>
#  include <DistributedDiagonalSOE.h>
#  include <DistributedDiagonalSolver.h>

#  include <MPIDiagonalSOE.h>
#  include <MPIDiagonalSolver.h>

#  define MPIPP_H
#  include <DistributedSuperLU.h>
#  include <DistributedProfileSPDLinSOE.h>
#endif

// TODO: remove
// extern DirectIntegrationAnalysis *theTransientAnalysis;
// extern LinearSOE *theSOE;

// LinearSOE*
// G3Parse_newLinearSOE(G3_Runtime*, int, G3_Char **);
LinearSOE*
// G3Parse_newLinearSOE(G3_Runtime* rt, int argc, G3_Char ** const argv)
G3Parse_newLinearSOE(ClientData, Tcl_Interp* interp, int, G3_Char **const);

LinearSOE*
TclDispatch_newPetscSOE(ClientData, Tcl_Interp *interp, int, G3_Char **const);

#if 0 // TODO: implement AnalysisBuilder->getLinearSOE();
int
TclCommand_systemSize(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char ** const argv)
{
  assert(clientData != nullptr);
  LinearSOE *theSOE = ((BasicAnalysisBuilder *)clientData)->getLinearSOE();

  char buffer[20];

  if (theSOE == 0) {
    sprintf(buffer, "NO SYSTEM SET");
    return TCL_OK;
  }

  sprintf(buffer, "%d", theSOE->getNumEqn());
  Tcl_SetResult(interp, buffer, TCL_VOLATILE);

  return TCL_OK;
}
#endif

int
specifySysOfEqnTable(ClientData clientData, Tcl_Interp *interp, int argc, G3_Char ** const argv)
{
  // make sure at least one other argument to contain type of system
  if (argc < 2) {
    opserr << G3_ERROR_PROMPT
           << "need to specify a system type" << "\n";
    return TCL_ERROR;
  }

  LinearSOE* theSOE = G3Parse_newLinearSOE(clientData, interp, argc, argv);

  if (theSOE == nullptr)
    return TCL_ERROR;

  BasicAnalysisBuilder* builder = (BasicAnalysisBuilder*)clientData;

  builder->set(theSOE);
  return TCL_OK;

}

LinearSOE*
G3Parse_newLinearSOE(ClientData clientData, Tcl_Interp* interp, int argc, G3_Char ** const argv)
{
  G3_Runtime* rt = G3_getRuntime(interp); 

  // argc is checked by calling function;

  LinearSOE *theSOE = nullptr;
  std::string sys_name{argv[1]};
  transform(sys_name.begin(), sys_name.end(), sys_name.begin(), ::tolower);
  auto ctor = soe_table.find(sys_name);

  if (ctor != soe_table.end()) {
    return ctor->second.ss(rt, argc, argv);

  } else if (strcasecmp(argv[1], "Umfpack")==0) {
    // TODO: if "umfpack" is in solver.hpp, this wont be reached
    return TclDispatch_newUmfpackLinearSOE(clientData, interp, argc, argv);
  } 
#if 0
  else if (strcmp(argv[2],"Thread") == 0) {  
      int blockSize = 4;
      int numThreads = 1;
      if (argc == 5) {
	if (Tcl_GetInt(interp, argv[3], &blockSize) != TCL_OK)
	  return nullptr; //TCL_ERROR;
	if (Tcl_GetInt(interp, argv[4], &numThreads) != TCL_OK)
	  return nullptr; //TCL_ERROR;
      }
      return new ProfileSPDLinSOE(
          *new ProfileSPDLinDirectThreadSolver(numThreads,blockSize,1.0e-12)
      );
  } 
#endif

#if defined(OPS_PETSC)
  else if (strcmp(argv[1], "petsc")==0 ||
           strcmp(argv[1], "Petsc")==0) {
    theSOE = TclDispatch_newPetscSOE(clientData, interp, argc, argv);
  }
#endif

  // Diagonal SOE & SOLVER
#ifdef _PARALLEL_INTERPRETERS
  else if (strcmp(argv[1], "ParallelProfileSPD") == 0) {
    ProfileSPDLinSolver *theSolver = new ProfileSPDLinDirectSolver();
    DistributedProfileSPDLinSOE *theParallelSOE =
        new DistributedProfileSPDLinSOE(*theSolver);
    theSOE = theParallelSOE;
    theParallelSOE->setProcessID(OPS_rank);
    theParallelSOE->setChannels(numChannels, theChannels);
  }
#endif

#ifdef _PARALLEL_INTERPRETERS
  else if (strcmp(argv[1], "MPIDiagonal") == 0) {
    setMPIDSOEFlag = true;
  }
#endif

  else {
    opserr << G3_ERROR_PROMPT << " system '" << argv[1] << "' is unknown or not installed\n";
    return nullptr;
  }

  return theSOE;
}


LinearSOE*
specify_SparseSPD(G3_Runtime *rt, int argc, G3_Char ** const argv)
{
//if ((strcmp(argv[1], "SparseSPD") == 0) ||
//         (strcmp(argv[1], "SparseSYM") == 0)) {
    Tcl_Interp *interp = G3_getInterpreter(rt);

    // determine ordering scheme
    //   1 -- MMD
    //   2 -- ND
    //   3 -- RCM

    int lSparse = 1;
    if (argc == 3) {
      if (Tcl_GetInt(interp, argv[2], &lSparse) != TCL_OK)
        return nullptr;
    }
    SymSparseLinSolver *theSolver = new SymSparseLinSolver();
    return new SymSparseLinSOE(*theSolver, lSparse);
}


#ifdef _THREADS
#  include "contrib/sys_of_eqn/ThreadedSuperLU/ThreadedSuperLU.h"
#else
#  include <SuperLU.h>
#endif
// TODO(cmp): Threaded SuperLU?

LinearSOE*
specifySparseGen(G3_Runtime* rt, int argc, G3_Char ** const argv)
{
  // SPARSE GENERAL SOE * SOLVER
//if ((strcmp(argv[1], "SparseGeneral") == 0) ||
//         (strcmp(argv[1], "SuperLU") == 0) ||
//         (strcmp(argv[1], "SparseGEN") == 0))
    Tcl_Interp *interp = G3_getInterpreter(rt);

    SparseGenColLinSolver *theSolver = nullptr;
    int count = 2;
    double thresh = 0.0;
    int npRow = 1;
    int npCol = 1;
    int np = 1;

    // defaults for threaded SuperLU
    while (count < argc) {
      if ((strcmp(argv[count], "p")    == 0) ||
          (strcmp(argv[count], "piv")  == 0) ||
          (strcmp(argv[count], "-piv") == 0)) {
        thresh = 1.0;
      } else if ((strcmp(argv[count], "-np") == 0) ||
                 (strcmp(argv[count], "np")  == 0)) {
        count++;
        if (count < argc)
          if (Tcl_GetInt(interp, argv[count], &np) != TCL_OK)
            return nullptr;
      } else if ((strcmp(argv[count], "npRow")  == 0) ||
                 (strcmp(argv[count], "-npRow") == 0)) {
        count++;
        if (count < argc)
          if (Tcl_GetInt(interp, argv[count], &npRow) != TCL_OK)
            return nullptr;
      } else if ((strcmp(argv[count], "npCol")  == 0) ||
                 (strcmp(argv[count], "-npCol") == 0)) {
        count++;
        if (count < argc)
          if (Tcl_GetInt(interp, argv[count], &npCol) != TCL_OK)
            return nullptr;
      }
      count++;
    }

    int permSpec = 0;
    int panelSize = 6;
    int relax = 6;

#ifdef _THREADS
    if (np != 0)
      theSolver = new ThreadedSuperLU(np, permSpec, panelSize, relax, thresh);
    else
      return nullptr;
// #endif
// 
// #ifdef _PARALLEL_PROCESSING
//     if (theSolver != 0)
//       delete theSolver;
//     theSolver = 0;
// 
//     if (npRow != 0 && npCol != 0) {
//       theSolver = new DistributedSuperLU(npRow, npCol);
//     }
#else
    char symmetric = 'N';
    double drop_tol = 0.0;
    while (count < argc) {
      if (strcmp(argv[count], "s") == 0 || strcmp(argv[count], "symmetric") ||
          strcmp(argv[count], "-symm")) {
        symmetric = 'Y';
      }
      count++;
    }
    // TODO(cmp) : SuperLU
    theSolver = new SuperLU(permSpec, drop_tol, panelSize, relax, symmetric);
#endif

#ifdef _PARALLEL_PROCESSING
    return new DistributedSparseGenColLinSOE(*theSolver);
#else
    return new SparseGenColLinSOE(*theSolver);
#endif
}


#if 0 // Some misc solvers i play with

else if (strcmp(argv[2],"Block") == 0) {
  int blockSize = 4;
  if (argc == 4) {
    if (Tcl_GetInt(interp, argv[3], &blockSize) != TCL_OK)
      return TCL_ERROR;
  }
  theSolver = theSolver = new ProfileSPDLinDirectBlockSolver(1.0e-12,blockSize);
}


  int blockSize = 4;
  int numThreads = 1;
  if (argc == 5) {
    if (Tcl_GetInt(interp, argv[3], &blockSize) != TCL_OK)
      return TCL_ERROR;
    if (Tcl_GetInt(interp, argv[4], &numThreads) != TCL_OK)
      return TCL_ERROR;
  }
  theSolver = new ProfileSPDLinDirectThreadSolver(numThreads,blockSize,1.0e-12); 

  } else if (strcmp(argv[2],"Thread") == 0) { 
    int blockSize = 4; 
    int numThreads = 1; 
    if (argc == 5) { 
      if (Tcl_GetInt(interp, argv[3], &blockSize) != TCL_OK) 
        return TCL_ERROR;
      if (Tcl_GetInt(interp, argv[4], &numThreads) != TCL_OK)
        return TCL_ERROR;
  }
  theSolver = new ProfileSPDLinDirectThreadSolver(numThreads,blockSize,1.0e-12);
}

else if (strcmp(argv[2],"Skypack") == 0) {
  if (argc == 5) {
    int mCols, mRows;
    if (Tcl_GetInt(interp, argv[3], &mCols) != TCL_OK)
      return TCL_ERROR;
    if (Tcl_GetInt(interp, argv[4], &mRows) != TCL_OK)
      return TCL_ERROR;
    theSolver = new ProfileSPDLinDirectSkypackSolver(mCols, mRows);
  } else
    theSolver = new ProfileSPDLinDirectSkypackSolver();
}
else
  theSolver = new ProfileSPDLinDirectSolver();

#endif // misc solvers


#if defined(_CULAS4) || defined(_CULAS5)
// CULA SPARSE
else if ((strcmp(argv[1], "CulaSparse") == 0)) {
  double absTol = 1.0e-6;
  double relTol = 1e-6;

  int maxInteration = 100000;

  int preCond = 5; // fainv
#  ifdef _CULAS4
  preCond = 1;
#  endif
  int solver = 0; // cg
  int count = 2;
  int single = 0;
  int host = 0;

  while (count < argc) {

    if (strcmp(argv[count], "-rTol") == 0) {
      count++;
      if (count < argc)
        if (Tcl_GetDouble(interp, argv[count], &relTol) != TCL_OK)
          return TCL_ERROR;
    } else if ((strcmp(argv[count], "-mInt") == 0)) {
      count++;
      if (count < argc)
        if (Tcl_GetInt(interp, argv[count], &maxInteration) != TCL_OK)
          return TCL_ERROR;
    } else if ((strcmp(argv[count], "-pre") == 0)) {
      count++;
      if (count < argc)
        if ((strcmp(argv[count], "none") == 0))
          preCond = 0;
        else if ((strcmp(argv[count], "jacobi") == 0))
          preCond = 1;
        else if ((strcmp(argv[count], "blockjacobi") == 0))
          preCond = 2;
        else if ((strcmp(argv[count], "ilu0") == 0))
          preCond = 3;
        else if ((strcmp(argv[count], "ainv") == 0))
          preCond = 4;
        else if ((strcmp(argv[count], "fainv") == 0))
          preCond = 5;
        else
          return TCL_ERROR;
    } else if ((strcmp(argv[count], "-solver") == 0)) {
      count++;
      if (count < argc)
        if ((strcmp(argv[count], "cg") == 0))
          solver = 0;
        else if ((strcmp(argv[count], "bicg") == 0))
          solver = 1;
        else if ((strcmp(argv[count], "blockstab") == 0))
          solver = 2;
        else if ((strcmp(argv[count], "blockstabl") == 0))
          solver = 3;
        else if ((strcmp(argv[count], "gmres") == 0))
          solver = 4;
        else if ((strcmp(argv[count], "minres") == 0))
          solver = 5;
        else
          return TCL_ERROR;
    } else if ((strcmp(argv[count], "-single") == 0)) {
      single = 1;
    } else if ((strcmp(argv[count], "-host") == 0)) {
      host = 1;
    }
    count++;
  }

#  ifdef _CULAS5
  CulaSparseSolverS5 *theSolver = new CulaSparseSolverS5(
      relTol, maxInteration, preCond, solver, single, host);
#  else
  CulaSparseSolverS4 *theSolver =
      new CulaSparseSolverS4(relTol, maxInteration, preCond, solver);
#  endif
  theSOE = new SparseGenRowLinSOE(*theSolver);
}
#endif
