//===----------------------------------------------------------------------===//
//
//        OpenSees - Open System for Earthquake Engineering Simulation
//
//===----------------------------------------------------------------------===//
//
//
#include <tcl.h>
//
// NOTE: wipe is not added to the table on purpose, it cannot be added to
// and removed from the interpreter as simply as the other commands.
static Tcl_CmdProc wipeAnalysis;
//
static Tcl_CmdProc specifyAnalysis;
static Tcl_CmdProc eigenAnalysis;
static Tcl_CmdProc modalProperties;
static Tcl_CmdProc responseSpectrum;
static Tcl_CmdProc printA;
static Tcl_CmdProc printB;
static Tcl_CmdProc initializeAnalysis;
static Tcl_CmdProc resetModel;
static Tcl_CmdProc analyzeModel;
static Tcl_CmdProc specifyConstraintHandler;
static Tcl_CmdProc modalDamping;

// commands/analysis/integrator.cpp
extern Tcl_CmdProc specifyIntegrator;

// commands/analysis/solver.cpp
extern Tcl_CmdProc specifySOE;
extern Tcl_CmdProc specifySysOfEqnTable;

// commands/analysis/algorithm.cpp
extern Tcl_CmdProc TclCommand_specifyAlgorithm;
extern Tcl_CmdProc TclCommand_numIter;
extern Tcl_CmdProc TclCommand_accelCPU;
extern Tcl_CmdProc TclCommand_totalCPU;
extern Tcl_CmdProc TclCommand_solveCPU;
extern Tcl_CmdProc TclCommand_numFact;

// from commands/analysis/ctest.cpp
extern Tcl_CmdProc specifyCTest;
extern Tcl_CmdProc getCTestNorms;
extern Tcl_CmdProc getCTestIter;
extern Tcl_CmdProc TclCommand_algorithmRecorder;

struct char_cmd {
  const char* name;
  Tcl_CmdProc*  func;
}  const tcl_analysis_cmds[] =  {
    {"system",              &specifySysOfEqnTable},

    {"test",                &specifyCTest},
    {"testIter",            &getCTestIter},
    {"testNorms",           &getCTestNorms},
    {"integrator",          &specifyIntegrator},
    {"constraints",         &specifyConstraintHandler},

    {"eigen",               &eigenAnalysis},
    {"analysis",            &specifyAnalysis},

    {"analyze",             &analyzeModel},
    {"initialize",          &initializeAnalysis},
    {"modalProperties",     &modalProperties},
    {"modalDamping",        &modalDamping},
    {"modalDampingQ",       &modalDamping},
    {"responseSpectrum",    &responseSpectrum},
    {"printA",              &printA},
    {"printB",              &printB},
    {"reset",               &resetModel},

  // From algorithm.cpp
    {"algorithm",           &TclCommand_specifyAlgorithm},
    {"numIter",             &TclCommand_numIter},
    {"numFact",             &TclCommand_numFact},
    {"accelCPU",            &TclCommand_accelCPU},
    {"totalCPU",            &TclCommand_totalCPU},
    {"solveCPU",            &TclCommand_solveCPU},
  // recorder.cpp
    {"algorithmRecorder",   &TclCommand_algorithmRecorder},
};

