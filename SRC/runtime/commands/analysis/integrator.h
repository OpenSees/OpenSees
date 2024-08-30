//===----------------------------------------------------------------------===//
//
//        OpenSees - Open System for Earthquake Engineering Simulation
//
//===----------------------------------------------------------------------===//
//
#include <runtimeAPI.h>

OPS_Routine OPS_Newmark;
OPS_Routine OPS_StagedNewmark;
OPS_Routine OPS_GimmeMCK;
OPS_Routine OPS_AlphaOS;
OPS_Routine OPS_AlphaOS_TP;
OPS_Routine OPS_AlphaOSGeneralized;
OPS_Routine OPS_AlphaOSGeneralized_TP;
OPS_Routine OPS_ExplicitDifference;
OPS_Routine OPS_CentralDifference;
OPS_Routine OPS_CentralDifferenceAlternative;
OPS_Routine OPS_CentralDifferenceNoDamping;
OPS_Routine OPS_Collocation;
OPS_Routine OPS_CollocationHSFixedNumIter;
OPS_Routine OPS_CollocationHSIncrLimit;
OPS_Routine OPS_CollocationHSIncrReduct;
OPS_Routine OPS_GeneralizedAlpha;
OPS_Routine OPS_HHT;
OPS_Routine OPS_HHT_TP;
OPS_Routine OPS_HHTExplicit;
OPS_Routine OPS_HHTExplicit_TP;
OPS_Routine OPS_HHTGeneralized;
OPS_Routine OPS_HHTGeneralized_TP;
OPS_Routine OPS_HHTGeneralizedExplicit;
OPS_Routine OPS_HHTGeneralizedExplicit_TP;
OPS_Routine OPS_HHTHSFixedNumIter;
OPS_Routine OPS_HHTHSFixedNumIter_TP;
OPS_Routine OPS_HHTHSIncrLimit;
OPS_Routine OPS_HHTHSIncrLimit_TP;
OPS_Routine OPS_HHTHSIncrReduct;
OPS_Routine OPS_HHTHSIncrReduct_TP;
OPS_Routine OPS_KRAlphaExplicit;
OPS_Routine OPS_KRAlphaExplicit_TP;
OPS_Routine OPS_NewmarkExplicit;
OPS_Routine OPS_NewmarkHSFixedNumIter;
OPS_Routine OPS_NewmarkHSIncrLimit;
OPS_Routine OPS_NewmarkHSIncrReduct;
OPS_Routine OPS_WilsonTheta;

#include <tcl.h>
#include <string>
#include <unordered_map>
#include <runtimeAPI.h>
#include <BasicAnalysisBuilder.h>

#include <Domain.h>
#include <Node.h>

// integrators
#include <LoadControl.h>
// #include <StagedLoadControl.h>
#include <ArcLength1.h>
#include <DisplacementControl.h>

#include <TRBDF2.h>
#include <TRBDF3.h>
#include <Houbolt.h>
#include <ParkLMS3.h>
#include <BackwardEuler.h>


StaticIntegrator* G3Parse_newHSIntegrator(ClientData, Tcl_Interp*, int, TCL_Char ** const);
StaticIntegrator* G3Parse_newLoadControl(ClientData, Tcl_Interp*, int argc, TCL_Char *argv[]);
StaticIntegrator* G3Parse_newEQPathIntegrator(ClientData, Tcl_Interp*, int argc, TCL_Char ** const);
StaticIntegrator* G3Parse_newArcLengthIntegrator(ClientData, Tcl_Interp*, int argc, TCL_Char ** const);
StaticIntegrator* G3Parse_newMinUnbalDispNormIntegrator(ClientData, Tcl_Interp*, int, TCL_Char ** const);
StaticIntegrator* G3Parse_newDisplacementControlIntegrator(ClientData, Tcl_Interp*, int, TCL_Char** const);


TransientIntegrator* G3Parse_newNewmark1Integrator(ClientData, Tcl_Interp*, int, TCL_Char ** const);
TransientIntegrator* TclCommand_newNewmarkIntegrator(ClientData, Tcl_Interp*, int, TCL_Char** const);


//
// Helpers
//

// Type 1
template <typename Type, OPS_Routine fn>
static int
dispatch(ClientData clientData, Tcl_Interp* interp, int argc, G3_Char** const argv)
{
  BasicAnalysisBuilder *builder = static_cast<BasicAnalysisBuilder*>(clientData);
  Type* theIntegrator = (Type*)fn( G3_getRuntime(interp), argc, argv );

  if (theIntegrator == nullptr)
    return TCL_ERROR;


  opsdbg << G3_DEBUG_PROMPT << "Set integrator to \n";
  theIntegrator->Print(opsdbg);
  builder->set(*theIntegrator);
  return TCL_OK;
}

// Type 2
template <typename Type, Type*(*fn)(ClientData, Tcl_Interp*, int, TCL_Char** const)>
static int
dispatch(ClientData clientData, Tcl_Interp* interp, int argc, TCL_Char** const argv)
{
  assert(clientData != nullptr);


  Type* theIntegrator = fn( clientData, interp, argc, argv );

  if (theIntegrator == nullptr)
    return TCL_ERROR;

  BasicAnalysisBuilder *builder = static_cast<BasicAnalysisBuilder*>(clientData);

  opsdbg << G3_DEBUG_PROMPT << "Set integrator to \n";
  theIntegrator->Print(opsdbg);
  builder->set(*theIntegrator);
  return TCL_OK;
}

// Type 3
template <int (*fn)(ClientData clientData, Tcl_Interp* interp, int, G3_Char** const)> 
static int
dispatch(ClientData clientData, Tcl_Interp* interp, int argc, G3_Char** const argv)
{
  assert(clientData != nullptr);
  return fn( clientData, interp, argc, argv );
}

#define DISPATCH(Type, Class)                                         \
  (Tcl_CmdProc*)[](ClientData clientData, Tcl_Interp*, int, G3_Char**const)->int{ \
    BasicAnalysisBuilder *builder = static_cast<BasicAnalysisBuilder*>(clientData); \
    Type* theIntegrator = new Class();                                \
    opsdbg << G3_DEBUG_PROMPT << "Set integrator to \n";              \
    theIntegrator->Print(opsdbg);                                     \
    builder->set(*theIntegrator);                                     \
    return TCL_OK;                                                    \
  }

//
//
//
#if 0
std::unordered_map<std::string, StaticIntegrator* (*)(ClientData, Tcl_Interp*, int argc, TCL_Char ** const)> 
StaticIntegratorLibrary = {
  {"LoadControl",                  G3Parse_newLoadControl},
//{"StagedLoadControl",            G3Parse_newStagedLoadControlIntegrator},
  {"EQPath",                       G3Parse_newEQPathIntegrator},
  {"ArcLength",                    G3Parse_newArcLengthIntegrator},
  {"MinUnbalDispNorm",             G3Parse_newMinUnbalDispNormIntegrator},
  {"DisplacementControl",          G3Parse_newDisplacementControlIntegrator},
};
#endif


std::unordered_map<std::string, Tcl_CmdProc*> 
StaticIntegratorLibrary = {
  {"LoadControl",         dispatch<StaticIntegrator, G3Parse_newLoadControl>},

  {"EQPath",              dispatch<StaticIntegrator, G3Parse_newEQPathIntegrator>},

  {"MinUnbalDispNorm",    dispatch<StaticIntegrator, G3Parse_newMinUnbalDispNormIntegrator>},

  {"DisplacementControl", dispatch<StaticIntegrator, G3Parse_newDisplacementControlIntegrator>},

  {"ArcLength",           dispatch<StaticIntegrator, G3Parse_newArcLengthIntegrator>},
};

std::unordered_map<std::string, Tcl_CmdProc*> 
TransientIntegratorLibrary = {
  // TRBDF2
  {"TRBDF2",          DISPATCH(TransientIntegrator, TRBDF2)},
  {"TRBDF2",          DISPATCH(TransientIntegrator, TRBDF2)},
  {"Bathe",           DISPATCH(TransientIntegrator, TRBDF2)},

  // TRBDF3
  {"TRBDF3",          DISPATCH(TransientIntegrator, TRBDF3)},
  {"Bathe3",          DISPATCH(TransientIntegrator, TRBDF3)},

  {"Houbolt",         DISPATCH(TransientIntegrator, Houbolt)},

  {"ParkLMS3",        DISPATCH(TransientIntegrator, ParkLMS3)},

  // MCK
  {"GimmeMCK",                dispatch<TransientIntegrator, OPS_GimmeMCK> },
  {"MCK",                     dispatch<TransientIntegrator, OPS_GimmeMCK> },
  {"ZZTop",                   dispatch<TransientIntegrator, OPS_GimmeMCK> },

  {"Newmark",                 dispatch<TransientIntegrator, TclCommand_newNewmarkIntegrator> },

  {"NewmarkExplicit",         dispatch<TransientIntegrator, OPS_NewmarkExplicit>},

  {"Newmark1",                dispatch<TransientIntegrator, G3Parse_newNewmark1Integrator>},

  {"NewmarkHSIncrReduct",     dispatch<TransientIntegrator, OPS_NewmarkHSIncrReduct>},

  {"NewmarkHSIncrLimit",      dispatch<TransientIntegrator, OPS_NewmarkHSIncrLimit>},

  {"NewmarkHSFixedNumIter",   dispatch<TransientIntegrator, OPS_NewmarkHSFixedNumIter>},

  {"HHT",                     dispatch<TransientIntegrator, OPS_HHT>},

  {"HHT_TP",                  dispatch<TransientIntegrator, OPS_HHT_TP>},

  {"HHTGeneralized",          dispatch<TransientIntegrator, OPS_HHTGeneralized>},

  {"HHTGeneralized_TP",       dispatch<TransientIntegrator, OPS_HHTGeneralized_TP>},

  {"HHTExplicit",             dispatch<TransientIntegrator, OPS_HHTExplicit>},

  {"HHTExplicit_TP",          dispatch<TransientIntegrator, OPS_HHTExplicit_TP>},

  {"HHTGeneralizedExplicit",  dispatch<TransientIntegrator, OPS_HHTGeneralizedExplicit>},

  {"HHTGeneralizedExplicit_TP", dispatch<TransientIntegrator, OPS_HHTGeneralizedExplicit_TP>},

  {"HHTHSIncrLimit",            dispatch<TransientIntegrator, OPS_HHTHSIncrLimit>},

  {"HHTHSIncrLimit_TP",         dispatch<TransientIntegrator, OPS_HHTHSIncrLimit_TP>},

  {"HHTHSIncrReduct",         dispatch<TransientIntegrator, OPS_HHTHSIncrReduct>},

  {"HHTHSIncrReduct_TP",      dispatch<TransientIntegrator, OPS_HHTHSIncrReduct_TP>},

  {"HHTHSFixedNumIter",       dispatch<TransientIntegrator, OPS_HHTHSFixedNumIter>},

  {"HHTHSFixedNumIter_TP",    dispatch<TransientIntegrator, OPS_HHTHSFixedNumIter_TP>},

  {"GeneralizedAlpha",        dispatch<TransientIntegrator, OPS_GeneralizedAlpha>},

  {"KRAlphaExplicit",         dispatch<TransientIntegrator, OPS_KRAlphaExplicit>},

  {"KRAlphaExplicit_TP",      dispatch<TransientIntegrator, OPS_KRAlphaExplicit_TP>},

  {"AlphaOS",                 dispatch<TransientIntegrator, OPS_AlphaOS>},

  {"AlphaOS_TP",              dispatch<TransientIntegrator, OPS_AlphaOS_TP>},

  {"AlphaOSGeneralized",      dispatch<TransientIntegrator, OPS_AlphaOSGeneralized>},

  {"AlphaOSGeneralized_TP",   dispatch<TransientIntegrator, OPS_AlphaOSGeneralized_TP>},

  {"Collocation",             dispatch<TransientIntegrator, OPS_Collocation>},

  {"CollocationHSIncrReduct", dispatch<TransientIntegrator, OPS_CollocationHSIncrReduct>},

  {"CollocationHSIncrLimit",  dispatch<TransientIntegrator, OPS_CollocationHSIncrLimit>},

  {"CollocationHSFixedNumIter",    dispatch<TransientIntegrator, OPS_CollocationHSFixedNumIter>},


  {"WilsonTheta",                  dispatch<TransientIntegrator, OPS_WilsonTheta>},

  {"ExplicitDifference",           dispatch<TransientIntegrator, OPS_ExplicitDifference>},

  {"CentralDifference",            dispatch<TransientIntegrator, OPS_CentralDifference>},

  {"CentralDifferenceAlternative", dispatch<TransientIntegrator, OPS_CentralDifferenceAlternative>},

  {"CentralDifferenceNoDamping",   dispatch<TransientIntegrator, OPS_CentralDifferenceNoDamping>},

};

