//===----------------------------------------------------------------------===//
//
//                                   xara
//
//===----------------------------------------------------------------------===//
//                              https://xara.so
//===----------------------------------------------------------------------===//
#include <tcl.h>
#include <Logging.h>
#include <Parsing.h>
#include <Domain.h>
#include <LoadPattern.h>
#include <Parameter.h>
#include <ParameterIter.h>
#include <BasicAnalysisBuilder.h>


int 
computeGradients(ClientData clientData, Tcl_Interp* interp, int argc, TCL_Char**const argv)
{
  BasicAnalysisBuilder *builder = static_cast<BasicAnalysisBuilder*>(clientData);

  if (builder->analyzeGradient() < 0) {
    opserr << OpenSees::PromptValueError << "failed to compute sensitivities\n";
    return TCL_ERROR;
  }
  
  return TCL_OK;
}


int
TclCommand_sensLambda(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char ** const argv)
{
  BasicAnalysisBuilder* builder = static_cast<BasicAnalysisBuilder*>(clientData);

  if (argc < 3) {
    opserr << OpenSees::PromptValueError << "insufficient arguments\n";
    return TCL_ERROR;
  }

  int pattern, paramTag;
  if (Tcl_GetInt(interp, argv[1], &pattern) != TCL_OK) {
    opserr << "ERROR reading load pattern tag\n";
    return TCL_ERROR;
  }

  LoadPattern *thePattern = builder->getDomain()->getLoadPattern(pattern);
  if (thePattern == nullptr) {
    opserr << "ERROR load pattern with tag " << pattern
           << " not found in domain\n";
    return TCL_ERROR;
  }

  if (Tcl_GetInt(interp, argv[2], &paramTag) != TCL_OK) {
    opserr << OpenSees::PromptValueError 
           << "sensLambda patternTag?  paramTag? - could not read "
              "paramTag? ";
    return TCL_ERROR;
  }

  Parameter *theParam = builder->getDomain()->getParameter(paramTag);
  if (theParam == nullptr) {
    opserr << OpenSees::PromptValueError 
           << "sensLambda: parameter " << paramTag << " not found" << "\n";
    return TCL_ERROR;
  }

  int gradIndex = theParam->getGradIndex();
  double value = thePattern->getLoadFactorSensitivity(gradIndex);

  Tcl_SetObjResult(interp, Tcl_NewDoubleObj(value));

  return TCL_OK;
}


int
TclCommand_sensitivityAlgorithm(ClientData clientData, Tcl_Interp* interp, int argc, TCL_Char**const argv)
{
  BasicAnalysisBuilder *builder = static_cast<BasicAnalysisBuilder*>(clientData);

  if (argc < 2) {
    opserr << "ERROR: Wrong number of parameters to sensitivity algorithm." << "\n";
    return TCL_ERROR;
  }

  // 1: compute at each step (default); 
  // 2: compute by command; 

  int analysisTypeTag = 1;
  if (strcmp(argv[1],"-computeAtEachStep") == 0)
      analysisTypeTag = 1;

  else if (strcmp(argv[1],"-computeByCommand") == 0)
      analysisTypeTag = 2;

  else {
      opserr << "Unknown sensitivity algorithm option: " << argv[1] << "\n";
      return TCL_ERROR;
  }

  if (builder->setGradientType(analysisTypeTag) < 0) {
    return TCL_ERROR;
  }

  return TCL_OK;
}

