//===----------------------------------------------------------------------===//
//
//        OpenSees - Open System for Earthquake Engineering Simulation    
//
//===----------------------------------------------------------------------===//
//
// Written: cmp and cc
//
#include <tcl.h>
#include <runtimeAPI.h>
#include <G3_Logging.h>

#include <Newmark1.h>
#include <Newmark.h>
#include <GeneralizedNewmark.h>

TransientIntegrator*
TclCommand_newNewmarkIntegrator(ClientData clientData, Tcl_Interp* interp, int argc, G3_Char ** const argv)
{

  if (argc < 4) {
    opserr << G3_ERROR_PROMPT << " incorrect number of args want Newmark $gamma $beta "
              "<-form $typeUnknown>\n";
    opserr << "        got ";
    for (int i=0; i<argc; i++)
      opserr << argv[i] << " ";
    opserr << "\n";
    return 0;
  }

  int dispFlag =  1;
  int initFlag =  1;
  double gamma, beta;
  double alphaM = 1.0,
         alphaF = 1.0;
  bool useGeneralized = false;

  // Keep track of required arguments
  bool gotBeta  = false,
       gotGamma = false;


  for (int argi=2; argi < argc; argi++) {
    const char *nextString = argv[argi];
    //
    // Keyword arguments
    //
    if ((strcmp(argv[argi], "-form") == 0) ||
        (strcmp(argv[argi], "-solve") == 0)) {
      if (argi == argc-1) {
        opserr << G3_ERROR_PROMPT << "form option requires an argument\n";
        return nullptr;
      } else {
        nextString = argv[++argi];
      }

      if ((nextString[0] == 'D') || (nextString[0] == 'd') ||
          (nextString[0] == 'U') || (nextString[0] == 'u'))
        dispFlag = 1;
      else if ((nextString[0] == 'V') || (nextString[0] == 'v'))
        dispFlag = 2;
      else if ((nextString[0] == 'A') || (nextString[0] == 'a'))
        dispFlag = 3;
      else {
        opserr << G3_ERROR_PROMPT << "invalid argument for parameter 'form'\n";
        return nullptr;
      }
    }
    else if ((strcmp(argv[argi], "-init") == 0)) {
      if (argi == argc-1) {
        opserr << G3_ERROR_PROMPT << "init option requires an argument\n";
        return nullptr;
      } else {
        nextString = argv[++argi];
      }
      if ((nextString[0] == 'D') || (nextString[0] == 'd') ||
          (nextString[0] == 'U') || (nextString[0] == 'u'))
        initFlag = 1;
      else if ((nextString[0] == 'V') || (nextString[0] == 'v'))
        initFlag = 2;
      else if ((nextString[0] == 'A') || (nextString[0] == 'a'))
        initFlag = 3;
      else {
        opserr << G3_ERROR_PROMPT << "invalid argument for parameter 'init'\n";
        return nullptr;
      }
    }
    else if ((strcmp(argv[argi], "-alpha") == 0)) {
      if (argi == argc-1) {
        opserr << G3_ERROR_PROMPT << "alpha option requires an argument\n";
        return nullptr;
      } else {
        nextString = argv[++argi];
      }
      useGeneralized = true;
      if (Tcl_GetDouble(interp, argv[argi], &alphaF) != TCL_OK) {
        opserr << G3_ERROR_PROMPT << "invalid arg at position '" << argi << "\n";
        return nullptr;
      }
    }
    //
    // Positional arguments
    //
    else if ((!gotGamma && strcmp(argv[argi], "-beta") != 0) || 
             strcmp(argv[argi], "-gamma") == 0) {
      // if given with flag (not by position), increment argi
      if (strcmp(argv[argi], "-gamma") == 0)
        argi++;

      gotGamma = true;

      if (Tcl_GetDouble(interp, argv[argi], &gamma) != TCL_OK) {
        opserr << G3_ERROR_PROMPT << "invalid value for gamma. Expected a number ";
        opserr << "but got '" << argv[2] << "'.\n";
        return nullptr;
      }
    }

    else if (!gotBeta || strcmp(argv[argi], "-beta") == 0) {
      // if given with flag (not by position), increment argi
      if (strcmp(argv[argi], "-beta") == 0)
        argi++;

      gotBeta = true;

      if (Tcl_GetDouble(interp, argv[argi], &beta) != TCL_OK) {
        opserr << G3_ERROR_PROMPT << "invalid value for beta. Expected a number ";
        opserr << "but got '" << argv[argi] << "'.\n";
        return nullptr;
      }
    }
    else {
      opserr << G3_ERROR_PROMPT << "unrecognized argument " << nextString << "\n";
      return nullptr;
    }
  }

  // Check that all required arguments were supplied
  if (!gotGamma) { 
    opserr << G3_ERROR_PROMPT << "missing required positional argument gamma\n";
    return nullptr;
  }

  if (!gotBeta) { 
    opserr << G3_ERROR_PROMPT << "missing required positional argument beta\n";
    return nullptr;
  }

  if (useGeneralized)
    return new GeneralizedNewmark(gamma, beta, alphaF, alphaM, dispFlag, initFlag);
  else
    return new Newmark(gamma, beta, dispFlag, initFlag);

  return nullptr;
}

TransientIntegrator*
G3Parse_newNewmark1Integrator(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char ** const argv)
{
    double gamma;
    double beta;
    double alphaM, betaK, betaKi, betaKc;
    if (argc != 4 && argc != 8) {
      opserr << "WARNING integrator Newmark1 gamma beta <alphaM> "
                "<betaKcurrent> <betaKi> <betaKlastCommitted>\n";
      return nullptr;
    }
    if (Tcl_GetDouble(interp, argv[2], &gamma) != TCL_OK) {
      opserr << "WARNING integrator Newmark1 gamma beta - undefined gamma\n";
      return nullptr;
    }
    if (Tcl_GetDouble(interp, argv[3], &beta) != TCL_OK) {
      opserr << "WARNING integrator Newmark1 gamma beta - undefined beta\n";
      return nullptr;
    }

    if (argc == 8 || argc == 7) {
      if (Tcl_GetDouble(interp, argv[4], &alphaM) != TCL_OK) {
        opserr << "WARNING integrator Newmark1 gamma beta alphaM betaK betaKi "
                  "betaKc - alphaM\n";
        return nullptr;
      }
      if (Tcl_GetDouble(interp, argv[5], &betaK) != TCL_OK) {
        opserr << "WARNING integrator Newmark1 gamma beta alphaM betaK betaKi "
                  "betaKc - betaK\n";
        return nullptr;
      }
      if (Tcl_GetDouble(interp, argv[6], &betaKi) != TCL_OK) {
        opserr << "WARNING integrator Newmark1 gamma beta alphaM betaK betaKi "
                  "betaKc - betaKi\n";
        return nullptr;
      }
      if (Tcl_GetDouble(interp, argv[7], &betaKc) != TCL_OK) {
        opserr << "WARNING integrator Newmark1 gamma beta alphaM betaK betaKi "
                  "betaKc - betaKc\n";
        return nullptr;
      }
    }
    if (argc == 4)
      return new Newmark1(gamma, beta);
    else
      return new Newmark1(gamma, beta, alphaM, betaK, betaKi, betaKc);
}

