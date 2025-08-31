//===----------------------------------------------------------------------===//
//
//                                   xara
//
//===----------------------------------------------------------------------===//
//                              https://xara.so
//===----------------------------------------------------------------------===//
// Description: This file implements commands that configure a 
// `ModelBuider`, including "model"
//
// Author: cmp
//
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <tcl.h>
#include <Logging.h>
#include <Parsing.h>
#include <runtimeAPI.h>
#include <Domain.h>
#include <FE_Datastore.h>

#include "BasicModelBuilder.h"

#ifdef _PARALLEL_PROCESSING
#  include <PartitionedDomain.h>
   extern PartitionedDomain theDomain;
#endif


bool builtModel = false;

FE_Datastore *theDatabase = nullptr;

extern int G3_AddTclAnalysisAPI(Tcl_Interp *, BasicModelBuilder&);
extern int G3_AddTclDomainCommands(Tcl_Interp *, Domain*);


// 
int
TclCommand_specifyModel(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char *argv[])
{
  G3_Runtime *rt = G3_getRuntime(interp);
  Domain *theNewDomain = (Domain*)clientData;

  BasicModelBuilder *theNewBuilder = nullptr;

  //
  //
  //
  if (clientData == nullptr) {
    theNewDomain = new Domain();

    // TODO: remove ops_TheActiveDomain
    ops_TheActiveDomain = theNewDomain;

    Tcl_CreateCommand(interp, "model", &TclCommand_specifyModel, theNewDomain, nullptr);

    G3_AddTclDomainCommands(interp, theNewDomain);
  }


  // make sure at least one other argument to contain model builder type given
  if (argc < 2) {
    opserr << OpenSees::PromptValueError << "need to specify a model type, valid types:\n";
    opserr << "\tBasicBuilder\n";
    return TCL_ERROR;
  }

  // check argv[1] for type of ModelBuilder and create the object
  if ((strcmp(argv[1], "basic") == 0)        ||
      (strcmp(argv[1], "Basic") == 0)        ||
      (strcmp(argv[1], "-ndm") == 0)         ||
      (strcmp(argv[1], "BasicBuilder") == 0) ||
      (strcmp(argv[1], "basicBuilder") == 0)) {

    if (argc < 3) {
      opserr << OpenSees::PromptValueError 
             << "incorrect number of arguments\n";
      return TCL_ERROR;
    }
    int ndm = 0;
    int ndf = 0;

    int posArg = 1; // track positional argument
    int argPos = 2;
    while (argPos < argc) {
      if (strcmp(argv[argPos], "-ndm") == 0 ||
          strcmp(argv[argPos], "-NDM") == 0) {
        argPos++;
        if (argPos < argc) {
          if (Tcl_GetInt(interp, argv[argPos], &ndm) != TCL_OK) {
            opserr << OpenSees::PromptValueError 
                   << "error reading ndm, got '" << argv[argPos] << "'\n";
            return TCL_ERROR;
          }
        }
        argPos++;
        posArg++;

      } else if (strcmp(argv[argPos], "-ndf") == 0 ||
                 strcmp(argv[argPos], "-NDF") == 0) {
        argPos++;
        if (argPos < argc)
          if (Tcl_GetInt(interp, argv[argPos], &ndf) != TCL_OK) {
            opserr << OpenSees::PromptValueError << "invalid parameter ndf";
            return TCL_ERROR;
          }
        argPos++;
        posArg++;

      } else if (posArg == 1) {
        if (Tcl_GetInt(interp, argv[argPos], &ndm) != TCL_OK) {
          opserr << OpenSees::PromptValueError << "invalid parameter ndm";
          return TCL_ERROR;
        }
        argPos++;
        posArg++;

      } else if (posArg == 2) {
        if (Tcl_GetInt(interp, argv[argPos], &ndf) != TCL_OK) {
          opserr << OpenSees::PromptValueError << "error reading ndf: " << argv[argPos] << "\n";
          return TCL_ERROR;
        }
        argPos++;
        posArg++;

      } else {
        // no matches, advance to next argument
        argPos++;
      }
    }

    // check that ndm was specified
    if (ndm == 0) {
      opserr << OpenSees::PromptValueError << "missing required argument ndm\n";
      return TCL_ERROR;
    }

    // check for ndf, if not assume one
    if (ndf == 0) {
      if (ndm == 1)
        ndf = 1;
      else if (ndm == 2)
        ndf = 3;
      else if (ndm == 3)
        ndf = 6;
      else {
        opserr << OpenSees::PromptValueError << "specified ndm, " << ndm << ", will not work\n";
        opserr << "        with any elements in BasicBuilder\n";
        return TCL_ERROR;
      }
    }

    // TODO: remove this
    int G3_setDomain(G3_Runtime*, Domain*);
    G3_setDomain(rt, theNewDomain);
    // create the model builder
    theNewBuilder = new BasicModelBuilder(*theNewDomain, interp, ndm, ndf);
    G3_setModelBuilder(rt, theNewBuilder);

    const char* analysis_option;
    if (!(analysis_option = Tcl_GetVar(interp,"opensees::pragma::analysis",TCL_GLOBAL_ONLY)) ||
         (strcmp(analysis_option, "off") != 0)) {
      G3_AddTclAnalysisAPI(interp, *theNewBuilder);
    }
  }
  else {
    opserr << OpenSees::PromptValueError 
           << "unknown model builder type '" << argv[1] << "' not supported"
           << OpenSees::SignalMessageEnd;
    return TCL_ERROR;
  }

  return TCL_OK;
}

int
TclCommand_wipeModel(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char *argv[])
{
  Tcl_Eval(interp, "_clearAnalysis");

  BasicModelBuilder *builder = static_cast<BasicModelBuilder*>(clientData);

  if (theDatabase != nullptr)
    delete theDatabase;

  if (builder != nullptr) {
    Domain* theDomain = builder->getDomain();
    theDomain->clearAll();
    ops_TheActiveDomain = nullptr;
    delete theDomain;
    delete builder;
    builtModel = false;
  }
  Tcl_CreateCommand(interp, "model", &TclCommand_specifyModel, nullptr, nullptr);
  Tcl_CreateCommand(interp, "wipe",  &TclCommand_wipeModel,    nullptr, nullptr);

  ops_Dt = 0.0;

#ifdef _PARALLEL_PROCESSING
  OPS_PARTITIONED = false;
#endif

  theDatabase = nullptr;

  // the domain deletes the record objects,
  // just have to delete the private array
  return TCL_OK;
}

// command invoked to build the model, i.e. to invoke buildFE_Model()
// on the ModelBuilder
int
buildModel(ClientData clientData, Tcl_Interp *interp, Tcl_Size argc, TCL_Char *argv[])
{
  G3_Runtime *rt = G3_getRuntime(interp);
  BasicModelBuilder* builder = (BasicModelBuilder*)G3_getModelBuilder(rt);

  // TODO: Remove `builtModel` var.
  // to build the model make sure the ModelBuilder has been constructed
  // and that the model has not already been constructed
  if (builder != 0 && builtModel == false) {
    builtModel = true;
    return builder->buildFE_Model();

  } else if (builder != nullptr && builtModel == true) {
    opserr << OpenSees::PromptValueError << "Model has already been built - not built again \n";
    return TCL_ERROR;

  } else {
    opserr << OpenSees::PromptValueError << "No ModelBuilder type has been specified \n";
    return TCL_ERROR;
  }
}

