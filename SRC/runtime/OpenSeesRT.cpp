//===----------------------------------------------------------------------===//
//
//        OpenSees - Open System for Earthquake Engineering Simulation
//
//===----------------------------------------------------------------------===//
//
// This file contains functions that are required by Tcl to load the
// OpenSeesRT library.
//
#ifndef OPENSEESRT_VERSION
#  define OPENSEESRT_VERSION "0.0.0"
#endif
//
#include <runtimeAPI.h>
#include "G3_Runtime.h"
#include <logging/G3_Logging.h>
#include <handler/OPS_Stream.h>
#include <StandardStream.h>
#include "commands/strings.cpp"
#include <stdio.h>
#include <stdlib.h>

// Determine when stdout is a TTY
#ifdef _WIN32
#  include <io.h>
#  define isatty _isatty
#  define STDERR_FILENO _fileno(stderr)
#else
#  include <unistd.h>               
#endif
//
extern int OpenSeesAppInit(Tcl_Interp *interp);
extern void G3_InitTclSequentialAPI(Tcl_Interp* interp);
extern int init_g3_tcl_utils(Tcl_Interp*);

//
// Tcl Command that returns the current OpenSees version
//
static int
version(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char **argv)
{
  char buffer[20];

  sprintf(buffer, "%s", OPENSEESRT_VERSION);
  Tcl_SetResult(interp, buffer, TCL_VOLATILE);

  return TCL_OK;
}


//
// Called when the library is loaded as a Tcl extension.
//
extern "C" int 
#ifdef _WIN32
__declspec(dllexport) // DLLEXPORT
#endif
Openseesrt_Init(Tcl_Interp *interp)
{
  if (Tcl_InitStubs(interp, TCL_VERSION, 0) == NULL)
    return TCL_ERROR;

  if (Tcl_PkgProvide(interp, "OpenSeesRT", OPENSEESRT_VERSION) == TCL_ERROR)
    return TCL_ERROR;

  // Create a runtime instance, and store it with the interpreter
  G3_Runtime *rt = new G3_Runtime{interp};
  Tcl_SetAssocData(interp, "G3_Runtime", NULL, (ClientData)rt);

  // Initialize OpenSees
  OpenSeesAppInit(interp);
  G3_InitTclSequentialAPI(interp); // Add sequential API
  init_g3_tcl_utils(interp);       // Add utility commands (linspace, range, etc.)

  char* verbosity = getenv("OPENSEESRT_VERBOSITY");
  if (verbosity != nullptr) {
    if (strcmp(verbosity, "DEBUG") == 0) {
      G3_SetStreamLevel(G3_LevelDebug, true);
    }
  }

  // Prevent coloring output when stderr is not a TTY
  if (isatty(STDERR_FILENO))
    G3_SetStreamColor(nullptr, G3_LevelWarn, 1);


  // Set some variables with package information
  Tcl_SetVar(interp, "opensees::copyright", copyright,      TCL_LEAVE_ERR_MSG);
  Tcl_SetVar(interp, "opensees::license",   license,        TCL_LEAVE_ERR_MSG);
  Tcl_SetVar(interp, "opensees::banner",    unicode_banner, TCL_LEAVE_ERR_MSG);
  Tcl_CreateCommand(interp, "version",      version,      nullptr, nullptr);
  return TCL_OK;
}

