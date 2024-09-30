//===----------------------------------------------------------------------===//
//
//        OpenSees - Open System for Earthquake Engineering Simulation
//
//===----------------------------------------------------------------------===//
//
//
#include <tcl.h>
#include <string.h>
#include <OPS_Globals.h>

int
TclObjCommand_pragma([[maybe_unused]] ClientData clientData, 
                     Tcl_Interp *interp, int objc, Tcl_Obj *const objv[])
{
  if (objc == 1)
    return TCL_OK;

  int argi = 1;
  const char* pragma = Tcl_GetString(objv[argi++]);

  if (strcmp(pragma, "analysis") == 0) {
    if (argi < objc && strcmp(Tcl_GetString(objv[argi]), "off") == 0) {
      Tcl_Eval(interp,
        "proc loadConst {args} {}\n"
        "proc wipeAnalysis	{args} {}\n"
        "proc constraints {args} {}\n"
        "proc numberer {args} {}\n"
        "proc system {args} {}\n"
        "proc test {args} {}\n"
        "proc algorithm {args} {}\n"
        "proc integrator {args} {}\n"
        "proc analysis {args} {}\n"
        "proc analyze {args} {return 0}\n"
        "proc eigen {args} {return 1}\n"
        "namespace eval opensees::pragma {set analysis off}\n"
      );
      return TCL_OK;
    } 
  }
  else if (strcmp(pragma, "openseespy") == 0) {
    if (argi < objc && strcmp(Tcl_GetString(objv[argi]), "off") == 0) {
      Tcl_Eval(interp, "namespace eval opensees::pragma {set openseespy 0}");

    } else if (argi < objc && strcmp(Tcl_GetString(objv[argi]), "on") == 0) {
      Tcl_Eval(interp, "namespace eval opensees::pragma {set openseespy 1}");

    } else {
      Tcl_Eval(interp, "namespace eval opensees::pragma {set openseespy 1}");
    }
  }
  return TCL_OK;
}

