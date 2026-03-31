//===----------------------------------------------------------------------===//
//
//                                   xara
//                              https://xara.so
//
//===----------------------------------------------------------------------===//
//
// Copyright (c) 2025, OpenSees/Xara Developers
// All rights reserved.  No warranty, explicit or implicit, is provided.
//
// This source code is licensed under the BSD 2-Clause License.
// See LICENSE file or https://opensource.org/licenses/BSD-2-Clause
//
//===----------------------------------------------------------------------===//
//
#include <tcl.h>
#include <string.h>
#include <Parsing.h>

class Domain;

int
TclObjCommand_pragma([[maybe_unused]] ClientData clientData, 
                     Tcl_Interp *interp, Tcl_Size objc, Tcl_Obj *const objv[])
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
        "proc sensitivityAlgorithm {args} {}\n"
        "proc integrator {args} {}\n"
        "proc analysis {args} {}\n"
        "proc analyze {args} {return 0}\n"
        "proc eigen {args} {return 1}\n"
        "namespace eval opensees::pragma {set analysis off}\n"
      );
      return TCL_OK;
    } else if (argi < objc && strcmp(Tcl_GetString(objv[argi]), "on") == 0) {
      Tcl_Eval(interp,
        "proc loadConst {args} {return 0}\n"
        "proc wipeAnalysis	{args} {return 0}\n"
        "proc constraints {args} {return 0}\n"
        "proc numberer {args} {return 0}\n"
        "proc system {args} {return 0}\n"
        "proc test {args} {return 0}\n"
        "proc algorithm {args} {return 0}\n"
        "proc sensitivityAlgorithm {args} {return 0}\n"
        "proc integrator {args} {return 0}\n"
        "proc analysis {args} {return 0}\n"
        "proc analyze {args} {}\n"
        "proc eigen {args} {}\n"
        "namespace eval opensees::pragma {set analysis on}\n"
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

