//===----------------------------------------------------------------------===//
//
//        OpenSees - Open System for Earthquake Engineering Simulation
//
//===----------------------------------------------------------------------===//
//
// Description: This file provides progress bar functionality to the
// interpreter.
//
// Author: Claudio Perez
//
#include <string>
#include <string.h>
#include <tcl.h>
#include <ProgressBar.hpp>


ProgressBar* progress_bar_ptr;

int
TclObjCommand_progress(ClientData clientData, Tcl_Interp *interp, int argc, Tcl_Obj* const*objv)
{
  if (strcmp(Tcl_GetString(objv[1]), "update") == 0) {
    if (clientData == nullptr || *(ProgressBar**)clientData == nullptr) {
      return TCL_ERROR;
    }

    std::string message = "";
    if (argc > 2)
      message = Tcl_GetString(objv[2]);
      
    (*(ProgressBar**)clientData)->update(message);
    return TCL_OK;

  } else if (strcmp(Tcl_GetString(objv[1]), "create") == 0) {
    int steps = 100;

    if (argc > 2 && Tcl_GetIntFromObj(interp, objv[2], &steps) == TCL_ERROR) {
      // Failed to read number of steps
    }

    if (*(ProgressBar**)clientData != nullptr) {
      delete *(ProgressBar**)clientData;
      *(ProgressBar**)clientData = nullptr; 
    }

    ProgressBar *bar = new ProgressBar(steps);
    bar->set_todo_char(L" "); //  ━━━━━━━━━━╸━━━━━━━━━━━━━━━╸━━━━━━━━
    bar->set_done_char(L"█");
    bar->set_opening_char(L"|");
    bar->set_closing_char(L"|");

    *(ProgressBar**)clientData = bar;

    return TCL_OK;
  }
  return TCL_ERROR;
}

