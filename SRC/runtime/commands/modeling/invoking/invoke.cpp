//===----------------------------------------------------------------------===//
//
//        OpenSees - Open System for Earthquake Engineering Simulation
//
//===----------------------------------------------------------------------===//
//
// Written: cmp
// Created: Spring 2023
//
#include <tcl.h>
#include <string>
#include <unordered_map>

#include <G3_Logging.h>

Tcl_CmdProc TclCommand_useUniaxialMaterial;
Tcl_CmdProc TclCommand_useCrossSection;
Tcl_CmdProc TclCommand_usePlaneStress;

const std::unordered_map<std::string, Tcl_CmdProc*> invoke_commands 
{
  {"UniaxialMaterial",  &TclCommand_useUniaxialMaterial       },

  {"FrameSection",      &TclCommand_useCrossSection           },
  {"section",           &TclCommand_useCrossSection           },

  {"PlaneStress",       &TclCommand_usePlaneStress            }
};

int
TclCommand_invoke(ClientData clientData, Tcl_Interp* interp, int argc, char const** const argv)
{
  // check number of arguments in command line
  if (argc < 4) {
    opserr << G3_ERROR_PROMPT << "bad arguments - want: using <obj-type> <obj-tag> {<operations>...}";
    return TCL_ERROR;
  }

  auto tcl_cmd = invoke_commands.find(std::string(argv[1]));

  if (tcl_cmd != invoke_commands.end()) {

    return (*tcl_cmd->second)(clientData, interp, argc, &argv[0]);

  } else {
    return TCL_ERROR;
  }

  return TCL_OK;

}

