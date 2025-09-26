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
#include <tcl.h>

Tcl_CmdProc TclCommand_wipeModel;
Tcl_CmdProc TclCommand_clearAnalysis;
Tcl_CmdProc TclCommand_specifyModel;


// formats.cpp
Tcl_CmdProc convertBinaryToText;
Tcl_CmdProc convertTextToBinary;
Tcl_CmdProc stripOpenSeesXML;

// domain/peri/commands.cpp
Tcl_CmdProc Tcl_Peri;

struct char_cmd {
  const char* name; Tcl_CmdProc*  func;
}
const InterpreterCommands[] =  {

  {"peri",                 Tcl_Peri},

  {"stripXML",             stripOpenSeesXML    },
  {"convertBinaryToText",  convertBinaryToText },
  {"convertTextToBinary",  convertTextToBinary },
};
