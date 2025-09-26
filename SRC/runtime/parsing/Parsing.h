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
#pragma once
#include <string.h>
#ifndef TCL_Char
typedef const char TCL_Char;
#endif
#ifndef G3_Char
typedef const char G3_Char;
#endif

#include "InputAPI.h"

// As suggested by https://core.tcl-lang.org/tcl/wiki?name=Migrating+C+extensions+to+Tcl+9
#ifndef TCL_SIZE_MAX
typedef int Tcl_Size;
#define TCL_SIZE_MAX INT_MAX
#endif

class Parameter;
namespace OpenSees {
namespace Parsing {
int
GetDoubleParam(Tcl_Interp *, Domain& , const char* arg, double* value, Parameter*&);
}
}