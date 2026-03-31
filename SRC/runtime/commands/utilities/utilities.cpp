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
#include "utilities/linalg.h"

int init_g3_tcl_utils(Tcl_Interp* interp)
{ 
  const int libsize = sizeof(linalg)/sizeof(char*);
  for (int i=0; i < libsize; i++)
    Tcl_Eval(interp, linalg[i]);
  return 0;
}
