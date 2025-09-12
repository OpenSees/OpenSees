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
#ifndef G3PARSE_H
#define G3PARSE_H
#ifndef G3_RUNTIME_H
#include <runtime/runtime/G3_Runtime.h>
#define G3_RUNTIME_H
#endif
#include <tcl.h>

#define G3_Char TCL_Char

enum SuccessFlag {
  G3_OK    = TCL_OK, 
  G3_ERROR = TCL_ERROR
};

typedef enum SuccessFlag SuccessFlag;

#endif // G3PARSE_H
