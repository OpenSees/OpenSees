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
// Written: fmk
//
#ifndef OPS_PACKAGES_H
#define OPS_PACKAGES_H
int getLibraryFunction(const char *libName, const char *functName,
                       void **libHandle, void **funcHandle);
#endif // OPS_PACKAGES_H
