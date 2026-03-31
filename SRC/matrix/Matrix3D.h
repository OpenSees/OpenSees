//===----------------------------------------------------------------------===//
//
//                                   xara
//                              https://xara.so
//
//===----------------------------------------------------------------------===//
//
// Copyright (c) 2025, Claudio M. Perez
// All rights reserved.  No warranty, explicit or implicit, is provided.
//
// This source code is licensed under the BSD 2-Clause License.
// See LICENSE file or https://opensource.org/licenses/BSD-2-Clause
//
//===----------------------------------------------------------------------===//

#ifndef Matrix3D_H
#define Matrix3D_H

#include "MatrixND.h"
#include <type_traits>

namespace OpenSees {

  using Matrix3D = MatrixND<3,3,double>;

  static_assert(std::is_trivially_copyable<Matrix3D>::value, "Matrix3D is not trivially copyable.");
  static_assert(std::is_trivial<Matrix3D>::value, "Matrix3D is not trivial.");
  static_assert(std::is_standard_layout<Matrix3D>::value, "Matrix3D is not standard layout.");
  static_assert(std::is_aggregate<Matrix3D>::value, "Matrix3D is not an aggregate type.");
}

#endif // Matrix3D_H
