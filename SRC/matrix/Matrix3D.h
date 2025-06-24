//===----------------------------------------------------------------------===//
//
//        OpenSees - Open System for Earthquake Engineering Simulation
//
//===----------------------------------------------------------------------===//
//
// Claudio Perez
//
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
