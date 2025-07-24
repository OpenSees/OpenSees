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
//
// Description: Versor implements a unit-normalized quaternion.
//
#pragma once
#include <cmath>
#include <Vector3D.h>
#include <Matrix3D.h>


namespace OpenSees {

struct Versor {
  Vector3D vector;
  double   scalar;

  template<typename Vec3T>
  static inline Versor from_vector(const Vec3T  &);

  static inline Versor from_matrix(const Matrix3D  &);

  template <typename Vec3T>
  inline Vec3T 
  rotate(const Vec3T& u) const {
    return u + 2.0 * vector.cross( scalar*u + vector.cross(u) );
  }

  inline void 
  normalize() {
    static constexpr double eps = 1e-12;

    double n2 = scalar*scalar + vector.dot(vector);

    // bad numbers; reset to identity
    if (!std::isfinite(n2) || n2 < eps) {
      scalar    = 1.0;
      vector[0] = vector[1] = vector[2] = 0.0;
      return;
    }

    // already close enough to unit
    constexpr double tol = 1e-6; 
    if (std::abs(n2 - 1.0) < tol)
      return;

    // normalize
    double inv_n = 1.0 / std::sqrt(n2);
    scalar      *= inv_n;
    vector[0]   *= inv_n;
    vector[1]   *= inv_n;
    vector[2]   *= inv_n;
  }

  inline constexpr Versor 
  conjugate() const noexcept {
    return Versor{{-vector[0], -vector[1], -vector[2]}, scalar};
  }

  inline Versor 
  conj_mult(const Versor& other) const noexcept
  {
    // Equivalent to R_IJ = R_I.T @ R_J, 
    // i.e.  q_IJ = conj(q_I)*q_J.
    Versor qij;
    qij.scalar = scalar * other.scalar + vector.dot(other.vector);
    qij.vector = (other.vector * scalar) - (vector * other.scalar) - vector.cross(other.vector);
    return qij;
  }

  inline Versor 
  mult_conj(const Versor& other) const noexcept
  {
    // Equivalent to R_IJ = R_I @ R_J^T, 
    // i.e.  q_IJ = q_I * conj(q_J).

    Versor out;
    out.scalar = scalar * other.scalar + vector.dot(other.vector);
    out.vector =      (vector    * other.scalar)
                    - (other.vector * scalar)
                    - vector.cross(other.vector);
    return out;
  }


  // Unary plus
  Versor operator+() const {
    return *this;
  }

  // Unary minus.
  Versor operator-() const {
    return {{-vector[0], -vector[1], -vector[2]}, -scalar};
  }
};

} //  namespace OpenSees

#include "Versor.tpp"

