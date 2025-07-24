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
#include "Versor.h"
#include <type_traits>

namespace OpenSees {

static_assert(std::is_trivially_copyable<Versor>::value, "Versor must be trivially copyable");

inline Versor
Versor::from_matrix(const Matrix3D &R)
{
  //===--------------------------------------------------------------------===//
  // Form a normalized quaternion (Versor) from a proper orthogonal matrix
  // using Spurrier's algorithm
  //===--------------------------------------------------------------------===//
  Versor q;

  // Trace of the rotation R
  const double trR = R(0,0) + R(1,1) + R(2,2);

  // a = max([trR R(0,0) R(1,1) R(2,2)]);
  double a = trR;
  for (int i = 0; i < 3; i++)
    if (R(i,i) > a)
      a = R(i,i);

  if (a == trR) {
    q.scalar = std::sqrt(1.0 + a)*0.5;

    for (int i = 0; i < 3; i++) {
      int j = (i+1)%3;
      int k = (i+2)%3;
      q.vector[i] = (R(k,j) - R(j,k))/(4.0*q.scalar);
    }
  }
  else {
    for (int i = 0; i < 3; i++)
      if (a == R(i,i)) {
        int j = (i+1)%3;
        int k = (i+2)%3;

        q.vector[i] = std::sqrt(std::max(a*0.5 + (1.0 - trR)/4.0, 0.0));
        q.scalar    = (R(k,j) - R(j,k))/(4.0*q.vector[i]);
        q.vector[j] = (R(j,i) + R(i,j))/(4.0*q.vector[i]);
        q.vector[k] = (R(k,i) + R(i,k))/(4.0*q.vector[i]);
      }
  }
  return q;
}

template<typename Vec3T>
inline Versor
Versor::from_vector(const Vec3T  &theta)
{
  double angle2 = 0.0;
  for (int i=0; i<3; i++)
    angle2 += theta[i]*theta[i];

  double angle = std::sqrt(angle2);

  Versor q;
  double sc, cs;
  if (angle2 < 1e-12) {
    sc = 0.5 - angle2 / 48.0; // + angle2*angle2 / 3840.0 - angle2*angle2*angle2 / 362880.0;
    cs = 1.0 - angle2 /  8.0; // + angle2*angle2 / 384.0 - angle2*angle2*angle2 / 40320.0;
  }
  else {
    sc = std::sin(angle*0.5) / angle;
    cs = std::cos(angle*0.5);
  }

  for (int i = 0; i < 3; i++)
    q.vector[i] = theta[i] * sc;

  q.scalar = cs;
  return q;
}
} // namespace OpenSees


inline OpenSees::Versor
operator*(const OpenSees::Versor &qa, const OpenSees::Versor &qb)
{
  OpenSees::Versor q12;
  q12.scalar = qa.scalar * qb.scalar - qa.vector.dot(qb.vector);
  q12.vector = (qb.vector * qa.scalar) + (qa.vector * qb.scalar) + qa.vector.cross(qb.vector);
  return q12;
}