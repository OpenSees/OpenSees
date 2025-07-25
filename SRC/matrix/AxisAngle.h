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
#pragma once
#include "Vector3D.h"
#include "VectorND.h"
#include "MatrixND.h"
#include "Matrix3D.h"
#include "GroupSO3.h"

namespace {

class AxisAngle {
public:
  Vector3D vector;
  double   angle, eta, mu;

  AxisAngle() 
  : vector{},
    angle(0.0),
    eta(0.0),
    mu(0.0)
  {
  }

  explicit
  AxisAngle(const Matrix3D& R)
  : AxisAngle(LogSO3(R))
  {
  }

  explicit
  AxisAngle(const Vector3D& v)
   : angle(v.norm()),
     vector(v)
  {
    Utility::dLogConst(angle, eta, mu);
  }

  inline Matrix3D
  dLog() const
  {
    //
    // d_R LogSO3(v) = Eye3 - 0.5*Sv + eta*Sv*Sv;
    //
    Matrix3D dH = Eye3;
    dH.addSpin(vector, -0.5);
    dH.addSpinSquare(vector, eta);
    return dH;
  }

  inline Vector3D
  dLog(const Vector3D&p) const
  {
    //
    // d_R LogSO3(v) = Eye3 - 0.5*Sv + eta*Sv*Sv;
    //
    Vector3D u = p;
    u -= vector.cross(p)*0.5;
    u += eta*vector.cross(vector.cross(p));
    return u;
  }

  inline Matrix3D
  ddLog(const Vector3D& p) const
  {
    //
    // -0.5*Hat(v) + eta*(Eye3*th.dot(v) + th.bun(v) - 2.*v.bun(th)) + mu*St2*v.bun(th);
    //
    const Vector3D& v = this->vector;

    Matrix3D St2 = Hat(vector);
    St2 = St2*St2;
    Matrix3D dH{};
    dH.addSpin(p, -0.5);  
    dH.addDiagonal( eta*vector.dot(p));
    dH.addTensorProduct(v, p, eta);
    dH.addTensorProduct(p, v, -2.0*eta);
    dH.addMatrixProduct(St2, p.bun(v), mu);

    return dH*this->dLog();
  }
};

} // namespace