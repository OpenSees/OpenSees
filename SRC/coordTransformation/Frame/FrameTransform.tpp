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
// Please cite the following resource in any derivative works:
//
// Perez, C.M., and Filippou F.C.. "On Nonlinear Geometric Transformations
//     of Finite Elements" Int. J. Numer. Meth. Engrg. 2024; 
//     https://doi.org/10.1002/nme.7506
//
//===----------------------------------------------------------------------===//
//
#pragma once
#include "FrameTransform.h"
#include <MatrixND.h>
#include <Vector3D.h>

namespace OpenSees {

template <int nn, int ndf>
VectorND<nn*ndf> 
FrameTransform<nn,ndf>::pushConstant(const VectorND<nn*ndf>& pl) 
{
  Matrix3D R;
  Vector3D x, y, z;
  getLocalAxes(x, y, z);
  for (int i=0; i<3; i++) {
    R(i,0) = x[i];
    R(i,1) = y[i];
    R(i,2) = z[i];
  }

  constexpr int N = nn * ndf;

  //
  // Initialize
  //
  VectorND<N> pg = pl;

  // (A) First pass: just do the direct transformations
  for (int i=0; i<nn; i++) {
    const int base = i * ndf;
    pg.insert(base,   R*Vector3D{pg[base  ], pg[base+1], pg[base+2]}, 1.0);
    pg.insert(base+3, R*Vector3D{pg[base+3], pg[base+4], pg[base+5]}, 1.0);
  }

  // (B) Second pass: add offset cross product ( r x F ) into the
  //     rotational DOFs [3..5]
  const std::array<Vector3D,nn> *offset = this->getRigidOffsets();
  if (offset) {
    const std::array<Vector3D, nn>& offsets = *offset;
    for (int i=0; i<nn; i++) {

      const int base = i * ndf;
      const Vector3D ni {
          pg[base+0],
          pg[base+1],
          pg[base+2]
      };

      // Add M = r x F
      pg.assemble(base+3, offsets[i].cross(ni), 1.0);
    }
  }

  return pg;
}


template<int nn, int ndf>
MatrixND<nn*ndf,nn*ndf>
FrameTransform<nn,ndf>::pushConstant(const MatrixND<nn*ndf,nn*ndf>& kl)
{
  //
  // Do diag(R)*M*diag(R)'
  //
  static MatrixND<nn*ndf,nn*ndf> Kg;
  Kg = kl;

  Matrix3D R;
  Vector3D x, y, z;
  getLocalAxes(x, y, z);
  for (int i=0; i<3; i++) {
    R(i,0) = x[i];
    R(i,1) = y[i];
    R(i,2) = z[i];
  }
  
  this->pushRotation(Kg, R);

  const std::array<Vector3D,nn> *offset = this->getRigidOffsets();
  if (offset) [[unlikely]] {
    this->pushOffsets(Kg, *offset);
  }
  return Kg;
}

template<int nn, int ndf>
constexpr void
FrameTransform<nn,ndf>::pushRotation(MatrixND<nn*ndf,nn*ndf>& Kg, const Matrix3D& R)
{
  //
  // Do diag(R)*M*diag(R)'
  //
  const Matrix3D RT = R.transpose();
  for (int i=0; i<nn; i++) {
    for (int j=0; j<nn; j++) {
      // Loop over 3x3 blocks for n and m
      for (int k=0; k<2; k++) {
        for (int l=0; l<2; l++) {
          Matrix3D Kab {{
            Kg(i*ndf+3*k+0, j*ndf+3*l  ),
            Kg(i*ndf+3*k+1, j*ndf+3*l  ),
            Kg(i*ndf+3*k+2, j*ndf+3*l  ),

            Kg(i*ndf+3*k+0, j*ndf+3*l+1),
            Kg(i*ndf+3*k+1, j*ndf+3*l+1),
            Kg(i*ndf+3*k+2, j*ndf+3*l+1),

            Kg(i*ndf+3*k+0, j*ndf+3*l+2),
            Kg(i*ndf+3*k+1, j*ndf+3*l+2),
            Kg(i*ndf+3*k+2, j*ndf+3*l+2)
          }};
          Kab = Kab*RT;
          Kab = R*Kab;
          Kg.insert(Kab, i*ndf+3*k, j*ndf+3*l, 1.0);
        }
      }
    }
  }
}

template<int nn, int ndf>
constexpr void
FrameTransform<nn,ndf>::pushOffsets(MatrixND<nn*ndf,nn*ndf>& Kg, 
                                     const std::array<Vector3D,nn>& offsets)
{
  for (int i=0; i<nn; i++) {
    for (int j=0; j<nn; j++) {
      {
        Matrix3D KmnW{};
        KmnW.addMatrixSpinProduct(Kg.template extract<3,3>(i*ndf+3, j*ndf), offsets[j], -1.0);
        Kg.assemble(KmnW, i*ndf+3, j*ndf+3, 1.0);
      }
      {
        Matrix3D WKnm{};
        WKnm.addSpinMatrixProduct(offsets[i], Kg.template extract<3,3>(i*ndf, j*ndf+3), 1.0);
        Kg.assemble(WKnm, i*ndf+3, j*ndf+3, 1.0);
      }
      {
        Matrix3D WKnn{};
        {
          Matrix3D Knn = Kg.template extract<3,3>(i*ndf, j*ndf);
          {
            Matrix3D KnnW{};
            KnnW.addMatrixSpinProduct(Knn, offsets[j], -1.0);
            Kg.assemble(KnnW, i*ndf, j*ndf+3, 1.0);
          }
          WKnn.addSpinMatrixProduct(offsets[i], Knn, 1.0);
          Kg.assemble(WKnn, i*ndf+3, j*ndf, 1.0);
        }
        Matrix3D WKmmW{};
        WKmmW.addMatrixSpinProduct(WKnn, offsets[j], -1.0);
        Kg.assemble(WKmmW, i*ndf+3, j*ndf+3, 1.0);
      }
    }
  }
}

} // namespace OpenSees