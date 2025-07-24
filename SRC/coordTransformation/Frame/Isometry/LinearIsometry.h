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
// Written: Claudio M. Perez, 
//          Filip C. Filippou
//          University of California, Berkeley
//
// Developed with FEDEASLab [2].
//
// References:
//
// [1] Perez, C.M., and Filippou F.C. (2024)
//     "On Nonlinear Geometric Transformations of Finite Elements" 
//     Int. J. Numer. Meth. Engrg.; https://doi.org/10.1002/nme.7506
//
// [2] Filippou, F.C. (1998)
//     "FEDEASLab: Finite Elements for Design Evaluation and Analysis of Structures"
//
// [3] Nour-Omid, B. and Rankin, C.C. (1991) "Finite rotation analysis and 
//     consistent linearization using projectors", 
//     Computer Methods in Applied Mechanics and Engineering, 93(3), pp. 353–384. 
//     Available at: https://doi.org/10.1016/0045-7825(91)90248-5.
//
#pragma once
#include <Vector3D.h>
#include <Matrix3D.h>
#include <MatrixND.h>
#include "EuclidIsometry.h"

class Node;
#define TRIAD C2

namespace OpenSees {

template <int nn>
class LinearIsometry : public AlignedIsometry<nn>
{
public:
  LinearIsometry(const Vector3D& vecxz)
  : AlignedIsometry<nn>(vecxz)
  {
  }

  Matrix3D
  update_basis(const Matrix3D&, const Matrix3D&, const Vector3D&) final 
  {
    return this->AlignedIsometry<nn>::R[init];
  }

  MatrixND<3,6> 
  getRotationGradient(int node) final {
    MatrixND<3,6> Gb{};

    constexpr Matrix3D ix = Hat(Vector3D {1, 0, 0});

    const double Ln = this->getLength();

    if (node == 0) {
      Gb.template insert<0,0>( ix, -1.0/Ln);
      Gb(0,3) =   0.5;
    }
    else if (node == nn-1) {
      Gb.template insert<0,0>( ix,  1.0/Ln);
      Gb(0,3) =  0.5;
    }
    return Gb;
  }
private:
  using AlignedIsometry<nn>::init;
  double L;
};
} // namespace OpenSees