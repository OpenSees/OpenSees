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
// [2] Filippou, F.C. (1998)
//     "FEDEASLab: Finite Elements for Design Evaluation and Analysis of Structures"
//
// [3] Battini, J.-M. and Pacoste, C. (2002) 
//     "Co-rotational beam elements with warping effects in instability problems", 
//     Computer Methods in Applied Mechanics and Engineering, 191(17–18), 
//     pp. 1755–1789. Available at: https://doi.org/10.1016/S0045-7825(01)00352-8.
//
#pragma once
#include <array>
#include <Vector3D.h>
#include <Matrix3D.h>
#include <MatrixND.h>
#include "EuclidIsometry.h"

class Node;

namespace OpenSees {

template <int nn>
class BattiniIsometry : public AlignedIsometry<nn>
{
public:
  BattiniIsometry(const Vector3D& vecxz)
  : AlignedIsometry<nn>(vecxz), 
    n(0)
  {
  }

  using AlignedIsometry<nn>::init;
  using AlignedIsometry<nn>::pres;

  Matrix3D
  update_basis(const Matrix3D& RI, const Matrix3D& RJ, const Vector3D& dx) final
  {
    Matrix3D R;
    {
      Vector3D e1 = dx;
      e1 /= e1.norm();

      constexpr static Vector3D D2 {0,1,0};
      const Vector3D E2 = this->AlignedIsometry<nn>::R[init]*D2;
      q = RI*E2; //*R[init];
      q.addVector(0.5, RJ*E2, 0.5);


      Vector3D e3 = e1.cross(q);
      e3 /= e3.norm();

      Vector3D e2 = e3.cross(e1);

      for (int i = 0; i < 3; i++) {
        R(i,0) = e1[i];
        R(i,1) = e2[i];
        R(i,2) = e3[i];
      }
    
      Vector3D Q = R^q;
      n = Q[0]/Q[1];

      Vector3D QI = R^(RI*E2);
      Vector3D QJ = R^(RJ*E2);
      n11 = QI[0]/Q[1];
      n12 = QI[1]/Q[1];
      n21 = QJ[0]/Q[1];
      n22 = QJ[1]/Q[1];
    }
    return R;
  }


  MatrixND<3,6> 
  getRotationGradient(int node) final {
    MatrixND<3,6> Gb{};

    constexpr Vector3D axis{1, 0, 0};
    constexpr Matrix3D ix = Hat(axis);

    double Ln = this->getLength();

    if (node == 0) {
      Gb.template insert<0,0>( ix, -1.0/Ln);
      Gb(0,2) =    n/Ln;
      Gb(0,3) =  n12/2.0; // - n;
      Gb(0,4) = -n11/2.0;
    }
    else if (node == nn-1) {
      Gb.template insert<0,0>( ix,   1.0/Ln);
      Gb(0,2) = -  n/Ln;
      Gb(0,3) =  n22/2.0;
      Gb(0,4) = -n12/2.0;
    }
    
    return Gb;
  }

private:
  Vector3D q {0,1,0};
  double n   = 0,
         n11 = 0,
         n12 = 1,
         n21 = 0,
         n22 = 1;
};
} // namespace OpenSees