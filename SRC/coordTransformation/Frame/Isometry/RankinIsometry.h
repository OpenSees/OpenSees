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

namespace OpenSees {

template <int nn>
class RankinIsometry : public AlignedIsometry<nn>
{
public:
  RankinIsometry(const Vector3D& vecxz)
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
      q = RI*E2;
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
    }
    return R;
  }


  MatrixND<12,12>
  getRotationJacobian(const VectorND<12>&pwx) final 
  {
    MatrixND<3,12> NWL{};
    const double Ln = this->getLength();

    constexpr static Matrix3D ex = Hat(Vector3D {1,0,0});

    for (int i=0; i<nn; i++)
      NWL.assemble(Hat(&pwx[i*6]), 0, i*6,  -1.0);

#if __cplusplus >= 202000L
    static constinit MatrixND<12,3> Gamma = MakeGamma();
    static constinit MatrixND<12,3> Psi0  = MakePsi();
    MatrixND<12,3> Psi = Psi0;
    Psi.template insert<6,0>(ex,  -Ln);
#else
    MatrixND<12,3> Gamma{};
    Gamma.template insert<0,0>(ex,  1.0);
    Gamma(3,0) = -1.0;
    Gamma.template insert<6,0>(ex, -1.0);

    MatrixND<12,3> Psi{};
    Psi.template insert<3,0>(Eye3, 1.0);
    Psi.template insert<6,0>(ex,   -Ln);
    Psi.template insert<9,0>(Eye3, 1.0);
#endif
    Matrix3D B = Gamma^Psi;
    Matrix3D A;
    B.invert(A);
    return Gamma*A.transpose()*NWL;
  }

  MatrixND<3,6>
  getBasisVariation(int ie, int node)
  {
    MatrixND<3,6> dei{};
    if (ie == 1) {
      Matrix3D A{};
      A(1,1) = A(2,2) = 1.0/this->getLength();
      if (node == 0) {
        dei.template insert<0,0>(A, -1.0);
      }
      else if (node == nn-1) {
        dei.template insert<0,0>(A,  1.0);
      }
    }

    else if (ie == 3) {
      Matrix3D A{};
      Vector3D q = this->getRotation()^this->q;
      Vector3D v = q.cross(Vector3D{1,0,0});

      A(0,0) = A(1,1) = 1.0/v.norm();
      Matrix3D Q = Hat(q);
      dei = A*(
        Hat(Vector3D{1,0,0})*Q
       +Q*this->getBasisVariation(1, node)
      );
    }
    return dei;
  }

  MatrixND<3,6> 
  getRotationGradient(int node) final {
    MatrixND<3,6> Gb{};

    constexpr Vector3D axis{1, 0, 0};
    constexpr Matrix3D ix = Hat(axis);

    const double Ln = this->getLength();

    if (node == 0) {
      Gb.template insert<0,0>( ix, -1.0/Ln);
      Gb(0,2) =  n/Ln;
      Gb(0,3) =   1.0;
      Gb(0,4) =    -n;
    }
    else if (node == nn-1) {
      Gb.template insert<0,0>( ix,  1.0/Ln);
      Gb(0,2) = -n/Ln;
      Gb(0,3) =  0.0;
    }
    return Gb;
  }

private:
  Vector3D q;
  double n   = 0;


#if __cplusplus >= 202000L
  static inline consteval MatrixND<12,3> 
  MakePsi()
  {
    MatrixND<12,3> Psi{};
    Psi.template insert<3,0>(Eye3, 1.0);
    Psi.template insert<9,0>(Eye3, 1.0);
    return Psi;
  }

  static inline consteval MatrixND<12,3> 
  MakeGamma()
  {
    MatrixND<12,3> Gamma{};
    constexpr Matrix3D ex = Hat(Vector3D {1,0,0});
    Gamma.template insert<0,0>(ex,  1.0);
    Gamma(3,0) = -1.0;
    Gamma.template insert<6,0>(ex, -1.0);
    return Gamma;
  }
#endif

};
} // namespace OpenSees