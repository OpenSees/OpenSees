//===----------------------------------------------------------------------===//
//
//                                   xara
//                              https://xara.so
//
//----------------------------------------------------------------------------//
//
// Please cite the following resource in any derivative works:
//
// [1] Perez, C.M., and Filippou F.C.. (2024) 
//     "On Nonlinear Geometric Transformations of Finite Elements" 
//     Int. J. Numer. Meth. Engrg.; https://doi.org/10.1002/nme.7506
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
// [3] Nour-Omid, B. and Rankin, C.C. (1991) "Finite rotation analysis and 
//     consistent linearization using projectors", 
//     Computer Methods in Applied Mechanics and Engineering, 93(3), pp. 353–384. 
//     Available at: https://doi.org/10.1016/0045-7825(91)90248-5.
//
#pragma once
#include <array>
#include <Vector3D.h>
#include <Matrix3D.h>
#include <MatrixND.h>
#include "EuclidIsometry.h"

class Node;
#define TRIAD C2

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

  virtual 
  MatrixND<12,12>
  getRotationJacobian(const VectorND<12>&pwx) final 
  {
    MatrixND<3,12> NWL{};
    double Ln = this->getLength();

    constexpr static Vector3D e1 {1,0,0};
    constexpr static Matrix3D ex = Hat(e1);

    for (int i=0; i<nn; i++)
      NWL.assemble(Hat(&pwx[i*6]), 0, i*6,  -1.0);

    MatrixND<12,3> Gamma{};
    Gamma.template insert<0,0>(ex,  1.0);
    Gamma(3,0) = -1.0;
    Gamma.template insert<6,0>(ex, -1.0);

    MatrixND<12,3> Psi{};
    Psi.template insert<3,0>(Eye3, 1.0);
    Psi.template insert<6,0>(ex,    Ln);
    Psi.template insert<9,0>(Eye3, 1.0);

    Matrix3D B = Gamma^Psi;
    B.invert();
    return Gamma*B.transpose()*NWL;
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

    double Ln = this->getLength();

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
};
} // namespace OpenSees