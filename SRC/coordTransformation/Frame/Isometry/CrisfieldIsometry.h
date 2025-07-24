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
// [3] Crisfield, M.A. (1990) "A consistent co-rotational formulation for
//     non-linear, three-dimensional, beam-elements", Computer Methods in Applied
//     Mechanics and Engineering, 81(2), pp. 131–150. Available at:
//     https://doi.org/10.1016/0045-7825(90)90106-V.
//
#pragma once
#include <limits>
#include <cmath>
#include <MatrixND.h>
#include <Matrix3D.h>
#include <GroupSO3.h>
#include "EuclidIsometry.h"

class Node;

namespace OpenSees {

namespace {
struct Triad {

  Triad(const OpenSees::Matrix3D &E)
    : e{{E(0,0),E(1,0),E(2,0)}, // e1
        {E(0,1),E(1,1),E(2,1)}, // e2
        {E(0,2),E(1,2),E(2,2)}} // e3
  {  
  }

  constexpr inline
  const Vector3D& operator[](int i) const {
    return e[i-1];
  }

  const Vector3D e[3];
};
}

template <int nn, bool orthogonal>
class CrisfieldIsometry : public AlignedIsometry<nn> {
public:
  CrisfieldIsometry(const Vector3D& vecxz)
   : AlignedIsometry<nn>{vecxz}, 
  //  Lr2{}, Lr3{}, 
   v{}
  {

  }

  Matrix3D
  update_basis(const Matrix3D& RI, const Matrix3D& RJ, const Vector3D& dx)
  {
    // Ln = dx.norm();
    this->AlignedIsometry<nn>::Ln = dx.norm();
    double Ln = this->getLength();
    {
      const Triad TrI{RI}, TrJ{RJ};
      rI[0] = TrI[1];
      rI[1] = TrI[2];
      rI[2] = TrI[3];
      rJ[0] = TrJ[1];
      rJ[1] = TrJ[2];
      rJ[2] = TrJ[3];
    }

    {
      // const Versor qI = Versor::from_matrix(RI);
      // const Versor qJ = Versor::from_matrix(RJ);
      // Versor qij = qJ.mult_conj(qI);
      // qij.normalize();
      Versor qij = Versor::from_matrix(RJ*RI.transpose());
      Vector3D gw = CayleyFromVersor(qij);

      gw *= 0.5;

      Rbar = CaySO3(gw)*RI;

      const Triad r{Rbar};
      r1 = r[1];
      r2 = r[2];
      r3 = r[3];
    }

    e[0]  = dx;
    e[0] /= Ln;

    //
    // Compute the base vectors e2, e3
    //
    // 'rotate' the mean rotation matrix Rbar on to e1 to
    // obtain e2 and e3 

    Matrix3D E;
    if constexpr (orthogonal)
    {
      constexpr double ktol = 1.0*std::numeric_limits<double>::epsilon();
      Vector3D r1 { Rbar(0,0), Rbar(1,0), Rbar(2,0) };

      // Clamp to avoid NaNs from acos
      double dot = std::max(-1.0, std::min(1.0, r1.dot(e[0])));

      if (std::fabs(dot - 1.0) < ktol) {                         // Rbar already aligned
        v.zero();
        for (int i=0; i<3; i++) {
          E(i,0) = e[0][i];
          E(i,1) = Rbar(i,1);
          E(i,2) = Rbar(i,2);
        }
      }
      else if (std::fabs(dot + 1.0) < ktol) {                    // opposite direction
        static constexpr double pi = 3.14159265358979323846;
        // choose any axis with numerical separation from r1
        v = r1.cross(Vector3D{1.0, 0.0, 0.0});
        if (v.dot(v) < ktol)
          v = r1.cross(Vector3D{0.0, 1.0, 0.0});
        v *= pi / v.norm();
        E = ExpSO3(v)*Rbar;
      }
      else {                                                     // general case
        v     = r1.cross(e[0]);
        double angle = std::atan2(v.norm(), dot);
        v    *= angle / v.norm();
        E = ExpSO3(v)*Rbar;
      }

      for (int i = 0; i < 3; i++)
        for (int j = 0; j < 3; j++)
          e[j][i] = E(i,j);
    }
    else
    {
      // Not working with frame-1007?
      //
      // Use the 'mid-point' procedure

      // e2 = r2 - (e1 + r1)*((r2^e1)*0.5);
  
      Vector3D tmp;
      tmp  = e[0];
      tmp += r1;

      e[1] = tmp;
      {
        e[1] *= 0.5*r2.dot(e[0]);
        e[1].addVector(-1.0,  r2, 1.0);
      }
  
      // e3 = r3 - (e1 + r1)*((r3^e1)*0.5);
      e[2] = tmp;
      {
        e[2] *= r3.dot(e[0])*0.5;
        e[2].addVector(-1.0,  r3, 1.0);
      }
      for (int i = 0; i < 3; i++)
        for (int j = 0; j < 3; j++)
          E(i,j) = e[j][i];
    }


    //
    //
    //
    for (int i = 0; i < 3; i++)
      for (int j = 0; j < 3; j++)
        A(i,j) = (double(i==j) - e[0][i]*e[0][j])/Ln;

    return E;
  }

  MatrixND<12,12>
  getRotationJacobian(const VectorND<12>&pl) final {
    MatrixND<12,12> dG{};

    return dG;
  }

  // MatrixND<3,6> getRotationGradient(int node) final;

  MatrixND<3,6> 
  getRotationGradient(int node) final {
    constexpr static Vector3D e1{1, 0, 0};

    double Ln = this->getLength();

    static constexpr
    MatrixND<1,3> E3 {{0.0, 0.0, 1.0}},
                  E2 {{0.0, 1.0, 0.0}};

    const Matrix3D& Tr = this->getRotation();
    Triad r{Tr^Rbar};
    Vector3D r1 = r[1], r3 = r[3];


    Matrix3D A{};
    for (int i = 0; i < 3; i++)
      for (int j = 0; j < 3; j++)
        A(i,j) = (double(i==j) - e1[i]*e1[j])/Ln;

    MatrixND<3,12> de3;
    if constexpr (!orthogonal)
      de3 = this->getLMatrix(r3, r1, e1, A).transpose();
    else
      de3 = this->getBasisVariation(r3, r1, e1, v, A);

    MatrixND<3,12> de1{};
    de1.template insert<0,0>(A, -1.0);
    de1.template insert<0,6>(A,  1.0);

    MatrixND<3,12> G{};
    G.template insert<0,0, 1,12>(E2*de3, -1.0);
    G.template insert<1,0, 1,12>(E3*de1, -1.0);
    G.template insert<2,0, 1,12>(E2*de1,  1.0);

    if (node == 0)
      return G.template extract<0,3,  0, 6>();
    else if (node == nn-1)
      return G.template extract<0,3,  6,12>();
    else
      return MatrixND<3,6>{};

  }

  //
  //
  //
  inline MatrixND<12,12>
  compute_tangent(const VectorND<12> &ul)
  {
    // extract columns of rotation matrices
    const Vector3D 
      &e1  =  e[0],
      &e2  =  e[1],
      &e3  =  e[2],

      &rI1 = rI[0],
      &rI2 = rI[1],
      &rI3 = rI[2],
      &rJ1 = rJ[0],
      &rJ2 = rJ[1],
      &rJ3 = rJ[2];
    MatrixND<12,3> Lr2,Lr3;
    Lr2 = this->getLMatrix(r2, r1, e[0], A);
    Lr3 = this->getLMatrix(r3, r1, e[0], A);
    //
    //
    //
    MatrixND<12,12> T{};

    //   T1 = [      O', (-S(rI3)*e2 + S(rI2)*e3)',        O', O']';

    // (-S(rI3)*e2 + S(rI2)*e3)
    Vector3D rxe  = rI2.cross(e3);
    rxe -= rI3.cross(e2);
    for (int i = 0; i < 3; i++)
      T(imx,i+3) =  rxe[i];

    //   T2 = [(A*rI2)', (-S(rI2)*e1 + S(rI1)*e2)', -(A*rI2)', O']';

    Vector3D At = A*rI2;

    // (-S(rI2)*e1 + S(rI1)*e2)'
    rxe  = rI1.cross(e2);
    rxe -= rI2.cross(e1);
    for (int i = 0; i < 3; i++) {
      T(imz,i  ) =  At[i];
      T(imz,i+3) =  rxe[i];
      T(imz,i+6) = -At[i];
    }

    //   T3 = [(A*rI3)', (-S(rI3)*e1 + S(rI1)*e3)', -(A*rI3)', O']';

    At = A*rI3;
    
    // -S(rI3)*e1 + S(rI1)*e3
    rxe  = rI1.cross(e3);
    rxe -= rI3.cross(e1);
    for (int i = 0; i < 3; i++) {
      T(imy,i  ) = - At[i];
      T(imy,i+3) = -rxe[i];
      T(imy,i+6) =   At[i];
    }

    //   T4 = [      O', O',        O', (-S(rJ3)*e2 + S(rJ2)*e3)']';
    rxe  = rJ2.cross(e3);
    rxe -= rJ3.cross(e2);
    for (int i = 0; i < 3; i++)
      T(jmx, i+9) =  rxe[i];   // S(rJ2)*e3 - S(rJ3)*e2

    // T5 = [(A*rJ2)', O', -(A*rJ2)', (-S(rJ2)*e1 + S(rJ1)*e2)']';
    At = A*rJ2;
    rxe  = rJ1.cross(e2); 
    rxe -= rJ2.cross(e1);
    for (int i = 0; i < 3; i++) {
        T(jmz, i  ) =  At[i];
        T(jmz, i+6) = -At[i];
        T(jmz, i+9) = rxe[i]; // (-S(rJ2)*e1 + S(rJ1)*e2)
    }

    // T6 = [(A*rJ3)', O', -(A*rJ3)', (-S(rJ3)*e1 + S(rJ1)*e3)']'
    At   = A*rJ3;
    rxe  = rJ1.cross(e3);  // (-S(rJ3)*e1 + S(rJ1)*e3)
    rxe -= rJ3.cross(e1);
    for (int i = 0; i < 3; i++) {
        T(jmy,i  ) =  -At[i];
        T(jmy,i+6) =   At[i];
        T(jmy,i+9) = -rxe[i];
    }

    //
    // Second part
    //

    // T(:,1) += Lr3*rI2 - Lr2*rI3;
    // T(:,2) +=           Lr2*rI1; z
    // T(:,3) += Lr3*rI1          ; y

    // T(:,4) += Lr3*rJ2 - Lr2*rJ3;
    // T(:,5) += Lr2*rJ1          ; z    // TODO ?????? check sign
    // T(:,6) += Lr3*rJ1          ; y    // TODO ?????? check sign

    // Bending Z
    for (int i = 0; i < 12; i++) {
      double T1i = 0;
      for (int k=0; k<3; k++)
        T1i += Lr2(i,k)*rI1[k];
      T(imz,i) += T1i;
    }

    for (int i = 0; i < 12; i++) {
      double T4i = 0;
      for (int k=0; k<3; k++)
        T4i += Lr2(i,k)*rJ1[k]; // Lr[i];
      T(jmz,i) += T4i;
    }

    // Torsion
    for (int i = 0; i < 12; i++) {
      double T0i = 0;
      for (int k=0; k<3; k++)
        T0i += Lr3(i,k)*rI2[k] - Lr2(i,k)*rI3[k];
      // T(jmx,i) += -T0i;
      T(imx,i) += T0i;
    }
    for (int i = 0; i < 12; i++) {
      double T3i = 0;
      for (int k=0; k<3; k++)
        T3i += Lr3(i,k)*rJ2[k] - Lr2(i,k)*rJ3[k];
      T(jmx,i) += T3i;
    }
    // Bending Y
    for (int i = 0; i < 12; i++) {
      double T2i = 0;
      for (int k=0; k<3; k++)
        T2i += Lr3(i,k)*rI1[k]; // Lr[i];
      T(imy,i) += T2i*-1;
    }
    for (int i = 0; i < 12; i++) {
      double T5i = 0;
      for (int k=0; k<3; k++)
        T5i += Lr3(i,k)*rJ1[k]; // Lr[i];
      T(jmy,i) += T5i*-1;
    }

    //
    //
    //
    for (int node=0; node < 2; node++)
      for (int j = 0; j < 3; j++) {
        const double c = 0.5 / std::cos(ul[(node? jmx : imx) + j]);
        for (int i = 0; i < 12; i++)
          T((node? jmx : imx) + j, i) *= c;
      }

    // Axial
    // T(:,7) = [-e1' O' e1' O']';
    for (int i = 0; i < 3; i++) {
        T(jnx,i  ) = -e1[i];
        T(jnx,i+6) =  e1[i];
    }

    return T;
  }


  //
  // Add geometric part of the transformation tangent
  //
  //  kg += T'*kl*T + ks1 + T * diag(m.*tan(thetal))*T'
  //      +  m(4)*(ks2r2t3_u3 + ks2r3u2_t2)
  //      +  m(2)*ks2r2t1 + m(3)*ks2r3t1
  //      +  m(5)*ks2r2u1 + m(6)*ks2r3u1
  //      +  ks3 + ks3' + ks4 + ks5;
  //
  int addTangent(MatrixND<12,12>& kg,
                                      const VectorND<12>& pl,
                                      const VectorND<12>& ul)
  {
    const Vector3D 
      &e1  =  e[0],
      &e2  =  e[1],
      &e3  =  e[2],

      &rI1 = rI[0],
      &rI2 = rI[1],
      &rI3 = rI[2],
      &rJ1 = rJ[0],
      &rJ2 = rJ[1],
      &rJ3 = rJ[2];

    double Ln = this->getLength();

    //
    // Ksigma1
    //
    {
      const double N = -pl[0]; // Axial force
      this->getKs1Matrix(kg, N);
    }

    //
    // Ksigma3
    //
    VectorND<6> m;
    m[0] =  0.5*pl[imx]/std::cos(ul(imx));
    m[2] = -0.5*pl[imy]/std::cos(ul(imy));
    m[1] =  0.5*pl[imz]/std::cos(ul(imz));

    m[3] =  0.5*pl[jmx]/std::cos(ul(jmx));
    m[5] = -0.5*pl[jmy]/std::cos(ul(jmy));
    m[4] =  0.5*pl[jmz]/std::cos(ul(jmz));

    this->getKs3Matrix(kg, m);

    //
    // Ksigma4
    //
    {
      Matrix3D ks33;
      ks33.zero();
      ks33.addSpinProduct(e2, rI3,  m[3]);
      ks33.addSpinProduct(e3, rI2, -m[3]);
      ks33.addSpinProduct(e2, rI1,  m[1]);
      ks33.addSpinProduct(e1, rI2, -m[1]);
      ks33.addSpinProduct(e3, rI1,  m[2]);
      ks33.addSpinProduct(e1, rI3, -m[2]);
      kg.assemble(ks33, 3, 3, 1.0);
    }

    {
      Matrix3D ks99;
      ks99.zero();
      ks99.addSpinProduct(e2, rJ3, -m[3]);
      ks99.addSpinProduct(e3, rJ2,  m[3]);
      ks99.addSpinProduct(e2, rJ1,  m[4]);
      ks99.addSpinProduct(e1, rJ2, -m[4]);
      ks99.addSpinProduct(e3, rJ1,  m[5]);
      ks99.addSpinProduct(e1, rJ3, -m[5]);

      kg.assemble(ks99, 9, 9, 1.0);
    }

    //
    // Ksigma5
    //
    //  Ks5 = [ Ks5_11   Ks5_12 | -Ks5_11   Ks5_14;
    //          Ks5_12'    O    | -Ks5_12'   O;
    //         -Ks5_11  -Ks5_12 |  Ks5_11  -Ks5_14;
    //          Ks5_14t     O   | -Ks5_14'   O];
    //
    //
    // v = (1/Ln)*(m(2)*rI2 + m(3)*rI3 + m(5)*rJ2 + m(6)*rJ3);
    //   = 1/Ln * (m[1]*rI2 + m[2]*rI3)
    //   + 1/Ln * (m[4]*rJ2 + m[5]*rJ3);
    //   = vi + vj
    //
    {
      Vector3D v;
      v.addVector(0.0, rI2, m[1]);
      v.addVector(1.0, rI3, m[2]);
      v.addVector(1.0, rJ2, m[4]);
      v.addVector(1.0, rJ3, m[5]);
      v /= Ln;

      // Ks5_11 = A*v*e1' + e1*v'*A + (e1'*v)*A;
      //        = A*vi*e1' + e1*vi'*A + (e1'*vi)*A
      //        + A*vj*e1' + e1*vj'*A + (e1'*vj)*A;
      //
      Matrix3D ks33{};
      ks33.addMatrix(A, e1.dot(v));

      Matrix3D m33{};
      m33.addTensorProduct(v, e1, 1.0);

      ks33.addMatrixProduct(A, m33, 1.0);

      for (int i = 0; i < 3; i++)
        for (int j = 0; j < 3; j++)
          m33(i,j) = e1[i]*v[j];

      ks33.addMatrixProduct(m33, A, 1.0);

      kg.assemble(ks33, 0, 0,  1.0);
      kg.assemble(ks33, 0, 6, -1.0);
      kg.assemble(ks33, 6, 0, -1.0);
      kg.assemble(ks33, 6, 6,  1.0);
    }

    // Ks5_12 = -(m(2)*A*S(rI2) + m(3)*A*S(rI3));
    {
      Matrix3D ks33;
      ks33.zero();
      ks33.addMatrixSpinProduct(A, rI2, -m[1]);
      ks33.addMatrixSpinProduct(A, rI3, -m[2]);

      kg.assemble(ks33, 0, 3,  1.0);
      kg.assemble(ks33, 6, 3, -1.0);
      kg.assembleTranspose(ks33, 3, 0,  1.0);
      kg.assembleTranspose(ks33, 3, 6, -1.0);

      //  Ks5_14 = -(m(5)*A*S(rJ2) + m(6)*A*S(rJ3));

      ks33.zero();
      ks33.addMatrixSpinProduct(A, rJ2, -m[4]);
      ks33.addMatrixSpinProduct(A, rJ3, -m[5]);

      kg.assemble(ks33, 0, 9,  1.0);
      kg.assemble(ks33, 6, 9, -1.0);

      kg.assembleTranspose(ks33, 9, 0,  1.0);
      kg.assembleTranspose(ks33, 9, 6, -1.0);
    }

    // Ksigma

    this->getKs2Matrix(kg, r2, rI3-rJ3, m[3]);
    this->getKs2Matrix(kg, r3, rJ2-rI2, m[3]);
    this->getKs2Matrix(kg, r2,     rI1, m[1]);
    this->getKs2Matrix(kg, r3,     rI1, m[2]);
    this->getKs2Matrix(kg, r2,     rJ1, m[4]);
    this->getKs2Matrix(kg, r3,     rJ1, m[5]);

    return 0;
  }


private:

  inline MatrixND<3,12>
  getBasisVariation(const Vector3D &ri, 
                    const Vector3D& r1, 
                    const Vector3D& e1, 
                    const Vector3D& v,
                    const Matrix3D& A) const noexcept
  {
    double nu = v.norm();
    Matrix3D L1 = ExpSO3(v)^dExpSO3(v, ri)*Hat(r1)*A; // nu
    Matrix3D L2 = ExpSO3(v)^dExpSO3(v, ri)*Hat(e1)*Hat(r1); // nu/2.0
    Matrix3D L3 = Hat(ri);

    MatrixND<3,12> L{};
    // i
    L.assemble(L1, 0, 0,  -nu);
    L.assemble(L2, 0, 3,   nu/2.0);
    L.assemble(L3, 0, 3,     -0.5);
    // j
    L.assemble(L1, 0, 6,   nu);
    L.assemble(L2, 0, 9,   nu/2.0);
    L.assemble(L3, 0, 9,     -0.5);

    return L;
  }

  [[nodiscard]] inline MatrixND<12,3>
  getLMatrix(const Vector3D &ri, 
             const Vector3D& r1,
             const Vector3D& e1,
             const Matrix3D& A) const noexcept
  {

    static Matrix3D rie1r1;
    static Matrix3D e1e1r1;

    const double rie1 = ri.dot(e1);

    for (int k = 0; k < 3; k++) {
      const double e1r1k = (e1[k] + r1[k]);
      for (int j = 0; j < 3; j++) {
        rie1r1(j,k) = ri[j]*e1r1k;
        e1e1r1(j,k) = e1[j]*e1r1k;
      }
    }

    static Matrix3D L1, L2;
    // L1  = ri'*e1 * A/2 + A*ri*(e1 + r1)'/2;
    L1.zero();
    L1.addMatrix(A, rie1*0.5);
    L1.addMatrixProduct(A, rie1r1, 0.5);

    // L2  = Sri/2 - ri'*e1*S(r1)/4 - Sri*e1*(e1 + r1)'/4;
    L2.zero();
    L2.addSpin(ri, 0.5);
    L2.addSpin(r1, -rie1/4.0);
    L2.addSpinMatrixProduct(ri, e1e1r1, -0.25);

    // L = [L1
    //      L2
    //     -L1
    //      L2];

    MatrixND<12,3> L{};
    L.assemble(L1, 0, 0,  1.0);
    L.assemble(L2, 3, 0,  1.0);
    L.assemble(L1, 6, 0, -1.0);
    L.assemble(L2, 9, 0,  1.0);

    return L;
  }

  inline void 
  getKs1Matrix(MatrixND<12,12> &Kg, double N) const noexcept
  {
    // a=0
    Kg.assemble(A, 0, 0,  N);
    Kg.assemble(A, 0, 6, -N);
    // a=1
    Kg.assemble(A, 6, 0, -N);
    Kg.assemble(A, 6, 6,  N);
  }

  inline void
  getKs3Matrix(MatrixND<12,12> &Kg, const VectorND<6>& m) const noexcept
  {
    //
    // Ksigma3
    //
    //  ks3 = [o kbar2  |  o kbar4];
    //
    //  where
    //
    //    kbar2 = -Lr2*(m(3)*S(rI3) + m(1)*S(rI1)) + Lr3*(m(3)*S(rI2) - m(2)*S(rI1)) ;
    //
    //    kbar4 =  Lr2*(m(3)*S(rJ3) - m(4)*S(rJ1)) - Lr3*(m(3)*S(rJ2) + m(5)*S(rJ1));
    //
    // or
    //
    //  ks3 = [o ka+kb  |  o kc+kd];
    //      = [o ka     |  o kc] + [o kb  |  o kd];
    //
    //  where
    //
    //    ka = -Lr2*S(rI3)*m(3)  
    //         +Lr2*S(rI1)*m(1);
    //    kb =  Lr3*S(rI2)*m(3)  
    //         -Lr3*S(rI1)*m(2);
    //
    //    kc =  Lr2*S(rJ3)*m(3)
    //         -Lr2*S(rJ1)*m(4);
    //    kd = -Lr3*S(rJ2)*m(3)  
    //         +Lr3*S(rJ1)*m(5);

    const Vector3D &rI1 = rI[0],
                   &rI2 = rI[1],
                   &rI3 = rI[2],
                   &rJ1 = rJ[0],
                   &rJ2 = rJ[1],
                   &rJ3 = rJ[2];

    
    MatrixND<12,3> Lr2,Lr3;
    Lr2 = this->getLMatrix(r2, r1, e[0], A);
    Lr3 = this->getLMatrix(r3, r1, e[0], A);

    Matrix3D Sm{};
    Sm.addSpin(rI3,  m[3]);
    Sm.addSpin(rI1,  m[1]);
    static MatrixND<12,3> kbar;
    kbar.zero();
    kbar.addMatrixProduct(Lr2, Sm, -1.0);

    Sm.zero();
    Sm.addSpin(rI2,  m[3]);
    Sm.addSpin(rI1, -m[2]);
    kbar.addMatrixProduct(Lr3, Sm,  1.0);

    Kg.assemble(kbar, 0, 3, 1.0);
    Kg.assembleTranspose(kbar, 3, 0, 1.0);

    Sm.zero();
    Sm.addSpin(rJ3,  m[3]);
    Sm.addSpin(rJ1, -m[4]);
    kbar.zero();
    kbar.addMatrixProduct(Lr2, Sm, 1.0);

    Sm.zero();
    Sm.addSpin(rJ2, m[3]);
    Sm.addSpin(rJ1, m[5]);
    kbar.addMatrixProduct(Lr3, Sm,  -1.0);

    Kg.assemble(kbar, 0, 9, 1.0);
    Kg.assembleTranspose(kbar, 9, 0, 1.0);
  }

  inline void
  getKs2Matrix(MatrixND<12,12> &Kg, 
               const Vector3D &ri, 
               const Vector3D &z,
               double scale) const noexcept
  {
    const Vector3D &e1 = e[0];

    const double Ln = this->getLength();

    //  Ksigma2 = [ K11   K12 -K11   K12
    //              K12'  K22 -K12'  K22
    //             -K11  -K12  K11  -K12
    //              K12'  K22 -K12'  K22];

    // U = (-1/2)*A*z*ri'*A + ri'*e1*A*z*e1'/(2*Ln)+...
    //      z'*(e1+r1)*A*ri*e1'/(2*Ln);

    const double rite1 = ri.dot(e1);
    const double zte1  =  z.dot(e1);
    const double ztr1  =  z.dot(r1);

    static Matrix3D zrit, ze1t;
    static Matrix rizt(3,3), rie1t(3,3);
    static Matrix3D e1zt;

    //  const Matrix3D e1zt = e1.bun(z);

    // Chrystal's looping order
    for (int j = 0; j < 3; j++) {
      for (int i = 0; i < 3; i++) {
        zrit(i,j)  = z[i]*ri[j];
        rizt(i,j)  = ri[i]*z[j];
        ze1t(i,j)  = z[i]*e1[j];
        e1zt(i,j)  = e1[i]*z[j];
        rie1t(i,j) = ri[i]*e1[j];
      }
    }

    Matrix3D U;
    U.addMatrixTripleProduct(0.0, A, zrit, -0.5);
    U.addMatrixProduct(A, ze1t,   rite1/(2*Ln));
    U.addMatrixProduct(A, rie1t, (zte1 + ztr1)/(2*Ln));


    // K11 = U + U' + ri'*e1*(2*(e1'*z)+z'*r1)*A/(2*Ln);

    Matrix3D ks{};
    ks.addMatrix(U, 1.0);

    // Add matrix U transpose
    for (int i = 0; i < 3; i++)
      for (int j = 0; j < 3; j++)
        ks(i,j) += U(j,i);

    ks.addMatrix(A, rite1*(2*zte1 + ztr1)/(2*Ln));

    Kg.assemble(ks, 0, 0,  scale);
    Kg.assemble(ks, 0, 6, -scale);
    Kg.assemble(ks, 6, 0, -scale);
    Kg.assemble(ks, 6, 6,  scale);

    // K12 = (1/4)*(-A*z*e1'*Sri - A*ri*z'*Sr1 - z'*(e1+r1)*A*Sri);
    Matrix3D m1{};
    m1.addMatrixProduct(A, ze1t, -1.0);
    ks.zero();
    ks.addMatrixSpinProduct(m1, ri, 0.25);

    m1.zero();
    m1.addMatrixProduct(A, rizt, -1.0);
    ks.addMatrixSpinProduct(m1, r1, 0.25);
    ks.addMatrixSpinProduct(A, ri, -0.25*(zte1+ztr1));

    Kg.assemble(ks, 0, 3,  scale);
    Kg.assemble(ks, 0, 9,  scale);
    Kg.assemble(ks, 6, 3, -scale);
    Kg.assemble(ks, 6, 9, -scale);

    Kg.assembleTranspose(ks, 3, 0,  scale);
    Kg.assembleTranspose(ks, 3, 6, -scale);
    Kg.assembleTranspose(ks, 9, 0,  scale);
    Kg.assembleTranspose(ks, 9, 6, -scale);

    // K22 = (1/8)*((-ri'*e1)*Sz*Sr1 + Sr1*z*e1'*Sri + ...
    //       Sri*e1*z'*Sr1 - (e1+r1)'*z*S(e1)*Sri + 2*Sz*Sri);

    ks.zero();
    ks.addSpinProduct(z, r1, -0.125*(rite1));

    m1.zero();
    m1.addSpinMatrixProduct( r1, ze1t, 1.0);
    ks.addMatrixSpinProduct( m1, ri, 0.125);

    m1.zero();
    m1.addSpinMatrixProduct(ri, e1zt, 1.0);
    ks.addMatrixSpinProduct(m1, r1, 0.125);

    ks.addSpinProduct(e1, ri, -0.125*(zte1 + ztr1));
    ks.addSpinProduct( z, ri, 0.25);

    // Ksigma2 = [ K11   K12 -K11   K12;
    //             K12t  K22 -K12t  K22;
    //            -K11  -K12  K11  -K12;
    //             K12t  K22 -K12t  K22];

    Kg.assemble(ks, 3, 3, scale);
    Kg.assemble(ks, 3, 9, scale);
    Kg.assemble(ks, 9, 3, scale);
    Kg.assemble(ks, 9, 9, scale);
  }


private:

  enum {
    inx= 0, // axial
    iny= 1, // Vy
    inz= 2, // Vz
    imx= 3, // torsion
    imy= 4, // rot y I
    imz= 5, // rot z I

    jnx= 6, // axial
    jny= 7,
    jnz= 8,
    jmx= 9, // torsion
    jmy=10, // rot y J
    jmz=11, // rot z J
  };

  Matrix3D A;
  Matrix3D Rbar;
  Vector3D v;
  Vector3D r1, r2, r3;
  Vector3D e[3], rI[3], rJ[3];
  // MatrixND<12,3> Lr2, Lr3;
  // double Ln;
};
}
