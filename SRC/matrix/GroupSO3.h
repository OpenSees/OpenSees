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
// References:
//
// [1] Perez, C.M., and Filippou F.C.. "On Nonlinear Geometric Transformations
//     of Finite Elements" Int. J. Numer. Meth. Engrg. 2024; 
//     https://doi.org/10.1002/nme.7506
//
// [2] Ritto-Corrêa, M. and Camotim, D. (2002) "On the differentiation of the
//     Rodrigues formula and its significance for the vector-like parameterization
//     of Reissner-Simo beam theory"
//     Int. J. Numer. Meth. Engrg., 55(9), pp.
//     1005–1032. Available at: https://doi.org/10.1002/nme.532.
//
// [3] Ibrahimbegović, A. and Mikdad, M.A. (1998) ‘Finite rotations in dynamics of
//     beams and implicit time‐stepping’, 41, pp. 781–814.
//
// [4] Pfister, F. (1998) ‘Bernoulli Numbers and Rotational Kinematics’,
//     Journal of Applied Mechanics, 65(3), pp. 758–763. 
//     Available at: https://doi.org/10.1115/1.2789120.
//
// [5] Nurlanov Z (2021) Exploring SO(3) logarithmic map: degeneracies and 
//     derivatives.
//
//===----------------------------------------------------------------------===//

//
// Functions related to exponential coordinates on SO(3). The implementation
// follows value semantics with the expectation that RVO will make this
// performant.
//
// Written: Claudio M. Perez
// 2023
//
#pragma once

#include <cmath>
#include <stdlib.h>
#include "Versor.h"
#include "Matrix3D.h"
#include "Vector3D.h"
using OpenSees::Matrix3D;
using OpenSees::Versor;

static constexpr double SmallAngle = 1e-12; // threshold for small angles
static constexpr double SmallAngle2 = SmallAngle * SmallAngle; // threshold for small angles squared

static constexpr Matrix3D Eye3 {{
  1, 0, 0,
  0, 1, 0,
  0, 0, 1
}};



static inline constexpr Vector3D
Vee(const Matrix3D &X) noexcept
{
//
// Return the axial vector x of the given skew-symmetric 3x3 matrix X.
//
  return {X(2,1), X(0,2), X(1,0)};
}

template <typename Vec3Type>
static inline constexpr Matrix3D
Hat(const Vec3Type &u) noexcept
{
  return Matrix3D {{  0  ,  u[2], -u[1],
                    -u[2],   0  ,  u[0],
                     u[1], -u[0],   0  }};
}


static inline double
GibSO3(const Vector3D &vec, double *a, double *b=nullptr, double *c=nullptr) noexcept
{
//
// Compute coefficients of the Rodrigues formula and their differentials.
//
//                       [1,2]     | [3]
//                       ----------+-------
//                       a1        |
//                       a2        | c4
//                       a3        | c5
//
//                       b1 = -c0  | c1
//                       b2        | c2 
//                       b3        | c3
//
//
// [2] Ritto-Corrêa, M. and Camotim, D. (2002) "On the differentiation of the
//     Rodrigues formula and its significance for the vector-like parameterization
//     of Reissner-Simo beam theory"
//     Int. J. Numer. Meth. Engrg., 55(9), pp.
//     1005–1032. Available at: https://doi.org/10.1002/nme.532.
//
// [3] Ibrahimbegović, A. and Mikdad, M.A. (1998) ‘Finite rotations in dynamics of
//     beams and implicit time‐stepping’, 41, pp. 781–814.
//
// [4] Pfister, F. (1998) ‘Bernoulli Numbers and Rotational Kinematics’,
//     Journal of Applied Mechanics, 65(3), pp. 758–763. 
//     Available at: https://doi.org/10.1115/1.2789120.
//
// ---------------------------------------------------------------------------
// Written: cmp                                                           2023
//===----------------------------------------------------------------------===//
//
  double angle2 =  vec.dot(vec);

  if (angle2  <= 1e-12) {
    //
    // Maclaurin series in Horner form
    //
    if (a != nullptr) {
      a[0]   =   1.0        - angle2*(1.0/2.0   - angle2*(1.0/6.0    - angle2/24.0));
      a[1]   =   1.0        - angle2*(1.0/6.0   - angle2*(1.0/120.0  - angle2/5040.0));
      a[2]   =   0.5        - angle2*(1.0/24.0  - angle2*(1.0/720.0  - angle2/40320.0));
      // a[3]   =   1.0/6.0    - angle2/(1.0/120.0 - angle2/(1.0/5040.0 - angle2/362880.0));
      a[3]   =   1.0/6.0    - angle2*(1.0/120.0 - angle2*(1.0/5040.0 - angle2/362880.0));
    }
    if (b != nullptr) {
      b[0]   = - a[1];
      b[1]   = - 1.0/3.0    + angle2*(1.0/30.0   - angle2*(1.0/840.0  - angle2/45360.0));
      b[2]   = - 1.0/12.0   + angle2*(1.0/180.0  - angle2*(1.0/6720.0 - angle2/453600.0));
      b[3]   = - 1.0/60.0   + angle2*(1.0/1260.0 - angle2*(1.0/60480 - angle2/4989600.0));
    }
    if (c != nullptr) {
      c[0]   = - b[1];
      c[1]   = b[3] - b[2];
      c[2]   = 1.0/90.0     - angle2*(1.0/1680.0  - angle2*(1.0/75600.0  - angle2/5987520.0)); 
      c[3]   = 1.0/630.0    - angle2*(1.0/15120.0 - angle2*(1.0/831600.0 - angle2/77837760.0));
    }

  } else {
    //
    // Use formulas of Equation (A8) from [1]
    //
    double angle  = vec.norm();
    double sn     = std::sin(angle);
    double cs     = std::cos(angle);
    double angle3 = angle*angle2;
    double angle4 = angle*angle3;
    double angle5 = angle*angle4;
                                                    
    if (a != nullptr) {
      a[0] = cs;
      a[1] = sn / angle;   
      a[2] = ( 1.0 - a[0] ) / angle2;    
  //  a[3] = ( 1.0 - a[1] ) / angle2;
      a[3] = (angle - sn)/(angle3);
    }

    if (b != nullptr) {
      b[0] = -a[1];
      b[1] = ( angle*cs - sn)/angle3;  
      b[2] = ( angle*sn - 2.0 + 2.0*cs)/angle4;
      b[3] = ( 3.0*sn - 2.0*angle -  angle*cs )/angle5;    
    }

    if (c != nullptr) {
      c[0] = -b[1];
      c[1] = (3.*sn - angle2*sn - 3.*angle*cs)/(angle5);
      c[2] = (8. - 8.*cs - 5.*angle*sn + angle2*cs)/(angle5*angle);
      c[3] = (8.*angle + 7.*angle*cs + angle2*sn - 15.*sn)/(angle5*angle2);
    }
  }

  return std::sqrt(angle2);
}

//
// Rotation Conversions
//

static inline Vector3D
CayleyFromVersor(const Versor &q)
{
  const double q0 = q.scalar;

  Vector3D w;
  for (int i = 0; i < 3; i++)
    w[i] = 2.0 * q.vector[i]/q0;

  return w;
}



// R = (q0^2 - q' * q)*I + 2 * q * q' + 2*q0*S(q);
static inline Matrix3D
MatrixFromVersor(const Versor &q)
{
  Matrix3D R{};

  const double factor = q.scalar*q.scalar - q.vector.dot(q.vector);

  for (int i = 0; i < 3; i++)
    R(i,i) = factor;

  R.addTensorProduct(q.vector, q.vector, 2.0);
  R.addSpin(q.vector, 2.0*q.scalar);

  return R;
}


// Form rotation axis vector from a quaternion
static inline Vector3D
VectorFromVersor(const Versor& q)
{
  // Initialize to zero
  Vector3D theta{};

  double q0 = q.scalar;

  //  Norm of the vector part
  double qn = q.vector.norm();

  const double* const qv = &q.vector[0];

  // Return zero-vector
  if (qn == 0)
     return theta;

  double factor = 1.0;
  if (q0 < 0) {
     q0     = -q0;
     factor = -1.0;
  }

  if (qn < q0)
     factor *= 2.0*std::asin(qn)/qn;
  else
     factor *= 2.0*std::acos(q0)/qn;


  for (int i=0; i<3; i++)
    theta[i] = qv[i]*factor;
  
  return theta;

}


// 
// Exponential Map
//

static inline Matrix3D
ExpSO3(const Vector3D &theta)
{
  // Form the first Gib coefficients
  double a[4];

  Matrix3D R = Eye3;
  if (GibSO3(theta, a) > 1e-16) {
    R.addSpin(theta, a[1]);
    R.addSpinSquare(theta, a[2]);
  }
  return R;
}


static inline Matrix3D
CaySO3(const Vector3D &cayley)
{
  //
  // Cayley map for a rotation matrix given the "tangent-scaled pseudo-vector"
  //
  // R = I + (S + S*S/2)/(1 + w' * w / 4);
  //
  //===--------------------------------------------------------------------===//

  const double c = 1.0/(1.0 + cayley.dot(cayley)/4.0);

  Matrix3D R;
  R.zero();
  R.addDiagonal(1.0);
  R.addSpin(cayley, c);
  R.addSpinSquare(cayley, 0.5*c);

  return R;
}



//
// Exponential Differentials
//

inline Matrix3D
TanSO3(const Vector3D &vec, char repr='L')
{
  //
  //  Compute right or left differential of the exponential.
  //
  //===--------------------------------------------------------------------===//
  double a[4];
  GibSO3(vec, a);

  Matrix3D T;

  // Assemble differential
  switch (repr) {
    case 'R':
      T(0,0)  =         a[1] + a[3]*vec[0]*vec[0];
      T(0,1)  = -vec[2]*a[2] + a[3]*vec[0]*vec[1];
      T(0,2)  =  vec[1]*a[2] + a[3]*vec[0]*vec[2];
      T(1,0)  =  vec[2]*a[2] + a[3]*vec[1]*vec[0];
      T(1,1)  =         a[1] + a[3]*vec[1]*vec[1];
      T(1,2)  = -vec[0]*a[2] + a[3]*vec[1]*vec[2];
      T(2,0)  = -vec[1]*a[2] + a[3]*vec[2]*vec[0];
      T(2,1)  =  vec[0]*a[2] + a[3]*vec[2]*vec[1];
      T(2,2)  =         a[1] + a[3]*vec[2]*vec[2];
      return T;

    case 'L':
    default:
      T(0,0)  =         a[1] + a[3]*vec[0]*vec[0];
      T(0,1)  =  vec[2]*a[2] + a[3]*vec[1]*vec[0];
      T(0,2)  = -vec[1]*a[2] + a[3]*vec[2]*vec[0];
      T(1,0)  = -vec[2]*a[2] + a[3]*vec[0]*vec[1];
      T(1,1)  =         a[1] + a[3]*vec[1]*vec[1];
      T(1,2)  =  vec[0]*a[2] + a[3]*vec[2]*vec[1];
      T(2,0)  =  vec[1]*a[2] + a[3]*vec[0]*vec[2];
      T(2,1)  = -vec[0]*a[2] + a[3]*vec[1]*vec[2];
      T(2,2)  =         a[1] + a[3]*vec[2]*vec[2];
      return T;

  }
}


inline Matrix3D
TExpSO3(const Vector3D &v)
{
// Compute the right differential of the exponential on SO(3) using the formula in
// Equation (A11) in [1]:
//
//        T = a[1]*Eye3 + a[2]*v.hat() + a[3]*v.bun(v);
//
// -----------------------------------------------------------------------------------------

  // Form first Gib coefficients
  double a[4];
  GibSO3(v, a);

  Matrix3D T{};
  T.addDiagonal(a[1])
   .addSpin(v, a[2])
   .addTensorProduct(v, v, a[3]);

  return T;
}


inline Matrix3D
ddTanSO3(const Vector3D &v, const Vector3D &p, const Vector3D &q)
{
  //
  //    return a[3]*psq + b[1]*p.dot(q)*Eye3
  //         + b[2]*(pxq.bun(v) + v.bun(pxq) + vxp.dot(q)*Eye3)
  //         + b[3]*( v.dot(p)*(q.bun(v) + v.bun(q)) 
  //              + v.dot(q)*(p.bun(v) + v.bun(p)) 
  //              + v.dot(p)*v.dot(q)*Eye3) 
  //         + vov*(c[1]*p.dot(q) + c[2]*(vxp.dot(q)) + c[3]*v.dot(p)*v.dot(q));
  //    
  // =========================================================================================
  double a[4], b[4], c[4];
  GibSO3(v, a, b, c);

  const Vector3D pxq = p.cross(q);
  const Vector3D vxp = v.cross(p);

  Matrix3D dT{};

  dT.addTensorProduct(p, q, a[3])
    .addTensorProduct(q, p, a[3])
    .addDiagonal(b[1]*p.dot(q))
    .addTensorProduct(pxq, v, b[2])
    .addTensorProduct(v, pxq, b[2])
    .addDiagonal(vxp.dot(q)*b[2])

    .addTensorProduct(q, v, b[3]*v.dot(p))
    .addTensorProduct(v, q, b[3]*v.dot(p))
    .addTensorProduct(p, v, b[3]*v.dot(q))
    .addTensorProduct(v, p, b[3]*v.dot(q))
    .addDiagonal(v.dot(p)*v.dot(q)*b[3])

    .addTensorProduct(v, v, c[1]*p.dot(q) + c[2]*(vxp.dot(q)) + c[3]*v.dot(p)*v.dot(q));

  return dT;
}



inline Matrix3D
dTanSO3(const Vector3D &v, const Vector3D &p, char repr='L')
{
  // Tangents of the right and left exponential differentials using Equations (A13) and (A14),
  // respectively from [1].
  //
  // repr     'L' or 'R' indicating left or right representation, 
  //          respectively, for the tangent space of SO(3)


  double a[4], b[4];
  GibSO3(v, a, b);

  Matrix3D Xi{};
  Xi.addSpin(p, a[2]*(repr == 'R' ? -1.0 : 1.0))
    .addDiagonal(a[3]*v.dot(p))
    .addTensorProduct(v, p, a[3])
    .addTensorProduct(p, v, b[1])
    .addTensorProduct(v.cross(p), v, b[2]*(repr == 'R' ? 1.0 : -1.0))
    .addTensorProduct(v, v, v.dot(p)*b[3]);
  return Xi;
}


inline Matrix3D
dExpSO3(const Vector3D &v, const Vector3D &p)
{
  double a[4], b[4];
  GibSO3(v, a, b);

  Matrix3D Xi{};
  Xi.addSpin(p, -a[1])
    .addDiagonal(a[2]*v.dot(p))
    .addTensorProduct(v, p, a[2])
    .addTensorProduct(p, v, b[0])
    .addTensorProduct(v.cross(p), v, b[1])
    .addTensorProduct(v, v, v.dot(p)*b[2]);

  return Xi;
}


inline Vector3D
LogSO3(const Matrix3D &R)
{
  //===--------------------------------------------------------------------===//
  //
  // Inverse of the exponential map on SO(3).
  //
  // Returns the axial parameters associated with the rotation `R`. The result
  // should satisfy the following equality for any 3-vector, `v`:
  //
  //       LogSO3(expm(Hat(v))) == v
  //
  // where `expm` is matrix exponential, and `Hat` is a function
  // which produces the skew-symmetric 3x3 matrix associated with vector `v`.
  //
  // Parameters
  //   R       (3x3)   Rotation (proper orthogonal) matrix.
  //
  //
  // References
  // 1. Nurlanov Z (2021) Exploring SO(3) logarithmic map: degeneracies and 
  //    derivatives.
  //
  //===--------------------------------------------------------------------===//

  return VectorFromVersor(Versor::from_matrix(R));
}


static inline Vector3D
LogC90(const Matrix3D &R)
{
  //===--------------------------------------------------------------------===//
  //
  // Crisfield's approximation to the logarithm on SO(3)
  //
  // - Crisfield M (1990) A consistent co-rotational formulation for non-linear,
  //   three-dimensional, beam-elements
  //
  //===--------------------------------------------------------------------===//

  return Vector3D {
    -std::asin(0.5*(R(1,2) - R(2,1))),
     std::asin(0.5*(R(0,2) - R(2,0))),
    -std::asin(0.5*(R(0,1) - R(1,0))),
  };
}

static inline Vector3D
LogC90_Skew(const Matrix3D &R)
{
  //===--------------------------------------------------------------------===//
  //
  //   Crisfield M (1990) A consistent co-rotational formulation for non-linear,
  //   three-dimensional, beam-elements
  //
  //===--------------------------------------------------------------------===//
  // Crisfield's approximation to the logarithm on SO(3)

  return Vector3D {
    std::asin(0.5*(R(1,2) - R(2,1))),
    std::asin(0.5*(R(0,1) - R(1,0))),
    std::asin(0.5*(R(0,2) - R(2,0))),
  };
}

namespace Utility {

template<typename T, std::size_t N>
constexpr T horner(T x2, const std::array<T,N>& c) noexcept
{
  // Evaluate  c0 + c1*x2 + c2*x2² + ...  with Horner's rule (and FMA if available).
  T y = c.back();
  for (std::size_t i = N-1; i-- > 0; )
      y = std::fma(x2, y, c[i]);
  return y;
}

static inline double 
dLogConst(double angle, double& a3, double& b3)
{
  //
  // a3 == eta
  //
  static constexpr double tol = 1.0/20.0;

  double ahalf  = angle/2.0;
  double angle2 = angle*angle;
  double angle3 = angle*angle2;
  double angle4 = angle*angle3;
  double angle5 = angle*angle4;
  double angle6 = angle*angle5;

  double a1;
  if (std::fabs(angle) > tol) [[unlikely]] {
    // a0 = 1.0 + angle2/2.0;
    a1 = ahalf/std::tan(ahalf);

    a3 = (1.0 - a1)/angle2;
    // a3 = (2.0*sn - angle*(1.0+cs))/(2.0*angle2*sn);
    // a3 = (1.0 - 0.5*angle*cot(0.5*angle))/angle2;

    // b3 = (angle*(angle + 2.0*hsn*hcs) - 8.0*hsn*hsn)/(4.0*angle4*hsn*hsn); // wrong?
    // b3 = (angle*(angle + std::sin(angle)) - 8.0*hsn*hsn)/(4.0*angle4*hsn*hsn);
    // b3 = std::fma(angle, angle + sn, 4.0*(cs - 1.0))/(2.0*angle4*(1.0 - cs));
    // b3 = std::fma(angle, angle + 2.0*hsn*hcs, -8.0*hsn*hsn)/(4.0*angle4*hsn*hsn);
    b3 = std::fma(a1, a1 + 1.0, 0.25 * angle2 - 2.0)/angle4;
  }
  else {
    // Maclaurin series
    a1 = 1.0 - angle2/12.0 + angle4/720.0 - angle6/30240.0;
    a3 = 1.0/12. + angle2/720. + angle4/30240. + angle6/1209600.;
    b3 = 1.0/360. + angle2/7560. + angle4/201600. + angle6/5987520.;
  }
  return 0.0;
}

static inline Vector3D
RescaleVector(const Vector3D &v, double &angle)
{
  static constexpr double pi     = 3.14159265358979323846;
  static constexpr double two_pi = 2.0 * pi;
  static constexpr double eps    = 1e-8;     // pole-avoidance band

  Vector3D u = v;

  // NOTE: Currently we are skipping the rescaling as it is not generally
  // needed by the formulations that have been ported to OpenSees from
  // FEDEASLab.
  return u;

  // Fast exit when |angle| already <= pi (no poles inside)
  if (std::fabs(angle) <= pi)
      return u;

  // Wrap once modulo 2pi
  double wrapped = std::remainder(angle, two_pi);   // (-pi, pi]

  // Nudge away from the sole pole at 0
  if (std::fabs(wrapped) < eps)
    wrapped = (wrapped < 0.0 ? -eps : eps);

// #define PRINCIPAL_WRAP
#ifndef PRINCIPAL_WRAP
  // Keep original magnitude unless it was near a pole
  if (std::fabs(angle) > pi && std::fabs(wrapped) >= eps)
    wrapped = angle;
#endif

  // Rescale u
  if (std::fabs(angle) > eps) {  // safe denominator
    double scale = wrapped / angle;
    for (int i = 0; i < 3; ++i)
      u[i] *= scale;
  }

  angle = wrapped;
  return u;
}
}

inline Matrix3D
dLogSO3(const Vector3D &v, double* a=nullptr)
{
  //===--------------------------------------------------------------------===//
  //
  // d_R LogSO3(v) = Eye3 - 0.5*Sv + eta*Sv*Sv;
  //
  //===--------------------------------------------------------------------===//
  //

  double angle = v.norm();
  const Vector3D u = Utility::RescaleVector(v, angle);

  double eta, mu;
  Utility::dLogConst(angle, eta, mu);
  
  Matrix3D dH = Eye3;
  dH.addSpin(u, -0.5);
  dH.addSpinSquare(u, eta);
  return dH;
}

inline Matrix3D 
ddLogSO3(const Vector3D& u, const Vector3D& p)
{
  //
  // -0.5*Hat(v) + eta*(Eye3*th.dot(v) + th.bun(v) - 2.*v.bun(th)) + mu*St2*v.bun(th);
  //
  double angle = u.norm();
  const Vector3D th = Utility::RescaleVector(u, angle);

  double eta, mu;
  Utility::dLogConst(angle, eta, mu);

  Matrix3D St2 = Hat(th);
  St2 = St2*St2;
  Matrix3D dH{};
  dH.addSpin(p, -0.5);  
  dH.addDiagonal( eta*th.dot(p));
  dH.addTensorProduct(th, p, eta);
  dH.addTensorProduct(p, th, -2.0*eta);
  dH.addMatrixProduct(St2, p.bun(th), mu);

  return dH*dLogSO3(th);
}
