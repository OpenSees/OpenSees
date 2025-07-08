//===----------------------------------------------------------------------===//
//
//                                   xara
//
//===----------------------------------------------------------------------===//
//                              https://xara.so
//===----------------------------------------------------------------------===//
//
// Description: Versor implements a unit-normalized quaternion.
//
#pragma once
#include <cmath>
#include <type_traits>
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


} //  namespace OpenSees

inline OpenSees::Versor
operator*(const OpenSees::Versor &qa, const OpenSees::Versor &qb)
{
  OpenSees::Versor q12;
  q12.scalar = qa.scalar * qb.scalar - qa.vector.dot(qb.vector);
  q12.vector = (qb.vector * qa.scalar) + (qa.vector * qb.scalar) + qa.vector.cross(qb.vector);
  return q12;
}
