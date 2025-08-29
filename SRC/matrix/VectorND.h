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
// Claudio Perez
//
#ifndef VectorND_H
#define VectorND_H
#include <cmath>
#include <assert.h>
#include <Vector.h>
#include <Matrix.h>


namespace OpenSees {

typedef int index_t;

template<int n, int m, typename T> struct MatrixND;

template <index_t N, typename T=double> 
struct VectorND {
  T values[N];

  operator Vector() { return Vector(values, N);}

  template<int n, int m, typename> friend struct MatrixND;

  //
  template <int ir, int nr> inline void
  assemble(const VectorND<nr> &v, double fact=1);

  template<int nr> void
  assemble(int a, const VectorND<nr>& v, double scale) noexcept;

  template <int ir, int nr> inline void
  insert(const VectorND<nr> &v, double fact=1);

  template<int nr> void
  insert(int a, const VectorND<nr>& v, double scale) noexcept;

  template<int nr>
  VectorND<nr> 
  extract(int a) noexcept;

  int
  addVector(const T thisFact, const VectorND<N> &other, const T otherFact) noexcept;

  VectorND<N> &addCross(const VectorND<N>& a, const VectorND<N> &b, double fact = 1.0) {
    static_assert(N == 3);
    values[0] += fact * (a.values[1] * b.values[2] - a.values[2] * b.values[1]);
    values[1] += fact * (a.values[2] * b.values[0] - a.values[0] * b.values[2]);
    values[2] += fact * (a.values[0] * b.values[1] - a.values[1] * b.values[0]);
    return *this;
  }

#ifdef XARA_VECTOR_FRIENDS
  template <int NC>
  inline int
  addMatrixVector(double thisFact, const MatrixND<N, NC, double> &m, const Vector& v, double otherFact);

  template <int NR>
  inline int
  addMatrixTransposeVector(double thisFact, const MatrixND<NR, N, double> &m, 
                           const Vector &v, double otherFact);

  inline int
  addMatrixVector(const double thisFact, const Matrix &m, const Vector &v, const double otherFact);

  int
  addVector(const T thisFact, const Vector &other, const T otherFact) noexcept;
#endif

  constexpr int
  size() const {
    return N;
  }

  constexpr inline void
  fill(double value) {
    for (T& item : values)
      item = value;
  }

  constexpr inline void
  zero() {
    for (T& item : values )
      item = 0.0;
  }

  template<typename VecT>
  constexpr T
  dot(const VecT &other) const noexcept {
    T sum = 0.0;
    for (index_t i = 0; i < N; ++i)
      sum += values[i] * other[i];

    return sum;
  }

  // Tensor product, also known as the "bun" product
  template <int nc>
  constexpr inline MatrixND<N,nc,double>
  bun(const VectorND<nc> &other) const noexcept {
    if constexpr (N == 3 && nc == 3)
      return MatrixND<N,nc,double> {{
         values[0]*other[0], values[1]*other[0], values[2]*other[0] ,
         values[0]*other[1], values[1]*other[1], values[2]*other[1] ,
         values[0]*other[2], values[1]*other[2], values[2]*other[2] 
      }};

    else {
      MatrixND<N,nc,double> prod;

      for (index_t j = 0; j < other.size(); ++j)
        for (index_t i = 0; i < this->size(); ++i)
          prod(i,j) = values[i] * other.values[j];

      return prod;
    }
  }

  // Return the cross product this vector with another vector, b.
  template <class VecB, class VecC>
  constexpr void 
  cross(const VecB& b, VecC& c) const noexcept {
    static_assert(N == 3, "Cross product is only defined for 3D vectors.");
    c[0] = values[1] * b[2] - values[2] * b[1];
    c[1] = values[2] * b[0] - values[0] * b[2];
    c[2] = values[0] * b[1] - values[1] * b[0];
    return c;
  }


  template <class Vec3T>
  constexpr VectorND<N> 
  cross(const Vec3T& b) const noexcept {
    static_assert(N == 3, "Cross product is only defined for 3D vectors.");
    // Return a new vector that is the cross product of this vector and b.
    VectorND<3> c;
    c[0] = values[1] * b[2] - values[2] * b[1];
    c[1] = values[2] * b[0] - values[0] * b[2];
    c[2] = values[0] * b[1] - values[1] * b[0];
    return c;
  }


  constexpr double
  norm() const noexcept {
    return std::sqrt(std::fabs(this->dot(*this)));
  }
  
  inline double 
  normalize() {
    double n = norm();

    if (n != 0.0)
      for (index_t i=0; i<N; i++)
        values[i] /= n;

    return n;
  }

  //
  inline constexpr T&
  operator[](index_t index) noexcept {
    assert(index >= 0 && index < N);
    return values[index];
  }

  inline constexpr const T&
  operator[](index_t index) const noexcept {
    assert(index >= 0 && index < N);
    return values[index];
  }

  inline constexpr T&
  operator()(index_t index) noexcept {
    return values[index];
  }

  inline constexpr const T&
  operator()(index_t index) const noexcept {
    return values[index];
  }

  template<typename VecT> inline
  VectorND<N> &operator=(const VecT &right) noexcept {
    for (index_t i=0; i< N; i++)
      values[i] = right[i];
    return *this;
  }

  inline
  VectorND<N> &
  operator/=(const double &right) noexcept {
    for (index_t i=0; i< N; i++)
      values[i] /= right;
    return *this;
  }

  inline
  VectorND<N>  
  operator/(const double &right) const noexcept {
    VectorND<N> res(*this);
    res /= right;
    return res;
  }

  inline constexpr
  VectorND<N> &
  operator*=(const double &right) noexcept{
    for (int i=0; i< N; i++)
      values[i] *= right;
    return *this;
  }

  inline constexpr
  VectorND<N>  
  operator*(const double &right) const noexcept {
    VectorND<N> res(*this);
    res *= right;
    return res;
  }

  inline constexpr
  VectorND<N> &
  operator+=(const VectorND<N> &right) noexcept {
    for (int i=0; i< N; i++)
      values[i] += right[i];
    return *this;
  }

  constexpr inline
  VectorND<N> &
  operator+=(const Vector &right) noexcept {
    assert(right.Size() == N);
    for (int i=0; i< N; i++)
      values[i] += right[i];
    return *this;
  }

  template <class VecT>
  VectorND<N> operator+(const VecT &right) const noexcept {
    VectorND<N> res {*this};
    res += right;
    return res;
  }

  template <class VecT> 
  constexpr inline
  VectorND<N> &operator-=(const VecT &right) noexcept {
    for (int i=0; i< N; i++)
      values[i] -= right[i];
    return *this;
  }

  template <class VecT>
  constexpr inline
  VectorND<N> operator-(const VecT &right) const noexcept {
    VectorND<N> res {*this};
    res -= right;
    return res;
  }

};
} // namespace OpenSees

#include "VectorND.tpp"

template<int N>
constexpr inline OpenSees::VectorND<N>
operator * (double a, const OpenSees::VectorND<N>& b) noexcept {
  return b * a;
}

#endif

