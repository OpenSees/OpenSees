//===----------------------------------------------------------------------===//
//
//        OpenSees - Open System for Earthquake Engineering Simulation
//
//===----------------------------------------------------------------------===//
//
// 
//  Objectives:
//  - little to no overhead above C-style arrays
//  - value semantics; objects do not decay to pointers;
//
//  This code is influenced by the following sources
//   list initialization:
//   - https://stackoverflow.com/questions/42068882/list-initialization-for-a-matrix-class
//
//   style/practices
//   - https://quuxplusone.github.io/blog/2021/04/03/static-constexpr-whittling-knife/
//  
//   Operator overloading / semantics
//   - https://stackoverflow.com/questions/9851188/does-it-make-sense-to-use-move-semantics-for-operator-and-or-operator/9851423#9851423
//
//   compile-time template restrictions/concepts:
//   - https://codereview.stackexchange.com/questions/259038/compile-time-matrix-class
//     (C++ 20)
//   - https://github.com/calebzulawski/cotila/
//     (C++ 17)
//
//===----------------------------------------------------------------------===//
//
// Claudio Perez
//
#ifndef VectorND_H
#define VectorND_H
#include <math.h>
#include <assert.h>
#include <array>
#include <stdexcept>
#include <functional>
#include <iostream> // overloading <<
#include <Vector.h>
#include <Matrix.h>
#include "blasdecl.h"

#if __cplusplus < 202000L
#  define consteval
#  define requires(X)
#endif

namespace OpenSees {

typedef int index_t;

template<int n, int m, typename T> struct MatrixND;

template <index_t N, typename T=double> 
requires(N > 0)
struct VectorND {
  T values[N];

  operator Vector() { return Vector(values, N);}

  template<int n, int m, typename> friend struct MatrixND;

  inline constexpr T&
  operator[](index_t index) {
    return values[index];
  }

  inline constexpr const T&
  operator[](index_t index) const {
    return values[index];
  }

  inline constexpr T&
  operator()(index_t index) {
    return values[index];
  }

  inline constexpr const T&
  operator()(index_t index) const {
    return values[index];
  }

  consteval int
  size() const {
    return N;
  }


  consteval inline void
  fill(double value) {
    for (T& item : values)
      item = value;
  }

  consteval inline void
  zero() {
    for (T& item : values )
      item = 0.0;
  }

  template<typename VecT> inline
  constexpr T
  dot(const VecT &other) const {
    T sum = 0.0;
    for (index_t i = 0; i < N; ++i) {
      sum += values[i] * other[i];
    }
    return sum;
  }

  // Tensor product, also known as the "bun" product
  template <int nc>
  inline OpenSees::MatrixND<N,nc,double>
  bun(const VectorND<nc> &other) const {
    OpenSees::MatrixND<N,nc,double> prod;

    for (int j = 0; j < other.size(); ++j)
      for (int i = 0; i < this->size(); ++i)
        prod(i,j) = values[i] * other.values[j];

    return prod;
  }

  constexpr T
  norm() const {
    return sqrt(this->dot(*this));
  }
  
  inline double normalize() {
    double n = norm();

    for (int i=0; i<N; i++)
      values[i] /= n;

    return n;
  }

  template<typename VecT> inline
  VectorND<N> &operator=(const VecT &right) {
    for (int i=0; i< N; i++)
      values[i] = right[i];
    return *this;
  }

  inline
  VectorND<N> &operator/=(const double &right) {
    for (int i=0; i< N; i++)
      values[i] /= right;
    return *this;
  }

  inline
  VectorND<N>  operator/(const double &right) const {
    VectorND<N> res(*this);
    res /= right;
    return res;
  }

  inline
  VectorND<N> &operator*=(const double &right) {
    for (int i=0; i< N; i++)
      values[i] *= right;
    return *this;
  }

  inline
  VectorND<N>  operator*(const double &right) const {
    VectorND<N> res(*this);
    res *= right;
    return res;
  }

  inline
  VectorND<N> &operator+=(const VectorND<N> &right) {
    for (int i=0; i< N; i++)
      values[i] += right[i];
    return *this;
  }

  VectorND<N> &operator+=(const Vector &right) {
    assert(right.Size() == N);
    for (int i=0; i< N; i++)
      values[i] += right[i];
    return *this;
  }

  template <class VecT>
  VectorND<N> operator+(const VecT &right) const {
    VectorND<N> res {*this};
    res += right;
    return res;
  }

  template <class VecT>
  VectorND<N> &operator-=(const VecT &right) {
    for (int i=0; i< N; i++)
      values[i] -= right[i];
    return *this;
  }

  template <class VecT>
  VectorND<N> operator-(const VecT &right) const {
    VectorND<N> res {*this};
    res -= right;
    return res;
  }

  // Return the cross product this vector with another vector, b.
  template <class VecB, class VecC> inline 
  void cross(const VecB& b, VecC& c) const requires(N==3) {
      c[0] = values[1] * b[2] - values[2] * b[1];
      c[1] = values[2] * b[0] - values[0] * b[2];
      c[2] = values[0] * b[1] - values[1] * b[0];
      return c;
  }


  template <class Vec3T> inline
  VectorND<N> cross(const Vec3T& b) const requires(N==3) {
      VectorND<3> c;
      c[0] = values[1] * b[2] - values[2] * b[1];
      c[1] = values[2] * b[0] - values[0] * b[2];
      c[2] = values[0] * b[1] - values[1] * b[0];
      return c;
  }

  friend std::ostream &
  operator<<(std::ostream &out, const VectorND &vec) {
    out << "{";
    for (int r=0; r<N; r++){
        out << vec[r] << ( r < N-1? ", ": "");
    }
    return out << "}\n";
  }

#include "VectorND.tpp"

};
} // namespace OpenSees

template<int N>
inline OpenSees::VectorND<N>
operator * (double a, const OpenSees::VectorND<N>& b) {
    return b * a;
}
#endif

