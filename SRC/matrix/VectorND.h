//===----------------------------------------------------------------------===//
//
//                                   xara
//
//===----------------------------------------------------------------------===//
//                              https://xara.so
//===----------------------------------------------------------------------===//
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
#include <Vector.h>
#include <Matrix.h>

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

  inline int
  addMatrixVector(const double thisFact, const Matrix &m, const Vector &v, const double otherFact);

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
  dot(const VecT &other) const noexcept {
    T sum = 0.0;
    for (index_t i = 0; i < N; ++i)
      sum += values[i] * other[i];

    return sum;
  }

  // Tensor product, also known as the "bun" product
  template <int nc>
  constexpr inline MatrixND<N,nc,double>
  bun(const VectorND<nc> &other) const {
    if constexpr (N == 3 && nc == 3)
      return MatrixND<N,nc,double> {{
        {values[0]*other[0], values[1]*other[0], values[2]*other[0]},
        {values[0]*other[1], values[1]*other[1], values[2]*other[1]},
        {values[0]*other[2], values[1]*other[2], values[2]*other[2]}
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
  template <class VecB, class VecC> inline constexpr
  void 
  cross(const VecB& b, VecC& c) requires(N==3) const noexcept {
      c[0] = values[1] * b[2] - values[2] * b[1];
      c[1] = values[2] * b[0] - values[0] * b[2];
      c[2] = values[0] * b[1] - values[1] * b[0];
      return c;
  }


  template <class Vec3T> inline constexpr
  VectorND<N> 
  cross(const Vec3T& b) const noexcept requires(N==3) {
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
  
  inline double normalize() {
    double n = norm();

    if (n != 0.0)
      for (index_t i=0; i<N; i++)
        values[i] /= n;

    return n;
  }

  //
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
  operator/=(const double &right) {
    for (index_t i=0; i< N; i++)
      values[i] /= right;
    return *this;
  }

  inline
  VectorND<N>  
  operator/(const double &right) const {
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

  inline
  VectorND<N> &
  operator+=(const Vector &right) {
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
inline OpenSees::VectorND<N>
operator * (double a, const OpenSees::VectorND<N>& b) {
  return b * a;
}

#endif

