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
// Desctiption: MatrixND is a fixed-size matrix class that is suitable for
// stack-allocation.
//
// Claudio Perez
//
#pragma once
#include <cmath>
#include <assert.h>
#include <array>

#include "VectorND.h"
#include "Matrix.h"
#include "Vector.h"

namespace OpenSees {

template <index_t NR, index_t NC, typename T=double>
struct alignas(64) MatrixND {
  std::array<T, NR*NC> values;

  //
  // Indexing
  //
  inline constexpr       T& operator()(int i, int j)       noexcept { return values[j*NR + i]; }
  inline constexpr const T& operator()(int i, int j) const noexcept { return values[j*NR + i]; }

  constexpr T* data() noexcept { return &(*this)(0,0); }
  inline constexpr const T* data() const noexcept { return values.data(); }

  // Convert to regular Matrix class
  operator Matrix() { return Matrix(&(*this)(0,0), NR, NC);}

  operator const Matrix() const { return Matrix(&(*this)(0,0), NR, NC);}

  constexpr void zero() noexcept;

  constexpr MatrixND<NC, NR> transpose() const noexcept;

  constexpr MatrixND<NR,NC,T>& addDiagonal(const double vol) noexcept;

  template <class MatT>
    void addMatrix(const MatT& A, const double scale);

  template <class VecA, class VecB>  constexpr MatrixND<NR,NC,T>& 
    addTensorProduct(const VecA& V, const VecB& W, const double scale) noexcept;

  template <class MatT, int nk> void 
    addMatrixProduct(const MatrixND<NR, nk, T> &, const MatT&, double scale);

  template <class MatT, int nk> void 
    addMatrixProduct(double, const MatrixND<NR, nk, T> &, const MatT&, double scale);

  // += A'B
  template <class MatT, int nk> void
    addMatrixTransposeProduct(double thisFact, const MatrixND<nk, NR, T> &, const MatT&, double scale);

  // += A'BA
  template <int nk> int 
    addMatrixTripleProduct(double thisFact, 
                           const MatrixND<nk, NR, T> &, 
                           const MatrixND<nk, nk, T>&, 
                           double scale);

  // += A'BC
  template <int nk, int nl> int 
    addMatrixTripleProduct(double thisFact, 
                           const MatrixND<nk,NR> &A, 
                           const MatrixND<nk,nl> &B, 
                           const MatrixND<nl,NC> &C, double otherFact);

  template <typename F> void map(F func) const;
  template <typename F> void map(F func, MatrixND<NR,NC,T>& destination);

  template<class VecT> constexpr MatrixND<NR,NC,T>& addSpin(const VecT& V) noexcept;
  template<class VecT> constexpr MatrixND<NR,NC,T>& addSpin(const VecT& V, double scale) noexcept;
  template<class VecT> constexpr MatrixND<NR,NC,T>& addSpinSquare(const VecT& V, double scale) noexcept;
  template<class VecT> constexpr void addSpinProduct(const VecT& a, const VectorND<NR,T>& b, double scale) noexcept;
  template<class VecT> constexpr void 
    addMatrixSpinProduct(const MatrixND<NR,NC,T>& A, const VecT& b, double scale) noexcept;
  template<class MatT> constexpr void 
    addSpinMatrixProduct(const VectorND<NR,T>& a, const MatT& B, double scale) noexcept;

  int invert(MatrixND<NR, NC, T> &) const;

  template <int row0, int row1, int col0, int col1>
  inline MatrixND<row1-row0,col1-col0>
  extract() const noexcept;

  template<int er, int ec>
  inline MatrixND<er,ec>
  extract(int row0, int col0) const noexcept;

  template <int init_row, int init_col, int nr, int nc> 
  inline constexpr void
  insert(const MatrixND<nr, nc, double> &M, double fact) noexcept;

  template <int nr, int nc> 
  inline constexpr void
  insert(const MatrixND<nr, nc, double> &M, int init_row, int init_col, double fact) noexcept;

  template <int nr, int nc> 
  constexpr inline void
  assemble(const MatrixND<nr, nc, double> &M, int init_row, int init_col, double fact) noexcept;

  template <int nr> inline void
  assemble(const VectorND<nr> &v, int init_row, int init_col, double fact) noexcept;
  
  template<int nr, int nc> inline void
  assembleTranspose(const MatrixND<nr, nc, double> &M, int init_row, int init_col, double fact) noexcept;

  template<int nr> void
  assembleTranspose(const VectorND<nr> &v, int init_row, int init_col, double scale) noexcept;


//
// Operators
//

  constexpr MatrixND &
  operator=(const Matrix &other) noexcept
  {
    for (index_t j = 0; j < NC; ++j) {
      for (index_t i = 0; i < NR; ++i) {
        (*this)(i,j) = other(i,j);
      }
    }
    return *this;
  }

  constexpr MatrixND &
  operator+=(const double value) noexcept {
    for (index_t j = 0; j < NC; ++j) {
      for (index_t i = 0; i < NR; ++i) {
        (*this)(i,j) += value;
      }
    }
    return *this;
  }

  constexpr MatrixND &
  operator+=(const MatrixND &other) noexcept {
    for (index_t j = 0; j < NC; ++j) {
      for (index_t i = 0; i < NR; ++i) {
        (*this)(i,j) += other(i,j);
      }
    }
    return *this;
  }
  
  constexpr MatrixND &
  operator-=(const MatrixND &other) noexcept
  {
    for (index_t j = 0; j < NC; ++j)
      for (index_t i = 0; i < NR; ++i)
        (*this)(i,j) -= other(i,j);

    return *this;
  }

  inline constexpr MatrixND &
  operator*=(T const scalar) noexcept
  {
    for (index_t j = 0; j < NC; ++j)
      for (index_t i = 0; i < NR; ++i)
        (*this)(i,j) *= scalar;

    return *this;
  }

  inline constexpr MatrixND &
  operator/=(T const scalar) noexcept
  {
    for (index_t j = 0; j < NC; ++j)
      for (index_t i = 0; i < NR; ++i)
        (*this)(i,j) /= scalar;

    return *this;
  }

  inline constexpr VectorND<NC>
  operator^(const VectorND<NR> &V) const noexcept
  {
    VectorND<NC> result;

    const double *dataPtr = &(*this)(0,0);
    for (int i=0; i<NC; i++) {
      result[i] = 0.0;
      for (int j=0; j<NR; j++)
        result[i] += *dataPtr++ * V[j];
    }

    return result;
  }

  friend constexpr MatrixND
  operator+(MatrixND left, const MatrixND &right) noexcept {
    left += right; 
    return left;
  }

  friend constexpr MatrixND
  operator-(MatrixND left, const MatrixND &right) noexcept {
    left -= right; 
    return left;
  }
  
  friend constexpr MatrixND // scalar * Matrix
  operator*(T scalar, MatrixND mat) noexcept {
    mat *= scalar;
    return mat;
  }

  friend constexpr MatrixND // Matrix * scalar
  operator*(MatrixND mat, T scalar) noexcept {
    mat *= scalar;
    return mat;
  }

  template <index_t J>
  inline constexpr friend MatrixND<NR, J>
  operator*(const MatrixND<NR, NC> &left, const MatrixND<NC, J> &right) noexcept {
    MatrixND<NR, J> prod;
    if constexpr (NR*NC > 16)
      prod.addMatrixProduct(0, left, right, 1);
    else
      for (index_t i = 0; i < NR; ++i) {
        for (index_t j = 0; j < J; ++j) {
          prod(i, j) = 0.0;
          for (index_t k = 0; k < NC; ++k) {
            prod(i, j) += left(i,k) * right(k,j);
          }
        }
      }
    return prod;
  }

  template <index_t J>
  friend  MatrixND<NR, J>
  operator*(const MatrixND<NR, NC> &left, const Matrix &right) noexcept {
    MatrixND<NR, J> prod;
    for (index_t i = 0; i < NR; ++i) {
      for (index_t j = 0; j < J; ++j) {
        prod(i, j) = 0.0;
        for (index_t k = 0; k < NC; ++k) {
          prod(i, j) += left(i,k) * right(k,j);
        }
      }
    }
    return prod;
  }

  // Matrix*Vector
  constexpr friend  VectorND<NR>
  operator*(const MatrixND<NR, NC> &left, const VectorND<NC> &right) noexcept {
    VectorND<NR> prod;
    for (index_t i = 0; i < NR; ++i) {
      prod[i] = 0.0;
      for (index_t k = 0; k < NC; ++k) {
        prod[i] += left(i,k) * right[k];
      }
    }
    return prod;
  }

  friend  VectorND<NR>
  operator*(const MatrixND<NR, NC> &left, const Vector &right) noexcept {
    VectorND<NR> prod;
    for (index_t i = 0; i < NR; ++i) {
        prod[i] = 0.0;
        for (index_t k = 0; k < NC; ++k) {
          prod[i] += left(i,k) * right(k);
        }
    }
    return prod;
  }

  template <index_t K>
  inline constexpr friend MatrixND<NC,K>
  operator^(const MatrixND<NR, NC> &left, const MatrixND<NR, K> &right) {
    MatrixND<NC, K> prod;
    if constexpr (NR*NC > 48)
      prod.addMatrixTransposeProduct(0.0, left, right, 1.0);
    else {
      for (index_t i = 0; i < NC; ++i) {
        for (index_t j = 0; j < K; ++j) {
          prod(i, j) = 0.0;
          for (index_t k = 0; k < NR; ++k) {
            prod(i, j) += left(k,i) * right(k,j);
          }
        }
      }
    }
    return prod;
  }

  friend constexpr MatrixND
  operator/(MatrixND mat, T scalar) {
    mat /= scalar; 
    return mat;
  }
};

} // namespace OpenSees
#include "MatrixND.tpp"
