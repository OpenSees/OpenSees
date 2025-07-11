//===----------------------------------------------------------------------===//
//
//                                   xara
//
//===----------------------------------------------------------------------===//
//                              https://xara.so
//===----------------------------------------------------------------------===//
// 
// Desctiption: MatrixND is a fixed-size matrix class that is suitable for
// stack-allocation.
//
//
// This code is influenced by the following sources
//  List initialization:
//  - https://stackoverflow.com/questions/42068882/list-initialization-for-a-matrix-class
//
//  Style/practices
//  - https://quuxplusone.github.io/blog/2021/04/03/static-constexpr-whittling-knife/
// 
//  Operator overloading / semantics
//  - https://stackoverflow.com/questions/9851188/does-it-make-sense-to-use-move-semantics-for-operator-and-or-operator/9851423#9851423
//
//  Compile-time template restrictions/concepts:
//  - https://codereview.stackexchange.com/questions/259038/compile-time-matrix-class
//    (C++ 20)
//  - https://github.com/calebzulawski/cotila/
//    (C++ 17)
//
//===----------------------------------------------------------------------===//
//
// Claudio Perez
//
#pragma once
#include <math.h>
#include <assert.h>
#include <array>

#include "VectorND.h"
#include "Matrix.h"
#include "Vector.h"

namespace OpenSees {

template <index_t NR, index_t NC, typename T=double>
struct alignas(64) MatrixND {
  // double values[NC][NR];
  std::array<T, NR*NC> values;

  //
  // Indexing
  //
  // (i,j) indexing
  // constexpr T &
  // operator()(index_t index_r, index_t index_c) noexcept {
  //   assert(index_r >= 0 && index_c >= 0);
  //   assert(index_r < NR && index_c < NC);
  //   return values[index_c][index_r];
  // }

  // inline constexpr const T & 
  // operator()(index_t index_r, index_t index_c) const noexcept {
  //   assert(index_r >= 0 && index_c >= 0);
  //   assert(index_r < NR && index_c < NC);
  //   return values[index_c][index_r];
  // }
  constexpr       T& operator()(int i, int j)       noexcept { return values[j*NR + i]; }
  constexpr const T& operator()(int i, int j) const noexcept { return values[j*NR + i]; }

  constexpr T* data() noexcept { return &(*this)(0,0); }


  // Convert to regular Matrix class
  operator Matrix() { return Matrix(&(*this)(0,0), NR, NC);}

  operator const Matrix() const { return Matrix(&(*this)(0,0), NR, NC);}

  int symeig(VectorND<NR>& vals);

  constexpr void zero() noexcept;

  constexpr double determinant() const ;

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

  //
  //
  //

  constexpr T
  trace() const noexcept
  {
    static_assert(NR == NC);
    T sum = 0.0;
    for (index_t i = 0; i < NR; ++i) {
      sum += (*this)(i,i);
    }
    return sum;
  }


  int
  solve(const VectorND<NR> &V, VectorND<NR> &res) const noexcept
  {
    static_assert(NR == NC);

    MatrixND<NR,NC> work = *this;
    int pivot_ind[NR];
    int nrhs = 1;
    int nr = NR;
    int nc = NC;
    int info = 0;
    res = V; // X will be overwritten with the solution
    DGESV(&nr, &nrhs, &work(0,0), &nr, pivot_ind, res.values, &nc, &info);
    return -abs(info);
  }


  template<index_t n>
  int solve(const MatrixND<n, n>& M, MatrixND<n, n>& X) const noexcept
  {
    static_assert(NR == NC, "Matrix must be square.");
    static_assert(n == NR, "RHS row-count must match A.");

    MatrixND<NR,NC,T> work = *this;               // copy of A to be factorised
    int ipiv[NR]{};

    int n_eq  = NR;               // order of the system
    int nrhs  = n;                // number of RHS columns
    int lda   = NR;               // leading dim of A
    int ldb   = NR;               // leading dim of X
    int info  = 0;

    X = M;                               // copy RHS, DGESV overwrites
    DGESV(&n_eq, &nrhs,
          &work(0,0), &lda,
          &ipiv[0],
          &X(0,0), &ldb,
          &info);

    return -std::abs(info);
  }

  // int solve(const Vector &V, Vector &res) const
  //   requires(NR == NC)
  // {

  //   MatrixND<NR,NC> work = *this;
  //   int pivot_ind[NR];
  //   int nrhs = 1;
  //   int nr = NR;
  //   int nc = NC;
  //   int info = 0;
  //   res = V; // X will be overwritten with the solution
  //   DGESV(&nr, &nrhs, &work.values[0][0], &nr, &pivot_ind[0], res.theData, &nc, &info);
  //   return -abs(info);
  // }


  // int
  // solve(const Matrix &M, Matrix &res)
  // {
  //   Matrix slver(*this);
  //   return slver.Solve(M, res);

  //   MatrixND<NR,NC> work = *this;
  //   int pivot_ind[NR];
  //   int nrhs = M.noCols();
  //   int nr = NR;
  //   int nc = NC;
  //   int info = 0;
  //   res = M; // M will be overwritten with the solution
  //   DGESV(&nr, &nrhs, &work(0,0), &nr, &pivot_ind[0], &res(0,0), &nc, &info);
  //   return -abs(info);
  // }
 

  template <int row0, int row1, int col0, int col1>
  inline MatrixND<row1-row0,col1-col0>
  extract() const noexcept
  {
    MatrixND<row1-row0,col1-col0> m;
    for (int i=0; i<row1-row0; i++)
      for (int j=0; j<col1-col0; j++)
        m(i,j) = (*this)(row0+i, col0+j);
    return m;
  }

  template<int er, int ec>
  inline MatrixND<er,ec>
  extract(int row0, int col0) const noexcept
  {
    MatrixND<er,ec> m;
    for (int i=0; i<er; i++)
      for (int j=0; j<ec; j++)
        m(i,j) = (*this)(row0+i, col0+j);
    return m;
  }

  template <int init_row, int init_col, int nr, int nc> 
  inline constexpr void
  insert(const MatrixND<nr, nc, double> &M, double fact) noexcept
  {
 
    constexpr int final_row = init_row + nr - 1;
    constexpr int final_col = init_col + nc - 1;
    static_assert((init_row >= 0) && (final_row < NR) && (init_col >= 0) && (final_col < NC), 
                  "MatrixND::insert: init_row, init_col, nr, nc out of bounds");

    for (int i=0; i<nc; i++) {
       int pos_Cols = init_col + i;
       for (int j=0; j<nr; j++) {
          int pos_Rows = init_row + j; 
          (*this)(pos_Rows,pos_Cols) = M(j,i)*fact;
       }
    }
  }

  template <int nr, int nc> 
  inline constexpr void
  insert(const MatrixND<nr, nc, double> &M, int init_row, int init_col, double fact) noexcept
  {
 
    [[maybe_unused]] int final_row = init_row + nr - 1;
    [[maybe_unused]] int final_col = init_col + nc - 1; 
    assert((init_row >= 0) && (final_row < NR) && (init_col >= 0) && (final_col < NC));

    for (int i=0; i<nc; i++) {
       int pos_Cols = init_col + i;
       for (int j=0; j<nr; j++) {
          int pos_Rows = init_row + j; 
          (*this)(pos_Rows,pos_Cols) = M(j,i)*fact;
       }
    }
  }

  template <int nr, int nc> 
  constexpr inline void
  assemble(const MatrixND<nr, nc, double> &M, int init_row, int init_col, double fact) noexcept
  {
 
    [[maybe_unused]] int final_row = init_row + nr - 1;
    [[maybe_unused]] int final_col = init_col + nc - 1; 
    assert((init_row >= 0) && (final_row < NR) && (init_col >= 0) && (final_col < NC));

    for (int i=0; i<nc; i++) {
       int pos_Cols = init_col + i;
       for (int j=0; j<nr; j++) {
          (*this)(init_row + j, pos_Cols) += M(j,i)*fact;
       }
    }
  }

  template <int nr> inline void
  assemble(const VectorND<nr> &v, int init_row, int init_col, double fact) noexcept;
  
  template<int nr, int nc> inline void
  assembleTranspose(const MatrixND<nr, nc, double> &M, int init_row, int init_col, double fact) noexcept;

  template<int nr> void
  assembleTranspose(const VectorND<nr> &v, int init_row, int init_col, double scale) noexcept;


//
// Operators
//
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
