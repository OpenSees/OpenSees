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
//----------------------------------------------------------------------------//
//
// Claudio Perez
//
#pragma once
#include <math.h>
#include <assert.h>
#include <array>
#include <stdexcept>

#include "VectorND.h"
#include "Matrix.h"
#include "Vector.h"

#if __cplusplus < 202000L
#  define consteval
#  define requires(X)
#endif

#define G23_STACK_MAX 10

namespace OpenSees {

template <index_t NR, index_t NC, typename T=double>
requires(NR > 0 && NC > 0)
struct MatrixND {
  double values[NC][NR];

  // Convert to dynamic Matrix class
  operator Matrix() { return Matrix(&values[0][0], NR, NC);}

  operator const Matrix() const { return Matrix(&values[0][0], NR, NC);}

  consteval void zero();

  constexpr MatrixND<NC, NR> transpose() const;

  MatrixND<NR,NC,T>& addDiagonal(const double vol) requires(NR == NC);

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

  template<class VecT> MatrixND<NR,NC,T>& addSpin(const VecT& V);
  template<class VecT> MatrixND<NR,NC,T>& addSpin(const VecT& V, double scale);
  template<class VecT> MatrixND<NR,NC,T>& addSpinSquare(const VecT& V, double scale);
  template<class VecT> void addSpinProduct(const VecT& a, const VectorND<NR,T>& b, double scale);
  template<class VecT> void addMatrixSpinProduct(const MatrixND<NR,NC,T>& A, const VecT& b, double scale);
  template<class MatT> void addSpinMatrixProduct(const VectorND<NR,T>& a, const MatT& B, double scale);

  int invert(MatrixND<NR, NC, T> &) const;
  int invert() {
    return Matrix(*this).Invert();
  }

  //
  // Indexing
  //
  constexpr std::array<T, NC> &
  operator[](index_t index) 
  {
    return values[index];
  }

  constexpr const std::array<T, NC> &
  operator[](index_t index) const 
  {
    return values[index];
  }
  
  // (i,j) indexing
  constexpr T &
  operator()(index_t index_r, index_t index_c) {
    assert(index_r >= 0 && index_c >= 0);
    assert(index_r < NR && index_c < NC);
    return values[index_c][index_r];
  }

  inline constexpr const T & 
  operator()(index_t index_r, index_t index_c) const {
    assert(index_r >= 0 && index_c >= 0);
    assert(index_r < NR && index_c < NC);
    return values[index_c][index_r];
  }


  constexpr MatrixND &
  operator=(const Matrix &other)
  {
    for (index_t j = 0; j < NC; ++j) {
      for (index_t i = 0; i < NR; ++i) {
        values[j][i] = other(i,j);
      }
    }
    return *this;
  }

  constexpr MatrixND &
  operator+=(const double value) {
    for (index_t j = 0; j < NC; ++j) {
      for (index_t i = 0; i < NR; ++i) {
        values[j][i] += value;
      }
    }
    return *this;
  }

  constexpr MatrixND &
  operator+=(const MatrixND &other) {
    for (index_t j = 0; j < NC; ++j) {
      for (index_t i = 0; i < NR; ++i) {
        values[j][i] += other.values[j][i];
      }
    }
    return *this;
  }
  
  constexpr MatrixND &
  operator-=(const MatrixND &other) 
  {
    for (index_t j = 0; j < NC; ++j)
      for (index_t i = 0; i < NR; ++i)
        values[j][i] -= other.values[j][i];

    return *this;
  }

  inline constexpr MatrixND &
  operator*=(T const scalar) 
  {
    for (index_t j = 0; j < NC; ++j)
      for (index_t i = 0; i < NR; ++i)
        values[j][i] *= scalar;

    return *this;
  }

  inline constexpr MatrixND &
  operator/=(T const scalar) 
  {
    for (index_t j = 0; j < NC; ++j)
      for (index_t i = 0; i < NR; ++i)
        values[j][i] /= scalar;

    return *this;
  }

  inline constexpr VectorND<NC>
  operator^(const VectorND<NR> &V) const
  {
      VectorND<NC> result;

      const double *dataPtr = &values[0][0];
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
  trace() const
  requires(NR == NC) 
  {
    T sum = 0.0;
    for (index_t i = 0; i < NR; ++i) {
      sum += values[i][i];
    }
    return sum;
  }


  int solve(const VectorND<NR> &V, VectorND<NR> &res) const
    requires(NR == NC)
  {
    MatrixND<NR,NC> work = *this;
    int pivot_ind[NR];
    int nrhs = 1;
    int nr = NR;
    int nc = NC;
    int info = 0;
    res = V; // X will be overwritten with the solution
    DGESV(&nr, &nrhs, &work.values[0][0], &nr, pivot_ind, res.values, &nc, &info);
    return -abs(info);
  }


  int solve(const Vector &V, Vector &res) const
    requires(NR == NC)
  {

    MatrixND<NR,NC> work = *this;
    int pivot_ind[NR];
    int nrhs = 1;
    int nr = NR;
    int nc = NC;
    int info = 0;
    res = V; // X will be overwritten with the solution
    DGESV(&nr, &nrhs, &work.values[0][0], &nr, &pivot_ind[0], res.theData, &nc, &info);
    return -abs(info);
  }


  int
  solve(const Matrix &M, Matrix &res)
  {
    Matrix slver(*this);
    return slver.Solve(M, res);

    MatrixND<NR,NC> work = *this;
    int pivot_ind[NR];
    int nrhs = M.noCols();
    int nr = NR;
    int nc = NC;
    int info = 0;
    res = M; // M will be overwritten with the solution
    DGESV(&nr, &nrhs, &work(0,0), &nr, &pivot_ind[0], &res(0,0), &nc, &info);
    return -abs(info);
  }
 

  template <int row0, int row1, int col0, int col1>
  inline MatrixND<row1-row0,col1-col0>
  extract() const
  {
    MatrixND<row1-row0,col1-col0> m;
    for (int i=0; i<row1-row0; i++)
      for (int j=0; j<col1-col0; j++)
        m(i,j) = (*this)(row0+i, col0+j);
    return m;
  }

  template<int er, int ec>
  inline MatrixND<er,ec>
  extract(int row0, int col0) const
  {
    MatrixND<er,ec> m;
    for (int i=0; i<er; i++)
      for (int j=0; j<ec; j++)
        m(i,j) = (*this)(row0+i, col0+j);
    return m;
  }

  template <int init_row, int init_col, int nr, int nc> inline void
  insert(const MatrixND<nr, nc, double> &M, double fact) 
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

  template <int nr, int nc> inline void
  insert(const MatrixND<nr, nc, double> &M, int init_row, int init_col, double fact) 
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

  template <int nr, int nc> inline void
  assemble(const MatrixND<nr, nc, double> &M, int init_row, int init_col, double fact) 
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
  operator+(MatrixND left, const MatrixND &right) {
    left += right; 
    return left;
  }

  friend constexpr MatrixND
  operator-(MatrixND left, const MatrixND &right) {
    left -= right; 
    return left;
  }
  
  friend constexpr MatrixND // scalar * Matrix
  operator*(T scalar, MatrixND mat) {
    mat *= scalar;
    return mat;
  }

  friend constexpr MatrixND // Matrix * scalar
  operator*(MatrixND mat, T scalar) {
    mat *= scalar;
    return mat;
  }

  template <index_t J>
  inline constexpr friend MatrixND<NR, J>
  operator*(const MatrixND<NR, NC> &left, const MatrixND<NC, J> &right) {
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

  template <index_t J>
  friend  MatrixND<NR, J>
  operator*(const MatrixND<NR, NC> &left, const Matrix &right) {
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

  template <index_t K>
  inline constexpr friend MatrixND<NC,K>
  operator^(const MatrixND<NR, NC> &left, const MatrixND<NR, K> &right) {
    MatrixND<NC, K> prod;
    if constexpr (0 && NR*NC > 16)
      prod.addMatrixTransposeProduct(0.0, left, right, 1.0);
    else
      for (index_t i = 0; i < NC; ++i) {
        for (index_t j = 0; j < K; ++j) {
          prod(i, j) = 0.0;
          for (index_t k = 0; k < NR; ++k) {
            prod(i, j) += left(k,i) * right(k,j);
          }
        }
      }
    return prod;
  }


  constexpr friend  VectorND<NR>
  operator*(const MatrixND<NR, NC> &left, const VectorND<NC> &right) {
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
  operator*(const MatrixND<NR, NC> &left, const Vector &right) {
    VectorND<NR> prod;
    for (index_t i = 0; i < NR; ++i) {
        prod[i] = 0.0;
        for (index_t k = 0; k < NC; ++k) {
          prod[i] += left(i,k) * right(k);
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
