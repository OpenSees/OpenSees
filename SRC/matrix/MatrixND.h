//===----------------------------------------------------------------------===//
//
//        OpenSees - Open System for Earthquake Engineering Simulation    
//    See https://opensees.berkeley.edu/OpenSees/copyright.php for license.
//
//===--- MatrixND.h - Matrix with fixed size ------------------------------===//
//
// Objectives:
// - little to no overhead above C-style arrays
// - value semantics; objects do not decay to pointers;
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
#ifndef MatrixND_H
#define MatrixND_H
#include <math.h>
#include <assert.h>
#include <array>
#include <stdexcept>
#include <functional>
#include <iostream> // overloading <<

#include "VectorND.h"
#include "Matrix.h"
#include "Vector.h"
#include "blasdecl.h"
#include "routines/cmx.h"
#include "routines/SY3.h"

#if __cplusplus < 202000L
#define consteval
#define requires(X)
#endif

#define G23_STACK_MAX 10

namespace OpenSees {

template <index_t NR, index_t NC, typename T=double>
requires(NR > 0 && NC > 0)
struct MatrixND {
  double values[NC][NR];

//MatrixND<NR, NC, T>(const MatrixND<NR, NC, T>&) = default;

  // Convert to regular Matrix class
  operator Matrix() { return Matrix(&values[0][0], NR, NC);}
  operator const Matrix() const { return Matrix(&values[0][0], NR, NC);}

  MatrixND<NR,NC,T>& addDiagonal(const double vol) requires(NR == NC);


  MatrixND<NR,NC,T> bun(const VectorND<NR,T>& a, const VectorND<NC,T> &b)
    requires(NR == NC == 3)
  {
    return MatrixND<NR,NC,T> {{{
      {a[0]*b[0], a[1]*b[0], a[2]*b[0]},
      {a[0]*b[1], a[1]*b[1], a[2]*b[1]},
      {a[0]*b[2], a[1]*b[2], a[2]*b[2]}
    }}};
  }

  int symeig(VectorND<NR>& vals) requires(NR == NC == 3) {
    double work[3][3];
    cmx_eigSY3(values, work, vals.values);
    return 0;
  }

  template <class MatT>
    void addMatrix(const MatT& A, const double scale);

  template <class VecA, class VecB>  MatrixND<NR,NC,T>& 
    addTensorProduct(const VecA& V, const VecB& W, const double scale);

  template <class MatT, int nk> void 
    addMatrixProduct(const MatrixND<NR, nk, T> &, const MatT&, double scale);

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

  int invert(MatrixND<NR, NC, T> &) const;
  int invert() {
    return Matrix(*this).Invert();
  }

  template<class VecT> MatrixND<NR,NC,T>& addSpin(const VecT& V);
  template<class VecT> MatrixND<NR,NC,T>& addSpin(const VecT& V, double scale);
  template<class VecT> MatrixND<NR,NC,T>& addSpinSquare(const VecT& V, double scale);
  template<class VecT> void addSpinProduct(const VecT& a, const VectorND<NR,T>& b, double scale);
  template<class VecT> void addMatrixSpinProduct(const MatrixND<NR,NC,T>& A, const VecT& b, double scale);
  template<class MatT> void addSpinMatrixProduct(const VectorND<NR,T>& a, const MatT& B, double scale);

//template<class VecT>
//void addSpinAtRow(const VecT& V, size_t row_index);
//template<class VecT>
//void addSpinAtRow(const VecT& V, size_t vector_index, size_t matrix_row_index);
//template<class VecT>
//MatrixND<NR,NC,T>& addSpin(const VecT& V, double mult) requires(NR == 3);
//template<class VecT>
//void addSpinAtRow(const VecT& V, double mult, size_t row_index);
//template<class VecT>
//void addSpinAtRow(const VecT& V, double mult, size_t vector_index, size_t matrix_row_index);

  //
  // Indexing
  //
  constexpr std::array<T, NC> &
  operator[](index_t index) {return values[index];}

  constexpr const std::array<T, NC> & // [i] indexing
  operator[](index_t index) const {return values[index];}
  
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

#if 0
  constexpr VectorND<NR>
  column(index_t index) const {
    assert(index >= 0);
    assert(index < NC);
    // TODO: ugly temporary implementation
    return *(VectorND<NR>*)(&(values[index]));
  }

  constexpr VectorND<NC>
  row(index_t index) const {
    assert(index >= 0);
    assert(index < NR);

    VectorND<NC,T> rw;
    for (index_t j = 0; j < NC; ++j) {
      rw[j] = values[j][index];
    }
    return rw;
  }

  // TODO: change to std::pair,array, or tuple.
  consteval VectorND<2,int>
  size() const {
    return {NR, NC};
  }
#endif

  //
  //
  //
  consteval void zero() {
    for (index_t j = 0; j < NC; ++j) {
      for (index_t i = 0; i < NR; ++i) {
        values[j][i] = 0.0;
      }
    }
  }

  constexpr MatrixND<NC, NR>
  transpose() const {
    MatrixND<NC, NR> result = {};
    for (index_t j = 0; j < NC; ++j) {
      for (index_t i = 0; i < NR; ++i) {
        result.values[i][j] = values[j][i];
      }
    }
    return result;
  }

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

  constexpr T
  determinant() const
  requires(NR == NC && NR == 2)
  {
    return values[0][0] * values[1][1] - values[0][1] * values[1][0];
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

//template <typename MatrixType> 
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

  template <int nr, int nc> inline void
  assemble(const MatrixND<nr, nc, double> &M, int init_row, int init_col, double fact) 
  {
 
    [[maybe_unused]] int final_row = init_row + nr - 1;
    [[maybe_unused]] int final_col = init_col + nc - 1; 
    assert((init_row >= 0) && (final_row < NR) && (init_col >= 0) && (final_col < NC));

    for (int i=0; i<nc; i++) {
       int pos_Cols = init_col + i;
       for (int j=0; j<nr; j++) {
          int pos_Rows = init_row + j; 
          (*this)(pos_Rows,pos_Cols) += M(j,i)*fact;
       }
    }
  }
  
  template<int nr, int nc> void
  assembleTranspose(const MatrixND<nr, nc, double> &M, int init_row, int init_col, double fact) 
  { 
    {
      [[maybe_unused]] int final_row = init_row + nc - 1;
      [[maybe_unused]] int final_col = init_col + nr - 1; 
      assert((init_row >= 0) && (final_row < NR) && (init_col >= 0) && (final_col < NC));
    }

    for (int i=0; i<nr; i++) {
       int pos_Cols = init_col + i;
       for (int j=0; j<nc; j++) {
          int pos_Rows = init_row + j; 
          (*this)(pos_Rows,pos_Cols) += M(i,j)*fact;
       }
    }
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

/*
  constexpr MatrixND &
  operator=(const MatrixND<NR,NC> &other)
  {  
    for (index_t j = 0; j < NC; ++j) {
      for (index_t i = 0; i < NR; ++i) {
        values[j][i] = other.values[j][i];
      }
    }
    return *this;
  }
*/

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
  operator-=(const MatrixND &other) {
    for (index_t j = 0; j < NC; ++j)
      for (index_t i = 0; i < NR; ++i)
        values[j][i] -= other.values[j][i];

    return *this;
  }

  constexpr MatrixND &
  operator*=(T const scalar) {
    for (index_t j = 0; j < NC; ++j)
      for (index_t i = 0; i < NR; ++i)
        values[j][i] *= scalar;

    return *this;
  }

  constexpr MatrixND &
  operator/=(T const scalar) {
    for (index_t j = 0; j < NC; ++j)
      for (index_t i = 0; i < NR; ++i)
        values[j][i] /= scalar;

    return *this;
  }


// Notes on operators:
// - define friend operators inside class for use as header-only library
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
  inline constexpr friend MatrixND<NR, J>
  operator^(const MatrixND<NR, NC> &left, const MatrixND<NC, J> &right) {
    MatrixND<NR, J> prod;
    for (index_t i = 0; i < NR; ++i) {
      for (index_t j = 0; j < J; ++j) {
        prod(i, j) = 0.0;
        for (index_t k = 0; k < NC; ++k) {
          prod(i, j) += left(k,i) * right(k,j);
        }
      }
    }
    return prod;
  }

  VectorND<NR>
  operator^(const VectorND<NR> &V) const
  {
      VectorND<NR> result;

      const double *dataPtr = &values[0][0];
      for (int i=0; i<NC; i++) {
        result[i] = 0.0;
        for (int j=0; j<NR; j++)
          result[i] += *dataPtr++ * V[j];
      }

      return result;
  }

#if 1
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
#endif

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

  friend constexpr MatrixND
  operator/(MatrixND mat, T scalar) {
    mat /= scalar; 
    return mat;
  }
  
  friend std::ostream &
  operator<<(std::ostream &out, MatrixND const &mat) {
    out << "{";
    for (int r=0; r<NR; r++){
      out << "{";
      for (int c=0; c<NC; c++)
        out << mat(r,c);
      out << "}, ";
    }
    return out << "}\n";
  }
}; // class MatrixND

#if 0
template<int n>
inline int
MatrixND<n,n>::invert(MatrixND<n,n>& M) const
//requires(n > 6)
{
  const Matrix A(*this);
  Matrix B(M);
  return A.Invert(B);
}
#endif

template <index_t nr, index_t nc, typename T>
template<typename F> inline void
MatrixND<nr, nc, T>::map(F func) const
{
  for (int i=0; i<nr; i++)
    for (int j = 0; j<nc; j++)
      func((*this)(i,j));
}


template <index_t nr, index_t nc, typename T>
template<typename F> inline void
MatrixND<nr, nc, T>::map(F func, MatrixND<nr,nc,T>& destination)
{
  for (int i=0; i<nr; i++)
    for (int j = 0; j<nc; j++)
      destination(i,j) = func((*this)(i,j));
}

// template<> inline int
template <index_t nr, index_t nc, typename T> inline int
MatrixND<nr,nc,T>::invert(MatrixND<nr,nc,T> &M) const
{
  int status;
  switch (nr) {
    case 2:
      cmx_inv2(&this->values[0][0], &M.values[0][0], &status);
      break;
    case 3:
      cmx_inv3(&this->values[0][0], &M.values[0][0], &status);
      break;
    case 4:
      cmx_inv4(&this->values[0][0], &M.values[0][0], &status);
      break;
    case 5:
      cmx_inv5(&this->values[0][0], &M.values[0][0], &status);
      break;
    case 6:
      cmx_inv6(&this->values[0][0], &M.values[0][0], &status);
      break;
    default:
      return -1;
  }
  return status;
}

template <index_t nr, index_t nc, typename T> inline
MatrixND<nr, nc, T>& 
MatrixND<nr, nc, T>::addDiagonal(const double diag)
{
  for (int i=0; i<nr; i++)
    (*this)(i,i) += diag;

  return *this;
}

template <index_t nr, index_t nc, typename T>
template <class MatT> inline
void MatrixND<nr, nc, T>::addMatrix(const MatT& A, const double scale)
{
  for (int i=0; i<nr; i++)
    for (int j=0; j<nc; j++)
      (*this)(i,j) += A(i,j)*scale;
}

template <index_t nr, index_t nc, typename T> 
template <class VecA, class VecB> inline
MatrixND<nr, nc, T>&
MatrixND<nr, nc, T>::addTensorProduct(const VecA& a, const VecB& b, const double scale)
{
  // Chrystal's bun order
  for (int j=0; j<nc; j++)
    for (int i=0; i<nr; i++)
      (*this)(i,j) += a[i]*b[j]*scale;
  return *this;
}

template <index_t nr, index_t nc, typename T> 
template <class MatT, int nk> inline
void
MatrixND<nr, nc, T>::addMatrixProduct(const MatrixND<nr, nk, T>& A, const MatT& B, const double scale)
{
  for (int i=0; i<nr; i++)
    for (int j=0; j<nc; j++)
      for (int k=0; k < nk; k++)
        (*this)(i,j) += A(i,k)*B(k,j)*scale;
}

// B'*C
template <index_t nr, index_t nc, typename T> 
template <class MatT, int nk> inline
void
MatrixND<nr, nc, T>::addMatrixTransposeProduct(double thisFact,
                                               const MatrixND<nk, nr, T>& B,
                                               const MatT& C,
                                               const double otherFact)
{
  if (thisFact == 1.0) {
    double *aijPtr = &values[0][0];
    for (int j=0; j<nc; j++) {
      for (int i=0; i<nr; i++) {
        const double *bkiPtr  = &(&B.values[0][0])[i*nk];
        const double *cjkPtr  = &(&C.values[0][0])[j*nk];
        double sum = 0.0;
        for (int k=0; k<nk; k++) {
          sum += *bkiPtr++ * *cjkPtr++;
        }
        *aijPtr++ += sum * otherFact;
      }
    } 
  } else if (thisFact == 0.0) {
    double *aijPtr = &values[0][0];
    for (int j=0; j<nc; j++) {
      for (int i=0; i<nr; i++) {
        const double *bkiPtr  = &(&B.values[0][0])[i*nk];
        const double *cjkPtr  = &(&C.values[0][0])[j*nk];
        double sum = 0.0;
        for (int k=0; k<nk; k++) {
          sum += *bkiPtr++ * *cjkPtr++;
        }
        *aijPtr++ = sum * otherFact;
      }
    } 
  } else {
    double *aijPtr = &values[0][0];
    for (int j=0; j<nc; j++) {
      for (int i=0; i<nr; i++) {
        const double *bkiPtr  = &(&B.values[0][0])[i*nk];
        const double *cjkPtr  = &(&C.values[0][0])[j*nk];
        double sum = 0.0;
        for (int k=0; k<nk; k++) {
          sum += *bkiPtr++ * *cjkPtr++;
        }
        *aijPtr = *aijPtr * thisFact + sum * otherFact;
        aijPtr++;
      }
    } 
  }
#if 0
  if (thisFact == 1.0) {
    for (int j=0; j<nc; j++)
      for (int i=0; i<nr; i++)
        for (int k=0; k < nk; k++)
          (*this)(i,j) += A(k,i)*B(k,j)*scale;
  } else if (thisFact == 0.0) {
    for (int j=0; j<nc; j++)
      for (int i=0; i<nr; i++) {
        double sum = 0.0;
        for (int k=0; k < nk; k++)
          sum  += A(k,i)*B(k,j);
        (*this)(i,j) = sum*scale;
      }
  } else {
    for (int j=0; j<nc; j++)
      for (int i=0; i<nr; i++)
        for (int k=0; k < nk; k++)
          (*this)(i,j)  = (*this)(i,j)*thisFact + A(k,i)*B(k,j)*scale;
  }
#endif
}

// A'BA
template <int nr, int nc, class scalar_t> 
template <int ncB> inline
int
MatrixND<nr,nc,scalar_t>::addMatrixTripleProduct( 
                               double thisFact,
                               const  MatrixND<ncB, nr, scalar_t> &T, 
                               const  MatrixND<ncB, ncB, scalar_t> &B, 
                               double otherFact)
  requires(nr == nc)
{
  if (otherFact == 0.0 && thisFact == 1.0)
    return 0;

  MatrixND<ncB, nc> BT;
  BT.zero();
  BT.addMatrixProduct(B, T, otherFact);
  this->addMatrixTransposeProduct(thisFact, T, BT, 1.0);

//{
//  int m   = B.numRows,
//      n   = T.numCols,
//      k   = B.numCols,
//      nrT = T.numRows;
//    //k = T.numRows;
//  double zero = 0.0,
//         one  = 1.0;

//  DGEMM ("N", "N", &m      , &n      , &k,&one      , B.data, &m, // m
//                                                      T.data, &nrT, // k
//                                          &zero,  matrixWork, &m);

//  DGEMM ("T", "N", &numRows, &numCols, &k,&otherFact, T.data, &nrT,
//                                                  matrixWork, &m, // k
//                                          &thisFact,    data, &numRows);
//  return 0;
//}

  return 0;
}


template <int nr, int nc, class scalar_t> 
template <int nk, int nl> inline int 
MatrixND<nr,nc,scalar_t>::addMatrixTripleProduct(double thisFact, 
                           const MatrixND<nk,nr> &A, 
                           const MatrixND<nk,nl> &B, 
                           const MatrixND<nl,nc> &C, double otherFact)
{

  if (otherFact == 0.0 && thisFact == 1.0)
    return 0;

  MatrixND<nk, nc> BC {0.0};
  BC.addMatrixProduct(B, C, otherFact);
  this->addMatrixTransposeProduct(thisFact, A, BC, 1.0);
  return 0;
}


template <int nr, int nc, typename T=double>
MatrixND(const T (&)[nc][nr])->MatrixND<nr, nc, T>;


//
// 3D
//

template<int NR, int NC, typename T>
template<class VecT>
inline MatrixND<NR,NC,T>& 
MatrixND<NR,NC,T>::addSpin(const VecT& v)  requires(NR == 3)
{
   const double v0 = v[0],
                v1 = v[1],
                v2 = v[2];

  (*this)(0, 0) += 0.0;   (*this)(0, 1) += -v2;     (*this)(0, 2) +=  v1;
  (*this)(1, 0) +=  v2;   (*this)(1, 1) +=  0.0;    (*this)(1, 2) += -v0;
  (*this)(2, 0) += -v1;   (*this)(2, 1) +=  v0;     (*this)(2, 2) += 0.0;

  return *this;
}


template<int NR, int NC, typename T>
template<class VecT>
inline MatrixND<NR,NC,T>&
MatrixND<NR,NC,T>::addSpin(const VecT& v, double mult)
{
   const double v0 = mult*v[0],
                v1 = mult*v[1],
                v2 = mult*v[2];

  (*this)(0, 0) += 0.0;   (*this)(0, 1) += -v2;     (*this)(0, 2) +=  v1;
  (*this)(1, 0) +=  v2;   (*this)(1, 1) += 0.00;    (*this)(1, 2) += -v0;
  (*this)(2, 0) += -v1;   (*this)(2, 1) +=  v0;     (*this)(2, 2) += 0.0;
  return *this;
}


template<int NR, int NC, typename T>
template <class VecT> inline
MatrixND<NR,NC,T>& 
MatrixND<NR,NC,T>::addSpinSquare(const VecT& v, const double scale)
  requires(NR == NC == 3)
{
  const double v1 = v[0],
               v2 = v[1],
               v3 = v[2];

  (*this)(0,0) += scale*( -v2*v2 - v3*v3 );
  (*this)(1,1) += scale*( -v1*v1 - v3*v3 );
  (*this)(2,2) += scale*( -v1*v1 - v2*v2 );

  (*this)(0,1) += scale*(  v1*v2 );
  (*this)(1,0) += scale*(  v1*v2 );
  (*this)(2,0) += scale*(  v1*v3 );
  (*this)(0,2) += scale*(  v1*v3 );
  (*this)(1,2) += scale*(  v2*v3 );
  (*this)(2,1) += scale*(  v2*v3 );
  return *this;
}


template<int NR, int NC, typename T>
template<class VecT> inline
void MatrixND<NR,NC,T>::addSpinProduct(const VecT& a, const VectorND<NR,T>& b, const double scale)
  requires(NR == NC == 3)
{
  // a^b^ = boa - a.b 1
  // where 'o' denotes the tensor product and '.' the dot product
  //
  this->addTensorProduct(b, a, scale);
  this->addDiagonal(-b.dot(a)*scale);
}

template<int NR, int NC, typename T>
template<class VecT> inline
void MatrixND<NR,NC,T>::addMatrixSpinProduct(const MatrixND<NR,NC,T>& A, const VecT& b, const double scale)
{
  // this += s*A*[b^]
  // where b^ is the skew-symmetric representation of the three-vector b, s is a scalar,
  // and A a 3x3 matrix.
  //
  (*this)(0, 0) += scale*( A(0,1)*b[2] - A(0,2)*b[1]);
  (*this)(0, 1) += scale*(-A(0,0)*b[2] + A(0,2)*b[0]);
  (*this)(0, 2) += scale*( A(0,0)*b[1] - A(0,1)*b[0]);
  (*this)(1, 0) += scale*( A(1,1)*b[2] - A(1,2)*b[1]);
  (*this)(1, 1) += scale*(-A(1,0)*b[2] + A(1,2)*b[0]);
  (*this)(1, 2) += scale*( A(1,0)*b[1] - A(1,1)*b[0]);
  (*this)(2, 0) += scale*( A(2,1)*b[2] - A(2,2)*b[1]);
  (*this)(2, 1) += scale*(-A(2,0)*b[2] + A(2,2)*b[0]);
  (*this)(2, 2) += scale*( A(2,0)*b[1] - A(2,1)*b[0]);
}

template<int NR, int NC, typename T>
template<class MatT> inline
void MatrixND<NR,NC,T>::addSpinMatrixProduct(const VectorND<NR,T>& a, const MatT& B, const double scale)
  requires(NR == NC == 3)
{
  // this += s*[a^]*B
  // where a^ is the skew-symmetric representation of the three-vector a, s is a scalar,
  // and B a 3x3 matrix.
  //
  (*this)(0, 0) += scale*( -B(1,0)*a[2] + B(2,0)*a[1]);
  (*this)(0, 1) += scale*( -B(1,1)*a[2] + B(2,1)*a[1]);
  (*this)(0, 2) += scale*( -B(1,2)*a[2] + B(2,2)*a[1]);
  (*this)(1, 0) += scale*(  B(0,0)*a[2] - B(2,0)*a[0]);
  (*this)(1, 1) += scale*(  B(0,1)*a[2] - B(2,1)*a[0]);
  (*this)(1, 2) += scale*(  B(0,2)*a[2] - B(2,2)*a[0]);
  (*this)(2, 0) += scale*( -B(0,0)*a[1] + B(1,0)*a[0]);
  (*this)(2, 1) += scale*( -B(0,1)*a[1] + B(1,1)*a[0]);
  (*this)(2, 2) += scale*( -B(0,2)*a[1] + B(1,2)*a[0]);
}

} // namespace OpenSees

#endif // MatrixND_H

