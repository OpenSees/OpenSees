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
#pragma once
#include "MatrixND.h"

#include "routines/xblas.h"
#include "routines/cmx.h"

namespace OpenSees {

template <int nr, int nc, typename T=double>
MatrixND(const T (&)[nc][nr])->MatrixND<nr, nc, T>;


template <index_t nr, index_t nc, typename T>
constexpr void
MatrixND<nr, nc, T>::zero() noexcept
{
  values.fill(T{});
}


template <index_t nr, index_t nc, typename T>
constexpr MatrixND<nc, nr>
MatrixND<nr, nc, T>::transpose() const noexcept
{
  MatrixND<nc, nr> result = {};
  for (index_t j = 0; j < nc; ++j) {
    for (index_t i = 0; i < nr; ++i) {
      result(j,i) = (*this)(i,j);
    }
  }
  return result;
}


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

//
//
//
template <index_t NR, index_t NC, typename T>
template <int nr> inline void
MatrixND<NR,NC,T>::assemble(const VectorND<nr> &v, int init_row, int init_col, double fact) noexcept
{

  [[maybe_unused]] int final_row = init_row + nr - 1;
  assert((init_row >= 0) && (final_row < NR));

  for (int j=0; j<nr; j++)
    (*this)(init_row+j, init_col) += v[j]*fact;
}


template <index_t NR, index_t NC, typename T>
template<int nr, int nc> void
MatrixND<NR,NC,T>::assembleTranspose(const MatrixND<nr, nc, double> &M, int init_row, int init_col, double fact) noexcept
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

template <index_t NR, index_t NC, typename T>
template<int nr> void
MatrixND<NR,NC,T>::assembleTranspose(const VectorND<nr> &v, int init_row, int init_col, double scale) noexcept
{ 
  {
    [[maybe_unused]] int final_col = init_col + nr - 1; 
    assert((init_row >= 0) && (init_col >= 0) && (final_col < NC));
  }

  for (int i=0; i<nr; i++)
    (*this)(init_row, init_col+i) += v[i]*scale;
}

//
// Solve
//
template <index_t nr, index_t nc, typename T> inline int
MatrixND<nr,nc,T>::invert(MatrixND<nr,nc,T> &M) const
{
  static_assert(nr == nc, "Matrix must be square");

  int status = -1;
  if constexpr (nr == 2) {
    cmx_inv2(&(*this)(0,0), &M(0,0), &status);
    return status;
  }
  if constexpr (nr == 3) {
    cmx_inv3(&(*this)(0,0), &M(0,0), &status);
    return status;
  }
  if constexpr (nr == 4) {
    cmx_inv4(&(*this)(0,0), &M(0,0), &status);
    return status;
  }
  if constexpr (nr == 5) {
    cmx_inv5(&(*this)(0,0), &M(0,0), &status);
    return status;
  }
  if constexpr (nr == 6) {
    cmx_inv6(&(*this)(0,0), &M(0,0), &status);
    return status;
  }

  // Use Lapack
  M.zero();
  M.addDiagonal(1.0); // Identity matrix 

  MatrixND<nr, nc, T> work = *this;

  int pivot_ind[nr];
  int nrhs = nr;
  int info = 0;
  int m = nr;
  int n = nc;
  DGESV(&m, 
        &nrhs,
        &work(0,0), &m,
        pivot_ind, 
        &M(0,0), 
        &n, 
        &info);

  if (info != 0)
    status = -std::abs(info);
  return status;
}

template <index_t nr, index_t nc, typename T>
constexpr MatrixND<nr, nc, T>& 
MatrixND<nr, nc, T>::addDiagonal(const double diag) noexcept
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
template <class VecA, class VecB> 
constexpr inline
MatrixND<nr, nc, T>&
MatrixND<nr, nc, T>::addTensorProduct(const VecA& a, const VecB& b, const double scale) noexcept
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
MatrixND<nr, nc, T>::addMatrixProduct(const MatrixND<nr, nk, T>& A, const MatT& B, double scale)
{
  if constexpr (nr*nc < 48)
    for (int i=0; i<nr; i++)
      for (int j=0; j<nc; j++)
        for (int k=0; k < nk; k++)
          (*this)(i,j) += scale*A(i,k)*B(k,j);
  else
  {
    int m = nr,
        n = nc,
        k = nk;
    double one = 1.0;
    DGEMM("N", "N", &m, &n, &k, &scale, 
                                const_cast<double*>(&A(0,0)), &m,
                                const_cast<double*>(&B(0,0)), &k,
                                &one,   &(*this)(0,0),        &m);
  }
}

#if 1
template <index_t nr, index_t nc, typename T> 
template <class MatT, int nk> inline
void
MatrixND<nr, nc, T>::addMatrixProduct(double scale_this, 
                                      const MatrixND<nr, nk, T>& A, 
                                      const MatT& B, double scale)
{
  int m = nr,
      n = nc,
      k = nk;
  DGEMM("N", "N", &m, &n, &k, &scale, 
                              const_cast<double*>(&A(0,0)), &m,
                              const_cast<double*>(&B(0,0)), &k,
                              &scale_this,   &(*this)(0,0), &m);
}
#endif

// B'*C
template <index_t nr, index_t nc, typename T> 
template <class MatT, int nk> inline
void
MatrixND<nr, nc, T>::addMatrixTransposeProduct(double thisFact,
                                               const MatrixND<nk, nr, T>& B,
                                               const MatT& C,
                                               double otherFact)
{
  if constexpr (nr*nc > 16) {
    int m = nr,
        n = nc,
        k = nk;
    DGEMM("T", "N", &m, &n, &k, &otherFact, 
                                const_cast<double*>(&B(0,0)), &k,
                                const_cast<double*>(&C(0,0)), &k,
                                &thisFact,   &(*this)(0,0),   &m);
    return;
  }

  if (thisFact == 1.0) {
    double *aijPtr =this->data();
    for (int j=0; j<nc; j++) {
      for (int i=0; i<nr; i++) {
        const double *bkiPtr  = &(&B(0,0))[i*nk];
        const double *cjkPtr  = &(&C(0,0))[j*nk];
        double sum = 0.0;
        for (int k=0; k<nk; k++) {
          sum += *bkiPtr++ * *cjkPtr++;
        }
        *aijPtr++ += sum * otherFact;
      }
    } 
  } else if (thisFact == 0.0) {
    double *aijPtr =this->data();
    for (int j=0; j<nc; j++) {
      for (int i=0; i<nr; i++) {
        const double *bkiPtr  = &(&B(0,0))[i*nk];
        const double *cjkPtr  = &(&C(0,0))[j*nk];
        double sum = 0.0;
        for (int k=0; k<nk; k++) {
          sum += *bkiPtr++ * *cjkPtr++;
        }
        *aijPtr++ = sum * otherFact;
      }
    } 
  } else {
    double *aijPtr =this->data();
    for (int j=0; j<nc; j++) {
      for (int i=0; i<nr; i++) {
        const double *bkiPtr  = &(&B(0,0))[i*nk];
        const double *cjkPtr  = &(&C(0,0))[j*nk];
        double sum = 0.0;
        for (int k=0; k<nk; k++) {
          sum += *bkiPtr++ * *cjkPtr++;
        }
        *aijPtr = *aijPtr * thisFact + sum * otherFact;
        aijPtr++;
      }
    } 
  }
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
{
  if (otherFact == 0.0 && thisFact == 1.0)
    return 0;

  MatrixND<ncB, nc> BT;
  BT.zero();
  BT.addMatrixProduct(B, T, otherFact);
  this->addMatrixTransposeProduct(thisFact, T, BT, 1.0);
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

  MatrixND<nk, nc> BC {};
  BC.addMatrixProduct(B, C, otherFact);
  this->addMatrixTransposeProduct(thisFact, A, BC, 1.0);
  return 0;
}


//
// 3D
//

template<int NR, int NC, typename T>
template<class VecT>
constexpr MatrixND<NR,NC,T>& 
MatrixND<NR,NC,T>::addSpin(const VecT& v) noexcept
{
  static_assert(NR == 3 && NC == 3, "addSpin requires a 3x3 matrix");

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
constexpr MatrixND<NR,NC,T>&
MatrixND<NR,NC,T>::addSpin(const VecT& v, double mult) noexcept
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
template <class VecT> 
constexpr
MatrixND<NR,NC,T>& 
MatrixND<NR,NC,T>::addSpinSquare(const VecT& v, const double scale) noexcept
{
  static_assert(NR == 3 && NC == 3, "addSpinSquare requires a 3x3 matrix");
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
template<class VecT> 
constexpr void 
MatrixND<NR,NC,T>::addSpinProduct(const VecT& a, const VectorND<NR,T>& b, const double scale) noexcept
{
  // a^b^ = boa - a.b 1
  // where 'o' denotes the tensor product and '.' the dot product
  //
  static_assert(NR == 3 && NC == 3, "Matrix must be 3x3");
  this->addTensorProduct(b, a, scale);
  this->addDiagonal(-b.dot(a)*scale);
}

template<int NR, int NC, typename T>
template<class VecT>
constexpr void 
MatrixND<NR,NC,T>::addMatrixSpinProduct(const MatrixND<NR,NC,T>& A, const VecT& b, const double scale) noexcept
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
template<class MatT>
constexpr void 
MatrixND<NR,NC,T>::addSpinMatrixProduct(const VectorND<NR,T>& a, const MatT& B, const double scale) noexcept
{
  static_assert(NR == 3 && NC == 3, "Matrix must be 3x3");
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
