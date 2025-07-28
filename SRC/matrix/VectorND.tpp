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
/* ****************************************************************** **
**    OpenSees - Open System for Earthquake Engineering Simulation    **
**          Pacific Earthquake Engineering Research Center            **
**                                                                    **
**                                                                    **
** (C) Copyright 1999, The Regents of the University of California    **
** All Rights Reserved.                                               **
**                                                                    **
** Commercial use of this program without express permission of the   **
** University of California, Berkeley, is strictly prohibited.  See   **
** file 'COPYRIGHT'  in main directory for information on usage and   **
** redistribution,  and for a DISCLAIMER OF ALL WARRANTIES.           **
**                                                                    **
** Developed by:                                                      **
**   Frank McKenna (fmckenna@ce.berkeley.edu)                         **
**   Gregory L. Fenves (fenves@ce.berkeley.edu)                       **
**   Filip C. Filippou (filippou@ce.berkeley.edu)                     **
**                                                                    **
** ****************************************************************** */
//
#include <assert.h>
#include "VectorND.h"
#include "routines/xblas.h"


namespace OpenSees {

template <index_t N, typename T>
template <int ir, int nr> inline void
VectorND<N,T>::assemble(const VectorND<nr> &v, double fact)
{
  static_assert((ir >= 0) && ((ir + nr - 1) < N));

  for (int j=0; j<nr; j++)
    (*this)(ir + j) += v[j]*fact;
}

template <index_t N, typename T>
template<int nr> void
VectorND<N,T>::assemble(int a, const VectorND<nr>& v, double scale) noexcept
{
  for (int i=0; i<nr; i++)
    (*this)[a+i] += v[i]*scale;
}

template <index_t N, typename T>
template <int ir, int nr> inline void
VectorND<N,T>::insert(const VectorND<nr> &v, double fact)
{
  static_assert((ir >= 0) && ((ir + nr - 1) < N));

  for (int j=0; j<nr; j++)
    (*this)(ir + j) = v[j]*fact;
}


template <index_t N, typename T>
template<int nr> void
VectorND<N,T>::insert(int a, const VectorND<nr>& v, double scale) noexcept
{
  for (int i=0; i<nr; i++)
    (*this)[a+i]  = v[i]*scale;
}


template <index_t N, typename T>
template<int nr>
VectorND<nr> 
VectorND<N,T>::extract(int a) noexcept
{
  VectorND<nr> v;
  for (int i=0; i<nr; i++)
    v[i] = (*this)[a+i];
  return v;
}

#ifdef XARA_VECTOR_FRIENDS
template <index_t N, typename T>
int
VectorND<N,T>::addVector(const T thisFact, const Vector &other, const T otherFact) noexcept
{
  if (otherFact == 0.0 && thisFact == 1.0)
    return 0; 

  else if (thisFact == 1.0) {
    // want: this += other * otherFact
    double *dataPtr = values;
    double *otherDataPtr = other.theData;
    if (otherFact == 1.0) { 
      // no point doing a multiplication if otherFact == 1.0
      for (int i=0; i<N; i++) 
        *dataPtr++ += *otherDataPtr++;
    } else if (otherFact == -1.0) { 
      // no point doing a multiplication if otherFact == 1.0
      for (int i=0; i<N; i++) 
        *dataPtr++ -= *otherDataPtr++;
    } else 
      for (int i=0; i<N; i++) 
        *dataPtr++ += *otherDataPtr++ * otherFact;

  } else if (thisFact == 0.0) {
      // want: this = other * otherFact
      double *dataPtr = values;
      double *otherDataPtr = other.theData;
      if (otherFact == 1.0) {
        for (int i=0; i<N; i++) 
          *dataPtr++ = *otherDataPtr++;
      } else if (otherFact == -1.0) {
        for (int i=0; i<N; i++) 
          *dataPtr++ = -(*otherDataPtr++);
      } else 
        for (int i=0; i<N; i++) 
          *dataPtr++ = *otherDataPtr++ * otherFact;
  } else {
    // want: this = this * thisFact + other * otherFact
    double *dataPtr = values;
    double *otherDataPtr = other.theData;
    if (otherFact == 1.0) {
      for (int i=0; i<N; i++) {
        double value = *dataPtr * thisFact + *otherDataPtr++;
        *dataPtr++ = value;
      }
    } else if (otherFact == -1.0) {
      for (int i=0; i<N; i++) {
        double value = *dataPtr * thisFact - *otherDataPtr++;
        *dataPtr++ = value;
      }
    } else 
      for (int i=0; i<N; i++) {
        double value = *dataPtr * thisFact + *otherDataPtr++ * otherFact;
        *dataPtr++ = value;
    }
  }

  return 0;
}


template <index_t N, typename T>
template <int NC>
inline int
VectorND<N,T>::addMatrixVector(double thisFact, const MatrixND<N, NC, double> &m, const Vector& v, double otherFact)
{
  // check the sizes are compatable
  assert(NC == v.sz);

  // see if quick return
  if (thisFact == 1.0 && otherFact == 0.0)
    return 0;

  else {
    int incr = 1,
            i = N,
            n = NC;
      DGEMV("N", &i, &n,
            &otherFact,
            &m(0,0), &i,
            v.theData, &incr,
            &thisFact,
            values,
            &incr);

    return 0;
  } 
}


template <index_t N, typename T>
template <int NR>
inline int
VectorND<N,T>::addMatrixTransposeVector(double thisFact, const MatrixND<NR, N, double> &m, const Vector &v, double otherFact)
{
  // check the sizes are compatable
  assert(NR == v.sz);


  if (thisFact == 1.0 && otherFact == 0.0)
    return 0;

  else {
    int incr = 1,
            i = NR,
            n = N;
    DGEMV("T", &i, &n,
          &otherFact,
          &m(0,0), &i,
          v.theData, &incr,
          &thisFact,
          values,   &incr);
    return 0;
  } 
}


template <index_t N, typename T>
inline int
VectorND<N,T>::addMatrixVector(const double thisFact, const Matrix &m, const Vector &v, const double otherFact)
{
  // check the sizes are compatable
  assert(N == m.noRows());
  assert(m.noCols() == v.sz);

  // see if quick return
  if (thisFact == 1.0 && otherFact == 0.0)
    return 0;

#ifdef VECTOR_BLAS
  else if (v.sz > 10) {
    int incr = 1,
            i = m.numRows,
            n = m.numCols;
    return
      DGEMV("N", &i, &n,
            &otherFact,
            m.data, &i,
            v.theData, &incr,
            &thisFact,
            values,   &incr);
  }
#endif

  else if (thisFact == 1.0) {

    // want: this += m * v * otherFact
    if (otherFact == 1.0) { // no point doing multiplication if otherFact = 1.0
      int otherSize = v.sz;
      double *matrixDataPtr = m.data;
      const double *otherDataPtr = v.theData;
      for (int i=0; i<otherSize; i++) {
        double otherData = *otherDataPtr++;
        for (int j=0; j < N; j++)
          values[j] += *matrixDataPtr++ * otherData;
      }
    } 
    else if (otherFact == -1.0) { // no point doing multiplication if otherFact = -1.0
      int otherSize = v.sz;
      double *matrixDataPtr = m.data;
      const double *otherDataPtr = v.theData;
      for (int i=0; i<otherSize; i++) {
        double otherData = *otherDataPtr++;
        for (int j=0; j < N; j++)
          values[j] -= *matrixDataPtr++ * otherData;
      }
    } 
    else { // have to do the multiplication
      int otherSize = v.sz;
      double *matrixDataPtr = m.data;
      double *otherDataPtr = v.theData;
      for (int i=0; i<otherSize; i++) {
        double otherData = *otherDataPtr++ * otherFact;
        for (int j=0; j < N; j++)
          values[j] += *matrixDataPtr++ * otherData;
      }
    }
  }

  else if (thisFact == 0.0) {
    
    // want: this = m * v * otherFact
    for (int i=0; i < N; i++)
      values[i] = 0.0;

    if (otherFact == 1.0) { 
      // avoid multiplication when otherFact = 1.0
      int otherSize = v.sz;
      const double *matrixDataPtr = m.data;
      const double *otherDataPtr = v.theData;
      for (int i=0; i<otherSize; i++) {
        const double otherData = *otherDataPtr++;
        for (int j=0; j < N; j++)
          values[j] += *matrixDataPtr++ * otherData;
      }
    } 

    else if (otherFact == -1.0) { 
      // avoid multiplication when otherFact = -1.0
      int otherSize = v.sz;
      double *matrixDataPtr = m.data;
      const double *otherDataPtr = v.theData;
      for (int i=0; i<otherSize; i++) {
        const double otherData = *otherDataPtr++;
        for (int j=0; j < N; j++)
          values[j] -= *matrixDataPtr++ * otherData;
      }

    } else {
      int otherSize = v.sz;
      double *matrixDataPtr = m.data;
      double *otherDataPtr = v.theData;
      for (int i=0; i<otherSize; i++) {
        double otherData = *otherDataPtr++ * otherFact;
        for (int j=0; j < N; j++)
          values[j] += *matrixDataPtr++ * otherData;
      }
    }
  }

  else {

    // want: this = this * thisFact + m * v * otherFact
    for (int i=0; i<N; i++)
      values[i] *= thisFact;

    if (otherFact == 1.0) { // no point doing multiplication if otherFact = 1.0
      int otherSize = v.sz;
      double *matrixDataPtr = m.data;
      double *otherDataPtr = v.theData;
      for (int i=0; i<otherSize; i++) {
        double otherData = *otherDataPtr++;
        for (int j=0; j < N; j++)
          values[j] += *matrixDataPtr++ * otherData;
      }
    } else if (otherFact == -1.0) { // no point doing multiplication if otherFact = 1.0
      int otherSize = v.sz;
      double *matrixDataPtr = m.data;
      double *otherDataPtr = v.theData;
      for (int i=0; i<otherSize; i++) {
        double otherData = *otherDataPtr++;
        for (int j=0; j < N; j++)
          values[j] -= *matrixDataPtr++ * otherData;
      }
    } else {
      int otherSize = v.sz;
      double *matrixDataPtr = m.data;
      double *otherDataPtr = v.theData;
      for (int i=0; i<otherSize; i++) {
        double otherData = *otherDataPtr++ * otherFact;
        for (int j=0; j < N; j++)
          values[j] += *matrixDataPtr++ * otherData;
      }
    }
  }

  // successfull
  return 0;
}

#endif // XARA_VECTOR_FRIENDS

template <index_t N, typename T>
int
VectorND<N,T>::addVector(const T thisFact, const VectorND<N> &other, const T otherFact) noexcept
{
  if (otherFact == 0.0 && thisFact == 1.0)
    return 0; 

  else if (thisFact == 1.0) {
    // want: this += other * otherFact
    double *dataPtr = values;
    const double * otherDataPtr = other.values;
    if (otherFact == 1.0) { // no point doing a multiplication if otherFact == 1.0
      for (int i=0; i<N; i++) 
        *dataPtr++ += *otherDataPtr++;
    } else if (otherFact == -1.0) { // no point doing a multiplication if otherFact == 1.0
      for (int i=0; i<N; i++) 
        *dataPtr++ -= *otherDataPtr++;
    } else 
      for (int i=0; i<N; i++) 
        *dataPtr++ += *otherDataPtr++ * otherFact;

  } else if (thisFact == 0.0) {
      // want: this = other * otherFact
      double *dataPtr = values;
      const double *otherDataPtr = other.values;
      if (otherFact == 1.0) {
        // no point doing a multiplication if otherFact == 1.0
        for (int i=0; i<N; i++) 
          *dataPtr++ = *otherDataPtr++;
      } else if (otherFact == -1.0) {
        // no point doing a multiplication if otherFact == 1.0
        for (int i=0; i<N; i++) 
          *dataPtr++ = -(*otherDataPtr++);
      } else 
        for (int i=0; i<N; i++) 
          *dataPtr++ = *otherDataPtr++ * otherFact;
  } else {
    // want: this = this * thisFact + other * otherFact
    double *dataPtr = values;
    const double *otherDataPtr = other.values;
    if (otherFact == 1.0) { // no point doing a multiplication if otherFact == 1.0
      for (int i=0; i<N; i++) {
        double value = *dataPtr * thisFact + *otherDataPtr++;
        *dataPtr++ = value;
      }
    } else if (otherFact == -1.0) { // no point doing a multiplication if otherFact ==-1.0
      for (int i=0; i<N; i++) {
        double value = *dataPtr * thisFact - *otherDataPtr++;
        *dataPtr++ = value;
      }
    } else 
      for (int i=0; i<N; i++) {
        double value = *dataPtr * thisFact + *otherDataPtr++ * otherFact;
        *dataPtr++ = value;
    }
  }

  return 0;
}


} //  namespace OpenSees
