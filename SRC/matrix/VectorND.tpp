//===----------------------------------------------------------------------===//
//
//                                   xara
//
//===----------------------------------------------------------------------===//
//                              https://xara.so
//===----------------------------------------------------------------------===//
//
#include <assert.h>
#include "VectorND.h"

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
    if (otherFact == 1.0)
      for (int i=0; i<N; i++) 
        *dataPtr++ += *otherDataPtr++;
    else if (otherFact == -1.0)
      for (int i=0; i<N; i++) 
        *dataPtr++ -= *otherDataPtr++;
    else 
      for (int i=0; i<N; i++) 
        *dataPtr++ += *otherDataPtr++ * otherFact;

  } else if (thisFact == 0.0) {
      // want: this = other * otherFact
      double *dataPtr = values;
      const double *otherDataPtr = other.values;
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
    const double *otherDataPtr = other.values;
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

} //  namespace OpenSees