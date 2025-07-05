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

} //  namespace OpenSees