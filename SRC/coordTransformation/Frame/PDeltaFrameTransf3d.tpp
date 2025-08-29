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
// Please cite the following resource in any derivative works:
//
// Perez, C.M., and Filippou F.C.. "On Nonlinear Geometric Transformations
//     of Finite Elements" Int. J. Numer. Meth. Engrg. 2024; 
//     https://doi.org/10.1002/nme.7506
//
//===----------------------------------------------------------------------===//
//
//        OpenSees - Open System for Earthquake Engineering Simulation
//
//===----------------------------------------------------------------------===//

//
// Adapted: Remo Magalhaes de Souza
//          04/2000
//
// Written: cmp
//
#include <Vector.h>
#include <Matrix.h>
#include <Matrix3D.h>
#include <Node.h>

#include "PDeltaFrameTransf3d.h"

namespace OpenSees {

template <int nn, int ndf>
PDeltaFrameTransf<nn,ndf>::PDeltaFrameTransf(int tag, 
                            const Vector3D &vecxz,
                            const std::array<Vector3D, nn> *offset,
                            int offset_flags)
  : FrameTransform<nn,ndf>(tag),
    offset_flags(offset_flags),
    linear(tag, vecxz, offset, offset_flags)
{

}


template <int nn, int ndf>
PDeltaFrameTransf<nn,ndf>::~PDeltaFrameTransf()
{

}

template <int nn, int ndf>
int
PDeltaFrameTransf<nn,ndf>::commit()
{
  return linear.commit();
}

template <int nn, int ndf>
int
PDeltaFrameTransf<nn,ndf>::revertToLastCommit()
{
  return linear.revertToLastCommit();
}

template <int nn, int ndf>
int
PDeltaFrameTransf<nn,ndf>::revertToStart()
{
  return linear.revertToStart();
}

template <int nn, int ndf>
int
PDeltaFrameTransf<nn,ndf>::initialize(std::array<Node*, nn>& new_nodes)
{
  return linear.initialize(new_nodes);
}

template <int nn, int ndf>
int
PDeltaFrameTransf<nn,ndf>::update()
{
  return linear.update();
}

template <int nn, int ndf>
int
PDeltaFrameTransf<nn,ndf>::getLocalAxes(Vector3D &XAxis, Vector3D &YAxis, Vector3D &ZAxis) const
{
  return linear.getLocalAxes(XAxis, YAxis, ZAxis);
}

template <int nn, int ndf>
double
PDeltaFrameTransf<nn,ndf>::getInitialLength()
{
  return linear.getInitialLength();
}

template <int nn, int ndf>
double
PDeltaFrameTransf<nn,ndf>::getDeformedLength()
{
  return linear.getDeformedLength();
}


template <int nn, int ndf>
// VectorND<nn*ndf>
int
PDeltaFrameTransf<nn,ndf>::push(VectorND<nn*ndf>&pl, Operation op)
{
  //
  // Include leaning column effects (P-Delta)
  //
  // Axial force
  const double N = pl[1*ndf+0];

  const Vector3D Du = linear.getDelta()/linear.getInitialLength();

  pl[0*ndf+1] -= Du[1] * N;
  pl[1*ndf+1] += Du[1] * N;

  pl[0*ndf+2] -= Du[2] * N;
  pl[1*ndf+2] += Du[2] * N;

  return linear.push(pl, op);
}


template <int nn, int ndf>
// MatrixND<nn*ndf,nn*ndf>
int
PDeltaFrameTransf<nn,ndf>::push(MatrixND<nn*ndf,nn*ndf>& kl, const VectorND<nn*ndf> &pl, Operation op)
{
  // Include geometric stiffness effects in local system;
  //
  // Kl += [ ]
  double NoverL = pl[6] / linear.getInitialLength();
  kl(1, 1) += NoverL;
  kl(2, 2) += NoverL;
  kl(7, 7) += NoverL;
  kl(8, 8) += NoverL;

  kl(1, 7) -= NoverL;
  kl(7, 1) -= NoverL;
  kl(2, 8) -= NoverL;
  kl(8, 2) -= NoverL;
  return linear.push(kl, pl, op);
}


template <int nn, int ndf>
VectorND<nn*ndf>
PDeltaFrameTransf<nn,ndf>::getStateVariation()
{
  return linear.getStateVariation();
}

template <int nn, int ndf>
Vector3D
PDeltaFrameTransf<nn,ndf>::getNodePosition(int tag)
{
  return linear.getNodePosition(tag);
}

template <int nn, int ndf>
Vector3D
PDeltaFrameTransf<nn,ndf>::getNodeRotationLogarithm(int tag)
{
  return linear.getNodeRotationLogarithm(tag);
}

template <int nn, int ndf>
FrameTransform<nn,ndf> *
PDeltaFrameTransf<nn,ndf>::getCopy() const
{
  Vector3D e1, e2, e3;
  linear.getLocalAxes(e1, e2, e3);

  return new PDeltaFrameTransf(this->getTag(), 
                               e3, 
                               linear.getRigidOffsets(),
                               offset_flags);
}


template <int nn, int ndf>
double
PDeltaFrameTransf<nn,ndf>::getLengthGrad()
{
  return linear.getLengthGrad();
}

template <int nn, int ndf>
void
PDeltaFrameTransf<nn,ndf>::Print(OPS_Stream &s, int flag)
{
  linear.Print(s, flag);
}

} // namespace OpenSees