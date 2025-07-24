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
// Written: Claudio M. Perez
// Created: 04/2025
//

#pragma once
#include <Vector.h>
#include <Matrix.h>
#include <Matrix3D.h>
#include <Node.h>
#include <OPS_Globals.h>
#include <AxisAngle.h>
#include <GroupSO3.h>
#include "EuclidFrameTransf.h"

namespace OpenSees {

template <int nn, int ndf, typename IsoT>
EuclidFrameTransf<nn,ndf,IsoT>::EuclidFrameTransf(int tag, 
                                           const Vector3D &vecxz, 
                                           const std::array<Vector3D, nn> *offset,
                                           int offset_flags)
  : FrameTransform<nn,ndf>(tag),
    L(0),
    nodes{},
    ur{},
    offsets{nullptr},
    offset_flags(offset_flags),
    basis{vecxz}
{
  R0.zero();
  R0.addDiagonal(1.0);
  double nz = vecxz.norm();
  for (int i=0; i<3; i++)
    vz[i] = vecxz[i]/nz;

  // Rigid joint offsets
  if (offset != nullptr) {
    offsets = new std::array<Vector3D, nn>{};
    *offsets = *offset;
    basis.setOffsets(offsets);
  }
}



template <int nn, int ndf, typename IsoT>
EuclidFrameTransf<nn,ndf,IsoT>::~EuclidFrameTransf()
{
  if (offsets != nullptr)
    delete offsets;
}

template <int nn, int ndf, typename IsoT>
int
EuclidFrameTransf<nn,ndf,IsoT>::commit()
{
  return 0;
}

template <int nn, int ndf, typename IsoT>
int
EuclidFrameTransf<nn,ndf,IsoT>::revertToLastCommit()
{
  return 0;
}

template <int nn, int ndf, typename IsoT>
int
EuclidFrameTransf<nn,ndf,IsoT>::revertToStart()
{
  return 0;
}


template <int nn, int ndf, typename IsoT>
int
EuclidFrameTransf<nn,ndf,IsoT>::initialize(std::array<Node*, nn>& new_nodes)
{
  for (int i=0; i<nn; i++) {
    nodes[i] = new_nodes[i];
    if (nodes[i] == nullptr) {
      opserr << "invalid pointers to element nodes\n";
      return -1;
    }
    // ensure the node is initialized
    nodes[i]->getTrialRotation();
  }
  
  const Vector &XI = nodes[   0]->getCrds();
  const Vector &XJ = nodes[nn-1]->getCrds();

  for (int i=0; i<3; i++) {
    xi[i] = XI[i];
    xj[i] = XJ[i];
  }
  
  Vector3D dx = xj - xi;

  if (offsets != nullptr) {
    for (int i=0; i<3; i++)
      dx(i) -= (*offsets)[   0][i];
    for (int i=0; i<3; i++)
      dx(i) += (*offsets)[nn-1][i];
  }

  // calculate the element length
  L = dx.norm();

  if (L == 0.0)
    return -2;

  int error = basis.initialize(nodes);

  R0 = basis.getRotation();
  return error;
}


template <int nn, int ndf, typename IsoT>
FrameTransform<nn,ndf> *
EuclidFrameTransf<nn,ndf,IsoT>::getCopy() const
{
  return new EuclidFrameTransf<nn,ndf,IsoT>(this->getTag(), vz, offsets);
}


template <int nn, int ndf, typename IsoT>
int
EuclidFrameTransf<nn,ndf,IsoT>::getLocalAxes(Vector3D &e1, Vector3D &e2, Vector3D &e3) const
{
  Matrix3D R = basis.getRotation();
  for (int i = 0; i < 3; i++) {
    e1[i] = R(i,0);
    e2[i] = R(i,1);
    e3[i] = R(i,2);
  }
  return 0;
}

template <int nn, int ndf, typename IsoT>
double
EuclidFrameTransf<nn,ndf,IsoT>::getInitialLength()
{
  return L;
}

template <int nn, int ndf, typename IsoT>
double
EuclidFrameTransf<nn,ndf,IsoT>::getDeformedLength()
{
  return basis.getLength();
}


//
// Pull
//
template <int nn, int ndf, typename IsoT>
int
EuclidFrameTransf<nn,ndf,IsoT>::update()
{
  if (basis.update(nodes) < 0) 
    return -1;

  const Matrix3D R = basis.getRotation();

  for (int i=0; i<nn; i++) {
    Versor q = nodes[i]->getTrialRotation();
    ur[i] = AxisAngle(R^(MatrixFromVersor(q)*R0));
  }

  return 0;
}

template <int nn, int ndf, typename IsoT>
Versor
EuclidFrameTransf<nn,ndf,IsoT>::getNodeRotation(int tag)
{
  return nodes[tag]->getTrialRotation();
}


template <int nn, int ndf, typename IsoT>
Vector3D 
EuclidFrameTransf<nn,ndf,IsoT>::getNodeLocation(int node)
{
  Vector3D xn = basis.getRotation()^nodes[node]->getCrds();

  xn += this->pullPosition<&Node::getTrialDisp>(node);
  
  return xn - basis.getLocation();
}


template <int nn, int ndf, typename IsoT>
Vector3D
EuclidFrameTransf<nn,ndf,IsoT>::getNodePosition(int node)
{
  Vector3D u = this->pullPosition<&Node::getTrialDisp>(node);
  u -= basis.getPosition();
  u += basis.getRotationDelta()^(nodes[node]->getCrds());
  return u;
}


template <int nn, int ndf, typename IsoT>
Vector3D
EuclidFrameTransf<nn,ndf,IsoT>::getNodeRotationLogarithm(int node)
{
  return ur[node].vector;
}


template <int nn, int ndf, typename IsoT>
VectorND<nn*ndf>
EuclidFrameTransf<nn,ndf,IsoT>::getStateVariation()
{

  static VectorND<nn*ndf> ul;
  for (int i=0; i<nn; i++) {
    const Vector &ddu = nodes[i]->getIncrDeltaDisp();
    for (int j = 0; j < ndf; j++) {
      ul[i*ndf+j] = ddu(j);
    }
  }

  const Matrix3D R = basis.getRotation();

  // (1) Global Offsets
  // Do ui -= ri x wi
  if constexpr (ndf >= 6)
    if (offsets && !(offset_flags&OffsetLocal)) [[unlikely]] {
      const std::array<Vector3D, nn>& offset = *offsets;
      for (int i=0; i<nn; i++) {

        const int j = i * ndf;
        Vector3D w {ul[j+3],ul[j+4],ul[j+5]};

        ul.assemble(j, (R*offset[i]).cross(w), -1.0);
      }
    }

  // (2) Rotations and translations

  for (int i=0; i<nn; i++) {
    const int j = i * ndf;
    ul.insert(j  , R^Vector3D{ul[j+0], ul[j+1], ul[j+2]}, 1.0);
    ul.insert(j+3, R^Vector3D{ul[j+3], ul[j+4], ul[j+5]}, 1.0);
  }

  // Projection
  {
    const Vector3D wr = basis.getRotationVariation(ndf, &ul[0]);
    const Vector3D dc = basis.getPositionVariation(ndf, &ul[0]);

    for (int i=0; i<nn; i++) {
      Vector3D ui = this->getNodePosition(i);
      ul.assemble(i*ndf+0, dc, -1.0);
      ul.assemble(i*ndf+0, ui.cross(wr), 1.0);
      ul.assemble(i*ndf+3, wr, -1.0);
    }
  }

  // 3) Offsets
  if constexpr (ndf >= 6)
    if (offsets && (offset_flags&OffsetLocal)) [[unlikely]] {
      const std::array<Vector3D, nn>& offset = *offsets;
      for (int i=0; i<nn; i++) {

        const int j = i * ndf;
        Vector3D w {ul[j+3],ul[j+4],ul[j+5]};

        ul.assemble(j, offset[i].cross(w), -1.0);
      }
    }

  // (5) Logarithm of rotations

  for (int i=0; i<nn; i++) {
    const int j = i * ndf+3;
    Vector3D v {ul[j+0], ul[j+1], ul[j+2]};
    ul.insert(i*ndf+3, ur[i].dLog(v), 1.0);
  }

  return ul;
}

//
// Push
//
template <int nn, int ndf, typename IsoT>
int
EuclidFrameTransf<nn,ndf,IsoT>::push(VectorND<nn*ndf>&p, Operation op)
{
  VectorND<nn*ndf>& pa = p;

  // 1) Logarithm
  for (int i=0; i<nn; i++) {
    Vector3D m {pa[i*ndf + 3], pa[i*ndf + 4], pa[i*ndf + 5]};
    pa.insert(i*ndf + 3, ur[i].dLog()^m, 1.0);
  }


  // 2.1) Sum of moments: m = sum_i mi + sum_i (xi x ni)
  Vector3D m{};
  for (int i=0; i<nn; i++) {
    // m += mi
    for (int j=0; j<3; j++)
      m[j] += pa[i*ndf+3+j];
    // m += xi x ni
    m += this->getNodeLocation(i).cross(Vector3D{pa[i*ndf+0], pa[i*ndf+1], pa[i*ndf+2]});
  }
  // 2.2) Adjust
  for (int i=0; i<nn; i++)
    pa.template assemble<6>(i*ndf, basis.getRotationGradient(i)^m, -1.0);


  // 3,4) Rotate and joint offsets
  Matrix3D R = basis.getRotation();
  for (int i=0; i<nn; i++) {
    const int base = i * ndf;
    pa.insert(base,   R*Vector3D{pa[base  ], pa[base+1], pa[base+2]}, 1.0);
    pa.insert(base+3, R*Vector3D{pa[base+3], pa[base+4], pa[base+5]}, 1.0);
  }

  // Offset
  return 0;
}


template <int nn, int ndf, typename IsoT>
int
EuclidFrameTransf<nn,ndf,IsoT>::push(MatrixND<nn*ndf,nn*ndf>&kb,
                                     const VectorND<nn*ndf>& pb, 
                                     Operation op)
{
  VectorND<nn*ndf> p = pb;

  MatrixND<nn*ndf,nn*ndf> Kb = kb;
  if (1) {//!(offset_flags & LogIter)) {
    for (int i=0; i<nn; i++) {
      Vector3D m{pb[i*ndf+3], pb[i*ndf+4], pb[i*ndf+5]};
      const Matrix3D Ai = ur[i].dLog();
      p.insert(i*ndf+3, Ai^m, 1.0);

      const Matrix3D kg = ur[i].ddLog(m);

      for (int j=0; j<nn; j++) {
        const Matrix3D Aj = ur[j].dLog();
        // loop over 3x3 blocks for n and m
        for (int k=0; k<2; k++) {
          for (int l=0; l<2; l++) {
            if (k == 0 && l == 0)
              continue;

            Matrix3D Kab {{
              Kb(i*ndf+3*k+0, j*ndf+3*l  ), Kb(i*ndf+3*k+1, j*ndf+3*l  ), Kb(i*ndf+3*k+2, j*ndf+3*l  ),
              Kb(i*ndf+3*k+0, j*ndf+3*l+1), Kb(i*ndf+3*k+1, j*ndf+3*l+1), Kb(i*ndf+3*k+2, j*ndf+3*l+1),
              Kb(i*ndf+3*k+0, j*ndf+3*l+2), Kb(i*ndf+3*k+1, j*ndf+3*l+2), Kb(i*ndf+3*k+2, j*ndf+3*l+2)
            }};
            if (k == 1)
              Kab = Ai^Kab; // row rotation block
            if (l == 1)
              Kab = Kab*Aj; // column rotation block

            Kb.insert(Kab, i*ndf+3*k, j*ndf+3*l, 1.0);
            if (i == j && k == 1 && l == 1)
              Kb.assemble(kg, i*ndf+3*k, j*ndf+3*l, 1.0);
          }
        }
      }
    }
  }

  // Kl = A ^ k * A

  MatrixND<nn*ndf,nn*ndf>& Kl = kb;
  const MatrixND<nn*ndf,nn*ndf> A = getProjection();
  Kl.addMatrixTripleProduct(0, A, Kb, 1);


  const VectorND<nn*ndf> Ap = A^p;

  // Kl += Kw * A
  Kb.zero();
  VectorND<12> qwx{};
  for (int i=0; i<nn; i++)
    for (int j=0; j<6; j++)
      qwx[i*6+j] = p[i*ndf+j] - Ap[i*ndf+j];

  if constexpr (ndf == 6) {
    Kl.addMatrixProduct(basis.getRotationJacobian(qwx), A, 1.0);
  }
  else {
    const MatrixND<12,12> Kw = basis.getRotationJacobian(qwx);
    Kb.assemble(Kw.template extract<0, 6,  0, 6>(),   0,   0, 1.0);
    Kb.assemble(Kw.template extract<0, 6,  6,12>(),   0, ndf, 1.0);
    Kb.assemble(Kw.template extract<6,12,  0, 6>(), ndf,   0, 1.0);
    Kb.assemble(Kw.template extract<6,12,  6,12>(), ndf, ndf, 1.0);
    Kl.addMatrixProduct(Kb, A, 1.0);
  }

  //
  // Kl += -W'*Pn'*A  - Pnm * W
  //
  Kb.zero();
  for (int j=0; j<nn; j++) {
    const MatrixND<3,6> Gj = basis.getRotationGradient(j);
    for (int i=0; i<nn; i++) {
      Kb.assemble(Hat(&p[i*ndf+0])*Gj,  i*ndf+0, j*ndf, -1.0);

      // Kl += -Pnm*W
      Kl.assemble(Hat(&Ap[i*ndf+0])*Gj, i*ndf+0, j*ndf, -1.0);
      Kl.assemble(Hat(&Ap[i*ndf+3])*Gj, i*ndf+3, j*ndf, -1.0);
    }
  }
  Kl.addMatrixTransposeProduct(1.0, Kb, A,  -1.0);

  // Kl = diag(R) * Kl * diag(R)^T
  FrameTransform<nn,ndf>::pushRotation(Kl, basis.getRotation());
  return 0;
}



//
// Sensitivity
//
template <int nn, int ndf, typename IsoT>
bool
EuclidFrameTransf<nn,ndf,IsoT>::isShapeSensitivity()
{
  int nodeParameterI = nodes[   0]->getCrdsSensitivity();
  int nodeParameterJ = nodes[nn-1]->getCrdsSensitivity();
  // TODO(sensitivity): implement dvz

  return (nodeParameterI != 0 || nodeParameterJ != 0);
}


template <int nn, int ndf, typename IsoT>
double
EuclidFrameTransf<nn,ndf,IsoT>::getLengthGrad()
{
  const int di = nodes[0]->getCrdsSensitivity();
  const int dj = nodes[1]->getCrdsSensitivity();

  Vector3D dxi{0.0};
  Vector3D dxj{0.0};

  if (di != 0)
    dxi(di-1) = 1.0;
  if (dj != 0)
    dxj(dj-1) = 1.0;

  return 1/L*(xj - xi).dot(dxj - dxi);
}

template <int nn, int ndf, typename IsoT>
double
EuclidFrameTransf<nn,ndf,IsoT>::getd1overLdh()
{
  return -getLengthGrad()/(L*L);
}


template <int nn, int ndf, typename IsoT>
void
EuclidFrameTransf<nn,ndf,IsoT>::Print(OPS_Stream &s, int flag)
{
  if (flag == OPS_PRINT_PRINTMODEL_JSON) {
    s << OPS_PRINT_JSON_MATE_INDENT << "{";
    s << "\"name\": " << this->getTag() << ", ";
    s << "\"type\": \"EuclidFrameTransf\"";
    s << ", \"vecxz\": [" 
      << vz[0] << ", " 
      << vz[1] << ", "
      << vz[2] << "]";
    if (offsets != nullptr) {
      s << ", \"offsets\": [";
      for (int i=0; i<nn; i++) {
        s << "["
          << (*offsets)[i][0] << ", " 
          << (*offsets)[i][1] << ", "
          << (*offsets)[i][2] << "]";
        if (i < nn-1)
          s << ", ";
      }
      s << "]";
    }
    s << "}";

    return;
  }
}

} // namespace OpenSees