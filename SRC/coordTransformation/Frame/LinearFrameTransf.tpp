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
// Description: LinearFrameTransf is a linear
// transformation for a space frame between the global
// and basic coordinate systems
//
// Written: Claudio M. Perez
//
#pragma once
#include <Vector.h>
#include <Versor.h>
#include <Matrix.h>
#include <Vector3D.h>
#include <Matrix3D.h>
#include <Node.h>
#include <OPS_Globals.h>
#include <GroupSO3.h>
#include "LinearFrameTransf.h"

namespace OpenSees {

static inline MatrixND<3,3>
FrameOrientationGradient(const Vector3D& xi,
                         const Vector3D& xj, 
                         const Vector3D& vz,
                         int di, int dj, int dv)
{
  Vector3D v1  = xj - xi;
  const double L     = v1.norm();
  Vector3D e1  = v1/L;

  Vector3D v2  = vz.cross(e1);

  Vector3D e2 = v2 / v2.norm();

  //
  Vector3D dvz{0.0};
  Vector3D dxi{0.0};
  Vector3D dxj{0.0};

  if (di != 0)
    dxi(di-1) = 1.0;
  if (dj != 0)
    dxj(dj-1) = 1.0;
  if (dv != 0)
    dvz(dv-1) = 1.0;

  double   dL  = 1/L*(xj - xi).dot(dxj - dxi);
  Vector3D dv1 = dxj - dxi;
  Vector3D de1 = 1/(L*L)*(dv1*L - v1*dL);

  double L2    = v2.norm();
  Vector3D dv2 = dvz.cross(e1) + vz.cross(de1);
  double dL2   = 1/L2*v2.dot(dv2);
  Vector3D de2 = 1/(L2*L2)*(dv2*L2 - v2*dL2);

  Vector3D de3 = de1.cross(e2) + e1.cross(de2);

  MatrixND<3,3> dR;
  dR(0,0) = de1(0);
  dR(1,0) = de1(1);
  dR(2,0) = de1(2);

  dR(0,1) = de2(0);
  dR(1,1) = de2(1);
  dR(2,1) = de2(2);

  dR(0,2) = de3(0);
  dR(1,2) = de3(1);
  dR(2,2) = de3(2);

  return dR;
}

template <int nn, int ndf>
LinearFrameTransf<nn,ndf>::LinearFrameTransf(int tag, 
                                             const Vector3D &vecxz, 
                                             const std::array<Vector3D, nn> *offset,
                                             int offset_flags)
  : FrameTransform<nn,ndf>(tag),
    Du{},
    ur{},
    L(0.0),
    offsets{nullptr},
    offset_flags(offset_flags),
    u_init{nullptr}, 
    initialDispChecked(false)
{
  R.zero();

  for (int i=0; i<3; i++)
    vz[i] = vecxz[i];

  R(0,2) = vz(0);
  R(1,2) = vz(1);
  R(2,2) = vz(2);

  // Rigid joint offsets
  if (offset != nullptr) {
    offsets = new std::array<Vector3D, nn>{};
    *offsets = *offset;
  }
}



template <int nn, int ndf>
LinearFrameTransf<nn,ndf>::~LinearFrameTransf()
{
  if (offsets != nullptr)
    delete offsets;

  for (int i=0; i<nn; i++)
    if (u_init[i] != nullptr)
      delete u_init[i];
}


template <int nn, int ndf>
int
LinearFrameTransf<nn,ndf>::commit()
{
  return 0;
}


template <int nn, int ndf>
int
LinearFrameTransf<nn,ndf>::revertToLastCommit()
{
  return 0;
}

template <int nn, int ndf>
int
LinearFrameTransf<nn,ndf>::revertToStart()
{
  return 0;
}


template <int nn, int ndf>
int
LinearFrameTransf<nn,ndf>::initialize(std::array<Node*, nn>& new_nodes)
{

  for (int i=0; i<nn; i++) {
    nodes[i] = new_nodes[i];
    if (nodes[i] == nullptr) {
      opserr << "invalid pointers to the element nodes\n";
      return -1;
    }
  }

  // Check for initial displacements at nodes
  if (initialDispChecked == false) {
    for (int i = 0; i<nn; i++) {
      const Vector &u = nodes[i]->getDisp();
      
      for (int j = 0; j < ndf; j++)
        if (u(j) != 0.0) {
          u_init[i] = new VectorND<ndf>{};
          for (int l = 0; l < ndf; l++)
            (*u_init[i])[l] = u(l);
          // break out
          j = ndf;
        }
    }
    initialDispChecked = true;
  }

  {
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


    if (u_init[0] != 0) {
      for (int i=0; i<3; i++)
        dx(i) -= (*u_init[0])[i];
    }

    if (u_init[nn-1] != 0) {
      for (int i=0; i<3; i++)
        dx(i) += (*u_init[nn-1])[i];
    }

    L = dx.norm();

    if (L == 0.0)
      return -2;

    return FrameTransform<nn,ndf>::Orient(dx, vz, R);
  }
}


template <int nn, int ndf>
int
LinearFrameTransf<nn,ndf>::getLocalAxes(Vector3D &e1, Vector3D &e2, Vector3D &e3) const
{
  for (int i = 0; i < 3; i++) {
    e1[i] = R(i,0);
    e2[i] = R(i,1);
    e3[i] = R(i,2);
  }
  return 0;
}

template <int nn, int ndf>
double
LinearFrameTransf<nn,ndf>::getInitialLength()
{
  return L;
}

template <int nn, int ndf>
double
LinearFrameTransf<nn,ndf>::getDeformedLength()
{
  return L+Du[0];
}


//
// Pull
//

template <int nn, int ndf>
int
LinearFrameTransf<nn,ndf>::update()
{
  constexpr Vector3D iv{1, 0, 0};

  Du  = this->pullPosition<&Node::getTrialDisp>(nn-1);
  Du -= this->pullPosition<&Node::getTrialDisp>(0);

  for (int node =0 ; node < nn; node++) {
    const Vector& u = nodes[node]->getTrialDisp();
    for (int i=0; i<3; i++) {
      ur[node][i] = 0;
      for (int j=0; j<3; j++)
        ur[node][i] += R(j,i) * u[3+j];
    }

    ur[node].addCrossProduct(iv, Du, -1.0/L);
  }
  return 0;
}


template <int nn, int ndf>
int
LinearFrameTransf<nn,ndf>::pull(VectorND<nn*ndf>& ul, 
             const Matrix3D& R,
             const std::array<Vector3D, nn> *offset,
             int offset_flags) 
{

  // (1)

  // (2) Global offsets
  //
  // Do ui -= ri x wi
  if constexpr (ndf >= 6)
    if (offset && !(offset_flags&OffsetLocal)) [[unlikely]] {
      const std::array<Vector3D, nn>& offsets = *offset;
      for (int i=0; i<nn; i++) {

        const int j = i * ndf;
        Vector3D w {ul[j+3],ul[j+4],ul[j+5]};

        ul.assemble(j, offsets[i].cross(w), -1.0);
      }
    }
  
  // (3) Rotation
  for (int i=0; i<nn; i++) {
    const int j = i * ndf;
    ul.insert(j  , R^Vector3D{ul[j+0], ul[j+1], ul[j+2]}, 1.0);
    ul.insert(j+3, R^Vector3D{ul[j+3], ul[j+4], ul[j+5]}, 1.0);
  }

  // 4) Offsets
  // Do ul -= R^ri x wi
  if constexpr (ndf >= 6)
    if (offset && (offset_flags&OffsetLocal)) [[unlikely]] {
      const std::array<Vector3D, nn>& offsets = *offset;
      for (int i=0; i<nn; i++) {

        const int j = i * ndf;
        Vector3D w {ul[j+3],ul[j+4],ul[j+5]};

        ul.assemble(j, offsets[i].cross(w), -1.0);
      }
    }

  // (5) Isometry
  // TODO (nn>2)
  constexpr static Vector3D iv {1,0,0};
  Vector3D uI = ul.template extract<3>(0);
  Vector3D Du = ul.template extract<3>((nn-1)*ndf) - uI;
  Vector3D ixDu = iv.cross(Du);
  for (int i=0; i<nn; i++) {
    // Translation
    ul.assemble(i*ndf, uI, -1.0);
    for (int j=1; j<3; j++)
      ul[i*ndf+j] -= double(i)/(nn-1.0)*Du[j];

    // Rotation
    ul.assemble(i*ndf+3, ixDu, -1.0/L);
  }

  return 0;
}

template <int nn, int ndf>
VectorND<nn*ndf>
LinearFrameTransf<nn,ndf>::getStateVariation()
{
  VectorND<nn*ndf> ug;
  for (int i=0; i<nn; i++) {
    const Vector &ddu = nodes[i]->getIncrDeltaDisp();
    for (int j = 0; j < ndf; j++) {
      ug[i*ndf+j] = ddu(j);
    }
  }
  LinearFrameTransf<nn,ndf>::pull(ug, R, offsets, offset_flags);
  return ug;
}


template <int nn, int ndf>
Vector3D
LinearFrameTransf<nn,ndf>::getNodePosition(int node)
{
  Vector3D v{};
  if (node == nn-1) {
    v[0] = Du[0];
  }
  // TODO(nn>2)
  return v;
}


template <int nn, int ndf>
Vector3D
LinearFrameTransf<nn,ndf>::getNodeRotationLogarithm(int node)
{
  return ur[node];
}


//
// Push
//
template <int nn, int ndf>
int
LinearFrameTransf<nn,ndf>::push(VectorND<nn*ndf>&p, Operation op)
{
  VectorND<nn*ndf> pa = p; // NOTE
  constexpr Vector3D iv{1, 0, 0};

  // 1.1) Sum of moments: m = sum_i mi + sum_i (xi x ni)
  Vector3D m{};
  for (int i=0; i<nn; i++) {
    // m += mi
    for (int j=0; j<3; j++)
      m[j] += p[i*ndf+3+j];

    const Vector3D n = Vector3D{p[i*ndf+0], p[i*ndf+1], p[i*ndf+2]};
    m.addVector(1, iv.cross(n), double(i)/double(nn-1)*L);
  }
  const Vector3D ixm = iv.cross(m);

  // 1.2) Adjust force part
  for (int i=0; i<nn; i++) {
    pa.assemble(i*ndf,  ixm,  (i? 1.0:-1.0)/L);
    pa[i*ndf+3] += m[0]*(i? -1.0:1.0)*0.5;
  }

  if (op == Operation::Isometry)
    return 0;

  // 2) Rotate and do joint offsets
  p = this->FrameTransform<nn,ndf>::pushConstant(pa);
  return 0;
}

template <int nn, int ndf>
int
LinearFrameTransf<nn,ndf>::push(MatrixND<nn*ndf,nn*ndf>&kb, 
                                const VectorND<nn*ndf>&, 
                                Operation op)
{

  static constexpr Vector3D axis{1, 0, 0};
  static constexpr Matrix3D ix  = Hat(Vector3D{1, 0, 0});
  static constexpr Matrix3D ioi = axis.bun(Vector3D{1, 0, 0});
  const Matrix3D RT = R.transpose();
  const Matrix3D ixRT = ix*RT;

  MatrixND<nn*ndf,nn*ndf> A{};
  if constexpr (ndf > 6)
    A.addDiagonal(1.0);

  {
    MatrixND<3,ndf> Gb{};
    Gb.template insert<0, 3>(ioi*RT, 0.5);
    Repeat<nn> ([&](auto a_) {
      constexpr static int a = a_.value;
      A.template insert<a*ndf  , a*ndf  >(RT);
      A.template insert<a*ndf+3, a*ndf+3>(RT);
      Repeat<nn> ([&](auto b_) {
        constexpr static int b = b_.value;

        // TODO(nn>2): Interpolate coordinate
        if constexpr (b == 0)
          Gb.template insert<0,0>(ixRT, -1/L);
        else if constexpr (b == nn-1)
          Gb.template insert<0,0>(ixRT,  1/L);

        // TODO(nn>2): Interpolate coordinate
        A.assemble(ix*Gb, a*ndf  , b*ndf,  double(a)/double(nn-1)*L);
        A.assemble(   Gb, a*ndf+3, b*ndf, -1.0);
      });
    });
  }

  // MatrixND<nn*ndf,nn*ndf> kl = kb;
  const MatrixND<nn*ndf,nn*ndf> KA = kb*A;
  kb.addMatrixTransposeProduct(0.0, A, KA, 1.0);

  if (offsets != nullptr) [[unlikely]]
    this->pushOffsets(kb, *offsets);

  return 0;
}


template <int nn, int ndf>
FrameTransform<nn,ndf> *
LinearFrameTransf<nn,ndf>::getCopy() const
{

  Vector3D xz;
  xz(0) = R(0,2);
  xz(1) = R(1,2);
  xz(2) = R(2,2);

  LinearFrameTransf *theCopy = new LinearFrameTransf<nn,ndf>(this->getTag(), xz, offsets);

  theCopy->nodes = nodes;
  theCopy->L     = L;
  for (int i = 0; i < 3; i++)
    for (int j = 0; j < 3; j++)
      theCopy->R(j,i) = R(j,i);

  return theCopy;
}


//
// Sensitivity
//
template <int nn, int ndf>
bool
LinearFrameTransf<nn,ndf>::isShapeSensitivity()
{
  int nodeParameterI = nodes[   0]->getCrdsSensitivity();
  int nodeParameterJ = nodes[nn-1]->getCrdsSensitivity();
  // TODO(sensitivity): implement dvz

  return (nodeParameterI != 0 || nodeParameterJ != 0);
}


template <int nn, int ndf>
double
LinearFrameTransf<nn,ndf>::getLengthGrad()
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


template <int nn, int ndf>
double
LinearFrameTransf<nn,ndf>::getd1overLdh()
{
  return -getLengthGrad()/(L*L);
}

template <int nn, int ndf>
void
LinearFrameTransf<nn,ndf>::pushGrad(VectorND<nn*ndf>& dp,
                                    VectorND<nn*ndf>& pl)
{  
  //
  // dp += T_{gl} dpl
  //
  double dL = this->getLengthGrad();
  double doneOverL = -dL/(L*L);

  constexpr Vector3D iv{1, 0, 0};

  // 1.1) Sum of moments: m = sum_i mi + sum_i (xi x ni)
  Vector3D m{};
  for (int i=0; i<nn; i++) {
    // m += mi
    for (int j=0; j<3; j++)
      m[j] += pl[i*ndf+3+j];

    const Vector3D n = Vector3D{pl[i*ndf+0], pl[i*ndf+1], pl[i*ndf+2]};
    m.addVector(1, iv.cross(n), double(i)/double(nn-1)*L);
  }
  const Vector3D ixm = iv.cross(m);

  // 1.2) Adjust force part
  for (int i=0; i<nn; i++)
    dp.assemble(i*ndf,  ixm,  (i? 1.0:-1.0)*doneOverL);


  // 2) Rotate and do joint offsets

  //
  // dp = T_{lg}' pl
  //
  int dv = 0; // TODO
  int di = nodes[0]->getCrdsSensitivity();
  int dj = nodes[1]->getCrdsSensitivity();

  for (int i=0; i<nn; i++) {
    const int base = i * ndf;
    dp.assemble(base,    R*Vector3D{dp[base  ], dp[base+1], dp[base+2]}, 1.0);
    dp.assemble(base+3,  R*Vector3D{dp[base+3], dp[base+4], dp[base+5]}, 1.0);
  }

  this->push(pl, Operation::Isometry);
  Matrix3D dR = FrameOrientationGradient(xi, xj, vz, di, dj, dv);
  
  for (int i=0; i<nn; i++) {
    const int base = i * ndf;
    dp.assemble(base,   dR*Vector3D{pl[base  ], pl[base+1], pl[base+2]}, 1.0);
    dp.assemble(base+3, dR*Vector3D{pl[base+3], pl[base+4], pl[base+5]}, 1.0);
  }


  return;
}


template <int nn, int ndf>
void
LinearFrameTransf<nn,ndf>::pullFixedGrad(VectorND<nn*ndf>& du)
{
  //
  // dub += (T_{bl}' T_{lg}  +   T_{bl} T_{lg}') * ug
  //

  //
  // Form ug
  //
  VectorND<nn*ndf> ug;
  for (int i = 0; i < nn; i++) {
    const Vector& u = nodes[i]->getTrialDisp();
    for (int j = 0; j < ndf; j++) {
      ug[i*ndf+j] = u(j);
    }
  }

  int dv = 0; // TODO(sensitivity)

  int di = nodes[0]->getCrdsSensitivity();
  int dj = nodes[1]->getCrdsSensitivity();

  //
  // du = Tbl dR^ug
  {
    Matrix3D dR = FrameOrientationGradient(xi, xj, vz, di, dj, dv);
  
    VectorND<nn*ndf> u1 = ug;
    LinearFrameTransf<nn,ndf>::pull(u1, dR, offsets, offset_flags);
    du += u1;
  }
  return;
}


template <int nn, int ndf>
void
LinearFrameTransf<nn,ndf>::pullTotalGrad(VectorND<nn*ndf>& du, int gradNumber)
{
  for (int n=0; n<nn; n++)
    for (int i = 0; i < ndf; i++) {
      du[i + n*ndf] = nodes[n]->getDispSensitivity((i + 1), gradNumber);
    }

  // dub = T_{bl} T_{lg} * ug'
  LinearFrameTransf<nn,ndf>::pull(du, R, offsets, offset_flags);

  // dub += (T_{bl}' T_{lg}  +   T_{bl} T_{lg}') * ug
  pullFixedGrad(du);

  return;
}


template <int nn, int ndf>
void
LinearFrameTransf<nn,ndf>::Print(OPS_Stream &s, int flag)
{
  if (flag == OPS_PRINT_PRINTMODEL_JSON) {
    s << OPS_PRINT_JSON_MATE_INDENT << "{";
    s << "\"name\": " << this->getTag() << ", ";
    s << "\"type\": \"LinearFrameTransf\"";
    s << ", \"vecxz\": [" 
      << R(0,2) << ", " 
      << R(1,2) << ", "
      << R(2,2) << "]";
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

  if (flag == OPS_PRINT_CURRENTSTATE) {
    s << "\nFrameTransform: " << this->getTag() << " Type: LinearFrameTransf\n";
    s << "\tOrientation: " << Matrix(&R(0,0), 3,3) << "\n";
  }
}

} // namespace OpenSees
