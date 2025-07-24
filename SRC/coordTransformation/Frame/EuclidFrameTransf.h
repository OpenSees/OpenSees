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
// Description: This file contains the implementation for the 
// EuclidFrameTransf class. EuclidFrameTransf is a euclidean transformation 
// of 3D space.
// When used with the RankinIsometry, it furnishes an improved corotational
// transformation for 3D frames.
//
//
// Written: Claudio M. Perez
// Created: 04/2025
//
#ifndef EuclidFrameTransf_hpp
#define EuclidFrameTransf_hpp

#include <array>
#include <AxisAngle.h>
#include <FrameTransform.h>
#include <Vector3D.h>
#include <MatrixND.h>

namespace OpenSees {

template <int nn, int ndf, typename IsoT>
class EuclidFrameTransf: public FrameTransform<nn,ndf>
{
public:
  constexpr static int n = nn*ndf;

  EuclidFrameTransf(int tag, 
                    const Vector3D &vecxz,
                    const std::array<Vector3D, nn> *offset=nullptr,
                    int offset_flags = 0);

  ~EuclidFrameTransf();

  using Operation = typename FrameTransform<nn,ndf>::Operation;

  const char *getClassType() const {return "EuclidFrameTransf";}
  
  virtual int getLocalAxes(Vector3D &x, Vector3D &y, Vector3D &z) const;
  
  FrameTransform<nn,ndf> *getCopy() const final;

  double getInitialLength() final;
  double getDeformedLength() final;
  const std::array<Vector3D,nn> *getRigidOffsets() const final {return offsets;}
  
  int initialize(std::array<Node*, nn>& new_nodes) final;
  int update() final;
  int commit() final;
  int revertToLastCommit() final;
  int revertToStart() final;

  VectorND<nn*ndf> getStateVariation() final;
  Vector3D getNodePosition(int tag) final;
  Versor   getNodeRotation(int tag) /* final */;
  Vector3D getNodeRotationLogarithm(int tag) final;

  int push(VectorND<nn*ndf>&pl, Operation) final;
  int push(MatrixND<nn*ndf,nn*ndf>& kl, const VectorND<nn*ndf>& pl, Operation) final;

  // Sensitivity
  bool isShapeSensitivity() final;
  double getLengthGrad() final;
  double getd1overLdh() final;

  // TaggedObject
  void Print(OPS_Stream &s, int flag) final;


private:

  Vector3D getNodeLocation(int tag);

  inline MatrixND<nn*ndf,nn*ndf> 
  getProjection() {

    MatrixND<nn*ndf,nn*ndf> A{};
    A.addDiagonal(1.0);

    MatrixND<3,ndf> Gb{};
    for (int a = 0; a<nn; a++) {
      for (int b = 0; b<nn; b++) {
        
        Gb.template insert<0,0>(basis.getRotationGradient(b), 1.0);
        Matrix3D Xa = Hat(this->getNodeLocation(a));
        A.assemble(Xa*Gb, a*ndf  , b*ndf,  1.0);
        A.assemble(   Gb, a*ndf+3, b*ndf, -1.0);
      }
    }

    return A;
  }

  template<const Vector& (Node::*Getter)()>
  const Vector3D
  pullPosition(int node)
  {
    const Vector &u = (nodes[node]->*Getter)();

    Vector3D v;
    for (int i=0; i<3; i++)
      v[i] = u[i];

    // 1) Offsets
    if (offsets) [[unlikely]] {
      if (!(offset_flags&OffsetLocal))  {
        Vector3D w {u[3], u[4], u[5]};
        v -= offsets->at(node).cross(w);
      }
    }

    // 2) Constant Rotation
    Matrix3D R = basis.getRotation();
    return R^v;
  }

  std::array<Node*, nn> nodes;
  std::array<AxisAngle, nn> ur; // rotation vector
  // std::array<Vector3D,  nn> ux; // displacement vector

  std::array<Vector3D, nn> *offsets;
  int offset_flags;
  Matrix3D R0;
  Vector3D xi, xj, vz;
  double L;           // undeformed element length

  IsoT basis;
};

} // namespace OpenSees

#include "EuclidFrameTransf.tpp"

#endif
