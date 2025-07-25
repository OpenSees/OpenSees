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
// Please cite the following resources in any derivative works:
//
// Perez, C.M., and Filippou F.C.. "On Nonlinear Geometric Transformations
//     of Finite Elements" Int. J. Numer. Meth. Engrg. 2024; 
//     https://doi.org/10.1002/nme.7506
//
//===----------------------------------------------------------------------===//

//
// Description: LinearFrameTransf implements a linearized euclidean transformation
// for a spatial frame
//
// Written: Claudio Perez
//
// References:
//
// Perez, C.M., and Filippou F.C.. "On Nonlinear Geometric Transformations
//   of Finite Elements" Int. J. Numer. Meth. Engrg. 2024; 
//   https://doi.org/10.1002/nme.7506
//
// Haukaas, Terje, and Michael H. Scott. 
//   “Shape Sensitivities in the Reliability Analysis of Nonlinear Frame Structures.” 
//   Computers & Structures 84, nos. 15–16 (2006): 964–77. 
//   https://doi.org/10.1016/j.compstruc.2006.02.014.
//
// Scott, Michael H., Paolo Franchin, Gregory L. Fenves, and Filip C. Filippou. 
//   “Response Sensitivity for Nonlinear Beam–Column Elements.” 
//   Journal of Structural Engineering 130, no. 9 (2004): 1281–88. 
//   https://doi.org/10.1061/(asce)0733-9445(2004)130:9(1281).
//
#ifndef LinearFrameTransf_hpp
#define LinearFrameTransf_hpp

#include <array>
#include <FrameTransform.h>
#include <Vector3D.h>
#include <MatrixND.h>

namespace OpenSees {

template <int nn, int ndf>
class LinearFrameTransf: public FrameTransform<nn,ndf>
{
public:
  constexpr static int n = nn*ndf;

  LinearFrameTransf(int tag,
                    const Vector3D &vecxz,
                    const std::array<Vector3D, nn> *offset=nullptr,
                    int offset_flags = 0);

  ~LinearFrameTransf();

  
  const char *getClassType() const {return "LinearFrameTransf";}
  
  virtual int getLocalAxes(Vector3D &x, Vector3D &y, Vector3D &z) const;
  
  virtual FrameTransform<nn,ndf> *getCopy() const;

  double getInitialLength() final;
  double getDeformedLength() final;
  const std::array<Vector3D,nn> *getRigidOffsets() const final {return offsets;}
  
  int initialize(std::array<Node*, nn>& new_nodes) final;
  int update() final;
  int commit() final;
  int revertToLastCommit() final;
  int revertToStart() final;

  Vector3D getNodePosition(int tag) final;
  Vector3D getNodeRotationLogarithm(int tag) final;
  VectorND<nn*ndf> getStateVariation() final;

  using Operation = typename FrameTransform<nn,ndf>::Operation;

  int push(VectorND<nn*ndf>&pl, Operation) final;
  int push(MatrixND<nn*ndf,nn*ndf>& kl, const VectorND<nn*ndf>& pl, Operation) final;


  // Sensitivity
  //
  void   pushGrad(VectorND<nn*ndf>& dp, VectorND<nn*ndf>& pl) override;
  void   pullFixedGrad(VectorND<nn*ndf>&) override;
  void   pullTotalGrad(VectorND<nn*ndf>&, int) override;
  bool   isShapeSensitivity() final;
  double getLengthGrad() final;
  double getd1overLdh() final;

  // TaggedObject
  void Print(OPS_Stream &s, int flag) final;

  // Personal
  Vector3D getDelta() {return Du;}
          
private:

  inline int
  pull(VectorND<nn*ndf>& ug, 
              const Matrix3D& R, 
              const std::array<Vector3D, nn> *offset = nullptr,
              int offset_flags = 0);


  template<const Vector& (Node::*Getter)()>
  const Vector3D
  pullPosition(int node)
  {
    const Vector &u = (nodes[node]->*Getter)();

    Vector3D v;
    for (int i=0; i<3; i++)
      v[i] = u[i];

    // 1) Offsets
    if (offsets) {
      if (!(offset_flags&OffsetLocal)) {
        Vector3D w {u[3], u[4], u[5]};
        v -= offsets->at(node).cross(w);
      }
    }

    // 2) Constant Rotation
    return R^v;
  }

  std::array<Node*, nn> nodes;
  Vector3D Du, ur[nn];

  std::array<Vector3D, nn> *offsets;
  int offset_flags;

  Vector3D xi, xj, vz;
  Matrix3D R;         // rotation matrix
  double L;           // undeformed element length

  std::array<VectorND<ndf>*, nn> u_init;
  bool initialDispChecked;
};
}

#include "LinearFrameTransf.tpp"
#endif

