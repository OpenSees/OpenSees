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
// Written: Claudio M. Perez, 
//          Filip C. Filippou
//          University of California, Berkeley
//
// Developed with FEDEASLab [2].
//
// References:
//
// [1] Perez, C.M., and Filippou F.C.. (2024) 
//     "On Nonlinear Geometric Transformations of Finite Elements" 
//     Int. J. Numer. Meth. Engrg.; https://doi.org/10.1002/nme.7506
//
// [2] Filippou, F.C. (1998)
//     "FEDEASLab: Finite Elements for Design Evaluation and Analysis of Structures"
//
#pragma once
#include <array>
#include <Vector3D.h>
#include <Matrix3D.h>
#include <MatrixND.h>
#include "EuclidIsometry.h"

template <int nn>
class SphericalIsometry : public Isometry {
public:
  SphericalIsometry(std::array<Node*,nn>& nodes, const Vector3D& vecxz);

  int initialize() final;
  int update() final;

  double getLength() const final;
  Matrix3D getRotation() const final;
  Vector3D getPosition() final;
  Vector3D getPositionVariation(int ndf, double* du) final;
  Vector3D getRotationVariation(int ndf, double* du) final;
  Matrix3D getRotationDelta() final;
  MatrixND<3,6> getRotationGradient(int node) final;
private:
  double L;
  double Ln;
  std::array<Node*,nn> nodes;
  Vector3D vz;   // vector in the x-z plane
  Vector3D dX;   // deformed length vector
  Vector3D c[2]; // current position of the center
  Matrix3D R[2]; // rotation matrices at the current and previous time step
  Vector3D Xc;   // center of the isometry
};


template <int nn>
SphericalIsometry<nn>::SphericalIsometry(std::array<Node*,nn>& nodes, const Vector3D& vecxz)
: nodes(nodes), vz(vecxz), Xc{},
  c{}, R{}, dX{}, Ln(0.0), L(0.0), 
  offsets(nullptr), offset_flags(0), ic(0)
{

}


template <int nn>
int
SphericalIsometry<nn>::update() {

  Vector3D e1 = dX;
  {
    //
    // Update state
    //
    const Vector& uI = nodes[   0]->getTrialDisp();
    const Vector& uJ = nodes[nn-1]->getTrialDisp();
    for (int k = 0; k < 3; k++)
      e1[k] += uJ(k) - uI(k);

    if (offsets != nullptr) [[unlikely]] {
      e1.addVector(1.0, (*offsets)[   0],  1.0);
      e1.addVector(1.0, nodes[0]->getTrialRotation().rotate((*offsets)[0]), -1.0);
      e1.addVector(1.0, (*offsets)[nn-1], -1.0);
      e1.addVector(1.0, nodes[nn-1]->getTrialRotation().rotate((*offsets)[nn-1]), 1.0);
    }

    // Calculate the deformed length
    Ln = e1.norm();

    if (Ln == 0.0) [[unlikely]] {
        opserr << "\nSouzaFrameTransf: deformed length is 0.0\n";
        return -2;
    }

    e1 /= Ln;
  }

  {
    constexpr static Vector3D D2 {0,1,0};
    const Vector3D E2 = R[init]*D2;
    Vector3D e2 = MatrixFromVersor(nodes[0]->getTrialRotation())*E2; //*R[init];
    e2.addVector(0.5, MatrixFromVersor(nodes[1]->getTrialRotation())*E2, 0.5);
    n = e2[0]/e2[1];
    Vector3D e3 = e1.cross(e2);
    e3 /= e3.norm();

    e2 = e3.cross(e1);

    for (int i = 0; i < 3; i++) {
      R[pres](i,0) = e1[i];
      R[pres](i,1) = e2[i];
      R[pres](i,2) = e3[i];
    }
    Vector3D e2 = vz.cross(e1);
  }

  Vector3D uc = nodes[ic]->getTrialDisp();
  if (offsets != nullptr) {
      uc.addVector(1.0, (*offsets)[ic], -1.0);
      uc.addVector(1.0, nodes[ic]->getTrialRotation().rotate((*offsets)[ic]), 1.0);
  }
  Vector3D X = nodes[ic]->getCrds();
  c[pres] = R[pres]^(X + uc);
  return 0;
}