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
// The purpose of this class is to wrap the more general FrameTransform<>
// templates to reproduce the legacy CrdTransf classes that were derived
// for elements in a "basic" coordinate system.
//
// Written: Claudio M. Perez
//
#ifndef BasicFrameTransf3d_h
#define BasicFrameTransf3d_h

#include <CrdTransf.h>
#include <FrameTransform.h>
#include <LinearFrameTransf.h>

class Vector;
class Matrix;

namespace OpenSees {

template<int ndf=6>
class BasicFrameTransf3d: public CrdTransf
{
public:
  explicit BasicFrameTransf3d(FrameTransform<2,ndf> *t);

  ~BasicFrameTransf3d() override;

  int getLocalAxes(Vector &x, Vector &y, Vector &z) final;

  CrdTransf *getCopy3d() final;

  double getInitialLength() final;
  double getDeformedLength() final;

  int initialize(Node *ni, Node *nj) final;
  int update() final;
  int commitState() final;
  int revertToLastCommit() final;
  int revertToStart() final;

  const Vector &getBasicTrialDisp()     final;
  const Vector &getBasicIncrDisp()      final;
  const Vector &getBasicIncrDeltaDisp() final;
  const Vector &getBasicTrialVel()      final;

  const Vector &getGlobalResistingForce(const Vector &basicForce, const Vector &p0) final;
  const Matrix &getGlobalStiffMatrix(const Matrix &basicStiff, const Vector &basicForce) final;
  const Matrix &getInitialGlobalStiffMatrix(const Matrix &basicStiff) final;

  // rotate consistent mass matrix
  const Matrix &getGlobalMatrixFromLocal(const Matrix &local) final;
 
  // methods used in post-processing only
  const Vector &getPointGlobalCoordFromLocal(const Vector &localCoords);
  const Vector &getPointGlobalDisplFromBasic(double xi, const Vector &basicDisps);
  const Vector &getPointLocalDisplFromBasic(double xi, const Vector &basicDisps);    

  //
  // Sensitivity
  //
  const Vector & getBasicTrialDispShapeSensitivity() final;
  const Vector & getBasicDisplSensitivity(int grad) final;
  const Vector &getGlobalResistingForceShapeSensitivity(const Vector &basicForce,
                                                        const Vector &p0, int grad) final;
  bool isShapeSensitivity() final;
  double getdLdh() final;
  double getd1overLdh() final;


  // MovableObject
  int sendSelf(int tag, Channel &) final;
  int recvSelf(int tag, Channel &, FEM_ObjectBroker &) final;
  const char *getClassType() const final {
    return "BasicFrameTransf3d";
  }
    
  // TaggedObject
  void Print(OPS_Stream &s, int flag) final;



private:
  using Operation = typename FrameTransform<2,ndf>::Operation;
  constexpr static int NBV = 6;
  constexpr static int NDF = ndf;
  enum : int {
      inx = -12, //  0
      iny = -12, //  1
      inz = -12, //  2
      imx = -12, //  3
      imy =   3, //  4
      imz =   1, //  5
      iwx =   6, //  6
      jnx =   0, //  6
      jny = -12, //  7
      jnz = -12, //  8
      jmx =   5, //  9
      jmy =   4, // 10
      jmz =   2, // 11
      jwx =   7,
  };


  static constexpr std::array<int, NDF*2> set_indices() {
    if constexpr (ndf-6 > 0) {
      return {
        inx, iny, inz, imx, imy, imz, iwx,
        jnx, jny, jnz, jmx, jmy, jmz, jwx
      };
    } else {
      return {
        inx, iny, inz, imx, imy, imz,
        jnx, jny, jnz, jmx, jmy, jmz
      };
    }
  }
  static constexpr auto iq = set_indices();

  FrameTransform<2,ndf> &t;
  LinearFrameTransf<2,ndf> linear;
};
} // namespace OpenSees

#include "BasicFrameTransf.tpp"
#endif
