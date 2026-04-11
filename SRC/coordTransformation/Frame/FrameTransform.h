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
// Written: cmp
//
#ifndef FrameTransform_h
#define FrameTransform_h

#include <VectorND.h>
#include <MatrixND.h>
#include <Matrix3D.h>
#include <TaggedObject.h>

using OpenSees::VectorND;
using OpenSees::MatrixND;
using OpenSees::Matrix3D;
class Information;
class Response;
class Node;

enum {
 CRDTR_TAG_CorotFrameTransfWarping3d,
 CRDTR_TAG_CorotFrameTransf3d,
 CRDTR_TAG_LinearFrameTransf3d,
 CRDTR_TAG_PDeltaFrameTransf3d
};

enum {
  OffsetGlobal     = 0, // 1<<0,
  OffsetLocal      = 1, // 1<<1,
  OffsetNormalized = 2, // 1<<2,

  LogIter          = 1<<3,
  LogIncr          = 1<<4,
  LogInit          = 1<<5,
  LogDefault       = 1<<6
};

namespace OpenSees {

template <int nn, int ndf>
class FrameTransform : public TaggedObject
{
public:
  explicit FrameTransform(int tag) : TaggedObject(tag) {}

  enum class Operation {
    Total       =    0,
    Logarithm   = 1<<0,
    LocalOffset = 1<<1,
    Isometry    = 1<<2,
    Rotation,
    GlobalOffset,
    Exponential
  };

  virtual FrameTransform<nn,ndf> *getCopy() const =0;


  virtual int initialize(std::array<Node*, nn>& nodes)=0;
  virtual int update() =0;
  virtual int commit() =0;
  virtual int revertToLastCommit() =0;
  virtual int revertToStart() =0;

  virtual Vector3D  getNodePosition(int tag) =0;
  virtual Vector3D  getNodeRotationLogarithm(int tag) =0;
  virtual VectorND<nn*ndf> getStateVariation() =0; // pull
  virtual int push(VectorND<nn*ndf>&pl, Operation=0) =0;
  virtual int push(MatrixND<nn*ndf,nn*ndf>& kl, const VectorND<nn*ndf>& pl, Operation=0) =0;

  virtual double getInitialLength() =0;
  virtual double getDeformedLength() =0;
  virtual const std::array<Vector3D,nn> *getRigidOffsets() const =0;

  //
  virtual int getLocalAxes(Vector3D &x, Vector3D &y, Vector3D &z) const =0;

  Vector3D getNormalVector() const {
    Vector3D x, y, z;
    if (getLocalAxes(x, y, z) < 0)
      return Vector3D{{0.0, 0.0, 1.0}};
    return z;
  }

  // Recorders
  virtual Response *setResponse(const char **argv, int argc, OPS_Stream &) {
    return nullptr;
  }
  virtual int getResponse(int responseID, Information &) {
    return -1;
  }

  // Sensitivity
  virtual void   pushGrad(VectorND<nn*ndf>& dp, VectorND<nn*ndf>& pl) {}
  virtual void   pullFixedGrad(VectorND<nn*ndf>&) {}
  virtual void   pullTotalGrad(VectorND<nn*ndf>&, int) {}
  virtual bool   isShapeSensitivity() {return false;}
  virtual double getLengthGrad() {return 0.0;}
  virtual double getd1overLdh() {return 0.0;}


  // deprecatred API
  virtual VectorND<nn*ndf>    pushResponse(VectorND<nn*ndf>&pl) final {
    VectorND<nn*ndf> pg{pl};
    push(pg, Operation::Total);
    return pg;
  }

  virtual MatrixND<nn*ndf,nn*ndf> pushResponse(MatrixND<nn*ndf,nn*ndf>& kl,
                                               const VectorND<nn*ndf>& pl) final {
    MatrixND<nn*ndf,nn*ndf> kg{kl};
    push(kg, pl, Operation::Total);
    return kg;              
  }



protected:
  constexpr static int ndm = 3;
  static inline constexpr void
  pushRotation(MatrixND<nn*ndf,nn*ndf>& Kg, const Matrix3D& R);

  static inline constexpr void
  pushOffsets(MatrixND<nn*ndf,nn*ndf>& Kg, const std::array<Vector3D,nn>& offsets);

  static int
  Orient(const Vector3D& dx, const Vector3D& vz, Matrix3D &R) {

    // calculate the element local x axis components wrt to the global coordinates

    Vector3D e1 = dx/dx.norm();

    //
    Vector3D e2 = vz.cross(e1);

    const double ynorm = e2.norm();

    if (ynorm == 0.0)
        return -1;

    e2 /= ynorm;

    Vector3D e3 = e1.cross(e2);

    for (int i = 0; i < 3; i++) {
      R(i,0) = e1[i];
      R(i,1) = e2[i];
      R(i,2) = e3[i];
    }
    return 0;
  }

  VectorND<nn*ndf>    pushConstant(const VectorND<nn*ndf>&pl);
  MatrixND<nn*ndf,nn*ndf> pushConstant(const MatrixND<nn*ndf,nn*ndf>& kl);

};
}
#include "FrameTransform.tpp"

#endif // include guard
