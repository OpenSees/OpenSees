//===----------------------------------------------------------------------===//
//
//                                   xara
//                              https://xara.so
//===----------------------------------------------------------------------===//
//
//        OpenSees - Open System for Earthquake Engineering Simulation
//
//===----------------------------------------------------------------------===//

//
// Description:  PDeltaFrameTransf generalizes the original PDeltaCrdTransf3d
// class to a more general FrameTransform interface.
//
// Adapted: Remo Magalhaes de Souza (rmsouza@ce.berkeley.edu)
// Created: 04/2000
//
#ifndef PDeltaFrameTransf_h
#define PDeltaFrameTransf_h

#include <array>
#include <FrameTransform.h>
#include "LinearFrameTransf.h"
#include <VectorND.h>
#include <MatrixND.h>

class Vector;
class Matrix;

namespace OpenSees {

template <int nn, int ndf>
class PDeltaFrameTransf: public FrameTransform<nn,ndf>
{
public:

  PDeltaFrameTransf(int tag, 
                    const Vector3D &vecxz,
                    const std::array<Vector3D, nn> *offset=nullptr,
                    int offset_flags = 0);

  ~PDeltaFrameTransf();
  
  const char *getClassType() const {return "PDeltaFrameTransf";}
  
  double getInitialLength()  final;
  double getDeformedLength() final;

  int initialize(std::array<Node*, nn>& new_nodes) final;
  int update() final;
  int commit() final;
  int revertToLastCommit() final;
  int revertToStart() final;

  VectorND<nn*ndf> getStateVariation() final;
  Vector3D getNodePosition(int tag) final;
  Vector3D getNodeRotationLogarithm(int tag) final;
  const std::array<Vector3D,nn> *getRigidOffsets() const final { return linear.getRigidOffsets();}
#if 0
  VectorND<nn*ndf>    pushResponse(VectorND<nn*ndf>&pl) final;
  MatrixND<nn*ndf,nn*ndf> pushResponse(MatrixND<nn*ndf,nn*ndf>& kl, const VectorND<nn*ndf>& pl) final;
#else
  using Operation = typename FrameTransform<nn,ndf>::Operation;
  int push(VectorND<nn*ndf>&pl, Operation) final;
  int push(MatrixND<nn*ndf,nn*ndf>& kl, const VectorND<nn*ndf>& pl, Operation) final;
#endif 

  FrameTransform<nn,ndf> *getCopy() const final;

  int getLocalAxes(Vector3D &x, Vector3D &y, Vector3D &z) const final;

  // Sensitivity
  double getLengthGrad() final;

  // Tagged Object
  void Print(OPS_Stream &s, int flag) final;

private:
  int offset_flags;
  LinearFrameTransf<nn,ndf> linear;

};

} // namespace OpenSees

#include "PDeltaFrameTransf3d.tpp"
#endif
