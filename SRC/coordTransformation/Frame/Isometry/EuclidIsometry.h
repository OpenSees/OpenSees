//===----------------------------------------------------------------------===//
//
//                                   xara
//                              https://xara.so
//----------------------------------------------------------------------------//
//
//                                 FEDEASLab
//       Finite Elements for Design Evaluation and Analysis of Structures
//
//----------------------------------------------------------------------------//
//
// Please cite the following resource in any derivative works:
//
// [1] Perez, C.M., and Filippou F.C.. "On Nonlinear Geometric Transformations
//     of Finite Elements" Int. J. Numer. Meth. Engrg. 2024; 
//     https://doi.org/10.1002/nme.7506
//
//===----------------------------------------------------------------------===//

#pragma once
#include <Node.h>
#include <Vector.h>
#include <Matrix.h>
#include <Versor.h>
#include <VectorND.h>
#include <Vector3D.h>
#include <MatrixND.h>
#include <Matrix3D.h>

#define TRIAD C2
namespace OpenSees {

class Isometry
{
public:
  virtual int       initialize() =0;
  virtual int       update() =0;

  virtual double    getLength() const =0;
  // x, \Lambda
  virtual Matrix3D  getRotation() const =0;
  virtual Vector3D  getPosition() =0;
  // \psi
  virtual Vector3D  getPositionVariation(int ndf, double* du) =0; 
  virtual Vector3D  getRotationVariation(int ndf, double* du) =0;
  virtual Matrix3D  getRotationDelta() =0;
  //
  virtual MatrixND<3,6> getRotationGradient(int node) =0;

};


template <int nn>
class AlignedIsometry : public Isometry
{
public:

  virtual Vector3D
  getRotationVariation(int ndf, double* du) {
    // psi_r = omega
    Vector3D w{};
    for (int i=0; i<nn; i++) {
      // const Vector &du = nodes[i]->getIncrDeltaDisp();
      auto Wi = this->getRotationGradient(i);
      for (int j=0; j<3; j++)
        for (int k=0; k<6; k++)
          w[j] += Wi(j,k) * du[ndf*i + k];
    }
    return w;
  }

protected:
};

} // namespace OpenSees
