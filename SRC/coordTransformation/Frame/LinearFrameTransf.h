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

//
// Description: This file contains the class definition for
// LinearFrameTransf.h. LinearFrameTransf provides the
// abstraction of a linear transformation for a spatial frame
// between the global and basic coordinate systems
//
// Adapted: Remo Magalhaes de Souza (rmsouza@ce.berkeley.edu)
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

  VectorND<nn*ndf> getStateVariation() final;
  Vector3D getNodePosition(int tag) final;
  Vector3D getNodeRotationLogarithm(int tag) final;

  VectorND<nn*ndf>        pushResponse(VectorND<nn*ndf>&pl) final;
  MatrixND<nn*ndf,nn*ndf> pushResponse(MatrixND<nn*ndf,nn*ndf>& kl, const VectorND<nn*ndf>& pl) final;

  // // method used to rotate consistent mass matrix
  // const Matrix &getGlobalMatrixFromLocal(const Matrix &local);
  

  // Sensitivity
  //
  const Vector & getBasicDisplFixedGrad();
  const Vector & getBasicDisplTotalGrad(int gradNumber);
  const Vector &getGlobalResistingForceShapeSensitivity (const Vector &basicForce, const Vector &p0, int grad);
  bool isShapeSensitivity() final;
  double getLengthGrad() final;
  double getd1overLdh() final;

  // TaggedObject
  void Print(OPS_Stream &s, int flag = 0) final;

  // Personal
  Vector3D getDelta() {return Du;}
          
private:

  inline VectorND<nn*ndf> 
  pullConstant(const VectorND<nn*ndf>& ug, 
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
  Vector3D Du;

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

