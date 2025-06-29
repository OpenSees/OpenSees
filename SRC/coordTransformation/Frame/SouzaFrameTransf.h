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
// SouzaFrameTransf. SouzaFrameTransf implements the formulation 
// of Crisfield (1990) with the objective of maintaining the
// original "Corotational" implementation by Remo Magalhaes de Souza, within
// the new framework proposed by Perez and Filippou (2024).
//
// Written by : cmp, March 2024
//
// Adapted from work by: Remo Magalhaes de Souza
//
#ifndef SouzaFrameTransf_hpp
#define SouzaFrameTransf_hpp

#include <array>
#include "FrameTransform.h"
#include <Vector.h>
#include <Matrix.h>
#include <Versor.h>
#include <Matrix3D.h>
#include <Vector3D.h>
#include "Isometry/CrisfieldTransform.h"

struct Triad;

namespace OpenSees {

template <int nn, int ndf>
class SouzaFrameTransf: public FrameTransform<nn,ndf>
{
public:
  SouzaFrameTransf(int tag, const Vector3D &vecxz,
                    const std::array<Vector3D, nn> *offset=nullptr,
                    int offset_flags = 0);

  ~SouzaFrameTransf();

  const char *getClassType() const {
    return "SouzaFrameTransf";
  }

  // NOTE: maybe add arg for rotation parameterization
  FrameTransform<nn,ndf> *getCopy() const final;

  int initialize(std::array<Node*, nn>& new_nodes) final;
  int update() final;
  int commit() final;
  int revertToLastCommit() final;        
  int revertToStart() final;
  int getLocalAxes(Vector3D &x, Vector3D &y, Vector3D &z) const final;
  const std::array<Vector3D,nn> *getRigidOffsets() const final { return offsets; }

  double getInitialLength();
  double getDeformedLength();

  VectorND<nn*ndf> getStateVariation() final;
  Vector3D getNodePosition(int tag) final;
  Vector3D getNodeRotationLogarithm(int tag) final;

  VectorND<nn*ndf>    pushResponse(VectorND<nn*ndf>&pl) final;
  MatrixND<nn*ndf,nn*ndf> pushResponse(MatrixND<nn*ndf,nn*ndf>& kl, const VectorND<nn*ndf>& pl) final;

  // Sensitivity
  double getLengthGrad() final;
  virtual const Vector &getBasicDisplTotalGrad(int grad);
  virtual const Vector &getBasicDisplFixedGrad();
  virtual const Vector &getGlobalResistingForceShapeSensitivity(const Vector &pb, const Vector &p0, int gradNumber);

  // Tagged Object
  void Print(OPS_Stream &s, int flag = 0) final;

private:
  constexpr static int n = nn*ndf;

  // compute the transformation matrix
  void compTransfMatrixBasicGlobal(const Versor&, const Versor*);

  int addTangent(MatrixND<12,12>& M, const VectorND<12>& pl);

  VectorND<6> pushResponse(const VectorND<6>& pa, int a, int b);

  enum {
    inx= 0, // axial
    iny= 1, // Vy
    inz= 2, // Vz
    imx= 3, // torsion
    imy= 4, // rot y I
    imz= 5, // rot z I

    jnx= 6, // axial
    jny= 7,
    jnz= 8,
    jmx= 9, // torsion
    jmy=10, // rot y J
    jmz=11, // rot z J
  };

  //
  // Member data
  //
  std::array<Node*, nn> nodes;

  Vector3D xAxis;                              // local x axis
  Vector3D vz;                                 // Vector that lies in local plane xz
  Vector3D dX;
  
  std::array<Vector3D, nn> *offsets;

  double *nodeIInitialDisp, *nodeJInitialDisp;
  bool  initialDispChecked;                    

  double L;                       // initial element length
  double Ln;                      // current element length (at trial state)

  Versor Q_past[nn];              // commited rotations
  Versor Q_pres[nn];              // trial rotations

  Vector3D alphaI;                // last trial rotations end i
  Vector3D alphaJ;                // last trial rotatations end j

  VectorND<n> ul;                 // local displacements (size=7)
  Vector3D    vr[nn];             //
  VectorND<n> ulcommit;           // commited local displacements
  VectorND<n> ulpr;               // previous local displacements

  OpenSees::MatrixND<n,n> T;     // transformation from local to global system

  OpenSees::Matrix3D R0;         // rotation from local to global coordinates
  CrisfieldTransform crs;

  // Static workspace variables
  Matrix3D A;
  MatrixND<12,3> Lr2, Lr3;   // auxiliary matrices
};
} // namespace OpenSees
#include "SouzaFrameTransf.tpp"
#endif
