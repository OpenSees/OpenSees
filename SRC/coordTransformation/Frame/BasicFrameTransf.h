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
// The purpose of this class is to wrap the more general FrameTransform<>
// templates to reproduce the legacy CrdTransf classes that were derived
// for elements in a "basic" coordinate system.
//
// cmp
//
#ifndef BasicFrameTransf3d_h
#define BasicFrameTransf3d_h

#include <array>
#include <CrdTransf.h>
#include <FrameTransform.h>
#include <Vector.h>
#include <Matrix.h>

namespace OpenSees {

template<int ndf=6>
class BasicFrameTransf3d: public CrdTransf
{
public:
  BasicFrameTransf3d(FrameTransform<2,ndf> *t);

  ~BasicFrameTransf3d();

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
  const Vector & getBasicDisplFixedGrad() final;
  const Vector & getBasicDisplTotalGrad(int grad) final;
  const Vector &getGlobalResistingForceShapeSensitivity(const Vector &basicForce, const Vector &p0, int grad) final;
  bool isShapeSensitivity() final;
  double getLengthGrad() final;
  double getd1overLdh() final;


  // MovableObject
  int sendSelf(int tag, Channel &) final;
  int recvSelf(int tag, Channel &, FEM_ObjectBroker &) final;
  const char *getClassType() const final {
      return "BasicFrameTransf3d";
  }
    
  // TaggedObject
  void Print(OPS_Stream &s, int flag = 0) final;


  FrameTransform<2,ndf> &t;

private:
  constexpr static int NBV = 6;
  constexpr static int NDF = ndf;
  enum : int {
      inx = -12, //  0
      iny = -12, //  1
      inz = -12, //  2
      imx = -12, //  3
      imy =   3, //  4
      imz =   1, //  5
      jnx =   0, //  6
      jny = -12, //  7
      jnz = -12, //  8
      jmx =   5, //  9
      jmy =   4, // 10
      jmz =   2, // 11
  };

  constexpr static int iq[] = {
      inx, iny, inz, imx, imy, imz,
      jnx, jny, jnz, jmx, jmy, jmz
  };

};
} // namespace OpenSees

#include "BasicFrameTransf.tpp"
#endif
