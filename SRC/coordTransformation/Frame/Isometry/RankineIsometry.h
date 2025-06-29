//===----------------------------------------------------------------------===//
//
//                                   xara
//                              https://xara.so
//----------------------------------------------------------------------------//
//
//                                 FEDEASLab
//       Finite Elements for Design Evaluation and Analysis of Structures
//                            https://fedeas.net
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
#include <array>
#include <Node.h>
#include <Vector3D.h>
#include <Matrix3D.h>
#include <MatrixND.h>
#include "EuclidIsometry.h"


namespace OpenSees {

template <int nn>
class RankineIsometry : public AlignedIsometry<nn>
{
public:
  RankineIsometry(std::array<Node*,nn>& nodes, const Vector3D& vecxz)
  : nodes{nodes}, vz(vecxz), Xc{}, c{}, R{}
  {
  }

  void
  setOffsets(std::array<Vector3D, nn>* offsets) {
    this->offsets = offsets;
  }

  virtual int 
  initialize() {

    for (int i=0; i<nn; i++)
      nodes[i]->getTrialRotation();

    const Vector &XI = nodes[   0]->getCrds();
    const Vector &XJ = nodes[nn-1]->getCrds();

    for (int i=0; i<3; i++)
      dX[i] = XJ[i] - XI[i];

    L = dX.norm();
    Ln = L;
    Vector3D e1 = dX/L;

    //
    Vector3D e2 = vz.cross(e1);

    const double ynorm = e2.norm();

    if (ynorm == 0.0)
      return -1;

    e2 /= ynorm;

    Vector3D e3 = e1.cross(e2);

    e2 = e3.cross(e1);

    for (int i = 0; i < 3; i++) {
      R[init](i,0) = e1[i];
      R[init](i,1) = e2[i];
      R[init](i,2) = e3[i];
    }

    Xc = nodes[ic]->getCrds();
    c[init] = R[init]^(Xc);

    this->update();
    return 0;
  }

  int
  update() final {

    Vector3D e1 = dX;
    {
      //
      // Update state
      //
      {
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
    }

    Vector3D uc = nodes[ic]->getTrialDisp();
    if (offsets != nullptr) {
      uc.addVector(1.0, (*offsets)[ic], -1.0);
      uc.addVector(1.0, nodes[ic]->getTrialRotation().rotate((*offsets)[ic]), 1.0);
    }
    Vector3D X = nodes[ic]->getCrds();
    c[pres] = R[pres]^(X + uc);
    return 0;
  };

  virtual MatrixND<3,6> 
  getRotationGradient(int node) {
    MatrixND<3,6> Gb{};

    constexpr Vector3D axis{1, 0, 0};
    constexpr Matrix3D ix = Hat(axis);
    constexpr Matrix3D ioi = axis.bun(axis);

    Gb.template insert<0, 3>(ioi, 0.5);
    if (node == 0)
      Gb.template insert<0,0>(ix, -1/Ln);
    else if (node == nn-1)
      Gb.template insert<0,0>(ix,  1/Ln);
    
    Gb(0,2) = (node == 0? 1.0 : -1.0)*n;
    return Gb;
  }

  virtual double
  getLength() const override {
    return Ln;
  }

  virtual Matrix3D
  getRotation() const override {
    return R[pres];
  }

  virtual Matrix3D 
  getRotationDelta() {
    return R[pres] - R[init];
  }

  Vector3D
  getLocation() {
    return c[pres];
  }

  virtual Vector3D
  getPosition() {
    // Return Delta c
    Vector3D X = nodes[ic]->getCrds();
    Vector3D Dc =  c[pres] - (R[init]^X) ; // (R[pres]^c[init]);
    return Dc;
  }

  virtual Vector3D 
  getPositionVariation(int ndf, double* du) {
    return Vector3D {du[ndf*ic+0], du[ndf*ic+1], du[ndf*ic+2]};
  }

private:
  constexpr static int ic = 0; // std::floor(0.5*(nn+1));
  enum { pres, prev, init};
  double L, Ln;
  Vector3D vz, dX, Xc;
  Matrix3D R[3];
  Vector3D c[3];
  Matrix3D dR;
  std::array<Node*,nn>& nodes;
  std::array<Vector3D, nn>* offsets = nullptr; // offsets

  double n   = 0,
         n11 = 1,
         n12 = 0,
         n21 = 0,
         n22 = 1;
};
} // namespace OpenSees