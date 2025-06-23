//===----------------------------------------------------------------------===//
//
//                                   xara
//
//===----------------------------------------------------------------------===//
//                              https://xara.so
//===----------------------------------------------------------------------===//
//
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

class FrameBasis
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
class RankineBasis : public FrameBasis
{
public:
  RankineBasis(std::array<Node*,nn>& nodes, const Vector3D& vecxz)
  : nodes{nodes}, vz(vecxz), Xc{}, c{}, R{} {
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

#if 1
    Xc = nodes[ic]->getCrds();
    c[init] = R[init]^(Xc);
#endif

    update();
    return 0;
  }

  virtual int
  update() {

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
#if 1 // TRIAD==R2
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

#elif 1 // TRIAD==C2
      Versor q0 = VersorFromMatrix(R[init]);
      Versor qI = nodes[0]->getTrialRotation()*q0;
      Versor qJ = nodes[nn-1]->getTrialRotation()*q0;
      Vector3D gammaw = CayleyFromVersor(qJ.mult_conj(qI));

      gammaw *= 0.5;

      //  Qbar = VersorProduct(VersorFromMatrix(CaySO3(gammaw)), qI);
      Matrix3D Rbar = CaySO3(gammaw)*MatrixFromVersor(qI); // *q0);
      Vector3D v { Rbar(0,0), Rbar(1,0), Rbar(2,0) };
      double dot = v.dot(e1);
      if (std::fabs(std::fabs(dot)-1.0) < 1.0e-10) {
        R[pres] = Rbar;
      } else {
        v  = v.cross(e1);
        double scale = std::acos(dot)/v.norm();
        v *= scale; // ::acos(r1.dot(e1));

        R[pres] = ExpSO3(v)*Rbar;

        Vector3D r1 { R[pres](0,0), R[pres](1,0), R[pres](2,0) };
        Vector3D r2 { R[pres](0,1), R[pres](1,1), R[pres](2,1) };
        Vector3D r3 { R[pres](0,2), R[pres](1,2), R[pres](2,2) };
        // opserr << Vector(r1-e1);
        // R[pres] = Rbar^ExpSO3(v);
      }
#else 
      Vector3D e2 = vz.cross(e1);
#endif

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

  virtual double 
  getLength() const override {
    return Ln;
  }


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
  double n   = 0,
         n11 = 1,
         n12 = 0,
         n21 = 0,
         n22 = 1;
  Vector3D vz, dX, Xc;
  Matrix3D R[3];
  Vector3D c[3];
  Matrix3D dR;
  std::array<Node*,nn>& nodes;
  std::array<Vector3D, nn>* offsets = nullptr; // offsets
};

} // namespace OpenSees
