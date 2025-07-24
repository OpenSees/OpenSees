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

#include <Vector.h>
#include <Matrix.h>
#include <Node.h>
#include <OPS_Globals.h>
#include <BasicFrameTransf.h>
#include "FrameTransform.h"
#include <utility/Unroll.h>

namespace OpenSees {

namespace {
  template<int ndf>
  const Vector&
  ShapeBasic(const VectorND<2*ndf>& ul)
  {
    static VectorND<6+(ndf-6)*2> ub;
    static Vector wrapper(ub);
    ub[0] = ul[1*ndf+0]; // Nj
    ub[1] = ul[0*ndf+5];
    ub[2] = ul[1*ndf+5];
    ub[3] = ul[0*ndf+4];
    ub[4] = ul[1*ndf+4];
    ub[5] = ul[1*ndf+3] - ul[0*ndf+3];
    return wrapper;
  }
}


template<int ndf>
BasicFrameTransf3d<ndf>::BasicFrameTransf3d(FrameTransform<2,ndf> *t)
: CrdTransf(t->getTag(), 0),
  t(*t),
  linear(t->getTag(), t->getNormalVector(), t->getRigidOffsets())
{

}

template<int ndf>
BasicFrameTransf3d<ndf>::~BasicFrameTransf3d()
{
  delete &t;
}

template<int ndf>
int
BasicFrameTransf3d<ndf>::commitState()
{
  return t.commit();
}

template<int ndf>
int
BasicFrameTransf3d<ndf>::revertToLastCommit()
{
  return t.revertToLastCommit();
}

template<int ndf>
int
BasicFrameTransf3d<ndf>::revertToStart()
{
  return t.revertToStart();
}

template<int ndf>
int
BasicFrameTransf3d<ndf>::update()
{
  return t.update();
}


template<int ndf>
int
BasicFrameTransf3d<ndf>::initialize(Node *i, Node *j)
{
  std::array<Node*, 2> nodes = {i, j};
  linear.initialize(nodes);
  return t.initialize(nodes);
}

template<int ndf>
int
BasicFrameTransf3d<ndf>::getLocalAxes(Vector &XAxis, Vector &YAxis, Vector &ZAxis)
{
  Vector3D x, y, z;
  int s = t.getLocalAxes(x, y, z);
  for (int i=0; i<3; i++) {
    XAxis(i) = x[i];
    YAxis(i) = y[i];
    ZAxis(i) = z[i];
  }
  return s;
}

template<int ndf>
double
BasicFrameTransf3d<ndf>::getInitialLength()
{
  return linear.getInitialLength();
}

template<int ndf>
double
BasicFrameTransf3d<ndf>::getDeformedLength()
{
  return t.getDeformedLength();
}


template<int ndf>
const Vector &
BasicFrameTransf3d<ndf>::getBasicTrialDisp()
{
  static VectorND<6> ub;
  static Vector wrapper(ub);
  Vector3D wi = t.getNodeRotationLogarithm(0),
           wj = t.getNodeRotationLogarithm(1);
  ub[0] = t.getDeformedLength() - t.getInitialLength(); // t.getNodePosition(1)[0]; //
  ub[1] = wi[2];
  ub[2] = wj[2];
  ub[3] = wi[1];
  ub[4] = wj[1];
  ub[5] = wj[0] - wi[0];
  return wrapper;
}


template<int ndf>
const Vector &
BasicFrameTransf3d<ndf>::getBasicIncrDeltaDisp()
{
  return ShapeBasic<ndf>(t.getStateVariation());
}

template<int ndf>
const Vector &
BasicFrameTransf3d<ndf>::getBasicIncrDisp()
{
  static VectorND<6> ub;
  static Vector wrapper(ub);
  opserr << "Unimplemented method\n";
  return wrapper;
}

template<int ndf>
const Vector &
BasicFrameTransf3d<ndf>::getBasicTrialVel()
{
  static VectorND<6> ub;
  static Vector wrapper(ub);
  opserr << "Unimplemented method\n";
  return wrapper;
}


template<int ndf>
const Vector &
BasicFrameTransf3d<ndf>::getGlobalResistingForce(const Vector &q_pres, const Vector &p0)
{
  // transform resisting forces from the basic system to local coordinates
  
  static constexpr int nwm = ndf - 6; // Number of warping DOFs

  static constexpr double c = 1.0;

  static VectorND<NDF*2> pl{};
  static Vector wrapper(pl); // to return reference
  pl.zero();
  pl[0*NDF+0]  = -q_pres[jnx] + p0[0]*c; // Ni
  pl[0*NDF+1]  =  p0[1]*c;               //
  pl[0*NDF+2]  =  p0[3]*c;               //
  pl[0*NDF+3]  = -q_pres[jmx];           // Ti
  pl[0*NDF+4]  =  q_pres[imy];
  pl[0*NDF+5]  =  q_pres[imz];
  pl[1*NDF+0]  =  q_pres[jnx];           // Nj
  pl[1*NDF+1]  =  p0[2]*c;
  pl[1*NDF+2]  =  p0[4]*c;
  pl[1*NDF+3]  =  q_pres[jmx];           // Tj
  pl[1*NDF+4]  =  q_pres[jmy];
  pl[1*NDF+5]  =  q_pres[jmz];

  if constexpr (nwm > 0) // Warping DOFs
    for (int i=0; i<nwm; i++) {
      // TODO
      pl[0*NDF+6+i] = -q_pres[6+i];
      pl[1*NDF+6+i] =  q_pres[6+i];
    }
  //

  t.push(pl, Operation::Total);

  return wrapper;
}


template<int ndf>
const Matrix &
BasicFrameTransf3d<ndf>::getGlobalStiffMatrix(const Matrix &kb, const Vector &q_pres)
{
  static constexpr int nwm = ndf - 6; // Number of warping DOFs

  static VectorND<NDF*2> pl{};
  pl.zero();
  pl[0*NDF+4]  =  q_pres[imy];
  pl[0*NDF+5]  =  q_pres[imz];
  pl[1*NDF+0]  =  q_pres[jnx];      // Nj
  pl[1*NDF+3]  =  q_pres[jmx];      // Tj
  pl[1*NDF+4]  =  q_pres[jmy];
  pl[1*NDF+5]  =  q_pres[jmz];
  //
  pl[0*NDF+0]  = -q_pres[jnx];      // Ni
  pl[0*NDF+3]  = -q_pres[jmx];      // Ti

  if constexpr (nwm > 0) // Warping DOFs
    for (int i=0; i<nwm; i++) {
      // TODO
      pl[0*NDF+6+i] = -q_pres[6+i];
      pl[1*NDF+6+i] =  q_pres[6+i];
    }

  //
  static MatrixND<2*NDF,2*NDF> kl;
  static Matrix Wrapper(kl);

  Repeat<NDF*2> ([&](auto j_) {
  constexpr static int j = j_.value;
    constexpr int jj = iq[j] >= 0 ? iq[j] : -iq[j];

    Repeat<NDF*2> ([&](auto i_) {
      constexpr static int i = i_.value;
      constexpr int ii = iq[i] >= 0 ? iq[i] : -iq[i];
      if constexpr (ii >= NBV || jj >= NBV) {
        kl(i,j) = 0.0;
        return;
      }

      kl(i,j) = kb(ii, jj);
    });
  });

  Repeat<NDF*2> ([&](auto i_) {
    constexpr static int i = i_.value;
    kl(0*NDF+0, i) = kl(i, 0*NDF+0) =  i==0? kl(NDF+0, NDF+0): (i==3? kl(NDF+0, NDF+3) : -kl( NDF+0, i));
    kl(0*NDF+3, i) = kl(i, 0*NDF+3) =  i==0? kl(NDF+3, NDF+0): (i==3? kl(NDF+3, NDF+3) : -kl( NDF+3, i));
  });

  t.push(kl, pl, Operation::Total);

  return Wrapper;
}


template<int ndf>
const Matrix &
BasicFrameTransf3d<ndf>::getInitialGlobalStiffMatrix(const Matrix &KB)
{
  static double kb[6][6];     // Basic stiffness

  static MatrixND<2*ndf,2*ndf> kl;  // Local stiffness

  double tmp[6][12]{};

  for (int i = 0; i < 6; i++)
    for (int j = 0; j < 6; j++)
      kb[i][j] = KB(i, j);

  // Transform basic stiffness to local system
  // First compute kb*T_{bl}
  for (int i = 0; i < 6; i++) {
    tmp[i][0]  = -kb[i][0];
    tmp[i][3]  = -kb[i][5];
    tmp[i][4]  =  kb[i][3];
    tmp[i][5]  =  kb[i][1];
    tmp[i][6]  =  kb[i][0];
    tmp[i][9]  =  kb[i][5];
    tmp[i][10] =  kb[i][4];
    tmp[i][11] =  kb[i][2];
  }

  kl.zero();
  // Now compute T'_{bl}*(kb*T_{bl})
  for (int i = 0; i < 12; i++) {
    kl( 0, i) = -tmp[0][i];
    kl( 3, i) = -tmp[5][i];
    kl( 4, i) = tmp[3][i];
    kl( 5, i) = tmp[1][i];

    kl( 6, i) = tmp[0][i];
    kl( 9, i) = tmp[5][i];
    kl(10, i) = tmp[4][i];
    kl(11, i) = tmp[2][i];
  }

  static MatrixND<ndf*2,ndf*2> kg;
  static Matrix M(kg);
  static constexpr VectorND<ndf*2> p0{};

  kg = linear.pushResponse(kl, p0);

  return M;
}


template<int ndf>
CrdTransf *
BasicFrameTransf3d<ndf>::getCopy3d()
{
  return new BasicFrameTransf3d(t.getCopy());
}


template<int ndf>
const Matrix &
BasicFrameTransf3d<ndf>::getGlobalMatrixFromLocal(const Matrix &M)
{
  //
  // Do diag(R)*M*diag(R)'
  //
  static MatrixND<ndf*2,ndf*2> Kout;
  static Matrix wrapper(Kout);
  static constexpr VectorND<ndf*2> p0{};
  wrapper = M;
  MatrixND<ndf*2,ndf*2> Kg = linear.pushResponse(Kout, p0);
  Kout = Kg;
  return wrapper;
}


template<int ndf>
const Vector &
BasicFrameTransf3d<ndf>::getPointGlobalCoordFromLocal(const Vector &xl)
{
  static Vector xg(3);
  return xg;
}


template<int ndf>
const Vector &
BasicFrameTransf3d<ndf>::getPointGlobalDisplFromBasic(double xi, const Vector &uxb)
{
  static Vector uxg(3);
  return uxg;
}

template<int ndf>
const Vector &
BasicFrameTransf3d<ndf>::getPointLocalDisplFromBasic(double xi, const Vector &uxb)
{
  static Vector uxl(3);
  return uxl;
}

//
// Sensitivity
//
template<int ndf>
bool
BasicFrameTransf3d<ndf>::isShapeSensitivity()
{
  return t.isShapeSensitivity();
}


template<int ndf>
double
BasicFrameTransf3d<ndf>::getdLdh()
{
  return t.getLengthGrad();
}

template<int ndf>
double
BasicFrameTransf3d<ndf>::getd1overLdh()
{
  double L = t.getInitialLength();
  return -getLengthGrad()/(L*L);
}

template<int ndf>
const Vector &
BasicFrameTransf3d<ndf>::getGlobalResistingForceShapeSensitivity(const Vector &q_pres,
                                                           const Vector &p0,
                                                           int gradNumber)
{
  // return t.getGlobalResistingForceShapeSensitivity(pb, p0, gradNumber);

  static constexpr int nwm = ndf - 6; // Number of warping DOFs

  static constexpr double c = 1.0;

  static VectorND<NDF*2> pl{};
  pl.zero();
  pl[0*NDF+0]  = -q_pres[jnx] + p0[0]*c; // Ni
  pl[0*NDF+1]  =  p0[1]*c;               //
  pl[0*NDF+2]  =  p0[3]*c;               //
  pl[0*NDF+3]  = -q_pres[jmx];           // Ti
  pl[0*NDF+4]  =  q_pres[imy];
  pl[0*NDF+5]  =  q_pres[imz];
  pl[1*NDF+0]  =  q_pres[jnx];           // Nj
  pl[1*NDF+1]  =  p0[2]*c;
  pl[1*NDF+2]  =  p0[4]*c;
  pl[1*NDF+3]  =  q_pres[jmx];           // Tj
  pl[1*NDF+4]  =  q_pres[jmy];
  pl[1*NDF+5]  =  q_pres[jmz];

  if constexpr (nwm > 0) // Warping DOFs
    for (int i=0; i<nwm; i++) {
      // TODO
      pl[0*NDF+6+i] = -q_pres[6+i];
      pl[1*NDF+6+i] =  q_pres[6+i];
    }
  //
  
  static VectorND<2*ndf> dp;
  static Vector wrapper(dp); // to return reference
  dp.zero();

  t.pushGrad(dp, pl);

  return wrapper;
}


template<int ndf>
const Vector &
BasicFrameTransf3d<ndf>::getBasicTrialDispShapeSensitivity()
{
  static VectorND<2*ndf> du;
  t.pullFixedGrad(du);
  return ShapeBasic<ndf>(du);
}

template<int ndf>
const Vector &
BasicFrameTransf3d<ndf>::getBasicDisplSensitivity(int gradNumber)
{
  static VectorND<2*ndf> du;
  t.pullTotalGrad(du, gradNumber);
  return ShapeBasic<ndf>(du);
}


template<int ndf>
void
BasicFrameTransf3d<ndf>::Print(OPS_Stream &s, int flag)
{
  t.Print(s, flag);
}


template<int ndf>
int
BasicFrameTransf3d<ndf>::sendSelf(int cTag, Channel &theChannel)
{
  return -1;
}


template<int ndf>
int
BasicFrameTransf3d<ndf>::recvSelf(int cTag, Channel &,
                            FEM_ObjectBroker &theBroker)
{
  return -1;
}

}
