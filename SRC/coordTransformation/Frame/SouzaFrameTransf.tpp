//===----------------------------------------------------------------------===//
//
//                                   xara
//
//===----------------------------------------------------------------------===//
//                              https://xara.so
//===----------------------------------------------------------------------===//
//
// Description: This file contains the implementation for the
// SouzaFrameTransf class. SouzaFrameTransf is a Corotational
// transformation for a spatial frame element between the global
// and basic coordinate systems.
//
// Written: Claudio Perez
// Created: 05/2024
//
// Adapted from: Remo Magalhaes de Souza (rmsouza@ce.berkeley.edu)
//
#include <math.h>
#include <Node.h>
#include <Channel.h>
#include <Logging.h>
#include <SouzaFrameTransf.hpp>

#include <Triad.h>
#include <Vector.h>
#include <Matrix.h>
#include <VectorND.h>
#include <Vector3D.h>
#include <MatrixND.h>
#include <Matrix3D.h>
#include <Rotations.hpp>
#include "Orient/CrisfieldTransform.h"

using namespace OpenSees;

template <int nn, int ndf>
SouzaFrameTransf<nn,ndf>::SouzaFrameTransf(int tag, const Vector3D &vz,
                                            const std::array<Vector3D, nn> *offset,
                                            int offset_flags)
  : FrameTransform<nn,ndf>(tag),
    vz(vz),
    offsets{nullptr},
    L(0), Ln(0),
    nodeIInitialDisp(0), nodeJInitialDisp(0),
    initialDispChecked(false)
{
  alphaI.zero();
  alphaJ.zero();

  // Rigid joint offsets
  if (offset != nullptr) {
    offsets = new std::array<Vector3D, nn>{};
    *offsets = *offset;
  }
}


template <int nn, int ndf>
SouzaFrameTransf<nn,ndf>::~SouzaFrameTransf()
{
  if (offsets != nullptr)
    delete offsets;
}


template <int nn, int ndf>
double
SouzaFrameTransf<nn,ndf>::getInitialLength()
{
  return L;
}


template <int nn, int ndf>
double
SouzaFrameTransf<nn,ndf>::getDeformedLength()
{
  return Ln;
}


template <int nn, int ndf>
FrameTransform<nn,ndf> *
SouzaFrameTransf<nn,ndf>::getCopy() const
{

  SouzaFrameTransf *theCopy =
    new SouzaFrameTransf<nn,ndf>(this->getTag(), vz, offsets);

  theCopy->nodes[0]  = nodes[0];
  theCopy->nodes[1]  = nodes[1];
  theCopy->xAxis     = xAxis;
  theCopy->L         = L;
  theCopy->Ln        = Ln;
  theCopy->R0        = R0;
  theCopy->Q_pres[0] = Q_pres[0];
  theCopy->Q_pres[1] = Q_pres[1];
  theCopy->Q_past[0] = Q_past[0];
  theCopy->Q_past[1] = Q_past[1];
  theCopy->ul = ul;
  theCopy->ulcommit = ulcommit;
  return theCopy;
}


template <int nn, int ndf>
int
SouzaFrameTransf<nn,ndf>::revertToStart()
{
  ul.zero();
  Q_pres[0] = VersorFromMatrix(R0);
  for (int i=1; i<nn; i++)
    Q_pres[i] = Q_pres[0];

  alphaI.zero();
  alphaJ.zero();

  this->update();
  return 0;
}


template <int nn, int ndf>
int
SouzaFrameTransf<nn,ndf>::commit()
{
  ulcommit  = ul;
  Q_past[0] = Q_pres[0];
  Q_past[1] = Q_pres[1];
  return 0;
}


template <int nn, int ndf>
int
SouzaFrameTransf<nn,ndf>::revertToLastCommit()
{
  // determine global displacement increments from last iteration
  const Vector &dispI = nodes[0]->getTrialDisp();
  const Vector &dispJ = nodes[1]->getTrialDisp();

  for (int k = 0; k < 3; k++) {
    alphaI(k) =  dispI(k+3);
    alphaJ(k) =  dispJ(k+3);
  }

  ul        = ulcommit;
  Q_pres[0] = Q_past[0];
  Q_pres[1] = Q_past[1];

  this->update();

  return 0;
}


template <int nn, int ndf>
int
SouzaFrameTransf<nn,ndf>::initialize(std::array<Node*, nn>& new_nodes)
{
  for (int i=0; i<nn; i++) {
    nodes[i] = new_nodes[i];
    if (nodes[i] == nullptr) {
      opserr << "invalid pointers to the element nodes\n";
      return -1;
    }
  }

  dX  = nodes[nn-1]->getCrds() - nodes[0]->getCrds();

  // Add initial displacements at nodes
  if (initialDispChecked == false) {
    const Vector &nodeIDisp = nodes[0]->getDisp();
    const Vector &nodeJDisp = nodes[1]->getDisp();
    for (int i = 0; i<6; i++)
      if (nodeIDisp[i] != 0.0) {
        nodeIInitialDisp = new double [6];
        for (int j = 0; j<6; j++)
          nodeIInitialDisp[j] = nodeIDisp[j];
        i = 6;
      }

    for (int j = 0; j<6; j++)
      if (nodeJDisp[j] != 0.0) {
        nodeJInitialDisp = new double [6];
        for (int i = 0; i<6; i++)
          nodeJInitialDisp[i] = nodeJDisp[i];
        j = 6;
      }
    initialDispChecked = true;
  }

  // if (nodeIInitialDisp != nullptr) {
  //   dX[0] -= nodeIInitialDisp[0];
  //   dX[1] -= nodeIInitialDisp[1];
  //   dX[2] -= nodeIInitialDisp[2];
  // }

  // if (nodeJInitialDisp != nullptr) {
  //   dX[0] += nodeJInitialDisp[0];
  //   dX[1] += nodeJInitialDisp[1];
  //   dX[2] += nodeJInitialDisp[2];
  // }


  //
  // Length and Orientation
  //
  Vector3D dx;

  dx = nodes[nn-1]->getCrds() - nodes[0]->getCrds();

  L = dx.norm();

  if (L == 0.0) {
      opserr << "\nSouzaFrameTransf::computeElemtLengthAndOrien: 0 length\n";
      return -2;
  }

  //
  // Set rotation matrix
  //
  int error = FrameTransform<nn,ndf>::Orient(dx, vz, R0);
  if (error)
    return error;

  // Compute initial pseudo-vectors for nodal triads
  Q_pres[0] = Q_pres[1] = VersorFromMatrix(R0);

  ul.zero();
  ulpr.zero();

  for (int i=0; i<nn; i++)
    vr[i].zero();

  this->commit();

  return 0;
}

template <int nn, int ndf>
VectorND<nn*ndf>
SouzaFrameTransf<nn,ndf>::getStateVariation()
{
  return ul - ulpr;
}

template <int nn, int ndf>
Vector3D
SouzaFrameTransf<nn,ndf>::getNodePosition(int tag)
{
  Vector3D u;
  for (int i=0; i<3; i++)
    u[i] = ul[tag*ndf+i];
  return u;
}

template <int nn, int ndf>
Vector3D
SouzaFrameTransf<nn,ndf>::getNodeRotationLogarithm(int tag)
{
  return vr[tag];
}

template <int nn, int ndf>
void inline
SouzaFrameTransf<nn,ndf>::compTransfMatrixBasicGlobal(
                                                const Versor& Qbar,
                                                const Versor* Q_pres)
{
  // extract columns of rotation matrices
  const Triad r {MatrixFromVersor(Qbar)},
              rI{MatrixFromVersor(Q_pres[0])},
              rJ{MatrixFromVersor(Q_pres[1])};
  const Vector3D 
    &e1  =  crs.getBasisE1(), // E[1],
    &e2  =  crs.getBasisE2(), // E[2],
    &e3  =  crs.getBasisE3(), // E[3],
    &r1  =  r[1],
    &r2  =  r[2],
    &r3  =  r[3],
    &rI1 = rI[1],
    &rI2 = rI[2],
    &rI3 = rI[3],
    &rJ1 = rJ[1],
    &rJ2 = rJ[2],
    &rJ3 = rJ[3];

  // Compute the transformation matrix from the basic to the
  // global system
  //   A = (1/Ln)*(I - e1*e1');
  // Matrix3D A;
  for (int i = 0; i < 3; i++)
    for (int j = 0; j < 3; j++)
      A(i,j) = (double(i==j) - e1[i]*e1[j])/Ln;

  // This must be called up here
  CrisfieldTransform::getLMatrix(A, e1, r1, r2, Lr2);
  CrisfieldTransform::getLMatrix(A, e1, r1, r3, Lr3);

  //               3 |             3            |     3    |           3              |
  //   T1 = [      O', (-S(rI3)*e2 + S(rI2)*e3)',        O',                        O']'; imx
  //   T2 = [(A*rI2)', (-S(rI2)*e1 + S(rI1)*e2)', -(A*rI2)',                        O']'; imz
  //   T3 = [(A*rI3)', (-S(rI3)*e1 + S(rI1)*e3)', -(A*rI3)',                        O']'; imy
  //
  //   T4 = [      O',                        O',        O', ( S(rJ2)*e3 - S(rJ3)*e2 )']'; jmx
  //   T5 = [(A*rJ2)',                        O', -(A*rJ2)', ( S(rJ1)*e2 - S(rJ2)*e1 )']'; jmz
  //   T6 = [(A*rJ3)',                        O', -(A*rJ3)', ( S(rJ1)*e3 - S(rJ3)*e1 )']'; jmy

  T.zero();

  //   T1 = [      O', (-S(rI3)*e2 + S(rI2)*e3)',        O', O']';

  // (-S(rI3)*e2 + S(rI2)*e3)
  Vector3D Se  = rI2.cross(e3);
  Se -= rI3.cross(e2);
  for (int i = 0; i < 3; i++)
    // T(jmx,i+3) =  -Se[i];
    T(imx,i+3) =  Se[i];

  //   T2 = [(A*rI2)', (-S(rI2)*e1 + S(rI1)*e2)', -(A*rI2)', O']';

  Vector3D At = A*rI2;

  // (-S(rI2)*e1 + S(rI1)*e2)'
  Se  = rI1.cross(e2);
  Se -= rI2.cross(e1);
  for (int i = 0; i < 3; i++) {
      T(imz,i  ) =  At[i];
      T(imz,i+3) =  Se[i];
      T(imz,i+6) = -At[i];
  }

  //   T3 = [(A*rI3)', (-S(rI3)*e1 + S(rI1)*e3)', -(A*rI3)', O']';

  At = A*rI3;
  
  // -S(rI3)*e1 + S(rI1)*e3
  Se  = rI1.cross(e3);
  Se -= rI3.cross(e1);
  for (int i = 0; i < 3; i++) {
      T(imy,i  ) =  At[i]*-1;
      T(imy,i+3) =  Se[i]*-1;
      T(imy,i+6) = -At[i]*-1;
  }

  //   T4 = [      O', O',        O', (-S(rJ3)*e2 + S(rJ2)*e3)']';
  Se  = rJ2.cross(e3);
  Se -= rJ3.cross(e2);
  for (int i = 0; i < 3; i++)
    T(jmx, i+9) =  Se[i];   // S(rJ2)*e3 - S(rJ3)*e2

  // T5 = [(A*rJ2)', O', -(A*rJ2)', (-S(rJ2)*e1 + S(rJ1)*e2)']';
  At = A*rJ2;
  Se  = rJ1.cross(e2); 
  Se -= rJ2.cross(e1);
  for (int i = 0; i < 3; i++) {
      T(jmz, i  ) =  At[i];
      T(jmz, i+6) = -At[i];
      T(jmz, i+9) =  Se[i]; // (-S(rJ2)*e1 + S(rJ1)*e2)
  }

  // T6 = [(A*rJ3)', O', -(A*rJ3)', (-S(rJ3)*e1 + S(rJ1)*e3)']'
  At  = A*rJ3;
  Se  = rJ1.cross(e3);  // (-S(rJ3)*e1 + S(rJ1)*e3)
  Se -= rJ3.cross(e1);
  for (int i = 0; i < 3; i++) {
      T(jmy,i  ) =  At[i]*-1;
      T(jmy,i+6) = -At[i]*-1;
      T(jmy,i+9) =  Se[i]*-1;
  }

  //
  // Second part
  //

  // T(:,1) += Lr3*rI2 - Lr2*rI3;
  // T(:,2) +=           Lr2*rI1; z
  // T(:,3) += Lr3*rI1          ; y

  // T(:,4) += Lr3*rJ2 - Lr2*rJ3;
  // T(:,5) += Lr2*rJ1          ; z    // ?????? check sign
  // T(:,6) += Lr3*rJ1          ; y    // ?????? check sign

  // Bending Z
  for (int i = 0; i < 12; i++) {
    double T1i = 0;
    for (int k=0; k<3; k++)
      T1i += Lr2(i,k)*rI1[k];
    T(imz,i) += T1i;
  }

  for (int i = 0; i < 12; i++) {
    double T4i = 0;
    for (int k=0; k<3; k++)
      T4i += Lr2(i,k)*rJ1[k]; // Lr[i];
    T(jmz,i) += T4i;
  }

  // Torsion
  for (int i = 0; i < 12; i++) {
    double T0i = 0;
    for (int k=0; k<3; k++)
      T0i += Lr3(i,k)*rI2[k] - Lr2(i,k)*rI3[k];
    // T(jmx,i) += -T0i;
    T(imx,i) += T0i;
  }
  for (int i = 0; i < 12; i++) {
    double T3i = 0;
    for (int k=0; k<3; k++)
      T3i += Lr3(i,k)*rJ2[k] - Lr2(i,k)*rJ3[k];
    T(jmx,i) += T3i;
  }
  // Bending Y
  for (int i = 0; i < 12; i++) {
    double T2i = 0;
    for (int k=0; k<3; k++)
      T2i += Lr3(i,k)*rI1[k]; // Lr[i];
    T(imy,i) += T2i*-1; // TODO: Check
  }
  for (int i = 0; i < 12; i++) {
    double T5i = 0;
    for (int k=0; k<3; k++)
      T5i += Lr3(i,k)*rJ1[k]; // Lr[i];
    T(jmy,i) += T5i*-1; // TODO: Check
  }

  //
  //
  //
  for (int node=0; node < 2; node++)
    for (int j = 0; j < 3; j++) {
      const double c = 0.5 / std::cos(ul[(node? jmx : imx) + j]);
      for (int i = 0; i < 12; i++)
        T((node? jmx : imx) + j, i) *= c;
    }

  // Axial
  // T(:,7) = [-e1' O' e1' O']';
  for (int i = 0; i < 3; i++) {
      T(jnx,i  ) = -e1[i];
      T(jnx,i+6) =  e1[i];
  }

  // Combine torsion
  for (int i=0; i<12; i++) {
    T(jmx,i) -= T(imx,i);
    T(imx,i) = 0;
  }
}

//
// Set RI,RJ,Rbar, Ln, e and ul
//
template <int nn, int ndf>
int
SouzaFrameTransf<nn,ndf>::update()
{
  // determine global displacement increments from last iteration

  const Vector& dispI = nodes[   0]->getTrialDisp();
  const Vector& dispJ = nodes[nn-1]->getTrialDisp();

  //
  // Update state
  //
  // 1.1 Relative translation
  Vector3D dx = dX;// dx = dX + dJI;
  {
    // Relative translation
    for (int k = 0; k < 3; k++)
      dx[k] += dispJ(k) - dispI(k);


    // Calculate the deformed length
    Ln = dx.norm();

    if (Ln == 0.0) {
      opserr << "\nSouzaFrameTransf: deformed length is 0.0\n";
      return -2;
    }
  }

  // 1.2 Rotational displacement increments
  {
    Vector3D dAlphaI, dAlphaJ;

    for (int k = 0; k < 3; k++) {
      dAlphaI[k] =  dispI(k+3) - alphaI[k];
      alphaI[k]  =  dispI(k+3);
      dAlphaJ[k] =  dispJ(k+3) - alphaJ[k];
      alphaJ[k]  =  dispJ(k+3);
    }

    // Update the nodal rotations
    Q_pres[0] = VersorProduct(Q_pres[0],  Versor::from_vector(dAlphaI));
    Q_pres[1] = VersorProduct(Q_pres[1],  Versor::from_vector(dAlphaJ));
  }

  //
  // 2) Form transformation
  //

  crs.update(Q_pres[0], Q_pres[1], dx);
  // Form the transformation tangent
  this->compTransfMatrixBasicGlobal(crs.getReference(), Q_pres);

  //
  // 3) Local deformations
  //

  // Save previous state
  ulpr = ul;

  // Rotations
  {
    Matrix3D e = crs.getRotation();
    vr[0] = LogC90(e^MatrixFromVersor(Q_pres[0]));
    for (int i=0; i<3; i++)
      ul[imx+i] = vr[0][i];

    vr[1] = LogC90(e^MatrixFromVersor(Q_pres[1]));
    for (int i=0; i<3; i++)
      ul[jmx+i] = vr[1][i];
  }

  // Axial
  ul(inx) = 0;
  ul(jnx) = Ln - L;

  return 0;
}


template <int nn, int ndf>
inline VectorND<nn*ndf>
SouzaFrameTransf<nn,ndf>::pushResponse(VectorND<nn*ndf>&pl)
{
  // return T^pl;
  VectorND<nn*ndf> pg{};
  for (int a = 0; a<nn; a++) {
    VectorND<6> pa {pl(a*ndf+0), pl(a*ndf+1), pl(a*ndf+2), 
                    pl(a*ndf+3), pl(a*ndf+4), pl(a*ndf+5)};

    for (int b = 0; b<2; b++) {
      VectorND<6> pab = pushResponse(pa, a, b);
      pg.assemble(b*6, pab, 1.0);
    }
  }
  return pg;
}

template <int nn, int ndf>
VectorND<6>
SouzaFrameTransf<nn,ndf>::pushResponse(const VectorND<6>&pa, int a, int b)
{
  VectorND<6> pg{};
  for (int i = 0; i < 6; i++)
    for (int j = 0; j < 6; j++)
      pg[j] += T(a*6 + i, b*6+j) * pa[i];

  return pg;
}

// do 
//    K = ag'*Km*ag + Kp
//
template <int nn, int ndf>
MatrixND<nn*ndf,nn*ndf>
SouzaFrameTransf<nn,ndf>::pushResponse(MatrixND<nn*ndf,nn*ndf>& kl, const VectorND<nn*ndf>& pl)
{    
  MatrixND<12,12> K;
  K.addMatrixTripleProduct(0.0, T, kl, 1.0);

  // Add geometric part kg
  this->addTangent(K, pl);

  return K;
}


//
// Add geometric part of the transformation tangent
//
//  kg += T'*kl*T + ks1 + T * diag(m.*tan(thetal))*T' + ...
//         m(4)*(ks2r2t3_u3 + ks2r3u2_t2) + ...
//         m(2)*ks2r2t1 + m(3)*ks2r3t1 + ...
//         m(5)*ks2r2u1 + m(6)*ks2r3u1 + ...
//         ks3 + ks3' + ks4 + ks5;
template <int nn, int ndf>
int
SouzaFrameTransf<nn,ndf>::addTangent(MatrixND<12,12>& kg, const VectorND<12>& pl)
{    
    const Triad r {MatrixFromVersor(crs.getReference())},
                rI{MatrixFromVersor(Q_pres[0])},
                rJ{MatrixFromVersor(Q_pres[1])};
    const Vector3D 
      &e1  =  crs.getBasisE1(), // E[1],
      &e2  =  crs.getBasisE2(), // E[2],
      &e3  =  crs.getBasisE3(), // E[3],
      &r1  =  r[1], // .rotate(E1), 
      &r2  =  r[2], // .rotate(E2), 
      &r3  =  r[3], // .rotate(E3),
      &rI1 = rI[1], // .rotate(E1), 
      &rI2 = rI[2], // .rotate(E2), 
      &rI3 = rI[3], // .rotate(E3),
      &rJ1 = rJ[1], // .rotate(E1), 
      &rJ2 = rJ[2], // .rotate(E2), 
      &rJ3 = rJ[3]; // .rotate(E3);

    // NOTE[cmp]
    // SouzaFrameTransf::compTransfMatrixBasicGlobal must be 
    // called first to set Lr1, Lr2 and T

    // Matrix3D A;
    // for (int i = 0; i < 3; i++)
    //   for (int j = 0; j < 3; j++)
    //     A(i,j) = (double(i==j) - e1[i]*e1[j])/Ln;
    // getLMatrix(A, e1, r1, r2, Lr2);
    // getLMatrix(A, e1, r1, r3, Lr3);

    //
    // Ksigma1
    //
    {
      const double N = -pl[0]; // Axial force
      // a=0
      kg.assemble(A, 0, 0,  N);
      kg.assemble(A, 0, 6, -N);
      // a=1
      kg.assemble(A, 6, 0, -N);
      kg.assemble(A, 6, 6,  N);
    }

    //
    // Ksigma3
    //
    //  ks3 = [o kbar2  |  o kbar4];
    //
    //  where
    //
    //    kbar2 = -Lr2*(m(3)*S(rI3) + m(1)*S(rI1)) + Lr3*(m(3)*S(rI2) - m(2)*S(rI1)) ;
    //
    //    kbar4 =  Lr2*(m(3)*S(rJ3) - m(4)*S(rJ1)) - Lr3*(m(3)*S(rJ2) + m(5)*S(rJ1));
    //
    // or
    //
    //  ks3 = [o ka+kb  |  o kc+kd];
    //      = [o ka     |  o kc] + [o kb  |  o kd];
    //
    //  where
    //
    //    ka = -Lr2*S(rI3)*m(3)  
    //         +Lr2*S(rI1)*m(1);
    //    kb =  Lr3*S(rI2)*m(3)  
    //         -Lr3*S(rI1)*m(2);
    //
    //    kc =  Lr2*S(rJ3)*m(3)
    //         -Lr2*S(rJ1)*m(4);
    //    kd = -Lr3*S(rJ2)*m(3)  
    //         +Lr3*S(rJ1)*m(5);

    VectorND<6> m;
    m[0] =  0.5*pl[imx]/std::cos(ul(imx));
    m[2] = -0.5*pl[imy]/std::cos(ul(imy));
    m[1] =  0.5*pl[imz]/std::cos(ul(imz));

    m[3] =  0.5*pl[jmx]/std::cos(ul(jmx));
    m[5] = -0.5*pl[jmy]/std::cos(ul(jmy));
    m[4] =  0.5*pl[jmz]/std::cos(ul(jmz));


    static Matrix3D Sm;
    Sm.zero();
    Sm.addSpin(rI3,  m[3]);
    Sm.addSpin(rI1,  m[1]);
    static MatrixND<12,3> kbar;
    kbar.zero();
    kbar.addMatrixProduct(Lr2, Sm, -1.0);

    Sm.zero();
    Sm.addSpin(rI2,  m[3]);
    Sm.addSpin(rI1, -m[2]);
    kbar.addMatrixProduct(Lr3, Sm,  1.0);

    kg.assemble(kbar, 0, 3, 1.0);
    kg.assembleTranspose(kbar, 3, 0, 1.0);

    Sm.zero();
    Sm.addSpin(rJ3,  m[3]);
    Sm.addSpin(rJ1, -m[4]);
    kbar.zero();
    kbar.addMatrixProduct(Lr2, Sm, 1.0);

    Sm.zero();
    Sm.addSpin(rJ2, m[3]);
    Sm.addSpin(rJ1, m[5]);
    kbar.addMatrixProduct(Lr3, Sm,  -1.0);

    kg.assemble(kbar, 0, 9, 1.0);
    kg.assembleTranspose(kbar, 9, 0, 1.0);


    //
    // Ksigma4
    //
    {
      Matrix3D ks33;
      ks33.zero();
      ks33.addSpinProduct(e2, rI3,  m[3]);
      ks33.addSpinProduct(e3, rI2, -m[3]);
      ks33.addSpinProduct(e2, rI1,  m[1]);
      ks33.addSpinProduct(e1, rI2, -m[1]);
      ks33.addSpinProduct(e3, rI1,  m[2]);
      ks33.addSpinProduct(e1, rI3, -m[2]);
      kg.assemble(ks33, 3, 3, 1.0);
    }

    //
    // Ksigma4
    //
    {
      Matrix3D ks33;
      ks33.zero();
      ks33.addSpinProduct(e2, rJ3, -m[3]);
      ks33.addSpinProduct(e3, rJ2,  m[3]);
      ks33.addSpinProduct(e2, rJ1,  m[4]);
      ks33.addSpinProduct(e1, rJ2, -m[4]);
      ks33.addSpinProduct(e3, rJ1,  m[5]);
      ks33.addSpinProduct(e1, rJ3, -m[5]);
  
      kg.assemble(ks33, 9, 9, 1.0);
    }


    //
    // Ksigma5
    //
    //  Ks5 = [ Ks5_11   Ks5_12 | -Ks5_11   Ks5_14;
    //          Ks5_12'    O    | -Ks5_12'   O;
    //         -Ks5_11  -Ks5_12 |  Ks5_11  -Ks5_14;
    //          Ks5_14t     O   | -Ks5_14'   O];
    //
    //
    // v = (1/Ln)*(m(2)*rI2 + m(3)*rI3 + m(5)*rJ2 + m(6)*rJ3);
    //   = 1/Ln * (m[1]*rI2 + m[2]*rI3)
    //   + 1/Ln * (m[4]*rJ2 + m[5]*rJ3);
    //   = vi + vj
    //
    {
      Vector3D v;
      v.addVector(0.0, rI2, m[1]);
      v.addVector(1.0, rI3, m[2]);
      v.addVector(1.0, rJ2, m[4]);
      v.addVector(1.0, rJ3, m[5]);
      v /= Ln;

      // Ks5_11 = A*v*e1' + e1*v'*A + (e1'*v)*A;
      //        = A*vi*e1' + e1*vi'*A + (e1'*vi)*A
      //        + A*vj*e1' + e1*vj'*A + (e1'*vj)*A;
      //
      Matrix3D ks33;
      ks33.zero();
      ks33.addMatrix(A, e1.dot(v));

      static Matrix3D m33;
      m33.zero();
      m33.addTensorProduct(v, e1, 1.0);

      ks33.addMatrixProduct(A, m33, 1.0);

      for (int i = 0; i < 3; i++)
        for (int j = 0; j < 3; j++)
          m33(i,j) = e1[i]*v[j];

      ks33.addMatrixProduct(m33, A, 1.0);

      kg.assemble(ks33, 0, 0,  1.0);
      kg.assemble(ks33, 0, 6, -1.0);
      kg.assemble(ks33, 6, 0, -1.0);
      kg.assemble(ks33, 6, 6,  1.0);
    }

    // Ks5_12 = -(m(2)*A*S(rI2) + m(3)*A*S(rI3));

    Matrix3D ks33;
    ks33.zero();
    ks33.addMatrixSpinProduct(A, rI2, -m[1]);
    ks33.addMatrixSpinProduct(A, rI3, -m[2]);

    kg.assemble(ks33, 0, 3,  1.0);
    kg.assemble(ks33, 6, 3, -1.0);
    kg.assembleTranspose(ks33, 3, 0,  1.0);
    kg.assembleTranspose(ks33, 3, 6, -1.0);

    //  Ks5_14 = -(m(5)*A*S(rJ2) + m(6)*A*S(rJ3));

    ks33.zero();
    ks33.addMatrixSpinProduct(A, rJ2, -m[4]);
    ks33.addMatrixSpinProduct(A, rJ3, -m[5]);

    kg.assemble(ks33, 0, 9,  1.0);
    kg.assemble(ks33, 6, 9, -1.0);

    kg.assembleTranspose(ks33, 9, 0,  1.0);
    kg.assembleTranspose(ks33, 9, 6, -1.0);

    // Ksigma -------------------------------
    Vector3D rm = rI3;

    rm.addVector(1.0, rJ3, -1.0);
    kg.addMatrix(CrisfieldTransform::getKs2Matrix(A, e1, r1, Ln, r2, rm), m[3]);

//  rm = rJ2;
    rm.addVector(0.0, rJ2, -1.0);
    rm.addVector(1.0, rI2, -1.0);
    kg.addMatrix(CrisfieldTransform::getKs2Matrix(A, e1, r1, Ln, r3,  rm), m[3]);
    kg.addMatrix(CrisfieldTransform::getKs2Matrix(A, e1, r1, Ln, r2, rI1), m[1]);
    kg.addMatrix(CrisfieldTransform::getKs2Matrix(A, e1, r1, Ln, r3, rI1), m[2]);
    //
    kg.addMatrix(CrisfieldTransform::getKs2Matrix(A, e1, r1, Ln, r2, rJ1), m[4]);
    kg.addMatrix(CrisfieldTransform::getKs2Matrix(A, e1, r1, Ln, r3, rJ1), m[5]);

    //
    //  T' * diag (M .* tan(thetal))*T
    //

    for (int node=0; node<2; node++) {
      for (int k = 0; k < 3; k++) {
        const double factor =  pl[6*node+3+k] * std::tan(ul[(node ? jmx : imx) + k]);
        for (int i = 0; i < 12; i++) {
          const double Tki = T((node ? jmx : imx) + k,i);
          for (int j = 0; j < 12; j++)
            kg(i,j) += Tki * factor * T((node ? jmx : imx) + k, j);
        }
      }
    }

    return 0;
}


template <int nn, int ndf>
int
SouzaFrameTransf<nn,ndf>::getLocalAxes(Vector3D &e1, Vector3D &e2, Vector3D &e3) const
{
  for (int i = 0; i < 3; i++) {
    e1[i] = R0(i,0);
    e2[i] = R0(i,1);
    e3[i] = R0(i,2);
  }
  return 0;
}


template <int nn, int ndf>
double
SouzaFrameTransf<nn,ndf>::getLengthGrad()
{
  const int di = nodes[0]->getCrdsSensitivity();
  const int dj = nodes[1]->getCrdsSensitivity();

  Vector3D dxi{0.0};
  Vector3D dxj{0.0};

  if (di != 0)
    dxi(di-1) = 1.0;
  if (dj != 0)
    dxj(dj-1) = 1.0;

  return 1/L * dX.dot(dxj - dxi);
}

template <int nn, int ndf>
const Vector &
SouzaFrameTransf<nn,ndf>::getBasicDisplTotalGrad(int gradNumber)
{
    opserr << "WARNING CrdTransf::getBasicDisplTotalGrad() - this method "
        << " should not be called.\n";

    static Vector dummy(1);
    return dummy;
}

template <int nn, int ndf>
const Vector &
SouzaFrameTransf<nn,ndf>::getBasicDisplFixedGrad()
{
    opserr << "ERROR CrdTransf::getBasicDisplFixedGrad() - has not been"
           << " implemented yet for the chosen transformation\n.";

    static Vector dummy(1);
    return dummy;
}


template <int nn, int ndf>
const Vector &
SouzaFrameTransf<nn,ndf>::getGlobalResistingForceShapeSensitivity(const Vector &pb,
						   const Vector &p0,
						   int gradNumber)
{
    opserr << "ERROR CrdTransf::getGlobalResistingForceSensitivity() - has not been"
        << " implemented yet for the chosen transformation." << endln;
    
    static Vector dummy(1);
    return dummy;
}

template <int nn, int ndf>
void
SouzaFrameTransf<nn,ndf>::Print(OPS_Stream &s, int flag)
{

  if (flag == OPS_PRINT_CURRENTSTATE) {
      s << "\nFrameTransform: " << this->getTag() << " Type: SouzaFrameTransf";
      s << "\tvxz: " << Vector(vz);
  }

  if (flag == OPS_PRINT_PRINTMODEL_JSON) {
    s << OPS_PRINT_JSON_MATE_INDENT << "{";
    s << "\"name\": " << this->getTag() << ", ";
    s << "\"type\": \"SouzaFrameTransf\"" << ", ";
    s << "\"vecxz\": [" << vz(0) << ", " << vz(1) << ", " << vz(2) << "]";

    if (offsets != nullptr) {
      s << ", \"offsets\": [";
      for (int i=0; i<nn; i++) {
        s << "["
          << (*offsets)[i][0] << ", " 
          << (*offsets)[i][1] << ", "
          << (*offsets)[i][2] << "]";
        if (i < nn-1)
          s << ", ";
      }
      s << "]";
    }
    s << "}";
    return;
  }
}

