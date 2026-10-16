/* ****************************************************************** **
**    OpenSees - Open System for Earthquake Engineering Simulation    **
**          Pacific Earthquake Engineering Research Center            **
**                                                                    **
**                                                                    **
** (C) Copyright 1999, The Regents of the University of California    **
** All Rights Reserved.                                               **
**                                                                    **
** Commercial use of this program without express permission of the   **
** University of California, Berkeley, is strictly prohibited.  See   **
** file 'COPYRIGHT'  in main directory for information on usage and   **
** redistribution,  and for a DISCLAIMER OF ALL WARRANTIES.           **
**                                                                    **
** Developed by:                                                      **
**   Frank McKenna (fmckenna@ce.berkeley.edu)                         **
**   Gregory L. Fenves (fenves@ce.berkeley.edu)                       **
**   Filip C. Filippou (filippou@ce.berkeley.edu)                     **
**                                                                    **
** ****************************************************************** */

// Created: derived from SSPbrickUP (C.McGann, UW, 10.2011) - 2026
//
// Description: This file contains the implementation of the SSPTethraUP class
//                Single-Point Tetrahedron element (4 nodes, 1 Gauss point at
//                the element centroid) with a u-p formulation for 3D analysis
//                of saturated porous media.
//
//                Formulation follows the SSPbrickUP element (Stabilized
//                Single-Point Brick with u-p formulation, Zienkiewicz-Shiomi
//                generalized Biot u-p approximation). For a linear 4-node
//                tetrahedron the strain field is constant over the element
//                and the one-point rule at the centroid is exact for linear
//                shape functions: the element is rank-sufficient (no
//                hourglass modes), so the membrane stabilization of the
//                SSP brick is not needed. The pressure stabilization
//                (parameter alpha) is retained.
//
//                Matrix layout (OpenSees u-p convention, same as SSPbrickUP):
//                  getTangentStiff: [ Ks   0  ]
//                                   [ 0    0  ]
//                  getDamp:         [ C    -Q ]
//                                   [-Q'  -H  ]
//                  getMass:         [ Ms   0  ]
//                                   [ 0 -(S+Ps)]
//
// Reference:   Zienkiewicz, O.C. and Shiomi, T. (1984). "Dynamic behavior of
//                saturated porous media; the generalized Biot formulation and
//                its numerical solution." International Journal for Numerical
//                and Analytical Methods in Geomechanics, 8, 71-96.

#include "SSPTethraUP.h"

#include <elementAPI.h>
#include <Information.h>
#include <ElementResponse.h>
#include <ElementalLoad.h>
#include <ID.h>
#include <Domain.h>
#include <Node.h>
#include <Channel.h>
#include <Message.h>
#include <FEM_ObjectBroker.h>
#include <Renderer.h>
#include <G3Globals.h>
#include <ErrorHandler.h>
#include <NDMaterial.h>
#include <Parameter.h>

#include <math.h>
#include <stdlib.h>
#include <stdio.h>

#define OPS_Export

static int num_SSPTethraUP = 0;

OPS_Export void *
OPS_SSPTethraUP(void)
{
  if (num_SSPTethraUP == 0) {
    num_SSPTethraUP++;
    opserr << "SSPTethraUP element - derived from SSPbrickUP: C.McGann, P.Arduino, P.Mackenzie-Helnwein, U.Washington\n";
  }

  // Pointer to an element that will be returned
  Element *theElement = 0;

  int numRemainingInputArgs = OPS_GetNumRemainingInputArgs();

  if (numRemainingInputArgs < 13) {
    opserr << "Invalid #args, want: element SSPTethraUP eleTag? iNode? jNode? kNode? lNode? matTag? fBulk? fDen? k1? k2? k3? e? alpha? <b1? b2? b3?>\n";
    return 0;
  }

  int iData[6];
  double dData[10];
  dData[7] = 0.0;
  dData[8] = 0.0;
  dData[9] = 0.0;

  int numData = 6;
  if (OPS_GetIntInput(&numData, iData) != 0) {
    opserr << "WARNING invalid integer data: element SSPTethraUP " << iData[0] << endln;
    return 0;
  }

  int matID = iData[5];
  NDMaterial *theMaterial = OPS_getNDMaterial(matID);
  if (theMaterial == 0) {
    opserr << "WARNING element SSPTethraUP " << iData[0] << endln;
    opserr << " Material: " << matID << "not found\n";
    return 0;
  }

  numData = 7;
  if (OPS_GetDoubleInput(&numData, dData) != 0) {
    opserr << "WARNING invalid double data: element SSPTethraUP " << iData[0] << endln;
    return 0;
  }

  int massType = 0;
  int currentParam = 7;
  while (OPS_GetNumRemainingInputArgs() > 0) {
    const char *arg = OPS_GetString();
    if (strcmp(arg, "-lumped") == 0) {
      massType = 1;
    } else {
      if (currentParam < 10) {
        dData[currentParam++] = atof(arg);
      }
    }
  }

  // parsing was successful, allocate the element
  theElement = new SSPTethraUP(iData[0], iData[1], iData[2], iData[3], iData[4],
                               *theMaterial,
                               dData[0], dData[1], dData[2], dData[3], dData[4], dData[5], dData[6],
                               dData[7], dData[8], dData[9], massType);

  if (theElement == 0) {
    opserr << "WARNING could not create element of type SSPTethraUP\n";
    return 0;
  }

  return theElement;
}

SSPTethraUP::SSPTethraUP(int tag, int Nd1, int Nd2, int Nd3, int Nd4,
                         NDMaterial &theMat, double Kf, double Rf,
                         double k1, double k2, double k3,
                         double eVoid, double alpha,
                         double b1, double b2, double b3, int theMassType)
  :Element(tag, ELE_TAG_SSPTethraUP),
  theMaterial(0),
  mExternalNodes(STUP_NUM_NODE),
  mTangentStiffness(STUP_NUM_DOF, STUP_NUM_DOF),
  mInternalForces(STUP_NUM_DOF),
  Q(STUP_NUM_DOF),
  mMass(STUP_NUM_DOF, STUP_NUM_DOF),
  mDamp(STUP_NUM_DOF, STUP_NUM_DOF),
  mVol(0),
  Bnot(6, 12),
  mSolidK(12, 12),
  mSolidM(12, 12),
  mPerm(4, 4),
  mPressStab(4, 4),
  dNmod(4, 3),
  fBulk(Kf),
  fDens(Rf),
  mPorosity(eVoid/(1.0+eVoid)),
  mAlpha(alpha),
  applyLoad(0),
  massType(theMassType)
{
  if (b1 != 0.0 || b2 != 0.0 || b3 != 0.0) {
    appliedB[0] = b1;
    appliedB[1] = b2;
    appliedB[2] = b3;
    applyLoad = 1;
  }

  mExternalNodes(0) = Nd1;
  mExternalNodes(1) = Nd2;
  mExternalNodes(2) = Nd3;
  mExternalNodes(3) = Nd4;

  fBulk  = Kf;
  fDens  = Rf;
  mAlpha = alpha;

  b[0] = b1;
  b[1] = b2;
  b[2] = b3;

  appliedB[0] = 0.0;
  appliedB[1] = 0.0;
  appliedB[2] = 0.0;

  perm[0] = k1;
  perm[1] = k2;
  perm[2] = k3;

  mPorosity = eVoid/(1.0 + eVoid);

  // get copy of the material object
  NDMaterial *theMatCopy = theMat.getCopy("ThreeDimensional");
  if (theMatCopy != 0) {
    theMaterial = (NDMaterial *)theMatCopy;
  } else {
    opserr << "SSPTethraUP::SSPTethraUP - failed to get copy of material model\n";
  }

  // check material
  if (theMaterial == 0) {
    opserr << "SSPTethraUP::SSPTethraUP - failed to allocate material model pointer\n";
    exit(-1);
  }
}

// null constructor
SSPTethraUP::SSPTethraUP()
  :Element(0, ELE_TAG_SSPTethraUP),
  theMaterial(0),
  mExternalNodes(STUP_NUM_NODE),
  mTangentStiffness(STUP_NUM_DOF, STUP_NUM_DOF),
  mInternalForces(STUP_NUM_DOF),
  Q(STUP_NUM_DOF),
  mMass(STUP_NUM_DOF, STUP_NUM_DOF),
  mDamp(STUP_NUM_DOF, STUP_NUM_DOF),
  mVol(0),
  Bnot(6, 12),
  mSolidK(12, 12),
  mSolidM(12, 12),
  mPerm(4, 4),
  mPressStab(4, 4),
  dNmod(4, 3),
  fBulk(0),
  fDens(0),
  mAlpha(0),
  mPorosity(0),
  applyLoad(0),
  massType(0)
{
  b[0] = 0.0;  b[1] = 0.0;  b[2] = 0.0;
  perm[0] = 0.0;  perm[1] = 0.0;  perm[2] = 0.0;
  appliedB[0] = 0.0;  appliedB[1] = 0.0;  appliedB[2] = 0.0;
}

// destructor
SSPTethraUP::~SSPTethraUP()
{
  if (theMaterial != 0) {
    delete theMaterial;
  }
}

int
SSPTethraUP::getNumExternalNodes(void) const
{
  return STUP_NUM_NODE;
}

const ID &
SSPTethraUP::getExternalNodes(void)
{
  return mExternalNodes;
}

Node **
SSPTethraUP::getNodePtrs(void)
{
  return theNodes;
}

int
SSPTethraUP::getNumDOF(void)
{
  return STUP_NUM_DOF;
}

void
SSPTethraUP::setDomain(Domain *theDomain)
{
  theNodes[0] = theDomain->getNode(mExternalNodes(0));
  theNodes[1] = theDomain->getNode(mExternalNodes(1));
  theNodes[2] = theDomain->getNode(mExternalNodes(2));
  theNodes[3] = theDomain->getNode(mExternalNodes(3));

  for (int i = 0; i < 4; i++) {
    if (theNodes[i] == 0) {
      return;  // don't go any further - otherwise segmentation fault
    }
  }

  // establish element geometry and kinematic quantities
  // (volume, shape function gradients, strain-displacement matrix)
  setTetrahedronData();

  // establish mass matrix for solid phase (constant, only need to compute once)
  GetSolidMass();

  // establish permeability matrix (constant, only need to compute once)
  GetPermeabilityMatrix();

  // call the base-class method
  this->DomainComponent::setDomain(theDomain);
}

int
SSPTethraUP::commitState(void)
{
  int retVal = 0;
  // call element commitState to do any base class stuff
  if ((retVal = this->Element::commitState()) != 0) {
    opserr << "SSPTethraUP::commitState() - failed in base class\n";
  }
  retVal = theMaterial->commitState();

  return retVal;
}

int
SSPTethraUP::revertToLastCommit(void)
{
  return theMaterial->revertToLastCommit();
}

int
SSPTethraUP::revertToStart(void)
{
  return theMaterial->revertToStart();
}

int
SSPTethraUP::update(void)
// this function updates variables for an incremental step n to n+1
{
  // if the element geometry is invalid, do nothing
  if (mVol <= 0.0) {
    return 0;
  }

  // get trial displacement
  const Vector &mDisp_1 = theNodes[0]->getTrialDisp();
  const Vector &mDisp_2 = theNodes[1]->getTrialDisp();
  const Vector &mDisp_3 = theNodes[2]->getTrialDisp();
  const Vector &mDisp_4 = theNodes[3]->getTrialDisp();

  // assemble displacement vector (displacement dof only)
  Vector u(12);
  u(0) =  mDisp_1(0);
  u(1) =  mDisp_1(1);
  u(2) =  mDisp_1(2);
  u(3) =  mDisp_2(0);
  u(4) =  mDisp_2(1);
  u(5) =  mDisp_2(2);
  u(6) =  mDisp_3(0);
  u(7) =  mDisp_3(1);
  u(8) =  mDisp_3(2);
  u(9) =  mDisp_4(0);
  u(10) = mDisp_4(1);
  u(11) = mDisp_4(2);

  // compute strain and send it to the material
  Vector strain(6);
  strain = Bnot*u;
  theMaterial->setTrialStrain(strain);

  return 0;
}

const Matrix &
SSPTethraUP::getTangentStiff(void)
// this function computes the tangent stiffness matrix for the element
{
  // solid phase stiffness matrix
  GetSolidStiffness();

  // assemble full element stiffness matrix  [ K  0 ]
  // comprised of K submatrix                [ 0  0 ]
  mTangentStiffness.Zero();
  for (int i = 0; i < 4; i++) {

    int I    = 3*i;
    int Ip1  = 3*i+1;
    int Ip2  = 3*i+2;
    int II   = 4*i;
    int IIp1 = 4*i+1;
    int IIp2 = 4*i+2;

    for (int j = 0; j < 4; j++) {

      int J    = 3*j;
      int Jp1  = 3*j+1;
      int Jp2  = 3*j+2;
      int JJ   = 4*j;
      int JJp1 = 4*j+1;
      int JJp2 = 4*j+2;

      mTangentStiffness(II,JJ)     = mSolidK(I,J);
      mTangentStiffness(IIp1,JJ)   = mSolidK(Ip1,J);
      mTangentStiffness(IIp2,JJ)   = mSolidK(Ip2,J);
      mTangentStiffness(IIp1,JJp1) = mSolidK(Ip1,Jp1);
      mTangentStiffness(IIp2,JJp1) = mSolidK(Ip2,Jp1);
      mTangentStiffness(IIp2,JJp2) = mSolidK(Ip2,Jp2);
      mTangentStiffness(II,JJp1)   = mSolidK(I,Jp1);
      mTangentStiffness(II,JJp2)   = mSolidK(I,Jp2);
      mTangentStiffness(IIp1,JJp2) = mSolidK(Ip1,Jp2);
    }
  }

  return mTangentStiffness;
}

const Matrix &
SSPTethraUP::getInitialStiff(void)
// this function computes the initial tangent stiffness matrix for the element
{
  return getTangentStiff();
}

const Matrix &
SSPTethraUP::getDamp(void)
{
  Matrix dampC(12,12);

  // if the element geometry is invalid, return zero damping
  if (mVol <= 0.0) {
    mDamp.Zero();
    return mDamp;
  }

  // solid phase stiffness matrix
  GetSolidStiffness();

  // contribution of mass matrix for Rayleigh damping
  if (alphaM != 0.0) {
    dampC.addMatrix(0.0, mSolidM, alphaM);
  }
  // contribution of stiffness matrix for Rayleigh damping
  if (betaK != 0.0) {
    dampC.addMatrix(1.0, mSolidK, betaK);
  } if (betaK0 != 0.0) {
    dampC.addMatrix(1.0, mSolidK, betaK0);
  } if (betaKc != 0.0) {
    dampC.addMatrix(1.0, mSolidK, betaKc);
  }

  // compute coupling matrix Q
  Matrix couple(12,4);
  Matrix INp(6,4);

  INp.Zero();
  INp(0,0) = 0.25; INp(0,1) = 0.25; INp(0,2) = 0.25; INp(0,3) = 0.25;
  INp(1,0) = 0.25; INp(1,1) = 0.25; INp(1,2) = 0.25; INp(1,3) = 0.25;
  INp(2,0) = 0.25; INp(2,1) = 0.25; INp(2,2) = 0.25; INp(2,3) = 0.25;

  couple.Zero();
  couple.addMatrixTransposeProduct(0.0, Bnot, INp, mVol);

  // assemble full element damping matrix   [  C  -Q ]
  // comprised of C, Q, and H submatrices   [ -Q' -H ]
  mDamp.Zero();
  for (int i = 0; i < 4; i++) {

    int I    = 3*i;
    int Ip1  = 3*i+1;
    int Ip2  = 3*i+2;
    int II   = 4*i;
    int IIp1 = 4*i+1;
    int IIp2 = 4*i+2;
    int IIp3 = 4*i+3;

    for (int j = 0; j < 4; j++) {

      int J    = 3*j;
      int Jp1  = 3*j+1;
      int Jp2  = 3*j+2;
      int JJ   = 4*j;
      int JJp1 = 4*j+1;
      int JJp2 = 4*j+2;
      int JJp3 = 4*j+3;

      // contribution of solid phase damping
      mDamp(II,JJ)     = dampC(I,J);
      mDamp(IIp1,JJ)   = dampC(Ip1,J);
      mDamp(IIp2,JJ)   = dampC(Ip2,J);
      mDamp(IIp1,JJp1) = dampC(Ip1,Jp1);
      mDamp(IIp2,JJp1) = dampC(Ip2,Jp1);
      mDamp(IIp2,JJp2) = dampC(Ip2,Jp2);
      mDamp(II,JJp1)   = dampC(I,Jp1);
      mDamp(II,JJp2)   = dampC(I,Jp2);
      mDamp(IIp1,JJp2) = dampC(Ip1,Jp2);

      // contribution of fluid-solid coupling
      mDamp(JJp3,II)   = -couple(I,j);
      mDamp(JJp3,IIp1) = -couple(Ip1,j);
      mDamp(JJp3,IIp2) = -couple(Ip2,j);
      mDamp(II,JJp3)   = -couple(I,j);
      mDamp(IIp1,JJp3) = -couple(Ip1,j);
      mDamp(IIp2,JJp3) = -couple(Ip2,j);

      // contribution of permeability matrix
      mDamp(IIp3,JJp3) = -mPerm(i,j);
    }
  }

  return mDamp;
}

const Matrix &
SSPTethraUP::getMass(void)
{
  mMass.Zero();

  // if the element geometry is invalid, return zero mass
  if (mVol <= 0.0) {
    return mMass;
  }

  // compute compressibilty matrix term S
  // (Np_i = Np_j = 1/4 at the single gauss point of the tetrahedron)
  double oneOverQ = -0.0625*mVol*mPorosity/fBulk;

  // full mass matrix for the element [ M    0   ]
  //  includes M and S submatrices    [ 0 -(S+Ps)]
  for (int i = 0; i < 4; i++) {

    int I    = 3*i;
    int Ip1  = 3*i+1;
    int Ip2  = 3*i+2;
    int II   = 4*i;
    int IIp1 = 4*i+1;
    int IIp2 = 4*i+2;
    int IIp3 = 4*i+3;

    for (int j = 0; j < 4; j++) {

      int J    = 3*j;
      int Jp1  = 3*j+1;
      int Jp2  = 3*j+2;
      int JJ   = 4*j;
      int JJp1 = 4*j+1;
      int JJp2 = 4*j+2;
      int JJp3 = 4*j+3;

      // contribution of solid phase mass
      mMass(II,JJ)     = mSolidM(I,J);
      mMass(IIp1,JJ)   = mSolidM(Ip1,J);
      mMass(IIp2,JJ)   = mSolidM(Ip2,J);
      mMass(IIp1,JJp1) = mSolidM(Ip1,Jp1);
      mMass(IIp2,JJp1) = mSolidM(Ip2,Jp1);
      mMass(IIp2,JJp2) = mSolidM(Ip2,Jp2);
      mMass(II,JJp1)   = mSolidM(I,Jp1);
      mMass(II,JJp2)   = mSolidM(I,Jp2);
      mMass(IIp1,JJp2) = mSolidM(Ip1,Jp2);

      // contribution of compressibility and stabilization matrices
      if (massType == 1) {
        if (i == j) {
          double rowSum = 0.0;
          for (int k = 0; k < 4; k++) {
            rowSum += mPressStab(i, k);
          }
          mMass(IIp3, JJp3) = 4.0 * oneOverQ - rowSum;
        } else {
          mMass(IIp3, JJp3) = 0.0;
        }
      } else {
        mMass(IIp3,JJp3) = oneOverQ - mPressStab(i,j);
      }
    }
  }

  return mMass;
}

void
SSPTethraUP::zeroLoad(void)
{
  applyLoad = 0;
  appliedB[0] = 0.0;
  appliedB[1] = 0.0;
  appliedB[2] = 0.0;
  Q.Zero();

  return;
}

int
SSPTethraUP::addLoad(ElementalLoad *theLoad, double loadFactor)
{
  // body forces can be applied in a load pattern
  int type;
  const Vector &data = theLoad->getData(type, loadFactor);

  if (type == LOAD_TAG_SelfWeight) {
    applyLoad = 1;
    appliedB[0] += loadFactor*data(0)*b[0];
    appliedB[1] += loadFactor*data(1)*b[1];
    appliedB[2] += loadFactor*data(2)*b[2];
    return 0;
  } else {
    opserr << "SSPTethraUP::addLoad - load type unknown for ele with tag: " << this->getTag() << endln;
    return -1;
  }

  return -1;
}

int
SSPTethraUP::addInertiaLoadToUnbalance(const Vector &accel)
{
  // get mass density from the material
  double density = theMaterial->getRho();

  // do nothing if density is zero
  if (density == 0.0) {
    return 0;
  }

  // Get R * accel from the nodes
  const Vector &Raccel1 = theNodes[0]->getRV(accel);
  const Vector &Raccel2 = theNodes[1]->getRV(accel);
  const Vector &Raccel3 = theNodes[2]->getRV(accel);
  const Vector &Raccel4 = theNodes[3]->getRV(accel);

  static double ra[16];
  ra[0] =  Raccel1(0);
  ra[1] =  Raccel1(1);
  ra[2] =  Raccel1(2);
  ra[3] =  0.0;
  ra[4] =  Raccel2(0);
  ra[5] =  Raccel2(1);
  ra[6] =  Raccel2(2);
  ra[7] =  0.0;
  ra[8] =  Raccel3(0);
  ra[9] =  Raccel3(1);
  ra[10] = Raccel3(2);
  ra[11] = 0.0;
  ra[12] = Raccel4(0);
  ra[13] = Raccel4(1);
  ra[14] = Raccel4(2);
  ra[15] = 0.0;

  // compute mass matrix
  this->getMass();

  for (int i = 0; i < 16; i++) {
    Q(i) -= mMass(i,i)*ra[i];
  }

  return 0;
}

const Vector &
SSPTethraUP::getResistingForce(void)
// this function computes the resisting force vector for the element
{
  Vector f1(12);
  f1.Zero();
  Vector f2(4);
  f2.Zero();

  // if the element geometry is invalid, return zero forces
  if (mVol <= 0.0) {
    mInternalForces.Zero();
    mInternalForces.addVector(1.0, Q, -1.0);
    return mInternalForces;
  }

  // get stress from the material
  Vector mStress = theMaterial->getStress();

  // get trial displacement
  const Vector &mDisp_1 = theNodes[0]->getTrialDisp();
  const Vector &mDisp_2 = theNodes[1]->getTrialDisp();
  const Vector &mDisp_3 = theNodes[2]->getTrialDisp();
  const Vector &mDisp_4 = theNodes[3]->getTrialDisp();

  // assemble displacement vector
  Vector d(12);
  d(0) =  mDisp_1(0);
  d(1) =  mDisp_1(1);
  d(2) =  mDisp_1(2);
  d(3) =  mDisp_2(0);
  d(4) =  mDisp_2(1);
  d(5) =  mDisp_2(2);
  d(6) =  mDisp_3(0);
  d(7) =  mDisp_3(1);
  d(8) =  mDisp_3(2);
  d(9) =  mDisp_4(0);
  d(10) = mDisp_4(1);
  d(11) = mDisp_4(2);

  // add internal force from the stress  ->  fint = Bnot'*stress*volume
  f1.addMatrixTransposeVector(1.0, Bnot, mStress, mVol);

  // subtract body forces from internal force vector
  // (for the linear tetrahedron the exact nodal volume weight is V/4)
  double density = theMaterial->getRho();

  if (applyLoad == 0) {
    for (int i = 0; i < 4; i++) {
      f1(3*i)   -= density*b[0]*mVol/4.0;
      f1(3*i+1) -= density*b[1]*mVol/4.0;
      f1(3*i+2) -= density*b[2]*mVol/4.0;
    }
  } else {
    for (int i = 0; i < 4; i++) {
      f1(3*i)   -= density*appliedB[0]*mVol/4.0;
      f1(3*i+1) -= density*appliedB[1]*mVol/4.0;
      f1(3*i+2) -= density*appliedB[2]*mVol/4.0;
    }
  }

  // account for fluid body forces
  Matrix k(3,3);
  Vector body(3);
  // permeability tensor
  k(0,0) = perm[0];
  k(1,1) = perm[1];
  k(2,2) = perm[2];
  // body force vector
  if (applyLoad == 0) {
    body(0) = b[0];
    body(1) = b[1];
    body(2) = b[2];
  } else {
    body(0) = appliedB[0];
    body(1) = appliedB[1];
    body(2) = appliedB[2];
  }
  f2 = mVol*fDens*dNmod*k*body;

  // assemble full internal force vector for the element
  mInternalForces(0)  = f1(0);
  mInternalForces(1)  = f1(1);
  mInternalForces(2)  = f1(2);
  mInternalForces(3)  = f2(0);
  mInternalForces(4)  = f1(3);
  mInternalForces(5)  = f1(4);
  mInternalForces(6)  = f1(5);
  mInternalForces(7)  = f2(1);
  mInternalForces(8)  = f1(6);
  mInternalForces(9)  = f1(7);
  mInternalForces(10) = f1(8);
  mInternalForces(11) = f2(2);
  mInternalForces(12) = f1(9);
  mInternalForces(13) = f1(10);
  mInternalForces(14) = f1(11);
  mInternalForces(15) = f2(3);

  // inertial unbalance load
  mInternalForces.addVector(1.0, Q, -1.0);

  return mInternalForces;
}

const Vector &
SSPTethraUP::getResistingForceIncInertia()
{
  // compute current resisting force, mass, and damping
  this->getResistingForce();
  this->getMass();
  this->getDamp();

  // terms stemming from acceleration
  const Vector &accel1 = theNodes[0]->getTrialAccel();
  const Vector &accel2 = theNodes[1]->getTrialAccel();
  const Vector &accel3 = theNodes[2]->getTrialAccel();
  const Vector &accel4 = theNodes[3]->getTrialAccel();

  Vector a(16);
  a(0) =  accel1(0);
  a(1) =  accel1(1);
  a(2) =  accel1(2);
  a(3) =  accel1(3);
  a(4) =  accel2(0);
  a(5) =  accel2(1);
  a(6) =  accel2(2);
  a(7) =  accel2(3);
  a(8) =  accel3(0);
  a(9) =  accel3(1);
  a(10) = accel3(2);
  a(11) = accel3(3);
  a(12) = accel4(0);
  a(13) = accel4(1);
  a(14) = accel4(2);
  a(15) = accel4(3);

  mInternalForces.addMatrixVector(1.0, mMass, a, 1.0);

  // terms stemming from velocity
  const Vector &vel1 = theNodes[0]->getTrialVel();
  const Vector &vel2 = theNodes[1]->getTrialVel();
  const Vector &vel3 = theNodes[2]->getTrialVel();
  const Vector &vel4 = theNodes[3]->getTrialVel();

  Vector v(16);
  v(0) =  vel1(0);
  v(1) =  vel1(1);
  v(2) =  vel1(2);
  v(3) =  vel1(3);
  v(4) =  vel2(0);
  v(5) =  vel2(1);
  v(6) =  vel2(2);
  v(7) =  vel2(3);
  v(8) =  vel3(0);
  v(9) =  vel3(1);
  v(10) = vel3(2);
  v(11) = vel3(3);
  v(12) = vel4(0);
  v(13) = vel4(1);
  v(14) = vel4(2);
  v(15) = vel4(3);

  mInternalForces.addMatrixVector(1.0, mDamp, v, 1.0);

  return mInternalForces;
}

int
SSPTethraUP::sendSelf(int commitTag, Channel &theChannel)
{
  int res = 0;

  // note: we don't check for dataTag == 0 for Element
  // objects as that is taken care of in a commit by the Domain
  // object - don't want to have to do the check if sending data
  int dataTag = this->getDbTag();

  // SSPTethraUP packs its data into a Vector and sends this to theChannel
  // along with its dbTag and the commitTag passed in the arguments
  static Vector data(13);
  data(0) = this->getTag();
  data(1) = fBulk;
  data(2) = fDens;
  data(3) = perm[0];
  data(4) = perm[1];
  data(5) = perm[2];
  data(6) = mPorosity;
  data(7) = mAlpha;
  data(8) = b[0];
  data(9) = b[1];
  data(10) = b[2];
  data(11) = theMaterial->getClassTag();

  // now the element sends the id of its material
  int matDbTag = theMaterial->getDbTag();

  // NOTE: we do have to ensure that the material has a database
  // tag if we are sending to a database channel.
  if (matDbTag == 0) {
    matDbTag = theChannel.getDbTag();
    if (matDbTag != 0)
      theMaterial->setDbTag(matDbTag);
  }
  data(12) = matDbTag;

  res += theChannel.sendVector(dataTag, commitTag, data);
  if (res < 0) {
    opserr << "WARNING SSPTethraUP::sendSelf() - " << this->getTag() << " failed to send Vector\n";
    return res;
  }

  // SSPTethraUP then sends the tags of its four nodes
  res += theChannel.sendID(dataTag, commitTag, mExternalNodes);
  if (res < 0) {
    opserr << "WARNING SSPTethraUP::sendSelf() - " << this->getTag() << " failed to send ID\n";
    return res;
  }

  // finally, SSPTethraUP asks its material object to send itself
  res = theMaterial->sendSelf(commitTag, theChannel);
  if (res < 0) {
    opserr << "WARNING SSPTethraUP::sendSelf() - " << this->getTag() << " failed to send its Material\n";
    return -3;
  }

  return 0;
}

int
SSPTethraUP::recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker)
{
  int res = 0;
  int dataTag = this->getDbTag();

  // SSPTethraUP creates a Vector, receives the Vector and then sets the
  // internal data with the data in the Vector
  static Vector data(13);
  res += theChannel.recvVector(dataTag, commitTag, data);
  if (res < 0) {
    opserr << "WARNING SSPTethraUP::recvSelf() - failed to receive Vector\n";
    return res;
  }

  this->setTag((int)data(0));
  fBulk = data(1);
  fDens = data(2);
  perm[0] = data(3);
  perm[1] = data(4);
  perm[2] = data(5);
  mPorosity = data(6);
  mAlpha = data(7);
  b[0] = data(8);
  b[1] = data(9);
  b[2] = data(10);

  // SSPTethraUP now receives the tags of its four external nodes
  res += theChannel.recvID(dataTag, commitTag, mExternalNodes);
  if (res < 0) {
    opserr << "WARNING SSPTethraUP::recvSelf() - " << this->getTag() << " failed to receive ID\n";
    return res;
  }

  // finally, SSPTethraUP creates a material object of the correct type, sets its
  // database tag, and asks this new object to receive itself
  int matClass = (int)data(11);
  int matDb    = (int)data(12);

  // check if material object exists and that it is the right type
  if ((theMaterial == 0) || (theMaterial->getClassTag() != matClass)) {

    // if old one, delete it
    if (theMaterial != 0)
      delete theMaterial;

    // create new material object
    NDMaterial *theMatCopy = theBroker.getNewNDMaterial(matClass);
    theMaterial = (NDMaterial *)theMatCopy;

    if (theMaterial == 0) {
      opserr << "WARNING SSPTethraUP::recvSelf() - " << this->getTag()
        << " failed to get a blank Material of type " << matClass << endln;
      return -3;
    }
  }

  // NOTE: we set the dbTag before we receive the material
  theMaterial->setDbTag(matDb);
  res = theMaterial->recvSelf(commitTag, theChannel, theBroker);
  if (res < 0) {
    opserr << "WARNING SSPTethraUP::recvSelf() - " << this->getTag() << " failed to receive its Material\n";
    return -3;
  }

  return 0;
}

int
SSPTethraUP::displaySelf(Renderer &theViewer, int displayMode, float fact, const char **modes, int numMode)
{
  // vertex display coordinate vectors
  static Vector v1(3);
  static Vector v2(3);
  static Vector v3(3);
  static Vector v4(3);
  theNodes[0]->getDisplayCrds(v1, fact, displayMode);
  theNodes[1]->getDisplayCrds(v2, fact, displayMode);
  theNodes[2]->getDisplayCrds(v3, fact, displayMode);
  theNodes[3]->getDisplayCrds(v4, fact, displayMode);

  // color vector
  static Vector values(3);
  values(0) = 0;
  values(1) = 0;
  values(2) = 0;

  // draw polygons for each tetrahedron face
  int res = 0;
  static Matrix coords(3, 3); // rows are face nodes

  // face 1 (1 3 2)
  for (int i = 0; i < 3; i++) {
    coords(0, i) = v1(i);
    coords(1, i) = v3(i);
    coords(2, i) = v2(i);
  }
  res += theViewer.drawPolygon(coords, values, this->getTag());

  // face 2 (1 2 4)
  for (int i = 0; i < 3; i++) {
    coords(0, i) = v1(i);
    coords(1, i) = v2(i);
    coords(2, i) = v4(i);
  }
  res += theViewer.drawPolygon(coords, values, this->getTag());

  // face 3 (1 4 3)
  for (int i = 0; i < 3; i++) {
    coords(0, i) = v1(i);
    coords(1, i) = v4(i);
    coords(2, i) = v3(i);
  }
  res += theViewer.drawPolygon(coords, values, this->getTag());

  // face 4 (2 3 4)
  for (int i = 0; i < 3; i++) {
    coords(0, i) = v2(i);
    coords(1, i) = v3(i);
    coords(2, i) = v4(i);
  }
  res += theViewer.drawPolygon(coords, values, this->getTag());

  return res;
}

void
SSPTethraUP::Print(OPS_Stream &s, int flag)
{
  if (flag == OPS_PRINT_CURRENTSTATE) {
    opserr << "SSPTethraUP, element id:  " << this->getTag() << endln;
    opserr << "   Connected external nodes:  ";
    for (int i = 0; i < STUP_NUM_NODE; i++) {
      opserr << mExternalNodes(i) << " ";
    }
    opserr << endln;
  }

  if (flag == OPS_PRINT_PRINTMODEL_JSON) {
    s << "\t\t\t{";
    s << "\"name\": " << this->getTag() << ", ";
    s << "\"type\": \"SSPTethraUP\", ";
    s << "\"nodes\": [" << mExternalNodes(0) << ", ";
    for (int i = 1; i < 3; i++)
      s << mExternalNodes(i) << ", ";
    s << mExternalNodes(3) << "], ";
    s << "\"bodyForces\": [" << b[0] << ", " << b[1] << ", " << b[2] << "], ";
    s << "\"material\": \"" << theMaterial->getTag() << "\"}";
  }
}

Response*
SSPTethraUP::setResponse(const char **argv, int argc, OPS_Stream &eleInfo)
{

  if (strcmp(argv[0],"stress3D6") == 0) {
    return new ElementResponse(this, 1, Vector(6));
  }
  if (strcmp(argv[0],"strain3D6") == 0) {
    return new ElementResponse(this, 2, Vector(6));
  } else if (strcmp(argv[0], "mass") == 0) {
    return new ElementResponse(this, 3, Matrix(STUP_NUM_DOF, STUP_NUM_DOF));
  } else if ((strcmp(argv[0], "pressure") == 0) ||
             (strcmp(argv[0], "porePressure") == 0)) {
    return new ElementResponse(this, 4, Vector(1));
  } else {
    // no special recorders for this element, call the method in the material class
    return theMaterial->setResponse(argv, argc, eleInfo);
  }
}

int
SSPTethraUP::getResponse(int responseID, Information &eleInfo)
{
  // no special recorders for this element, call the method in the material class
  if (responseID == 1) {
    // get stress from the material
    const Vector &mStress = theMaterial->getStress();
    return eleInfo.setVector(mStress);
  } else if (responseID == 2) {
    // get strain from the material
    const Vector &mStrain = theMaterial->getStrain();
    return eleInfo.setVector(mStrain);
  } else if (responseID == 3) {
    return eleInfo.setMatrix(this->getMass());
  } else if (responseID == 4) {
    // pore pressure at the single gauss point: eP = sum(Np_i * p_i), Np_i = 1/4
    double eP = 0.0;
    for (int i = 0; i < 4; i++) {
      const Vector &disp = theNodes[i]->getTrialDisp();
      eP += 0.25*disp(3);
    }
    Vector pressure(1);
    pressure(0) = eP;
    return eleInfo.setVector(pressure);
  } else {
    // no special recorders for this element, call the method in the material class
    return theMaterial->getResponse(responseID, eleInfo);
  }
}

int
SSPTethraUP::setParameter(const char **argv, int argc, Parameter &param)
{
  if (argc < 1) {
    return -1;
  }

  int res = -1;

  // check for element parameters first
  if (strcmp(argv[0],"xPerm") == 0) {
    // permeability in direction 1
    return param.addObject(3, this);
  } else if (strcmp(argv[0],"yPerm") == 0) {
    // permeability in direction 2
    return param.addObject(4, this);
  } else if (strcmp(argv[0],"zPerm") == 0) {
    // permeability in direction 3
    return param.addObject(6, this);
  } else {
    // default is to call setParameter in the material
    int matRes = res;
    matRes = theMaterial->setParameter(argv, argc, param);
    if (matRes != -1) {
      res = matRes;
    }
  }

  return res;
}

int
SSPTethraUP::updateParameter(int parameterID, Information &info)
{
  int res = -1;
  int matRes = res;

  if (parameterID == res) {
    return -1;
  } else if (parameterID == 3) {
    // update element permeability in direction 1
    perm[0] = info.theDouble;
    GetPermeabilityMatrix();
    return 0;
  } else if (parameterID == 4) {
    // update element permeability in direction 2
    perm[1] = info.theDouble;
    GetPermeabilityMatrix();
    return 0;
  } else if (parameterID == 6) {
    // update element permeability in direction 3
    perm[2] = info.theDouble;
    GetPermeabilityMatrix();
    return 0;
  } else {
    // update the material parameter
    matRes = theMaterial->updateParameter(parameterID, info);
    if (matRes != -1) {
      res = matRes;
    }
    return res;
  }
}

void
SSPTethraUP::setTetrahedronData(void)
// this function computes the element volume, the shape function gradients
// at the single gauss point (element centroid) and the strain-displacement
// matrix for the 4-node linear tetrahedron
{
  // natural derivatives of shape functions
  // parent element: node 1 at (0,0,0), node 2 at (1,0,0),
  //                 node 3 at (0,1,0), node 4 at (0,0,1)
  // N1 = 1-xi-eta-zeta, N2 = xi, N3 = eta, N4 = zeta
  static const double dNlocData[4][3] = {
    { -1.0, -1.0, -1.0 },
    {  1.0,  0.0,  0.0 },
    {  0.0,  1.0,  0.0 },
    {  0.0,  0.0,  1.0 }
  };
  Matrix dNloc(4,3);
  for (int i = 0; i < 4; i++) {
    for (int j = 0; j < 3; j++) {
      dNloc(i,j) = dNlocData[i][j];
    }
  }

  // assemble nodal coordinate matrix (3 x 4)
  const Vector &crd1 = theNodes[0]->getCrds();
  const Vector &crd2 = theNodes[1]->getCrds();
  const Vector &crd3 = theNodes[2]->getCrds();
  const Vector &crd4 = theNodes[3]->getCrds();

  Matrix coords(3,4);
  for (int a = 0; a < 3; a++) {
    coords(a,0) = crd1(a);
    coords(a,1) = crd2(a);
    coords(a,2) = crd3(a);
    coords(a,3) = crd4(a);
  }

  // jacobian matrix: columns are the edge vectors from node 1
  Matrix Jmat(3,3);
  Jmat = coords*dNloc;

  // jacobian determinant (detJ = 6*Volume)
  double detJ =  Jmat(0,0)*Jmat(1,1)*Jmat(2,2) + Jmat(0,1)*Jmat(1,2)*Jmat(2,0) + Jmat(0,2)*Jmat(1,0)*Jmat(2,1)
                -Jmat(0,2)*Jmat(1,1)*Jmat(2,0) - Jmat(0,1)*Jmat(1,0)*Jmat(2,2) - Jmat(0,0)*Jmat(1,2)*Jmat(2,1);

  // check element geometry: node ordering must define a positive volume
  // (nodes 1-2-3 counterclockwise when seen from node 4)
  if (detJ <= 0.0) {
    opserr << "SSPTethraUP::setTetrahedronData() - detJ <= 0 for element: " << this->getTag() << endln;
    opserr << " non-positive tetrahedron volume, detJ = " << detJ << endln;
    opserr << " check node ordering: nodes 1-2-3 must appear counterclockwise when viewed from node 4\n";
    mVol = 0.0;
    return;
  }

  mVol = detJ/6.0;

  // inverse of jacobian matrix
  Matrix Jinv(3,3);
  Jmat.Invert(Jinv);

  // shape function gradients at the single gauss point (element centroid)
  // for the linear tetrahedron they are constant over the element
  Matrix dN(4,3);
  dN = dNloc*Jinv;

  for (int i = 0; i < 4; i++) {
    dNmod(i,0) = dN(i,0);
    dNmod(i,1) = dN(i,1);
    dNmod(i,2) = dN(i,2);
  }

  // strain-displacement matrix (6 x 12)
  // rows: eps_xx, eps_yy, eps_zz, gamma_xy, gamma_yz, gamma_zx
  Bnot.Zero();
  for (int i = 0; i < 4; i++) {
    Bnot(0,3*i)   = dNmod(i,0);
    Bnot(1,3*i+1) = dNmod(i,1);
    Bnot(2,3*i+2) = dNmod(i,2);
    Bnot(3,3*i)   = dNmod(i,1);
    Bnot(3,3*i+1) = dNmod(i,0);
    Bnot(4,3*i+1) = dNmod(i,2);
    Bnot(4,3*i+2) = dNmod(i,1);
    Bnot(5,3*i)   = dNmod(i,2);
    Bnot(5,3*i+2) = dNmod(i,0);
  }

  return;
}

void
SSPTethraUP::GetSolidStiffness(void)
// this function computes the stiffness matrix for the solid phase
{
  // get material tangent
  const Matrix &Cmat = theMaterial->getTangent();

  // full element stiffness matrix
  // (no hourglass stabilization needed: the linear tetrahedron with
  //  one-point integration has constant strain and no zero-energy modes)
  mSolidK.Zero();
  mSolidK.addMatrixTripleProduct(1.0, Bnot, Cmat, mVol);

  return;
}

void
SSPTethraUP::GetSolidMass(void)
// this function computes the mass matrix for the solid phase
{
  mSolidM.Zero();

  // if the element geometry is invalid, return zero mass
  if (mVol <= 0.0) {
    return;
  }

  // get solid mass density from the material
  double density = theMaterial->getRho();

  // return zero matrix if density is zero
  if (density == 0.0) {
    return;
  }

  // use nodal volume weights to get lumped mass values
  // (exact for the linear tetrahedron: integral of N_i = V/4)
  double massTerm;
  for (int i = 0; i < 4; i++) {
    massTerm = density*mVol/4.0;
    mSolidM(3*i,3*i)   += massTerm;
    mSolidM(3*i+1,3*i+1) += massTerm;
    mSolidM(3*i+2,3*i+2) += massTerm;
  }

  return;
}

void
SSPTethraUP::GetPermeabilityMatrix(void)
// this function computes the permeability matrix for the element
{
  Matrix k(3,3);
  Matrix dNT(3,4);
  mPerm.Zero();
  mPressStab.Zero();

  // if the element geometry is invalid, return zero matrices
  if (mVol <= 0.0) {
    return;
  }

  // permeability tensor
  k(0,0) = perm[0];
  k(1,1) = perm[1];
  k(2,2) = perm[2];

  // the pressure shape function gradients are constant for the linear
  // tetrahedron, therefore the one-point rule at the centroid is exact
  for (int i = 0; i < 4; i++) {
    dNT(0,i) = dNmod(i,0);
    dNT(1,i) = dNmod(i,1);
    dNT(2,i) = dNmod(i,2);
  }

  // permeability matrix H = int (grad Np)^T k (grad Np) dV
  mPerm.addMatrixTripleProduct(1.0, dNT, k, mVol);

  // pressure stabilization matrix Ps = alpha * int (grad Np)^T (grad Np) dV
  mPressStab.addMatrixTransposeProduct(1.0, dNT, dNT, mVol*mAlpha);

  return;
}
