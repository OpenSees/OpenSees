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

// $Revision$
// $Date$
// $URL$

// File: ~/model/ModElasticBeam3d.C
//
// Written: fmk 11/95
// Revised:
//
// Purpose: This file contains the class definition for ModElasticBeam3d.
// ModElasticBeam3d is a 3d beam element. As such it can only
// connect to a node with 6-dof.

#include <Channel.h>
#include <Domain.h>
#include <FEM_ObjectBroker.h>
#include <ModElasticBeam3d.h>

#include <CrdTransf.h>
#include <ElementResponse.h>
#include <ElementalLoad.h>
#include <ID.h>
#include <Information.h>
#include <Parameter.h>
#include <Renderer.h>
#include <SectionForceDeformation.h>
#include <elementAPI.h>
#include <math.h>
#include <stdlib.h>
#include <string>

Matrix ModElasticBeam3d::K(12, 12);
Vector ModElasticBeam3d::P(12);
Matrix ModElasticBeam3d::kb(6, 6);

void *OPS_ModElasticBeam3d(void) {
  int numArgs = OPS_GetNumRemainingInputArgs();

  // Read the optional arguments first
  double mass = 0.0;
  int cMass = 0;
  int numData = 1;
  int numOptionalArgs = 0;
  while (OPS_GetNumRemainingInputArgs() > 0) {
    std::string theType = OPS_GetString();
    if (theType == "-mass") {
      numOptionalArgs++;
      if (OPS_GetNumRemainingInputArgs() > 0) {
        numData = 1;
        if (OPS_GetDoubleInput(&numData, &mass) < 0)
          return 0;
        numOptionalArgs++;
      }
    } else if (theType == "-cMass") {
      numOptionalArgs++;
      cMass = 1;
    }
  }

  if (numArgs > 0) {
    OPS_ResetCurrentInputArg(-numArgs);
  }
  numArgs = numArgs - numOptionalArgs;

  if (numArgs != 16) {
    opserr << "ERROR not enough args provided, want: element ModElasticBeam3d tag? iNode? jNode? A? E? G? J? Ix? Iy? K11x? K33x? K44x? K11y? K33y? K44y? transfType? <-rho rho?> <-cMass>" << endln;
    return 0;
  }

  int ndm = OPS_GetNDM();
  int ndf = OPS_GetNDF();
  if (ndm != 3 || ndf != 6) {
    opserr << "ndm must be 3 and ndf must be 6"  << endln;
    return 0;
  }

  // inputs:
  int iData[3];
  numData = 3;
  if (OPS_GetIntInput(&numData, &iData[0]) < 0)
    return 0;

  SectionForceDeformation *theSection = 0;
  CrdTransf *theTrans = 0;
  double data[12];
  int transfTag;

  numData = 12;
  if (OPS_GetDoubleInput(&numData, &data[0]) < 0)
	return 0;
  numData = 1;
  if (OPS_GetIntInput(&numData, &transfTag) < 0)
	return 0;
  theTrans = OPS_getCrdTransf(transfTag);
  if (theTrans == 0) {
	opserr << "no CrdTransf is found" << endln;
	return 0;
  }

  double a = data[0];
  double e = data[1];
  double g = data[2];
  double jx = data[3];
  double iy = data[4];
  double iz = data[5];
  double k11y = data[6];
  double k33y = data[7];
  double k44y = data[8];
  double k11z = data[9];
  double k33z = data[10];
  double k44z = data[11];
  int tag = iData[0];
  int Nd1 = iData[1];
  int Nd2 = iData[2];

  return new ModElasticBeam3d(tag, a, e, g, jx, iy, iz,
							  k11y, k33y, k44y, k11z, k33z, k44z,
                              Nd1, Nd2, *theTrans, mass, cMass);
}

ModElasticBeam3d::ModElasticBeam3d()
    : Element(0, ELE_TAG_ModElasticBeam3d), A(0.0), E(0.0), G(0.0), Jx(0.0),
      Iy(0.0), Iz(0.0),
	  K11y(0.0), K33y(0.0), K44y(0.00),
	  K11z(0.0), K33z(0.0), K44z(0.00),
      rho(0.0), cMass(0), Q(12), q(6), wx(0.0), wy(0.0),
      wz(0.0), connectedExternalNodes(2), theCoordTransf(0) {
  // does nothing
  q0[0] = 0.0;
  q0[1] = 0.0;
  q0[2] = 0.0;
  q0[3] = 0.0;
  q0[4] = 0.0;

  p0[0] = 0.0;
  p0[1] = 0.0;
  p0[2] = 0.0;
  p0[3] = 0.0;
  p0[4] = 0.0;

  // set node pointers to NULL
  for (int i = 0; i < 2; i++)
    theNodes[i] = 0;
}

ModElasticBeam3d::ModElasticBeam3d(int tag, double a, double e, double g,
                                   double jx, double iy, double iz, double k11y,
                                   double k33y, double k44y, double k11z,
                                   double k33z, double k44z, int Nd1, int Nd2,
                                   CrdTransf &coordTransf, double r, int cm)
    : Element(tag, ELE_TAG_ModElasticBeam3d), A(a), E(e), G(g), Jx(jx), Iy(iy),
      Iz(iz),
	  K11y(k11y), K33y(k33y), K44y(k44y),
	  K11z(k11z), K33z(k33z), K44z(k44z),
      rho(r), cMass(cm), Q(12), q(6), wx(0.0), wy(0.0), wz(0.0),
      connectedExternalNodes(2), theCoordTransf(0) {
  connectedExternalNodes(0) = Nd1;
  connectedExternalNodes(1) = Nd2;

  theCoordTransf = coordTransf.getCopy3d();

  if (!theCoordTransf) {
    opserr << "ModElasticBeam3d::ModElasticBeam3d -- failed to get copy of "
              "coordinate transformation" << endln;
    exit(-1);
  }

  q0[0] = 0.0;
  q0[1] = 0.0;
  q0[2] = 0.0;
  q0[3] = 0.0;
  q0[4] = 0.0;

  p0[0] = 0.0;
  p0[1] = 0.0;
  p0[2] = 0.0;
  p0[3] = 0.0;
  p0[4] = 0.0;

  // set node pointers to NULL
  for (int i = 0; i < 2; i++)
    theNodes[i] = 0;
}

ModElasticBeam3d::~ModElasticBeam3d() {
  if (theCoordTransf)
    delete theCoordTransf;
}

int ModElasticBeam3d::getNumExternalNodes(void) const { return 2; }

const ID &ModElasticBeam3d::getExternalNodes(void) {
  return connectedExternalNodes;
}

Node **ModElasticBeam3d::getNodePtrs(void) { return theNodes; }

int ModElasticBeam3d::getNumDOF(void) { return 12; }

void ModElasticBeam3d::setDomain(Domain *theDomain) {
  if (theDomain == 0) {
    opserr << "ModElasticBeam3d::setDomain -- Domain is null" << endln;
    exit(-1);
  }

  theNodes[0] = theDomain->getNode(connectedExternalNodes(0));
  theNodes[1] = theDomain->getNode(connectedExternalNodes(1));

  if (theNodes[0] == 0) {
    opserr << "ModElasticBeam3d::setDomain  tag: " << this->getTag()
           << " -- Node 1: " << connectedExternalNodes(0)
           << " does not exist" << endln;
    exit(-1);
  }

  if (theNodes[1] == 0) {
    opserr << "ModElasticBeam3d::setDomain  tag: " << this->getTag()
           << " -- Node 2: " << connectedExternalNodes(1)
           << " does not exist" << endln;
    exit(-1);
  }

  int dofNd1 = theNodes[0]->getNumberDOF();
  int dofNd2 = theNodes[1]->getNumberDOF();

  if (dofNd1 != 6) {
    opserr << "ModElasticBeam3d::setDomain  tag: " << this->getTag()
           << " -- Node 1: " << connectedExternalNodes(0)
           << " has incorrect number of DOF" << endln;
    exit(-1);
  }

  if (dofNd2 != 6) {
    opserr << "ModElasticBeam3d::setDomain  tag: " << this->getTag()
           << " -- Node 2: " << connectedExternalNodes(1)
           << " has incorrect number of DOF" << endln;
    exit(-1);
  }

  this->DomainComponent::setDomain(theDomain);

  if (theCoordTransf->initialize(theNodes[0], theNodes[1]) != 0) {
    opserr << "ModElasticBeam3d::setDomain  tag: " << this->getTag()
           << " -- Error initializing coordinate transformation" << endln;
    exit(-1);
  }

  double L = theCoordTransf->getInitialLength();

  if (L == 0.0) {
    opserr << "ModElasticBeam3d::setDomain  tag: " << this->getTag()
           << " -- Element has zero length" << endln;
    exit(-1);
  }
}

int ModElasticBeam3d::commitState() {
  int retVal = 0;
  // call element commitState to do any base class stuff
  if ((retVal = this->Element::commitState()) != 0) {
    opserr << "ModElasticBeam3d::commitState () - failed in base class";
  }
  retVal += theCoordTransf->commitState();
  return retVal;
}

int ModElasticBeam3d::revertToLastCommit() {
  int retVal = 0;
  retVal += theCoordTransf->revertToLastCommit();
  return retVal;
}

int ModElasticBeam3d::revertToStart() {
  int retVal = 0;
  retVal += theCoordTransf->revertToStart();
  return retVal;
}

int ModElasticBeam3d::update(void) { return theCoordTransf->update(); }

const Matrix &ModElasticBeam3d::getTangentStiff(void) {
  const Vector &v = theCoordTransf->getBasicTrialDisp();

  double L = theCoordTransf->getInitialLength();
  double oneOverL = 1.0 / L;
  double EoverL = E * oneOverL;
  double EAoverL = A * EoverL;        // EA/L
  double GJoverL = G * Jx * oneOverL; // GJ/L

  q(0) = EAoverL * v(0);
  q(5) = GJoverL * v(5);
  kb.Zero();
  kb(0, 0) = EAoverL;
  kb(5, 5) = GJoverL;

  double EIzoverL2 = K44z * Iz * EoverL;   // 2EIz/L, modified
  double EIzoverL4 = K11z * Iz * EoverL;   // 4EIz/L, modified
  double EIzoverL6 = K33z * Iz * EoverL;   // 4EIz/L, modified
  q(1) = EIzoverL4 * v(1) + EIzoverL2 * v(2);
  q(2) = EIzoverL2 * v(1) + EIzoverL6 * v(2);
  kb(1, 1) = EIzoverL4;
  kb(2, 2) = EIzoverL6;
  kb(2, 1) = kb(1, 2) = EIzoverL2;

  double EIyoverL2 = K44y * Iy * EoverL; // 2EIy/L, modified
  double EIyoverL4 = K11y * Iy * EoverL; // 4EIy/L, modified
  double EIyoverL6 = K33y * Iy * EoverL; // 4EIy/L, modified
  q(3) = EIyoverL4 * v(3) + EIyoverL2 * v(4);
  q(4) = EIyoverL2 * v(3) + EIyoverL6 * v(4);
  kb(3, 3) = EIyoverL4;
  kb(4, 4) = EIyoverL6;
  kb(4, 3) = kb(3, 4) = EIyoverL2;

  q(0) += q0[0];
  q(1) += q0[1];
  q(2) += q0[2];
  q(3) += q0[3];
  q(4) += q0[4];

  return theCoordTransf->getGlobalStiffMatrix(kb, q);
}

const Matrix &ModElasticBeam3d::getInitialStiff(void) {
  //  const Vector &v = theCoordTransf->getBasicTrialDisp();

  double L = theCoordTransf->getInitialLength();

  double oneOverL = 1.0 / L;
  double EoverL = E * oneOverL;
  double EAoverL = A * EoverL;        // EA/L
  double GJoverL = G * Jx * oneOverL; // GJ/L

  kb.Zero();
  kb(0, 0) = EAoverL;
  kb(5, 5) = GJoverL;

  double EIzoverL2 = K44z * Iz * EoverL;   // 2EIz/L, modified
  double EIzoverL4 = K11z * Iz * EoverL;   // 4EIz/L, modified
  double EIzoverL6 = K33z * Iz * EoverL;   // 4EIz/L, modified
  kb(1, 1) = EIzoverL4;
  kb(2, 2) = EIzoverL6;
  kb(2, 1) = kb(1, 2) = EIzoverL2;

  double EIyoverL2 = K44y * Iy * EoverL;   // 2EIy/L, modified
  double EIyoverL4 = K11y * Iy * EoverL;   // 4EIy/L, modified
  double EIyoverL6 = K33y * Iy * EoverL;   // 4EIy/L, modified
  kb(3, 3) = EIyoverL4;
  kb(4, 4) = EIyoverL6;
  kb(4, 3) = kb(3, 4) = EIyoverL2;

  return theCoordTransf->getInitialGlobalStiffMatrix(kb);
}

const Matrix &ModElasticBeam3d::getMass(void) {
  K.Zero();

  if (rho > 0.0) {
    // get initial element length
    double L = theCoordTransf->getInitialLength();
    if (cMass == 0) {
      // lumped mass matrix
      double m = 0.5 * rho * L;
      K(0, 0) = m;
      K(1, 1) = m;
      K(2, 2) = m;
      K(6, 6) = m;
      K(7, 7) = m;
      K(8, 8) = m;
    } else {
      // consistent mass matrix
      static Matrix ml(12, 12);
      double m = rho * L / 420.0;
      ml(0, 0) = ml(6, 6) = m * 140.0;
      ml(0, 6) = ml(6, 0) = m * 70.0;
      ml(3, 3) = ml(9, 9) = m * (Jx / A) * 140.0;
      ml(3, 9) = ml(9, 3) = m * (Jx / A) * 70.0;

      ml(2, 2) = ml(8, 8) = m * 156.0;
      ml(2, 8) = ml(8, 2) = m * 54.0;
      ml(4, 4) = ml(10, 10) = m * 4.0 * L * L;
      ml(4, 10) = ml(10, 4) = -m * 3.0 * L * L;
      ml(2, 4) = ml(4, 2) = -m * 22.0 * L;
      ml(8, 10) = ml(10, 8) = -ml(2, 4);
      ml(2, 10) = ml(10, 2) = m * 13.0 * L;
      ml(4, 8) = ml(8, 4) = -ml(2, 10);

      ml(1, 1) = ml(7, 7) = m * 156.0;
      ml(1, 7) = ml(7, 1) = m * 54.0;
      ml(5, 5) = ml(11, 11) = m * 4.0 * L * L;
      ml(5, 11) = ml(11, 5) = -m * 3.0 * L * L;
      ml(1, 5) = ml(5, 1) = m * 22.0 * L;
      ml(7, 11) = ml(11, 7) = -ml(1, 5);
      ml(1, 11) = ml(11, 1) = -m * 13.0 * L;
      ml(5, 7) = ml(7, 5) = -ml(1, 11);

      // transform local mass matrix to global system
      K = theCoordTransf->getGlobalMatrixFromLocal(ml);
    }
  }

  return K;
}

void ModElasticBeam3d::zeroLoad(void) {
  Q.Zero();

  q0[0] = 0.0;
  q0[1] = 0.0;
  q0[2] = 0.0;
  q0[3] = 0.0;
  q0[4] = 0.0;

  p0[0] = 0.0;
  p0[1] = 0.0;
  p0[2] = 0.0;
  p0[3] = 0.0;
  p0[4] = 0.0;

  wx = 0.0;
  wy = 0.0;
  wz = 0.0;

  return;
}

int ModElasticBeam3d::addLoad(ElementalLoad *theLoad, double loadFactor) {
  int type;
  const Vector &data = theLoad->getData(type, loadFactor);
  double L = theCoordTransf->getInitialLength();

  if (type == LOAD_TAG_Beam3dUniformLoad) {
    double wy = data(0) * loadFactor; // Transverse
    double wz = data(1) * loadFactor; // Transverse
    double wx = data(2) * loadFactor; // Axial (+ve from node I to J)

    this->wx += wx;
    this->wy += wy;
    this->wz += wz;

    double Vy = 0.5 * wy * L;
    double Mz = Vy * L / 6.0; // wy*L*L/12
    double Vz = 0.5 * wz * L;
    double My = Vz * L / 6.0; // wz*L*L/12
    double P = wx * L;

    // Reactions in basic system
    p0[0] -= P;
    p0[1] -= Vy;
    p0[2] -= Vy;
    p0[3] -= Vz;
    p0[4] -= Vz;

    // Fixed end forces in basic system
    q0[0] -= 0.5 * P;
    q0[1] -= Mz;
    q0[2] += Mz;
    q0[3] += My;
    q0[4] -= My;

  } else if (type == LOAD_TAG_Beam3dPartialUniformLoad) {
    double wa = data(2) * loadFactor; // Axial
    double wy = data(0) * loadFactor; // Transverse
    double wz = data(1) * loadFactor; // Transverse
    double a = data(3) * L;
    double b = data(4) * L;
    double c = 0.5 * (b + a);
    double cOverL = c / L;

    double P = wa * (b - a);
    double Fy = wy * (b - a);
    double Fz = wz * (b - a);

    // Reactions in basic system
    p0[0] -= P;
    double V1, V2;
    V1 = Fy * (1.0 - cOverL);
    V2 = Fy * cOverL;
    p0[1] -= V1;
    p0[2] -= V2;
    V1 = Fz * (1.0 - cOverL);
    V2 = Fz * cOverL;
    p0[3] -= V1;
    p0[4] -= V2;

    // Fixed end forces in basic system
    q0[0] -= P * cOverL;
    double M1, M2;
    double beta2 = (1 - cOverL) * (1 - cOverL);
    double alfa2 = (cOverL) * (cOverL);
    double gamma2 = (b - a) / L;
    gamma2 *= gamma2;

    M1 = -wy * (b - a) * (c * beta2 + gamma2 / 12.0 * (L - 3 * (L - c)));
    M2 = wy * (b - a) * ((L - c) * alfa2 + gamma2 / 12.0 * (L - 3 * c));
    q0[1] += M1;
    q0[2] += M2;
    M1 = -wz * (b - a) * (c * beta2 + gamma2 / 12.0 * (L - 3 * (L - c)));
    M2 = wz * (b - a) * ((L - c) * alfa2 + gamma2 / 12.0 * (L - 3 * c));
    q0[3] -= M1;
    q0[4] -= M2;
  } else if (type == LOAD_TAG_Beam3dPointLoad) {
    double Py = data(0) * loadFactor;
    double Pz = data(1) * loadFactor;
    double N = data(2) * loadFactor;
    double aOverL = data(3);

    if (aOverL < 0.0 || aOverL > 1.0)
      return 0;

    double a = aOverL * L;
    double b = L - a;

    // Reactions in basic system
    p0[0] -= N;
    double V1, V2;
    V1 = Py * (1.0 - aOverL);
    V2 = Py * aOverL;
    p0[1] -= V1;
    p0[2] -= V2;
    V1 = Pz * (1.0 - aOverL);
    V2 = Pz * aOverL;
    p0[3] -= V1;
    p0[4] -= V2;

    double L2 = 1.0 / (L * L);
    double a2 = a * a;
    double b2 = b * b;

    // Fixed end forces in basic system
    q0[0] -= N * aOverL;
    double M1, M2;
    M1 = -a * b2 * Py * L2;
    M2 = a2 * b * Py * L2;
    q0[1] += M1;
    q0[2] += M2;
    M1 = -a * b2 * Pz * L2;
    M2 = a2 * b * Pz * L2;
    q0[3] -= M1;
    q0[4] -= M2;
  } else {
    opserr << "ModElasticBeam3d::addLoad()  -- load type unknown for element "
              "with tag: "
           << this->getTag() << endln;
    return -1;
  }

  return 0;
}

int ModElasticBeam3d::addInertiaLoadToUnbalance(const Vector &accel) {
  if (rho == 0.0)
    return 0;

  // get R * accel from the nodes
  const Vector &Raccel1 = theNodes[0]->getRV(accel);
  const Vector &Raccel2 = theNodes[1]->getRV(accel);

  if (6 != Raccel1.Size() || 6 != Raccel2.Size()) {
    opserr << "ModElasticBeam3d::addInertiaLoadToUnbalance matrix and vector "
              "sizes are incompatible" << endln;
    return -1;
  }

  // want to add ( - fact * M R * accel ) to unbalance
  if (cMass == 0) {
    // take advantage of lumped mass matrix
    double L = theCoordTransf->getInitialLength();
    double m = 0.5 * rho * L;

    Q(0) -= m * Raccel1(0);
    Q(1) -= m * Raccel1(1);
    Q(2) -= m * Raccel1(2);

    Q(6) -= m * Raccel2(0);
    Q(7) -= m * Raccel2(1);
    Q(8) -= m * Raccel2(2);
  } else {
    // use matrix vector multip. for consistent mass matrix
    static Vector Raccel(12);
    for (int i = 0; i < 6; i++) {
      Raccel(i) = Raccel1(i);
      Raccel(i + 6) = Raccel2(i);
    }
    Q.addMatrixVector(1.0, this->getMass(), Raccel, -1.0);
  }

  return 0;
}

const Vector &ModElasticBeam3d::getResistingForceIncInertia() {
  P = this->getResistingForce();

  // add the damping forces if rayleigh damping
  if (alphaM != 0.0 || betaK != 0.0 || betaK0 != 0.0 || betaKc != 0.0)
    P.addVector(1.0, this->getRayleighDampingForces(), 1.0);

  if (rho == 0.0)
    return P;

  // add inertia forces from element mass
  const Vector &accel1 = theNodes[0]->getTrialAccel();
  const Vector &accel2 = theNodes[1]->getTrialAccel();

  if (cMass == 0) {
    // take advantage of lumped mass matrix
    double L = theCoordTransf->getInitialLength();
    double m = 0.5 * rho * L;

    P(0) += m * accel1(0);
    P(1) += m * accel1(1);
    P(2) += m * accel1(2);

    P(6) += m * accel2(0);
    P(7) += m * accel2(1);
    P(8) += m * accel2(2);
  } else {
    // use matrix vector multip. for consistent mass matrix
    static Vector accel(12);
    for (int i = 0; i < 6; i++) {
      accel(i) = accel1(i);
      accel(i + 6) = accel2(i);
    }
    P.addMatrixVector(1.0, this->getMass(), accel, 1.0);
  }

  return P;
}

const Vector &ModElasticBeam3d::getResistingForce() {
  const Vector &v = theCoordTransf->getBasicTrialDisp();

  double L = theCoordTransf->getInitialLength();
  double oneOverL = 1.0 / L;
  double EoverL = E * oneOverL;
  double EAoverL = A * EoverL;        // EA/L
  double GJoverL = G * Jx * oneOverL; // GJ/L

  q(0) = EAoverL * v(0);
  q(5) = GJoverL * v(5);

  double EIzoverL2 = K44z * Iz * EoverL;   // 2EIz/L, modified
  double EIzoverL4 = K11z * Iz * EoverL;   // 4EIz/L, modified
  double EIzoverL6 = K33z * Iz * EoverL;   // 4EIz/L, modified
  q(1) = EIzoverL4 * v(1) + EIzoverL2 * v(2);
  q(2) = EIzoverL2 * v(1) + EIzoverL6 * v(2);

  double EIyoverL2 = K44y * Iy * EoverL;   // 2EIy/L, modified
  double EIyoverL4 = K11y * Iy * EoverL;   // 4EIy/L, modified
  double EIyoverL6 = K33y * Iy * EoverL;   // 4EIy/L, modified
  q(3) = EIyoverL4 * v(3) + EIyoverL2 * v(4);
  q(4) = EIyoverL2 * v(3) + EIyoverL6 * v(4);

  q(0) += q0[0];
  q(1) += q0[1];
  q(2) += q0[2];
  q(3) += q0[3];
  q(4) += q0[4];

  Vector p0Vec(p0, 5);

  //  opserr << q;

  P = theCoordTransf->getGlobalResistingForce(q, p0Vec);

  // subtract external load P = P - Q
  if (rho != 0)
    P.addVector(1.0, Q, -1.0);

  return P;
}

int ModElasticBeam3d::sendSelf(int cTag, Channel &theChannel) {
  int res = 0;

  static Vector data(22);

  int indx = 0;

  data(indx++) = A;
  data(indx++) = E;
  data(indx++) = G;
  data(indx++) = Jx;
  data(indx++) = Iy;
  data(indx++) = Iz;
  data(indx++) = K11y;
  data(indx++) = K33y;
  data(indx++) = K44y;
  data(indx++) = K11z;
  data(indx++) = K33z;
  data(indx++) = K44z;
  data(indx++) = rho;
  data(indx++) = cMass;
  data(indx++) = this->getTag();
  data(indx++) = connectedExternalNodes(0);
  data(indx++) = connectedExternalNodes(1);

  int dbTag = theCoordTransf->getDbTag();

  if (dbTag == 0) {
    dbTag = theChannel.getDbTag();
    if (dbTag != 0)
      theCoordTransf->setDbTag(dbTag);
  }

  data(indx++) = dbTag;
  data(indx++) = alphaM;
  data(indx++) = betaK;
  data(indx++) = betaK0;
  data(indx++) = betaKc;

  // Send the data vector
  res += theChannel.sendVector(this->getDbTag(), cTag, data);
  if (res < 0) {
    opserr << "ModElasticBeam3d::sendSelf -- could not send data Vector" << endln;
    return res;
  }

  // Ask the CoordTransf to send itself
  res += theCoordTransf->sendSelf(cTag, theChannel);
  if (res < 0) {
    opserr << "ModElasticBeam3d::sendSelf -- could not send CoordTransf" << endln;
    return res;
  }

  return res;
}

int ModElasticBeam3d::recvSelf(int cTag, Channel &theChannel,
                               FEM_ObjectBroker &theBroker) {
  int res = 0;
  static Vector data(22);

  res += theChannel.recvVector(this->getDbTag(), cTag, data);
  if (res < 0) {
    opserr << "ModElasticBeam3d::recvSelf -- could not receive data Vector" << endln;
    return res;
  }

  int indx = 0;

  A = data(indx++);
  E = data(indx++);
  G = data(indx++);
  Jx = data(indx++);
  Iy = data(indx++);
  Iz = data(indx++);
  K11y = data(indx++);
  K33y = data(indx++);
  K44y = data(indx++);
  K11z = data(indx++);
  K33z = data(indx++);
  K44z = data(indx++);
  rho = data(indx++);
  cMass = (int)data(indx++);
  this->setTag((int)data(indx++));
  connectedExternalNodes(0) = (int)data(indx++);
  connectedExternalNodes(1) = (int)data(indx++);
  theCoordTransf->setDbTag((int)data(indx++));
  alphaM = data(indx++);
  betaK = data(indx++);
  betaK0 = data(indx++);
  betaKc = data(indx++);

  res += theCoordTransf->recvSelf(cTag, theChannel, theBroker);
  if (res < 0) {
    opserr << "ModElasticBeam3d::recvSelf -- could not receive CoordTransf" << endln;
    return res;
  }

  return res;
}

void ModElasticBeam3d::Print(OPS_Stream &s, int flag) {
  this->getResistingForce();

  if (flag == -1) {
    int eleTag = this->getTag();
    s << "EL_BEAM\t" << eleTag << "\t";
    s << "\t" << connectedExternalNodes(0) << "\t" << connectedExternalNodes(1);
    s << "\t0\t0.0000000" << endln;
  }

  else if (flag < -1) {
    int counter = (flag + 1) * -1;
    int eleTag = this->getTag();
    const Vector &force = this->getResistingForce();

    double P, MZ1, MZ2, VY, MY1, MY2, VZ, T;
    double L = theCoordTransf->getInitialLength();
    double oneOverL = 1.0 / L;

    P = q(0);
    MZ1 = q(1);
    MZ2 = q(2);
    VY = (MZ1 + MZ2) * oneOverL;
    MY1 = q(3);
    MY2 = q(4);
    VZ = (MY1 + MY2) * oneOverL;
    T = q(5);

    s << "FORCE\t" << eleTag << "\t" << counter << "\t0";
    s << "\t" << -P + p0[0] << "\t" << VY + p0[1] << "\t" << -VZ + p0[3]
      << endln;
    s << "FORCE\t" << eleTag << "\t" << counter << "\t1";
    s << "\t" << P << ' ' << -VY + p0[2] << ' ' << VZ + p0[4] << endln;
    s << "MOMENT\t" << eleTag << "\t" << counter << "\t0";
    s << "\t" << -T << "\t" << MY1 << "\t" << MZ1 << endln;
    s << "MOMENT\t" << eleTag << "\t" << counter << "\t1";
    s << "\t" << T << ' ' << MY2 << ' ' << MZ2 << endln;
  }

  else if (flag == 2) {
    this->getResistingForce(); // in case linear algo

    static Vector xAxis(3);
    static Vector yAxis(3);
    static Vector zAxis(3);

    theCoordTransf->getLocalAxes(xAxis, yAxis, zAxis);

    s << "#ElasticBeamColumn3D" << endln;
    s << "#LocalAxis " << xAxis(0) << " " << xAxis(1) << " " << xAxis(2);
    s << " " << yAxis(0) << " " << yAxis(1) << " " << yAxis(2) << " ";
    s << zAxis(0) << " " << zAxis(1) << " " << zAxis(2) << endln;

    const Vector &node1Crd = theNodes[0]->getCrds();
    const Vector &node2Crd = theNodes[1]->getCrds();
    const Vector &node1Disp = theNodes[0]->getDisp();
    const Vector &node2Disp = theNodes[1]->getDisp();

    s << "#NODE " << node1Crd(0) << " " << node1Crd(1) << " " << node1Crd(2)
      << " " << node1Disp(0) << " " << node1Disp(1) << " " << node1Disp(2)
      << " " << node1Disp(3) << " " << node1Disp(4) << " " << node1Disp(5)
      << endln;

    s << "#NODE " << node2Crd(0) << " " << node2Crd(1) << " " << node2Crd(2)
      << " " << node2Disp(0) << " " << node2Disp(1) << " " << node2Disp(2)
      << " " << node2Disp(3) << " " << node2Disp(4) << " " << node2Disp(5)
      << endln;

    double N, Mz1, Mz2, Vy, My1, My2, Vz, T;
    double L = theCoordTransf->getInitialLength();
    double oneOverL = 1.0 / L;

    N = q(0);
    Mz1 = q(1);
    Mz2 = q(2);
    Vy = (Mz1 + Mz2) * oneOverL;
    My1 = q(3);
    My2 = q(4);
    Vz = -(My1 + My2) * oneOverL;
    T = q(5);

    s << "#END_FORCES " << -N + p0[0] << ' ' << Vy + p0[1] << ' ' << Vz + p0[3]
      << ' ' << -T << ' ' << My1 << ' ' << Mz1 << endln;
    s << "#END_FORCES " << N << ' ' << -Vy + p0[2] << ' ' << -Vz + p0[4] << ' '
      << T << ' ' << My2 << ' ' << Mz2 << endln;
  }

  if (flag == OPS_PRINT_CURRENTSTATE) {

    this->getResistingForce(); // in case linear algo

    s << "\nModElasticBeam3d: " << this->getTag() << endln;
    s << "\tConnected Nodes: " << connectedExternalNodes;
    s << "\tCoordTransf: " << theCoordTransf->getTag() << endln;
    s << "\tmass density:  " << rho << ", cMass: " << cMass << endln;
    double N, Mz1, Mz2, Vy, My1, My2, Vz, T;
    double L = theCoordTransf->getInitialLength();
    double oneOverL = 1.0 / L;

    N = q(0);
    Mz1 = q(1);
    Mz2 = q(2);
    Vy = (Mz1 + Mz2) * oneOverL;
    My1 = q(3);
    My2 = q(4);
    Vz = -(My1 + My2) * oneOverL;
    T = q(5);

    s << "\tEnd 1 Forces (P Mz Vy My Vz T): " << -N + p0[0] << ' ' << Mz1 << ' '
      << Vy + p0[1] << ' ' << My1 << ' ' << Vz + p0[3] << ' ' << -T << endln;
    s << "\tEnd 2 Forces (P Mz Vy My Vz T): " << N << ' ' << Mz2 << ' '
      << -Vy + p0[2] << ' ' << My2 << ' ' << -Vz + p0[4] << ' ' << T << endln;
  }

  if (flag == OPS_PRINT_PRINTMODEL_JSON) {
    s << "\t\t\t{";
    s << "\"name\": " << this->getTag() << ", ";
    s << "\"type\": \"ModElasticBeam3d\", ";
    s << "\"nodes\": [" << connectedExternalNodes(0) << ", "
      << connectedExternalNodes(1) << "], ";
    s << "\"E\": " << E << ", ";
    s << "\"G\": " << G << ", ";
    s << "\"A\": " << A << ", ";
    s << "\"Jx\": " << Jx << ", ";
    s << "\"Iy\": " << Iy << ", ";
    s << "\"Iz\": " << Iz << ", ";
    s << "\"K11y\": " << K11y << ", ";
    s << "\"K33y\": " << K33y << ", ";
    s << "\"K44y\": " << K44y << ", ";
    s << "\"K11z\": " << K11z << ", ";
    s << "\"K33z\": " << K33z << ", ";
    s << "\"K44z\": " << K44z << ", ";
    s << "\"massperlength\": " << rho << ", ";
    s << "\"crdTransformation\": \"" << theCoordTransf->getTag() << "\"}";
  }
}

int ModElasticBeam3d::displaySelf(Renderer &theViewer, int displayMode,
                                  float fact, const char **modes, int numMode) {
  static Vector v1(3);
  static Vector v2(3);

  theNodes[0]->getDisplayCrds(v1, fact, displayMode);
  theNodes[1]->getDisplayCrds(v2, fact, displayMode);
  float d1 = 0.0;
  float d2 = 0.0;
  int res = 0;

  if (displayMode > 0 && numMode == 0)
    res += theViewer.drawLine(v1, v2, d1, d1, this->getTag(), 0);
  else if (displayMode < 0)
    return theViewer.drawLine(v1, v2, 0.0, 0.0, this->getTag(), 0);

  if (numMode > 0) {
    // calculate q for potential need below
    this->getResistingForce();
  }

  for (int i = 0; i < numMode; i++) {

    const char *theMode = modes[i];
    if (strcmp(theMode, "axialForce") == 0) {
      d1 = q(0);
      d2 = q(0);
      ;

      res += theViewer.drawLine(v1, v2, d1, d1, this->getTag(), i);

    } else if (strcmp(theMode, "endMoments") == 0) {
      d1 = q(1);
      d2 = q(2);
      static Vector delta(3);
      delta = v2 - v1;
      delta /= 10;
      res += theViewer.drawPoint(v1 + delta, d1, this->getTag(), i);
      res += theViewer.drawPoint(v2 - delta, d2, this->getTag(), i);
    }
  }

  return res;
}

Response *ModElasticBeam3d::setResponse(const char **argv, int argc,
                                        OPS_Stream &output) {

  Response *theResponse = 0;

  output.tag("ElementOutput");
  output.attr("eleType", "ModElasticBeam3d");
  output.attr("eleTag", this->getTag());
  output.attr("node1", connectedExternalNodes[0]);
  output.attr("node2", connectedExternalNodes[1]);

  // global forces
  if (strcmp(argv[0], "force") == 0 || strcmp(argv[0], "forces") == 0 ||
      strcmp(argv[0], "globalForce") == 0 ||
      strcmp(argv[0], "globalForces") == 0) {

    output.tag("ResponseType", "Px_1");
    output.tag("ResponseType", "Py_1");
    output.tag("ResponseType", "Pz_1");
    output.tag("ResponseType", "Mx_1");
    output.tag("ResponseType", "My_1");
    output.tag("ResponseType", "Mz_1");
    output.tag("ResponseType", "Px_2");
    output.tag("ResponseType", "Py_2");
    output.tag("ResponseType", "Pz_2");
    output.tag("ResponseType", "Mx_2");
    output.tag("ResponseType", "My_2");
    output.tag("ResponseType", "Mz_2");

    theResponse = new ElementResponse(this, 2, P);

    // local forces
  } else if (strcmp(argv[0], "localForce") == 0 ||
             strcmp(argv[0], "localForces") == 0) {

    output.tag("ResponseType", "N_1");
    output.tag("ResponseType", "Vy_1");
    output.tag("ResponseType", "Vz_1");
    output.tag("ResponseType", "T_1");
    output.tag("ResponseType", "My_1");
    output.tag("ResponseType", "Mz_1");
    output.tag("ResponseType", "N_2");
    output.tag("ResponseType", "Vy_2");
    output.tag("ResponseType", "Vz_2");
    output.tag("ResponseType", "T_2");
    output.tag("ResponseType", "My_2");
    output.tag("ResponseType", "Mz_2");

    theResponse = new ElementResponse(this, 3, P);

    // basic forces
  } else if (strcmp(argv[0], "basicForce") == 0 ||
             strcmp(argv[0], "basicForces") == 0) {

    output.tag("ResponseType", "N");
    output.tag("ResponseType", "Mz_1");
    output.tag("ResponseType", "Mz_2");
    output.tag("ResponseType", "My_1");
    output.tag("ResponseType", "My_2");
    output.tag("ResponseType", "T");

    theResponse = new ElementResponse(this, 4, Vector(6));
  }
  // basic stiffness -
  else if (strcmp(argv[0], "basicStiffness") == 0) {

    output.tag("ResponseType", "N");
    output.tag("ResponseType", "Mz_1");
    output.tag("ResponseType", "Mz_2");
    output.tag("ResponseType", "My_1");
    output.tag("ResponseType", "My_2");
    output.tag("ResponseType", "T");

    theResponse = new ElementResponse(this, 19, Matrix(6, 6));

  } else if (strcmp(argv[0], "deformations") == 0 ||
             strcmp(argv[0], "basicDeformations") == 0) {

    output.tag("ResponseType", "eps");
    output.tag("ResponseType", "theta11");
    output.tag("ResponseType", "theta12");
    output.tag("ResponseType", "theta21");
    output.tag("ResponseType", "theta22");
    output.tag("ResponseType", "phi");
    theResponse = new ElementResponse(this, 5, Vector(6));
  }

  else if (strcmp(argv[0], "sectionX") == 0) {
    if (argc > 2) {
      float xL = atof(argv[1]);
      if (xL < 0.0)
        xL = 0.0;
      if (xL > 1.0)
        xL = 1.0;
      if (strcmp(argv[2], "forces") == 0) {
        theResponse = new ElementResponse(this, 6, Vector(6));
        Information &info = theResponse->getInformation();
        info.theDouble = xL;
      }
    }
  }

  output.endTag(); // ElementOutput

  if (theResponse == 0)
    theResponse = theCoordTransf->setResponse(argv, argc, output);

  return theResponse;
}

int ModElasticBeam3d::getResponse(int responseID, Information &eleInfo) {
  double N, V, M1, M2, T;
  double L = theCoordTransf->getInitialLength();
  double oneOverL = 1.0 / L;
  static Vector Sd(3);
  static Vector Res(12);
  Res = this->getResistingForce();
  static Vector s(6);
  static Matrix kb(6, 6);

  switch (responseID) {
  case 1: // stiffness
    return eleInfo.setMatrix(this->getTangentStiff());

  case 2: // global forces
    return eleInfo.setVector(Res);

  case 3: // local forces
    // Axial
    N = q(0);
    P(6) = N;
    P(0) = -N + p0[0];

    // Torsion
    T = q(5);
    P(9) = T;
    P(3) = -T;

    // Moments about z and shears along y
    M1 = q(1);
    M2 = q(2);
    P(5) = M1;
    P(11) = M2;
    V = (M1 + M2) * oneOverL;
    P(1) = V + p0[1];
    P(7) = -V + p0[2];

    // Moments about y and shears along z
    M1 = q(3);
    M2 = q(4);
    P(4) = M1;
    P(10) = M2;
    V = (M1 + M2) * oneOverL;
    P(2) = -V + p0[3];
    P(8) = V + p0[4];

    return eleInfo.setVector(P);

  case 4: // basic forces

    return eleInfo.setVector(q);

  case 5:
    return eleInfo.setVector(theCoordTransf->getBasicTrialDisp());

  case 6: {
    double xL = eleInfo.theDouble;
    double x = xL * L;

    s(0) = q(0) + wx * (L - x);
    s(1) = q(1) * (xL - 1.0) + q(2) * xL + 0.5 * wy * x * (x - L);
    s(2) = (q(1) + q(2)) / L + wy * (x - 0.5 * L);
    s(3) = q(3) * (xL - 1.0) + q(4) * xL - 0.5 * wz * x * (x - L);
    s(4) = (q(3) + q(4)) / L - wz * (x - 0.5 * L);
    s(5) = q(5);

    return eleInfo.setVector(s);
  }

  case 19: // basic stiffness
    kb.Zero();
    kb(0, 0) = E * A / L;
    kb(5, 5) = G * Jx / L;
    kb(1, 1) = kb(2, 2) = 4 * E * Iz / L;
    kb(1, 2) = kb(2, 1) = 2 * E * Iz / L;

    kb(3, 3) = kb(4, 4) = 4 * E * Iy / L;
    kb(3, 4) = kb(4, 3) = 2 * E * Iy / L;
    return eleInfo.setMatrix(kb);

  default:
    break;
  }

  return -1;
}

int ModElasticBeam3d::setParameter(const char **argv, int argc,
                                   Parameter &param) {
  if (argc < 1)
    return -1;

  // E of the beam interior
  if (strcmp(argv[0], "E") == 0) {
    param.setValue(E);
    return param.addObject(1, this);
  }
  // A of the beam interior
  if (strcmp(argv[0], "A") == 0) {
    param.setValue(A);
    return param.addObject(2, this);
  }
  // Iz of the beam interior
  if (strcmp(argv[0], "Iz") == 0) {
    param.setValue(Iz);
    return param.addObject(3, this);
  }
  // Iy of the beam interior
  if (strcmp(argv[0], "Iy") == 0) {
    param.setValue(Iy);
    return param.addObject(4, this);
  }
  // K11y of the beam interior
  if (strcmp(argv[0], "K11y") == 0) {
    param.setValue(K11y);
    return param.addObject(5, this);
  }
  // K33y of the beam interior
  if (strcmp(argv[0], "K33y") == 0) {
    param.setValue(K33y);
    return param.addObject(6, this);
  }
  // K44y of the beam interior
  if (strcmp(argv[0], "K44y") == 0) {
    param.setValue(K44y);
    return param.addObject(7, this);
  }
  // K11z of the beam interior
  if (strcmp(argv[0], "K11z") == 0) {
    param.setValue(K11z);
    return param.addObject(8, this);
  }
  // K33z of the beam interior
  if (strcmp(argv[0], "K33z") == 0) {
    param.setValue(K33z);
    return param.addObject(9, this);
  }
  // K44z of the beam interior
  if (strcmp(argv[0], "K44z") == 0) {
    param.setValue(K44z);
    return param.addObject(10, this);
  }
  // G of the beam interior
  if (strcmp(argv[0], "G") == 0) {
    param.setValue(G);
    return param.addObject(11, this);
  }
  // J of the beam interior
  if (strcmp(argv[0], "J") == 0) {
    param.setValue(Jx);
    return param.addObject(12, this);
  }

  return -1;
}

int ModElasticBeam3d::updateParameter(int parameterID, Information &info) {
  switch (parameterID) {
  case -1:
    return -1;
  case 1:
    E = info.theDouble;
    return 0;
  case 2:
    A = info.theDouble;
    return 0;
  case 3:
    Iz = info.theDouble;
    return 0;
  case 4:
    Iy = info.theDouble;
    return 0;
  case 5:
    K11y = info.theDouble;
    return 0;
  case 6:
    K33y = info.theDouble;
    return 0;
  case 7:
    K44y = info.theDouble;
    return 0;
  case 8:
    K11z = info.theDouble;
    return 0;
  case 9:
    K33z = info.theDouble;
    return 0;
  case 10:
    K44z = info.theDouble;
    return 0;
  case 11:
    G = info.theDouble;
    return 0;
  case 12:
    Jx = info.theDouble;
    return 0;
  default:
    return -1;
  }
}
