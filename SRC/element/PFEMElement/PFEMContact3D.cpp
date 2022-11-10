/* ******************************************************************
***
**    OpenSees - Open System for Earthquake Engineering Simulation **
**          Pacific Earthquake Engineering Research Center **
** **
** **
** (C) Copyright 1999, The Regents of the University of California **
** All Rights Reserved. **
** **
** Commercial use of this program without express permission of the **
** University of California, Berkeley, is strictly prohibited.  See **
** file 'COPYRIGHT'  in main directory for information on usage and **
** redistribution,  and for a DISCLAIMER OF ALL WARRANTIES. **
** **
** Developed by: **
**   Frank McKenna (fmckenna@ce.berkeley.edu) **
**   Gregory L. Fenves (fenves@ce.berkeley.edu) **
**   Filip C. Filippou (filippou@ce.berkeley.edu) **
** **
** ******************************************************************
*/

// $Revision $
// $Date $

// Written: Minjie Zhu

#include "PFEMContact3D.h"

#include <Domain.h>
#include <Node.h>
#include <elementAPI.h>

#include <cmath>
#include <cstring>
#include <iostream>
#include <map>

Matrix PFEMContact3D::K;
Vector PFEMContact3D::P;

// for FEM_ObjectBroker, recvSelf must invoke
PFEMContact3D::PFEMContact3D()
    : Element(0, ELE_TAG_PFEMContact3D),
      ntags(8),
      nodes(8),
      Dc(0.0),
      E(0.0),
      rho(0.0),
      Ae(0.0),
      Fi(0.0),
      ndf(9),
      ndir(3) {}

// for object
PFEMContact3D::PFEMContact3D(int tag, int nd1, int nd2, int nd3,
                             int nd4, int nd5, int nd6, int nd7,
                             int nd8, double dc, double e, double r,
                             double nx, double ny, double nz,
                             double ae)
    : Element(tag, ELE_TAG_PFEMContact3D),
      ntags(8),
      nodes(8),
      Dc(dc),
      E(e),
      rho(r),
      Ae(ae),
      Fi(0.0),
      ndf(9),
      ndir(3) {
    ntags(0) = nd1;
    ntags(1) = nd2;
    ntags(2) = nd3;
    ntags(3) = nd4;
    ntags(4) = nd5;
    ntags(5) = nd6;
    ntags(6) = nd7;
    ntags(7) = nd8;

    ndir(0) = nx;
    ndir(1) = ny;
    ndir(2) = nz;
}

PFEMContact3D::~PFEMContact3D() {}

void PFEMContact3D::setDomain(Domain *theDomain) {
    DomainComponent::setDomain(theDomain);

    if (theDomain == 0) {
        return;
    }

    int ndm = OPS_GetNDM();

    int eletag = this->getTag();
    ndf[0] = 0;
    for (int i = 0; i < ntags.Size(); i++) {
        // get node
        nodes[i] = theDomain->getNode(ntags(i));
        if (nodes[i] == 0) {
            opserr << "WARNING: node " << ntags(i)
                   << " does not exist ";
            opserr << "in PFEMContact3D - setDomain() " << eletag
                   << "\n ";
            return;
        }
        if (nodes[i]->getNumberDOF() < ndm) {
            opserr << "WARNING: node " << ntags(i) << " ndf < ndm ";
            opserr << "in PFEMContact3D - setDomain() " << eletag
                   << "\n ";
            return;
        }
        ndf[i + 1] = ndf[i] + nodes[i]->getNumberDOF();
    }

    // impact velocity
    double vI = getVI();
    if (vI > 0) {
        Fi = vI * Ae * sqrt(E * rho);
    } else {
        Fi = 0;
    }
}

double PFEMContact3D::getP(double D) {
    if (inContact(D)) {
        return Fi;
    }
    return 0.0;
}

int PFEMContact3D::update() {
    // get D
    double D = getD();

    if (D < 1e-15) {
        opserr << "WARNING: D = " << D << "\n";
        return -1;
    }

    return 0;
}

const Matrix &PFEMContact3D::getMass() {
    int numdof = getNumDOF();

    K.resize(numdof, numdof);
    K.Zero();
    return K;
}

const Matrix &PFEMContact3D::getDamp() {
    int numdof = getNumDOF();

    K.resize(numdof, numdof);
    K.Zero();
    return K;
}

const Matrix &PFEMContact3D::getTangentStiff() {
    int numdof = getNumDOF();

    K.resize(numdof, numdof);
    K.Zero();

    return K;
}

const Matrix &PFEMContact3D::getInitialStiff() {
    return getTangentStiff();
}

int PFEMContact3D::addInertiaLoadToUnbalance(const Vector &accel) {
    return 0;
}

const Vector &PFEMContact3D::getResistingForce() {
    int numdof = getNumDOF();

    // get P
    P.resize(numdof);
    P.Zero();

    return P;
}

const Vector &PFEMContact3D::getResistingForceIncInertia() {
    getResistingForce();

    double D = getD();
    double pp = getP(D);

    for (int a = 0; a < ntags.Size(); ++a) {
        for (int i = 0; i < ndir.Size(); ++i) {
            if (a < 4) {
                P(ndf[a] + i) -= 0.25 * pp * ndir(i);
            } else {
                P(ndf[a] + i) += 0.25 * pp * ndir(i);
            }
        }
    }

    return P;
}

const char *PFEMContact3D::getClassType() const {
    return "PFEMContact3D";
}

int PFEMContact3D::sendSelf(int commitTag, Channel &theChannel) {
    return 0;
}

int PFEMContact3D::recvSelf(int commitTag, Channel &theChannel,
                            FEM_ObjectBroker &theBroker) {
    return 0;
}

void PFEMContact3D::Print(OPS_Stream &s, int flag) {
    s << "PFEMContact3D: " << this->getTag() << endln;
}

int PFEMContact3D::displaySelf(Renderer &theViewer, int displayMode,
                               float fact, const char **displayModes,
                               int numModes) {
    return 0;
}

double PFEMContact3D::getVI() {
    double vI = 0.0;

    for (int a = 0; a < ntags.Size(); ++a) {
        // get velocity
        const Vector &vel = nodes[a]->getTrialVel();
        double vn = 0.0;
        for (int i = 0; i < ndir.Size(); ++i) {
            vn += vel(i) * ndir(i);
        }
        if (a < 4) {
            vI -= vn;
        } else {
            vI += vn;
        }
    }

    return vI / 4.0;
}

bool PFEMContact3D::inContact(double D) {
    // distance > critical value
    if (Dc <= D) {
        return false;
    }

    return true;
}

double PFEMContact3D::getD() {
    // nodal coordinates projected to ndir
    double D = 0.0;

    // get projected coords
    for (int a = 0; a < ntags.Size(); ++a) {
        const Vector &coords = nodes[a]->getCrds();
        const Vector &disp = nodes[a]->getTrialDisp();

        double xn = 0.0;
        for (int i = 0; i < ndir.Size(); ++i) {
            xn += (coords(i) + disp(i)) * ndir(i);
        }

        if (a < 4) {
            D -= xn;
        } else {
            D += xn;
        }
    }

    if (D < 0) D = -D;

    return D / 4.0;
}