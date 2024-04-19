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

// $Revision $
// $Date $

// Written: Minjie Zhu


#include "PFEMContact2D.h"
#include <elementAPI.h>
#include <Domain.h>
#include <Node.h>
#include <cmath>
#include <cstring>
#include <map>
#include <iostream>

Matrix PFEMContact2D::K;
Vector PFEMContact2D::P;

// for FEM_ObjectBroker, recvSelf must invoke
PFEMContact2D::PFEMContact2D()
        : Element(0, ELE_TAG_PFEMContact2D), ntags(3),
          nodes(3), kdoverAd(0.0), mu(0.1), beta(0.3), Dc(0.0),
          alpha(0.5), E(0.0), rho(0.0), ndf(4), signvt0(0), F0(0) {
}

// for object
PFEMContact2D::PFEMContact2D(int tag, int nd1, int nd2, int nd3,
                             double k, double t, double m,
                             double b, double dc, double a,
                             double e, double r)
        : Element(tag, ELE_TAG_PFEMContact2D), ntags(3),
          nodes(3), kdoverAd(k), thk(t), mu(m), beta(b), Dc(dc),
          alpha(a), E(e), rho(r), ndf(4), signvt0(0), F0(0) {
    ntags(0) = nd1;
    ntags(1) = nd2;
    ntags(2) = nd3;
}


PFEMContact2D::~PFEMContact2D() {
}

void
PFEMContact2D::setDomain(Domain *theDomain) {
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
            opserr << "WARNING: node " << ntags(i) << " does not exist ";
            opserr << "in PFEMContact2D - setDomain() " << eletag << "\n ";
            return;
        }
        if (nodes[i]->getNumberDOF() < ndm) {
            opserr << "WARNING: node " << ntags(i) << " ndf < ndm ";
            opserr << "in PFEMContact2D - setDomain() " << eletag << "\n ";
            return;
        }
        ndf[i + 1] = ndf[i] + nodes[i]->getNumberDOF();
    }

    // impact velocityies
    Vector vn, vj;
    getV(vn, signvt0, vj);
    F0 = (vn(0) + vn(1)) / 2.0 - vn(2);

    // get kdoverAd
    double A, B, C, dx, dy, L;
    double x1, y1, x2, y2, x3, y3;
    getLine(A, B, C, dx, dy, x1, y1, x2, y2, x3, y3, L);

    kdoverAd *= L * thk;
    if (F0 > 0) {
        F0 *= L * thk * sqrt(E * rho);
    } else {
        F0 = 0;
    }
}

// return current distances between node 3 and edge 1-2
double
PFEMContact2D::getLine(double &A, double &B, double &C,
                       double &dx, double &dy,
                       double &x1, double &y1,
                       double &x2, double &y2,
                       double &x3, double &y3, double &L) {
    const Vector &coord1 = nodes[0]->getCrds();
    const Vector &disp1 = nodes[0]->getTrialDisp();
    const Vector &coord2 = nodes[1]->getCrds();
    const Vector &disp2 = nodes[1]->getTrialDisp();
    const Vector &coord3 = nodes[2]->getCrds();
    const Vector &disp3 = nodes[2]->getTrialDisp();

    x1 = coord1(0) + disp1(0);
    x2 = coord2(0) + disp2(0);
    x3 = coord3(0) + disp3(0);

    y1 = coord1(1) + disp1(1);
    y2 = coord2(1) + disp2(1);
    y3 = coord3(1) + disp3(1);

    dx = x2 - x1;
    dy = y2 - y1;
    L = sqrt(dx * dx + dy * dy);

    A = -dy / L;
    B = dx / L;
    C = (x1 * y2 - x2 * y1) / L;

    return A * x3 + B * y3 + C;
}

void
PFEMContact2D::getdL(double L, double dx, double dy, Vector &dL) {
    dL.resize(6);
    dL.Zero();
    dL(0) = -dx;
    dL(1) = -dy;
    dL(2) = dx;
    dL(3) = dy;
    dL /= L;
}

void
PFEMContact2D::getdA(double L, double A, const Vector &dL, Vector &dA) {
    dA = dL;
    dA *= -A / L;
    dA(1) += 1.0 / L;
    dA(3) += -1.0 / L;
}

void
PFEMContact2D::getdB(double L, double B, const Vector &dL, Vector &dB) {
    dB = dL;
    dB *= -B / L;
    dB(0) += -1.0 / L;
    dB(2) += 1.0 / L;
}

void
PFEMContact2D::getdC(double L, double C, double x1, double y1,
                     double x2, double y2, const Vector &dL,
                     Vector &dC) {
    dC = dL;
    dC *= -C / L;
    dC(0) += y2 / L;
    dC(1) += -x2 / L;
    dC(2) += -y1 / L;
    dC(3) += x1 / L;
}

void
PFEMContact2D::getdD(double A, double B, double x3, double y3,
                     const Vector &dA,
                     const Vector &dB,
                     const Vector &dC,
                     Vector &dD) {
    dD = dC;
    dD.addVector(1.0, dB, y3);
    dD.addVector(1.0, dA, x3);
    dD(4) += A;
    dD(5) += B;
}

double
PFEMContact2D::getP(double D) {

    if (inContact(D)) {
        double p = alpha * kdoverAd * (Dc - D);
        if (p < F0) {
            return p;
        } else {
            return F0;
        }
    }
    return 0.0;

}

void
PFEMContact2D::getdP(const Vector &dD, double D, Vector &dP) {
    dP = dD;
    if (inContact(D)) {
        double p = alpha * kdoverAd * (Dc - D);
        if (p < F0) {
            dP *= -kdoverAd * alpha;
            return;
        } else {
            dP.Zero();
            return;
        }
    } else {
        dP.Zero();
    }
}

int
PFEMContact2D::update() {


    // get D
    double A, B, C, dx, dy, L;
    double x1, y1, x2, y2, x3, y3;
    double D = getLine(A, B, C, dx, dy, x1, y1, x2, y2, x3, y3, L);

    if (fabs(L) < 1e-15 || D < 0) {
        opserr << "WARNING: L = " << L << " or D = " << D << "\n";
        return -1;
    }

    return 0;
}

const Matrix &
PFEMContact2D::getMass() {

    int numdof = getNumDOF();

    K.resize(numdof, numdof);
    K.Zero();
    return K;
}

const Matrix &
PFEMContact2D::getDamp() {

    int numdof = getNumDOF();

    K.resize(numdof, numdof);
    K.Zero();

    // get D
    double A, B, C, dx, dy, L;
    double x1, y1, x2, y2, x3, y3;
    getLine(A, B, C, dx, dy, x1, y1, x2, y2, x3, y3, L);

    for (int a = 0; a < ntags.Size(); ++a) {
        K(ndf[a], ndf[a]) += mu * A * A;
        K(ndf[a], ndf[a] + 1) += mu * A * B;
        K(ndf[a] + 1, ndf[a]) += mu * A * B;
        K(ndf[a] + 1, ndf[a] + 1) += mu * B * B;
    }

    return K;
}

const Matrix &
PFEMContact2D::getTangentStiff() {

    int numdof = getNumDOF();

    // get D
    double A, B, C, dx, dy, L;
    double x1, y1, x2, y2, x3, y3;
    double D = getLine(A, B, C, dx, dy, x1, y1, x2, y2, x3, y3, L);

    // get derivatives
    Vector dL, dA, dB, dC, dD, dP, vn, vj;
    getdL(L, dx, dy, dL);
    getdA(L, A, dL, dA);
    getdB(L, B, dL, dB);
    getdC(L, C, x1, y1, x2, y2, dL, dC);
    getdD(A, B, x3, y3, dA, dB, dC, dD);
    getdP(dD, D, dP);
    double vt;
    getV(vn, vt, vj);
    double pp = getP(D);

    // get K
    K.resize(numdof, numdof);
    K.Zero();

    double factor[] = {0.5, 0.5, 0.5, 0.5, -1.0, -1.0};
    double factor2[] = {-0.5, 0.5, -0.5, 0.5, 1, -1};
    for (int a = 0; a < ntags.Size(); ++a) {
        for (int b = 0; b < ntags.Size(); ++b) {

            // P
            K(ndf[a], ndf[b]) += factor[2 * a] * A * dP[2 * b];
            K(ndf[a], ndf[b] + 1) += factor[2 * a] * A * dP[2 * b + 1];
            K(ndf[a] + 1, ndf[b]) += factor[2 * a + 1] * B * dP[2 * b];
            K(ndf[a] + 1, ndf[b] + 1) += factor[2 * a + 1] * B * dP[2 * b + 1];

            K(ndf[a], ndf[b]) += pp * factor[2 * a] * dA(2 * b);
            K(ndf[a], ndf[b] + 1) += pp * factor[2 * a] * dA(2 * b + 1);
            K(ndf[a] + 1, ndf[b]) += pp * factor[2 * a + 1] * dB(2 * b);
            K(ndf[a] + 1, ndf[b] + 1) += pp * factor[2 * a + 1] * dB(2 * b + 1);

            // friction
            K(ndf[a], ndf[b]) += beta * signvt0 * factor2[2 * a] * B * dP[2 * b];
            K(ndf[a], ndf[b] + 1) += beta * signvt0 * factor2[2 * a] * B * dP[2 * b + 1];
            K(ndf[a] + 1, ndf[b]) += beta * signvt0 * factor2[2 * a + 1] * A * dP[2 * b];
            K(ndf[a] + 1, ndf[b] + 1) += beta * signvt0 * factor2[2 * a + 1] * A * dP[2 * b + 1];

            K(ndf[a], ndf[b]) += beta * pp * signvt0 * factor2[2 * a] * dB(2 * b);
            K(ndf[a], ndf[b] + 1) += beta * pp * signvt0 * factor2[2 * a] * dB(2 * b + 1);
            K(ndf[a] + 1, ndf[b]) += beta * pp * signvt0 * factor2[2 * a + 1] * dA(2 * b);
            K(ndf[a] + 1, ndf[b] + 1) += beta * pp * signvt0 * factor2[2 * a + 1] * dA(2 * b + 1);

            // related to damping
            K(ndf[a], ndf[b]) += mu * (2 * vj[2 * a] * A * dA(2 * b) +
                                       vj[2 * a + 1] * A * dB(2 * b) +
                                       vj[2 * a + 1] * B * dA(2 * b));
            K(ndf[a], ndf[b] + 1) += mu * (vj[2 * a] * A * dA(2 * b + 1) +
                                           vj[2 * a + 1] * A * dB(2 * b + 1) +
                                           vj[2 * a + 1] * B * dA(2 * b + 1));
            K(ndf[a] + 1, ndf[b]) += mu * (2 * vj[2 * a + 1] * B * dB(2 * b) +
                                           vj[2 * a] * A * dB(2 * b) +
                                           vj[2 * a] * B * dA(2 * b));
            K(ndf[a] + 1, ndf[b] + 1) += mu * (2 * vj[2 * a + 1] * B * dB(2 * b + 1) +
                                               vj[2 * a] * A * dB(2 * b + 1) +
                                               vj[2 * a] * B * dA(2 * b + 1));
        }
    }

    return K;
}


const Matrix &
PFEMContact2D::getInitialStiff() {
    return getTangentStiff();
}

int
PFEMContact2D::addInertiaLoadToUnbalance(const Vector &accel) {
    return 0;
}

const Vector &
PFEMContact2D::getResistingForce() {

    int numdof = getNumDOF();

    // get D
    double A, B, C, dx, dy, L;
    double x1, y1, x2, y2, x3, y3;
    double D = getLine(A, B, C, dx, dy, x1, y1, x2, y2, x3, y3, L);

    // get derivatives
    double pp = getP(D);

    // get P
    P.resize(numdof);
    P.Zero();

    double factor[] = {0.5, 0.5, 0.5, 0.5, -1.0, -1.0};

    for (int a = 0; a < ntags.Size(); ++a) {
        P(ndf[a]) += factor[2 * a] * A * pp;
        P(ndf[a] + 1) += factor[2 * a + 1] * B * pp;
    }

    return P;
}

const Vector &
PFEMContact2D::getResistingForceIncInertia() {

    getResistingForce();

    double A, B, C, dx, dy, L;
    double x1, y1, x2, y2, x3, y3;
    double D = getLine(A, B, C, dx, dy, x1, y1, x2, y2, x3, y3, L);

    Vector vn, vj;
    double vt;
    getV(vn, vt, vj);
    double pp = getP(D);

    double factor[] = {-0.5, 0.5, -0.5, 0.5, 1, -1};

    for (int a = 0; a < ntags.Size(); ++a) {
        P(ndf[a]) += mu * vn(a) * A;
        P(ndf[a] + 1) += mu * vn(a) * B;

        P(ndf[a]) += beta * pp * signvt0 * factor[2 * a] * B;
        P(ndf[a] + 1) += beta * pp * signvt0 * factor[2 * a + 1] * A;

    }

    return P;
}


const char *
PFEMContact2D::getClassType() const {
    return "PFEMContact2D";
}

int
PFEMContact2D::sendSelf(int commitTag, Channel &theChannel) {
    return 0;
}

int
PFEMContact2D::recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker) {
    return 0;
}

void
PFEMContact2D::Print(OPS_Stream &s, int flag) {
    s << "PFEMContact2D: " << this->getTag() << endln;
}

int
PFEMContact2D::displaySelf(Renderer &theViewer, int displayMode, float fact, const char **displayModes, int numModes) {
    return 0;
}

void
PFEMContact2D::getV(Vector &vn, double &vt, Vector &vj) {

    // get D
    double A, B, C, dx, dy, L;
    double x1, y1, x2, y2, x3, y3;
    getLine(A, B, C, dx, dy, x1, y1, x2, y2, x3, y3, L);

    vn.resize(ntags.Size());
    vn.Zero();
    Vector vit(ntags.Size());
    vj.resize(ntags.Size() * 2);
    vj.Zero();

    for (int a = 0; a < ntags.Size(); ++a) {

        // get velocity
        const Vector &vel = nodes[a]->getTrialVel();
        vn(a) = vel(0) * A + vel(1) * B;
        vit(a) = vel(0) * B - vel(1) * A;
        vj(2 * a) = vel(0);
        vj(2 * a + 1) = vel(1);
    }

    vt = vit(2) - (vit(0) + vit(1)) / 2.0;
}

bool
PFEMContact2D::inContact(double D) {

    // distance > critical value
    if (Dc <= D) {
        return false;
    }

//    // normal velocity
//    Vector vn, vt;
//    getV(vn, vt);
//
//    // average velocity of line
//    double vl = (vn(0) + vn(1)) / 2.0;
//    double vl0 = (vn0(0) + vn0(1)) / 2.0;
//
//    // if in impact
//    if (vl == 0 && vn(2) == 0) {
//        return false;
//    } else if (vl >= 0 && vn(2) <= 0) {
//        return true;
//    } else if (vl >= 0 && vn(2) >= 0 && vl > vn(2)) {
//        return true;
//    } else if (vl <= 0 && vn(2) <= 0 && vl < vn(2)) {
//        return true;
//    }
//
//    // if not in impact
//    if (fabs(vn(2)) <= mu * fabs(vn0(2)) ||
//        fabs(vl) <= mu * fabs(vl0)) {
//        return true;
//    }

    return true;
}