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

// Minjie
#include "CurvedPipe.h"

#include <CrdTransf.h>

#include <cmath>

// 20 gauss points in (wi, xi)
std::vector<double> CurvedPipe::gaussPts = {
    0.1527533871307258,  -0.0765265211334973, 0.1527533871307258,
    0.0765265211334973,  0.1491729864726037,  -0.2277858511416451,
    0.1491729864726037,  0.2277858511416451,  0.1420961093183820,
    -0.3737060887154195, 0.1420961093183820,  0.3737060887154195,
    0.1316886384491766,  -0.5108670019508271, 0.1316886384491766,
    0.5108670019508271,  0.1181945319615184,  -0.6360536807265150,
    0.1181945319615184,  0.6360536807265150,  0.1019301198172404,
    -0.7463319064601508, 0.1019301198172404,  0.7463319064601508,
    0.0832767415767048,  -0.8391169718222188, 0.0832767415767048,
    0.8391169718222188,  0.0626720483341091,  -0.9122344282513259,
    0.0626720483341091,  0.9122344282513259,  0.0406014298003869,
    -0.9639719272779138, 0.0406014298003869,  0.9639719272779138,
    0.0176140071391521,  -0.9931285991850949, 0.0176140071391521,
    0.9931285991850949};

// function object for

void *OPS_CurvedPipeElement() {
    // check inputs
    if (OPS_GetNumRemainingInputArgs() < 10) {
        opserr << "Invalid #args,  want: element CurvedPipe "
                  "tag? nd1? nd2? transfTag? pipeMatTag? pipeSecTag?"
                  "xC? yC? zC?"
                  "<-T0 T0? -p p? -cMass? -tolWall? tolWall?>\n";
        return 0;
    }

    // get tag
    int iData[6];
    int numData = 6;
    if (OPS_GetIntInput(&numData, iData) < 0) {
        opserr << "WARNING invalid integer input for curved pipe "
                  "element\n";
        return 0;
    }

    // get center
    Vector center(3);
    numData = 3;
    if (OPS_GetDoubleInput(&numData, &center(0)) < 0) {
        opserr << "WARNING invalid center or radius input for curved "
                  "pipe element\n";
        return 0;
    }

    // get data
    double T0 = 0.0, pressure = 0.0;
    double tolWall = 0.1;
    int cMass = 0;
    numData = 1;
    while (OPS_GetNumRemainingInputArgs() > 0) {
        const char *theType = OPS_GetString();
        if (strcmp(theType, "-T0") == 0) {
            if (OPS_GetNumRemainingInputArgs() > 0) {
                if (OPS_GetDoubleInput(&numData, &T0) < 0) {
                    opserr << "WARNING: failed to read T0\n";
                    return 0;
                }
            }
        } else if (strcmp(theType, "-p") == 0) {
            if (OPS_GetNumRemainingInputArgs() > 0) {
                if (OPS_GetDoubleInput(&numData, &pressure) < 0) {
                    opserr << "WARNING: failed to read internal "
                              "pressure\n";
                    return 0;
                }
            }
        } else if (strcmp(theType, "-tolWall") == 0) {
            if (OPS_GetNumRemainingInputArgs() > 0) {
                if (OPS_GetDoubleInput(&numData, &tolWall) < 0) {
                    opserr << "WARNING: failed to read fraction of "
                              "wall thickness "
                              "to be used for dimensional tolerance "
                              "tests\n";
                    return 0;
                }
            }
            if (tolWall < 0 || tolWall > 1) {
                opserr << "WARNING: tolWall < 0 or > 1\n";
                return 0;
            }
        } else if (strcmp(theType, "-cMass") == 0) {
            cMass = 1;
        }
    }

    auto *theSect = dynamic_cast<PipeSection *>(
        OPS_getSectionForceDeformation(iData[5]));
    if (theSect == 0) {
        opserr << "WARNING: section " << iData[5]
               << " is not found or not a curved pipe section\n";
        return 0;
    }

    auto *theMat = dynamic_cast<PipeMaterial *>(
        OPS_getUniaxialMaterial(iData[4]));
    if (theMat == 0) {
        opserr << "WARNING: uniaxialMaterial " << iData[4]
               << " is not found or not a curved pipe material\n";
        return 0;
    }

    auto *theTrans = OPS_getCrdTransf(iData[3]);
    if (theTrans == 0) {
        opserr << "WARNING: CrdTransf " << iData[3]
               << " is not found\n";
        return 0;
    }

    auto *ele = new CurvedPipe(iData[0], iData[1], iData[2],
                               *theTrans, *theMat, *theSect, center,
                               T0, pressure, cMass);

    return ele;
}

CurvedPipe::CurvedPipe()
    : Pipe(),
      center(3),
      radius(0.0),
      theta0(0.0),
      tolWall(0.1),
      wa(0.0),
      wy(0.0),
      wz(0.0) {}

CurvedPipe::CurvedPipe(int tag, int nd1, int nd2,
                       CrdTransf &theTransf, PipeMaterial &mat,
                       PipeSection &sect, const Vector &c, double to,
                       double pre, int cm, double tol)
    : Pipe(tag, ELE_TAG_CurvedPipe),
      center(3),
      radius(0.0),
      theta0(0.0),
      tolWall(tol),
      wa(0.0),
      wy(0.0),
      wz(0.0) {
    if (Pipe::createPipe(nd1, nd2, theTransf, mat, sect, cm, 0, 0) <
        0) {
        opserr << "WARNING: failed to create curved pipe element\n";
        exit(-1);
    }
    if (center.Size() != 3) {
        center.resize(3);
    }
    for (int i = 0; i < c.Size(); ++i) {
        if (i >= 3) {
            break;
        }
        center(i) = c(i);
    }
    if (radius <= 0) {
        opserr << "WARNING: radius <= 0\n";
        exit(-1);
    }
}

CurvedPipe::~CurvedPipe() {}

const char *CurvedPipe::getClassType(void) const {
    return "CurvedPipe";
};

void CurvedPipe::setDomain(Domain *theDomain) {
    // set domain
    this->Pipe::setDomain(theDomain);

    // compute theta0
    if (this->getTheta0() < 0) {
        opserr << "WARNING: failed to compute theta0\n";
        exit(-1);
    }
}

void CurvedPipe::zeroLoad(void) {
    // update section data
    if (updateSectionData() < 0) {
        opserr << "CurvedPipe::setDomain failed to update section "
                  "data\n";
        return;
    }

    // update material data
    if (updateMaterialData() < 0) {
        opserr << "CurvedPipe::setDomain failed to update material "
                  "data\n";
        return;
    }

    this->ElasticBeam3d::zeroLoad();

    // due to thermal
    double temp = aveTemp();
    if (temp > 0) {
        ElasticBeam3d::q0[0] -= E * A * alp * temp;
    }

    // due to internal pressure
    if (pressure != 0) {
        double dout = theSect->DOUT();
        double thk = theSect->WALL();

        ElasticBeam3d::q0[0] -=
            0.25 * pressure * (dout - thk) * (1. - 2 * nu) * A / thk;
    }
}

int CurvedPipe::getTheta0() {
    // compute radius
    const auto &crds1 = theNodes[0]->getCrds();
    const auto &crds2 = theNodes[1]->getCrds();

    double R1 = (center - crds1).Norm();
    double R2 = (center - crds2).Norm();
    radius = (R1 + R2) / 2.0;

    double thk = theSect->WALL();

    if (fabs(R1 - R2) > tolWall * thk) {
        opserr << "WARNING: the computed radius from node I is "
                  "different to "
                  "the one computed from node J more than "
               << tolWall << " * wall thickness\n";
    }

    // compute theta0
    double Lhalf = theCoordTransf->getInitialLength() / 2.0;
    theta0 = 2.0 * asin(radius / Lhalf);

    return 0;
}

void CurvedPipe::bx(double theta, Matrix &mat) {
    mat.resize(6, 6);
    mat.Zero();
    double c = cos(theta);
    double s = sin(theta);
    double L = theCoordTransf->getInitialLength();
    double R = radius;
    double H = R * (c - cos(theta0));
    double x = L / 2.0 + R * s;
    double invL = 1.0 / L;
    mat(0, 0) = c;
    mat(0, 1) = -s * invL;
    mat(0, 2) = -s * invL;
    mat(1, 0) = -H;
    mat(1, 1) = R * s * invL - 0.5;
    mat(1, 2) = R * s * invL + 0.5;
    mat(2, 3) = x * invL - 1;
    mat(2, 4) = x * invL;
    mat(3, 5) = 1;
    mat(4, 0) = s;
    mat(4, 1) = c * invL;
    mat(4, 2) = c * invL;
    mat(4, 3) = invL;
    mat(4, 4) = invL;
}

void CurvedPipe::Spx(double theta, Vector &vec) {
    vec.resize(6);
    vec.Zero();
    double ct = cos(theta);
    double ct0 = cos(theta0);
    double st = sin(theta);
    double s2t = sin(2 * theta);
    double st2 = st * st;
    double L = theCoordTransf->getInitialLength();
    double L2 = L * L;
    double R = radius;
    double R2 = R * R;
    vec(0) = L * wa * ct / 2.0 - R * wa * s2t / 2.0 - R * wy * st2;
    vec(1) = -L2 * wy / 8.0 - L * R * wa * ct / 2.0 +
             L * R * wa * ct0 / 2.0 - R2 * wa * st * ct0 +
             R2 * wa * s2t / 2.0 + R2 * wy * st2 / 2.0;
    vec(2) = wz * (-L2 + 4 * R2 * st2) / 8.0;
    vec(4) = st * (L * wa - 2 * R * wa * st + 2 * R * wy * ct) / 2.0;
    vec(5) = R * wz * st;
}

void CurvedPipe::fs(double theta, Matrix &mat) {
    mat.resize(6, 6);
    mat.Zero();

    mat(0, 0) = 1.0 / (E * A);
    mat(1, 1) = 1.0 / (E * Iz);
    mat(2, 2) = 1.0 / (E * Iy);
    mat(3, 3) = 1.0 / (G * Jx);

    double alphaV = theSect->ALFAV();
    if (alphaV > 0) {
        double Ay = A / alphaV;
        double Az = A / alphaV;
        mat(4, 4) = 1.0 / (G * Ay);
        mat(5, 5) = 1.0 / (G * Az);
    }
}

void CurvedPipe::fb(double theta, Matrix &mat) {
    mat.resize(6, 6);
    mat.Zero();

    Matrix bxmat, fsmat;
    bx(theta, bxmat);
    fs(theta, fsmat);
    mat.addMatrixTripleProduct(0.0, bxmat, fsmat, 1.0);
}

void CurvedPipe::ubno(double theta, Vector &vec) {
    vec.resize(6);
    vec.Zero();

    Matrix fsmat, bxmat;
    Vector spxvec;
    fs(theta, fsmat);
    bx(theta, bxmat);
    Spx(theta, spxvec);

    Vector temp(6);
    temp.addMatrixVector(0.0, fsmat, spxvec, 1.0);
    vec.addMatrixTransposeVector(0.0, bxmat, temp, 1.0);
}

int CurvedPipe::kb(double theta, Matrix &mat, Vector &vec) {
    Matrix fbmat;
    Vector ubnovec;

    integrateGauss(-theta0 / 2.0, theta0 / 2.0, fbmat, ubnovec);

    if (fbmat.Invert(mat) < 0) {
        return -1;
    }

    vec.resize(6);
    vec.addMatrixVector(0.0, mat, ubnovec, -1.0);

    return 0;
}

void CurvedPipe::integrateGauss(double a, double b, Matrix &resm,
                                Vector &resv) {
    double p = (b - a) / 2.0;
    double q = (b + a) / 2.0;
    resm.resize(6, 6);
    resm.Zero();
    resv.resize(6);
    resv.Zero();
    for (int i = 0; i < (int)gaussPts.size() / 2; ++i) {
        double xi = gaussPts[2 * i + 1] * p + q;
        double wi = gaussPts[2 * i] * p;
        Matrix mat;
        Vector vec;
        fb(xi, mat);
        resm.addMatrix(1.0, mat, wi);
        ubno(xi, vec);
        resv.addVector(1.0, vec, wi);
    }
}