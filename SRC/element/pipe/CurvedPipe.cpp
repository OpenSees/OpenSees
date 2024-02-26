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
    : Pipe(), center(3), radius(0.0), theta0(0.0), tolWall(0.1) {}

CurvedPipe::CurvedPipe(int tag, int nd1, int nd2,
                       CrdTransf &theTransf, PipeMaterial &mat,
                       PipeSection &sect, const Vector &c, double to,
                       double pre, int cm, double tol)
    : Pipe(tag, ELE_TAG_CurvedPipe),
      center(3),
      radius(0.0),
      theta0(0.0),
      tolWall(tol) {
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

int CurvedPipe::flexibility(Matrix &fb) {
    // update section data
    if (updateSectionData() < 0) {
        opserr << "CurvedPipe::flexibility failed to update section "
                  "data\n";
        return;
    }

    // update material data
    if (updateMaterialData() < 0) {
        opserr << "CurvedPipe::flexibility failed to update material "
                  "data\n";
        return;
    }

    // material and section properties
    double A = ElasticBeam3d::A;
    double alphaV = theSect->ALFAV();
    double G = ElasticBeam3d::G;
    double E = ElasticBeam3d::E;
    double Iz = ElasticBeam3d::Iz;
    double Iy = ElasticBeam3d::Iy;
    double R = radius;
    double c = cos(theta0);
    double s = sin(theta0);
    double s2 = sin(theta0 / 2.0);
    double t0 = theta0;
    double R2 = R * R;
    double Ay = A / 1e6;
    double Az = A / 1e6;

    if (alphaV > 0) {
        Ay = A / alphaV;
        Az = A / alphaV;
    }

    // fb
    fb.resize(6, 6);
    fb.Zero();

    if (alphaV > 0) {
        fb(0, 0) = A * Ay * G * R2;
        fb(0, 0) *= 2 * t0 * c * c + t0 - 8 * s2 * c + s;
        fb(0, 0) += A * E * Iz * (t0 - s);
        fb(0, 0) += Ay * G * Iz * (t0 + s);
        fb(0, 0) /= 2 * A * Ay * E * G * Iz;
    }

    // fb(1, 0) = ;

    return 0;
}
int CurvedPipe::basicDeform(Vector &ub0) { return 0; }

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