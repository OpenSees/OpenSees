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

// Written: Minjie
//
// Description: material for pipe element

#include <PipeMaterial.h>

void *OPS_PipeMaterial(void) {
    // line: 9233
    // check inputs
    if (OPS_GetNumRemainingInputArgs() < 6) {
        opserr << "Invalid #args,  want: uniaxialMaterial Pipe "
                  "tag? nt? T1? E1? xnu1? alpT1? <T2? E2? xnu2? alpT2? "
                  "...>\n";
        return 0;
    }

    // get tag and nt
    int iData[2];
    int numData = 2;
    if (OPS_GetIntInput(&numData, iData) < 0) {
        opserr << "WARNING invalid tag or nt for pipe material\n";
        return 0;
    }

    int tag = iData[0];
    int nt = iData[1];
    if (nt < 1) {
        nt = 1;
    }

    if (OPS_GetNumRemainingInputArgs() < 4 * nt) {
        opserr << "WARNING: not enough inputs for " << nt
               << " temperature points\n";
        return 0;
    }

    // create material
    auto *mat = new PipeMaterial(tag);

    // get temperature points
    double prevT = 0.0;
    for (int i = 0; i < nt; ++i) {
        numData = 4;
        double ddata[4];
        if (OPS_GetDoubleInput(&numData, ddata) < 0) {
            opserr << "WARNING: invalid input for " << i + 1
                   << " th temperature point\n";
            delete mat;
            return 0;
        }
        if (ddata[1] <= 0) {
            opserr << "WARNING: E cannot <= 0\n";
            delete mat;
            return 0;
        }
        if (ddata[2] <= 0) {
            opserr << "WARNING: poisson's ratio cannot <= 0\n";
            delete mat;
            return 0;
        }
        if (ddata[3] <= 0) {
            opserr << "WARNING: thermal expansion coefficient cannot "
                      "<= 0\n";
            delete mat;
            return 0;
        }
        if (i > 0) {
            if (ddata[0] <= prevT) {
                opserr << "WARNING: the temperature points should be "
                          "in ascending order\n";
                delete mat;
                return 0;
            }
        }
        mat->addPoint(ddata[0], ddata[1], ddata[2], ddata[3]);
        prevT = ddata[0];
    }
    return mat;
}

PipeMaterialTemperaturePoint::PipeMaterialTemperaturePoint(double t,
                                                           double e,
                                                           double nu,
                                                           double a)
    : T(t), E(e), xnu(nu), alp(a) {}
PipeMaterialTemperaturePoint::PipeMaterialTemperaturePoint()
    : T(0), E(0), xnu(0), alp(0) {}

PipeMaterial::PipeMaterial(int tag)
    : UniaxialMaterial(tag, MAT_TAG_PipeMaterial) {}
PipeMaterial::PipeMaterial()
    : UniaxialMaterial(0, MAT_TAG_PipeMaterial) {}
PipeMaterial::~PipeMaterial() {}

const char *PipeMaterial::getClassType(void) const {
    return "PipeMaterial";
}

int PipeMaterial::setTrialStrain(double strain, double strainRate) {
    // trialStrain = strain;
    // trialStrainRate = strainRate;
    return 0;
}
int PipeMaterial::setTrialStrain(double strain, double temperature,
                                 double strainRate) {
    return 0;
}
int PipeMaterial::setTrial(double strain, double &stress,
                           double &tangent, double strainRate) {
    return 0;
}
int PipeMaterial::setTrial(double strain, double temperature,
                           double &stress, double &tangent,
                           double &thermalElongation,
                           double strainRate) {
    return 0;
}

double PipeMaterial::getStrain(void) { return 0.0; }
double PipeMaterial::getStrainRate(void) { return 0.0; }
double PipeMaterial::getStress(void) { return 0.0; }
double PipeMaterial::getTangent(void) { return 0.0; }
double PipeMaterial::getDampTangent(void) { return 0.0; }
double PipeMaterial::getInitialTangent(void) { return 0.0; }

int PipeMaterial::commitState(void) { return 0; }
int PipeMaterial::revertToLastCommit(void) { return 0; }
int PipeMaterial::revertToStart(void) { return 0; }

UniaxialMaterial *PipeMaterial::getCopy(void) {
    PipeMaterial *mat = new PipeMaterial(this->getTag());
    mat->points = this->points;
    mat->parameterID = this->parameterID;
    return mat;
}

int PipeMaterial::sendSelf(int commitTag, Channel &theChannel) {
    return 0;
}
int PipeMaterial::recvSelf(int commitTag, Channel &theChannel,
                           FEM_ObjectBroker &theBroker) {
    return 0;
}

void PipeMaterial::Print(OPS_Stream &s, int flag) {}

int PipeMaterial::setParameter(const char **argv, int argc,
                               Parameter &param) {
    return -1;
}
int PipeMaterial::updateParameter(int parameterID,
                                  Information &info) {
    return -1;
}

// AddingSensitivity:BEGIN
// //////////////////////////////////////////
int PipeMaterial::activateParameter(int parameterID) {
    this->parameterID = parameterID;
    return 0;
}
double PipeMaterial::getStressSensitivity(int gradIndex,
                                          bool conditional) {
    return 0;
}
double PipeMaterial::getTangentSensitivity(int gradIndex) {
    return 0.0;
}
double PipeMaterial::getInitialTangentSensitivity(int gradIndex) {
    return 0.0;
}
int PipeMaterial::commitSensitivity(double strainGradient,
                                    int gradIndex, int numGrads) {
    return 0;
}
// AddingSensitivity:END
// ///////////////////////////////////////////

void PipeMaterial::addPoint(double t, double e, double xnu,
                            double a) {
    points.push_back(PipeMaterialTemperaturePoint(t, e, xnu, a));
}

PipeMaterialTemperaturePoint PipeMaterial::selectPoint(double T,
                                                       int &retVal) {
    // line: 14148
    PipeMaterialTemperaturePoint pt;

    // only 1 point
    if (points.size() == 1) {
        pt = points[0];
        retVal = 0;
        return pt;
    }

    // more than 1 point
    int n = -1;
    for (int k = 1; k < (int)points.size(); ++k) {
        if (T >= points[k - 1].T && T <= points[k].T) {
            n = k;
            break;
        }
    }

    // not found
    if (n < 0) {
        opserr << "WARNING: the average temperature " << T
               << " for pipe material " << this->getTag()
               << " is not found\n";
        retVal = -1;
        return pt;
    }

    // dt
    auto &pt1 = points[n - 1];
    auto &pt2 = points[n];
    double dt = pt2.T - pt1.T;
    if (dt < 1e-8) {
        opserr << "WARNING: zero or negative temperature "
                  "difference between points ("
               << pt1.T << ") and (" << pt2.T << ")\n";
        retVal = -1;
        return pt;
    }

    // interpolate
    double ratio = (T - pt1.T) / dt;
    pt.E = pt1.E + ratio * (pt2.E - pt1.E);
    pt.xnu = pt1.xnu + ratio * (pt2.xnu - pt1.xnu);
    pt.alp = pt1.alp + ratio * (pt2.alp - pt1.alp);

    return pt;
}
