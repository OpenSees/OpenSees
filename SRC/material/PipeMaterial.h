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

#ifndef PipeMaterial_h
#define PipeMaterial_h

// Written: Minjie
//
// Description: material for pipe element

#include <UniaxialMaterial.h>
#include <elementAPI.h>

#include <vector>

class PipeMaterialTemperaturePoint {
   public:
    PipeMaterialTemperaturePoint(double t, double e, double nu,
                                 double a)
        : T(t), E(e), xnu(nu), alp(a) {}

   public:
    double T;    // temperature
    double E;    // young's modulus
    double xnu;  // poisson's ratio
    double alp;  // termal expansion coefficient
};

class PipeMaterial : public UniaxialMaterial {
   public:
    explicit PipeMaterial(int tag)
        : UniaxialMaterial(tag, MAT_TAG_PipeMaterial) {}
    PipeMaterial() : UniaxialMaterial(0, MAT_TAG_PipeMaterial) {}
    ~PipeMaterial() {}

    const char *getClassType(void) const { return "PipeMaterial"; }

    int setTrialStrain(double strain, double strainRate = 0.0) {
        // trialStrain = strain;
        // trialStrainRate = strainRate;
    }
    int setTrialStrain(double strain, double temperature,
                       double strainRate) {}
    int setTrial(double strain, double &stress, double &tangent,
                 double strainRate = 0.0) {}
    int setTrial(double strain, double temperature, double &stress,
                 double &tangent, double &thermalElongation,
                 double strainRate = 0.0) {}

    double getStrain(void) { return 0.0; }
    double getStrainRate(void) { return 0.0; }
    double getStress(void) { return 0.0; }
    double getTangent(void) { return 0.0; }
    double getDampTangent(void) { return 0.0; }
    double getInitialTangent(void) { return 0.0; }

    int commitState(void) { return 0; }
    int revertToLastCommit(void) { return 0; }
    int revertToStart(void) { return 0; }

    UniaxialMaterial *getCopy(void) {
        PipeMaterial *mat = new PipeMaterial(this->getTag());
        mat->points = this->points;
        mat->parameterID = this->parameterID;
        return mat;
    }

    int sendSelf(int commitTag, Channel &theChannel) { return 0; }
    int recvSelf(int commitTag, Channel &theChannel,
                 FEM_ObjectBroker &theBroker) {
        return 0;
    }

    void Print(OPS_Stream &s, int flag = 0) {}

    int setParameter(const char **argv, int argc, Parameter &param) {
        return -1;
    }
    int updateParameter(int parameterID, Information &info) {
        return -1;
    }

    // AddingSensitivity:BEGIN
    // //////////////////////////////////////////
    int activateParameter(int parameterID) {
        this->parameterID = parameterID;
    }
    double getStressSensitivity(int gradIndex, bool conditional) {
        return 0;
    }
    double getTangentSensitivity(int gradIndex) { return 0.0; }
    double getInitialTangentSensitivity(int gradIndex) { return 0.0; }
    int commitSensitivity(double strainGradient, int gradIndex,
                          int numGrads) {
        return 0;
    }
    // AddingSensitivity:END
    // ///////////////////////////////////////////

    void addPoint(double t, double e, double xnu, double a) {
        points.push_back(PipeMaterialTemperaturePoint(t, e, xnu, a));
    }

    int selectPoint(double T, PipeMaterialTemperaturePoint &pt) {
        // line: 11371 - 11427

        // only 1 point
        if (points.size() == 1) {
            pt = points[0];
            return 0;
        }

        // more than 1 point
        int n = -1;
        for (int k = 2; k < (int)points.size(); ++k) {
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
            return -1;
        }

        // dt
        auto &pt1 = points[n - 1];
        auto &pt2 = points[n];
        double dt = pt2.T - pt1.T;
        if (dt < 1e8) {
            opserr << "WARNING: zero or negative temperature "
                      "difference between points ("
                   << pt1.T << ") and (" << pt2.T << "\n";
            return -1;
        }

        // interpolate
        double ratio = (T - pt1.T) / dt;
        pt.E = pt1.E + ratio * (pt2.E - pt1.E);
        pt.xnu = pt1.xnu + ratio * (pt2.xnu - pt1.xnu);
        pt.alp = pt1.alp + ratio * (pt2.alp - pt1.alp);

        return 0;
    }

   protected:
   private:
    // double trialStrain;
    // double trialStrainRate;
    // double trialTemperature;
    std::vector<PipeMaterialTemperaturePoint> points;

    // AddingSensitivity:BEGIN
    // //////////////////////////////////////////
    int parameterID;
    // AddingSensitivity:END
    // ///////////////////////////////////////////
};

void *OPS_PipeMaterial(void) {
    // line: 7607 - 7643
    // check inputs
    if (OPS_GetNumRemainingInputArgs() < 6) {
        opserr << "Invalid #args,  want: uniaxialMaterial Pipe "
                  "tag? nt? T1? E1? xnu1? alp1? <T2? E2? xnu2? alp2? "
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

#endif
