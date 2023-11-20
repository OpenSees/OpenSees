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

#ifndef PipeSection_h
#define PipeSection_h

#include <Matrix.h>
#include <SectionForceDeformation.h>
#include <Vector.h>
#include <classTags.h>
#include <elementAPI.h>

double pi = 3.141592653589793;

class Channel;
class FEM_ObjectBroker;
class Information;

class PipeSection : public SectionForceDeformation {
    // 8305 - 8339
   public:
    PipeSection(int tag, double d, double t, double a, double w,
                double r)
        : SectionForceDeformation(tag, SEC_TAG_PipeSection),
          dout(d),
          thk(t),
          alphaV(a),
          wgt(w),
          rho(r),
          e(),
          mat(),
          code(2) {
        if (rho <= 1e-12) {
            rho = wgt / 386.4;
        }
        rout = dout * 0.5;
        rin = rout - thk;
        double rout2 = rout * rout;
        double rin2 = rin * rin;
        Ax = pi * (rout2 - rin2);
        Iy = 0.25 * pi * (rout2 * rout2 - rin2 * rin2);
        if (alphaV < 1e-8) {
            double dum2 = 4.0 * (rout2 * rout - rin2 * rin) / 3.0;
            double dum3 = (rout2 + rin2) * thk;
            if (dum3 < 1e-8) {
                opserr << "WARNING: (ro^2+ri^2) * t < 1e-8. AlphaV "
                          "is set to 100.\n";
                alphaV = 100.0;
            } else {
                alphaV = dum2 / dum3;
            }
        }
    }
    ~PipeSection(void) {}

    int commitState(void) { return 0; }
    int revertToLastCommit(void) { return 0; }
    int revertToStart(void) { return 0; }

    const char *getClassType(void) const { return "PipeSection"; };

    int setTrialSectionDeformation(const Vector &) { return 0; }
    const Vector &getSectionDeformation(void) { return e; }

    const Vector &getStressResultant(void) { return e; }
    const Matrix &getSectionTangent(void) { return mat; }
    const Matrix &getInitialTangent(void) { return mat; }
    const Matrix &getSectionFlexibility(void) { return mat; }
    const Matrix &getInitialFlexibility(void) { return mat; }

    SectionForceDeformation *getCopy(void) {
        auto *theCopy =
            new PipeSection(this->getTag(), this->dout, this->thk,
                            this->alphaV, this->wgt, this->rho);
        theCopy->parameterID = this->parameterID;
        return theCopy;
    }
    const ID &getType(void) { return code; }
    int getOrder(void) const { return 2; }

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
    int activateParameter(int parameterID) {
        this->parameterID = parameterID;
    }
    const Vector &getStressResultantSensitivity(int gradIndex,
                                                bool conditional) {
        return e;
    }
    const Matrix &getSectionTangentSensitivity(int gradIndex) {
        return mat;
    }
    const Matrix &getInitialTangentSensitivity(int gradIndex) {
        return mat;
    }
    const Matrix &getSectionFlexibilitySensitivity(int gradIndex) {
        return mat;
    }
    const Matrix &getInitialFlexibilitySensitivity(int gradIndex) {
        return mat;
    }

   protected:
   private:
    double dout;
    double thk;
    double alphaV;
    double wgt;
    double rho;

    Vector e;
    Matrix mat;
    ID code;

    double rout;
    double rin;
    double Ax;
    double Iy, Iz;
    double Jx;
    double Ay, Az;

    int parameterID;
};

void *OPS_PipeSection(void) {
    // 7648 - 7669
    // check inputs
    if (OPS_GetNumRemainingInputArgs() < 5) {
        opserr << "Invalid #args,  want: section Pipe "
                  "tag? do? t? alphaV? w? <rho?>\n";
        return 0;
    }

    // get tag
    int iData[1];
    int numData = 1;
    if (OPS_GetIntInput(&numData, iData) < 0) {
        opserr << "WARNING invalid tag  for pipe section\n";
        return 0;
    }

    // get data: do, t, alphaV, w
    double data[4];
    numData = 4;
    if (OPS_GetDoubleInput(&numData, data) < 0) {
        opserr << "WARNING: invalid data for pipe section\n";
        return 0;
    }
    if (data[0] <= 1e-8) {
        opserr << "WARNING: outside diameter is zero or negative\n";
        return 0;
    }
    if (data[1] <= 1e-8) {
        opserr << "WARNING: thickness is zero or negative\n";
        return 0;
    }
    if (data[3] <= 0.0) {
        opserr << "WARNING: weight is zero or negative\n";
        return 0;
    }

    double rho = 0.0;
    numData = 1;
    if (OPS_GetNumRemainingInputArgs() > 0) {
        if (OPS_GetDoubleInput(&numData, data) < 0) {
            opserr << "WARNING: invalid rho for pipe section\n";
            return 0;
        }
    }

    auto *sect = new PipeSection(iData[0], data[0], data[1], data[2],
                                 data[3], rho);

    return sect;
}

#endif
