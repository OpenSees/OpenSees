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

#ifndef CoulombDamperMaterial_h
#define CoulombDamperMaterial_h

#include <UniaxialMaterial.h>

class CoulombDamperMaterial : public UniaxialMaterial {
   public:
    CoulombDamperMaterial(int tag, double k, double fc, double t,
                          double damp, int m, int n);
    CoulombDamperMaterial();
    ~CoulombDamperMaterial();

    const char *getClassType(void) const {
        return "CoulombDamperMaterial";
    }

    int setTrialStrain(double strain, double strainRate = 0.0);
    int setTrial(double strain, double &stress, double &tangent,
                 double strainRate = 0.0);
    double getStrain(void) { return trialStrain; }
    double getStrainRate(void) { return trialStrainRate; }
    double getStress(void);
    double getTangent(void);
    double getDampTangent(void);
    double getInitialTangent(void);

    int commitState(void);
    int revertToLastCommit(void);
    int revertToStart(void);

    UniaxialMaterial *getCopy(void);

    int sendSelf(int commitTag, Channel &theChannel);
    int recvSelf(int commitTag, Channel &theChannel,
                 FEM_ObjectBroker &theBroker);

    void Print(OPS_Stream &s, int flag = 0);

    int setParameter(const char **argv, int argc, Parameter &param);
    int updateParameter(int parameterID, Information &info);

    // AddingSensitivity:BEGIN
    // //////////////////////////////////////////
    int activateParameter(int parameterID);
    double getStressSensitivity(int gradIndex, bool conditional);
    double getTangentSensitivity(int gradIndex);
    double getDampTangentSensitivity(int gradIndex);
    double getInitialTangentSensitivity(int gradIndex);
    int commitSensitivity(double strainGradient, int gradIndex,
                          int numGrads);
    // AddingSensitivity:END
    // ///////////////////////////////////////////

   protected:
   private:
    double trialStrain;
    double trialStrainRate;
    double tangent;
    double friction;
    double commitTrialStrainRate;
    int flipped;
    double tol;
    double dampOutTangent;
    int method;
    int numFlipped;

    // AddingSensitivity:BEGIN
    // //////////////////////////////////////////
    int parameterID;
    // AddingSensitivity:END
    // ///////////////////////////////////////////

    double sign();
    double dsign();
    double factor();
};

#endif
