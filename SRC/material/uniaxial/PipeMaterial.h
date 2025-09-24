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
                                 double a);
    PipeMaterialTemperaturePoint();

   public:
    double T;    // temperature
    double E;    // young's modulus
    double xnu;  // poisson's ratio
    double alp;  // thermal expansion coefficient
};

class PipeMaterial : public UniaxialMaterial {
   public:
    explicit PipeMaterial(int tag);
    PipeMaterial();
    ~PipeMaterial();

    const char *getClassType(void) const;

    int setTrialStrain(double strain, double strainRate = 0.0);
    int setTrialStrain(double strain, double temperature,
                       double strainRate);
    int setTrial(double strain, double &stress, double &tangent,
                 double strainRate = 0.0);
    int setTrial(double strain, double temperature, double &stress,
                 double &tangent, double &thermalElongation,
                 double strainRate = 0.0);

    double getStrain(void);
    double getStrainRate(void);
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
    double getInitialTangentSensitivity(int gradIndex);
    int commitSensitivity(double strainGradient, int gradIndex,
                          int numGrads);
    // AddingSensitivity:END
    // ///////////////////////////////////////////

    void addPoint(double t, double e, double xnu, double a);

    PipeMaterialTemperaturePoint selectPoint(double T, int &retVal);

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

#endif
