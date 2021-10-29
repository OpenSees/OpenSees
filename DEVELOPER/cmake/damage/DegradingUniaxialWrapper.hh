#ifndef DegradingUniaxialWrapper_h
#define DegradingUniaxialWrapper_h

#include <UniaxialMaterial.h>

class DegradingUniaxialWrapper : public UniaxialMaterial {
public:
  DegradingUniaxialWrapper(int tag, UniaxialMaterial &material, double min, double max);
  DegradingUniaxialWrapper();
  ~DegradingUniaxialWrapper();

  const char *
  getClassType(void) const {return "DegradingUniaxialWrapper";}

  int    setTrialStrain(double strain, double strainRate = 0.0);
  int    setTrialStrain(double strain, double FiberTemperature, double strainRate);
  double getStrain(void);
  double getStrainRate(void);
  double getStress(void);
  double getTangent(void);
  double getDampTangent(void);
  double getInitialTangent(void){return theMaterial->getInitialTangent();}

  int commitState(void);
  int revertToLastCommit(void);
  int revertToStart(void);

  UniaxialMaterial *getCopy(void);

  int sendSelf(int commitTag, Channel &theChannel);
  int recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker);

  void Print(OPS_Stream &s, int flag = 0);
  bool
  hasFailed(void)
  {
    return Cfailed;
  }

  int setParameter(const char **argv, int argc, Parameter &param);
  int updateParameter(int parameterID, Information &info);

  // AddingSensitivity:BEGIN //////////////////////////////////////////
  double getStressSensitivity(int gradIndex, bool conditional);
  double getStrainSensitivity(int gradIndex);
  double getInitialTangentSensitivity(int gradIndex);
  double getDampTangentSensitivity(int gradIndex);
  double getRhoSensitivity(int gradIndex);
  int commitSensitivity(double strainGradient, int gradIndex, int numGrads);
  // AddingSensitivity:END ///////////////////////////////////////////

protected:
private:
  UniaxialMaterial *theMaterial;

  double minStrain;
  double maxStrain;

  bool Tfailed;
  bool Cfailed;
};

#endif // DegradingUniaxialWrapper_H
