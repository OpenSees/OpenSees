#ifndef DegradingUniaxialWrapper_h
#define DegradingUniaxialWrapper_h
#include <string>
#include <functional> // std::function
#include <tcl.h>
#include <UniaxialMaterial.h>
#include "FedeasAPI.h"

class DegradingUniaxialWrapper : public UniaxialMaterial {
public:
  DegradingUniaxialWrapper(int tag, UniaxialMaterial &material, StateOperator* damage);
  DegradingUniaxialWrapper();
  ~DegradingUniaxialWrapper();

//   static UniaxialMaterial* parseNew(Tcl_Interp*, void*, int, TCL_Char **);

  const char *
  getClassType(void) const {return "DegradingUniaxialWrapper";}

  int    setTrialStrain(double strain, double strainRate = 0.0);
  int    setTrialStrain(double strain, double temperature, double strainRate);
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
  int recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker);

  void Print(OPS_Stream &s, int flag = 0);

  int setCoupling(double);
  int setParameter(const char **argv, int argc, Parameter &param);
  int updateParameter(int parameterID, Information &info);

  double getStressSensitivity(int gradIndex, bool conditional);
  double getStrainSensitivity(int gradIndex);
  double getInitialTangentSensitivity(int gradIndex);
  double getDampTangentSensitivity(int gradIndex);
  double getRhoSensitivity(int gradIndex);
  int commitSensitivity(double strainGradient, int gradIndex, int numGrads);
  
//  int setDamageWrapper(Tcl_Interp*, std::string);

protected:
private:
  UniaxialMaterial *theMaterial;

  double m_stress,
         m_tangent, 
         m_rate_tol=1e-6;

  struct UniaxialState {
  /* This struct defines the interface between
   * the public wrapper and it's external 
   * implementation. */
   double  e, ep, De, se, kt, ke;
  };

  UniaxialState past,pres;
  // typedef std::function<int(void*, void*)> degrade_f;
  StateOperator* degrade = NULL;
};

#endif // DegradingUniaxialWrapper_H

