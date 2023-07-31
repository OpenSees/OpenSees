#ifndef CreepMaterial_h
#define CreepMaterial_h

#include <UniaxialMaterial.h>

class CreepMaterial : public UniaxialMaterial
{
 public:
  CreepMaterial(int tag, UniaxialMaterial &wrapped);
  CreepMaterial();

  ~CreepMaterial();

  const char *getClassType(void) const {return "CreepMaterial";}
  UniaxialMaterial *getCopy(void);

  int setTrialStrain(double strain, double strainRate = 0.0);
  double getStrain(void);
  double getStress(void);
  double getTangent(void);
  double getInitialTangent(void);

  int commitState(void);
  int revertToLastCommit(void);
  int revertToStart(void);

  int sendSelf(int commitTag, Channel &theChannel);
  int recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker);

  void Print(OPS_Stream &s, int flag = 0);

  Response *setResponse(const char **argv, int argc, OPS_Stream &theOutput);
  int getResponse(int responseID, Information &matInfo);
  
 private:
  UniaxialMaterial *wrappedMaterial;

  enum {maxNumSteps = 500};
  enum {growNumSteps = 100};
  
  double PHI_i[maxNumSteps];
  double E_i[maxNumSteps];
  double DSIG_i[maxNumSteps];
  double dsig_i[maxNumSteps];
  double TIME_i[maxNumSteps];
  double DTIME_i[maxNumSteps];
};

#endif
