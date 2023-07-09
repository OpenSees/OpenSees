 
// Written: Q. Gu and Z.Zeng
// Created: 2014
//Reference:Seismic Performance of SRC Columns with High Ratio of Encased Steel written by YIN Xiaowei,supervised by Prof. LU Xilin,Tongji University 2012 
// Description: This file contains the class definition for 
// ResilienceLow. 

#ifndef ResilienceLow_h
#define ResilienceLow_h

#include <UniaxialMaterial.h>

class ResilienceLow : public UniaxialMaterial
{
 public:
  ResilienceLow(int tag,double PY, double DPmax, double Pmax, double Ke,double Kd);
     
  ~ResilienceLow();

  const char *getClassType(void) const {return "ResilienceLow";};

  int setTrialStrain(double strain, double strainRate = 0.0); 
  double getStrain(void);
  double getStress(void);
  double getTangent(void);
  double getInitialTangent(void);
  int commitState(void);
  int revertToLastCommit(void);    
  int revertToStart(void);        
  UniaxialMaterial *getCopy(void);
  
  int sendSelf(int commitTag, Channel &theChannel);  
  int recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker);    
  void Print(OPS_Stream &s, int flag =0);
 int determineState();
 private:
  double strain;   // trial strain
  double stress;   // trial stress
  double tangent;  // trial tangent
  double Cstrain;   // commit strain
  double Cstress;   // commit stress
  double Ctangent;  // commit tangent
// ---------- envelop curve ------
  double  DY, PY, DPmax, Pmax,Kd;
  double  Ke, Kui,Kri,CKui,CKri,Di,CDi; 
  int mode, Cmode,Flag,CFlag;
  double   strainRFMode2,  stressRFMode2,  strainRFMode4,  stressRFMode4,  strainRFMode6,  stressRFMode6, strainRFMode8,  stressRFMode8,  strainRFMode10,  stressRFMode10,strainRFMode11,  stressRFMode11 ;
  double   CstrainRFMode2,  CstressRFMode2,  CstrainRFMode4,  CstressRFMode4,  CstrainRFMode6,CstressRFMode6, CstrainRFMode8,  CstressRFMode8, CstrainRFMode10,CstressRFMode10,CstrainRFMode11, CstressRFMode11  ;
};
#endif

