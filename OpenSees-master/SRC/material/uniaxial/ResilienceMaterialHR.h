 /* ********************************************************************
**    OpenSees - Open System for Earthquake Engineering Simulation    **
**          Pacific Earthquake Engineering Research Center	      **
**                                                                    **
** (C) Copyright by Quan Gu and Yongdou Liu @ Xiamen University	      **
**								      **
** ****************************************************************** */
// Created: 2014

// Referenced from "Seismic Performance of SRC Columns with High Ratio of Encased Steel ".Prof. LU Xilin.    and  YIN Xiaowei 
// Description: This file contains the class interface for ResilienceMaterialHR

#ifndef ResilienceMaterialHR_h
#define ResilienceMaterialHR_h

#include <UniaxialMaterial.h>

class ResilienceMaterialHR : public UniaxialMaterial
{
 public:
  ResilienceMaterialHR(int tag, double  DY, double PY, double DPmax, double Pmax, double Ke,double Kd, double k);
  ~ResilienceMaterialHR();
  const char *getClassType(void) const {return "ResilienceMaterialHR";};
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
  double  DY, PY, DPmax, Pmax;
  double  Ke, Ku,Kr,Kd,coefficient; 
  int mode, Cmode; 
  double   strainP8,  strainP9,  strainP11,  strainP12,  strainP10,  strainP7,  strainRFMode9,  stressRFMode9,  strainRFMode12,  stressRFMode12;
  double   strainRFMode2,  stressRFMode2,  strainRFMode4,  stressRFMode4,  strainRFMode6,  stressRFMode6,  stressP7,  stressP10,  strainRFMode13;
  double   stressRFMode13,CstressRFMode12,  CstrainRFMode2,  CstressRFMode2,  CstrainRFMode4,  CstressRFMode4,  CstrainRFMode6;
  double   CstrainP8,  CstrainP9,  CstrainP11,  CstrainP12,  CstrainP10,  CstrainP7,  CstrainRFMode9,  CstressRFMode9,  CstrainRFMode12;
  double    CstressRFMode6,  CstressP7,  CstressP10,  CstrainRFMode13,  CstressRFMode13;
};
#endif

