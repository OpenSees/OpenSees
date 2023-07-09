/* ****************************************************************** **
**    OpenSees - Open System for Earthquake Engineering Simulation    **
**          Pacific Earthquake Engineering Research Center            **
**                                                                    **
**                                                                    **
** (C) Copyright 1999, The Regents of the University of California    **
** All Rights Reserved.                                               **
**                                                                    **
** Commercial use of this program without express permission of the   **
** University of California, Berkeley, is strictly prohibited.  See   **
** file 'COPYRIGHT'  in main directory for information on usage and   **
** redistribution,  and for a DISCLAIMER OF ALL WARRANTIES.           **
**                                                                    **
** Developed by:                                                      **
**   Frank McKenna (fmckenna@ce.berkeley.edu)                         **
**   Gregory L. Fenves (fenves@ce.berkeley.edu)                       **
**   Filip C. Filippou (filippou@ce.berkeley.edu)                     **
**                                                                    **
** ****************************************************************** */
 
/* 
** Description: This file contains the class implementation for a       
** Constitutive (FEM) model for FRP and Steel-Confined Concrete for 
** Circular Concrete Sections  **Created in May 2011**                 
** written and developed by:                                        
** Konstantinos G. Megalooikonomou (C)                              
*/

#ifndef FRPConfinedConcrete_h
#define FRPConfinedConcrete_h

#include <UniaxialMaterial.h>

class FRPConfinedConcrete : public UniaxialMaterial
{
 public:
  FRPConfinedConcrete(int tag, double fpc1, double fpc2, double epsc0, double D, double c,double Ej, double Sj, double tj, double eju, double S, double fyl, double fyh, double dlong, double dtrans, double Es, double v0, double k, double useBuck);
  FRPConfinedConcrete();
  ~FRPConfinedConcrete();

  const char *getClassType(void) const {return "FRPConfinedConcrete";};

  int setTrialStrain(double strain, double strainRate = 0.0);
  int setTrial (double strain, double &stress, double &tangent, double strainRate);

  double ComputeTendStrain(void);
  double getStress(void);
  double getStrain(void);
  double getTangent(void);
  double getInitialTangent(void);
  double getLatstress (void);
  double getLatStrain (void);

  int commitState(void);
  int revertToLastCommit(void);    
  int revertToStart(void);        
  
  UniaxialMaterial *getCopy(void);
  
  int sendSelf(int commitTag, Channel &theChannel);  
  int recvSelf(int commitTag, Channel &theChannel, 
	       FEM_ObjectBroker &theBroker);    
  
  void Print(OPS_Stream &s, int flag =0);

  int setParameter(const char **argv, int argc, Parameter &param);
  int updateParameter(int parameterID, Information &info);
  int activateParameter(int passedParameterID);
  double getStressSensitivity(int gradNumber, bool conditional);
  int commitSensitivity(double TstrainSensitivity, int gradNumber, int numGrads);

  int getVariable(const char *varName, Information &theInfo);
  
 protected:
  
 private:
  double fpc1;
  double fpc2;
  double epsc0;
  double D;
  double c;
  double Ej;
  double Sj;
  double tj;
  double eju;
  double S;
  double fyl;
  double fyh;
  double dlong;
  double dtrans;
  double Es;
  double v0;
  double k;
  double useBuck; //practically boolean but declared as double to maintain input uniformity

// History variables from last converged state (Past)
  double CminStrain;
  double CunloadSlope;
  double CendStrain;
  double CbLatstress;
  bool   CConvFlag;
  double CConfRat;
  double CConfStrain;
  double CLBuck;


// State variables from last converged state (Past)
  double Cstrain;
  double Cstress;
  double Ctangent;
  double CLatStrain;
  double CaLatstress;

// History variables from last converged state (Present)
  double TminStrain; 
  double TunloadSlope;
  double TendStrain;
  double TbLatstress;
  bool   TConvFlag;
  double TConfRat;
  double TConfStrain;
  double TLBuck;

  // State variables from last converged state (Present)
  double Tstrain;
  double Tstress;
  double Ttangent;
  double TLatStrain;
  double TaLatstress;

  void determineTrialState (double dStrain);
  void reload (void);
  void   unload (void);
  void envelope (void);
  void flat (double flcover_n, double arrayLat[6]);

  // AddingSensitivity:BEGIN //////////////////////////////////////////
  int parameterID;
  Matrix *SHVs;
 // AddingSensitivity:END ///////////////////////////////////////////
	
  bool buckCrInit;
  //addon for help
  bool myRegulaFalsi(double Pcr, double EIred, double Es, double Ash, double Dcore, double S, int mBuck, double& xRes, bool& returnFlag);
  double PCriticalSolve(double n, double Pcr, double EIred, double Es, double Ash, double Dcore, double S, int mBuck);
};


#endif
