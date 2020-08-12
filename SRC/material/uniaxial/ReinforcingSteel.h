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
                                                                        
// $Revision: 1.4 $
// $Date: 2006-08-03 23:42:19 $
// $Source: /usr/local/cvs/OpenSees/SRC/material/uniaxial/ReinforcingSteel.h,v $

/* ****************************************************************** **
** THIS FILE WAS DEVELOPED AT UC DAVIS                                **
**                                                                    **
** Programmed by: Jon Mohle (jfmohle@ucdavis.edu)                     **
** Supervisor: Sashi Kunnath (skkunnath@ucdavis.edu)                  **
**                                                                    **
********************************************************************* */
// Written: Jon Mohle
// Created: October 2003
// Updated: September 2005
// Description: This file contains the class definition for 
// ReinforcingSteel.

#ifndef ReinforcingSteel_h
#define ReinforcingSteel_h
#define HelpDebugMat        //Debugging info, JFM

#include "UniaxialMaterial.h"

// LastRule_RS must be greater than or equal to 12.  Add branches 4 at a time
// eq 12, 16, 20, 24 etc... failure to add branches 4 at a time will cause problems.
// A higher number will result in better tracking of the loop memory effects at the cost of additional memory usage
// each pair of rules requires 11 additional double variables.
const int LastRule_RS=20;  // must be divisable by 4!!!!!!!!!!!

class ReinforcingSteel : public UniaxialMaterial
{
 public:
  ReinforcingSteel(int tag, double fyield, double fultimate, double youngs, double youngs_hard, 
		   double estrainhard, double eultimate, int buckModel, double slenderness, double alpha, double r, 
		   double gama, double Fatigue1, double Fatigue2, double Degrade1, 
		   double RC1, double RC2, double RC3, double A1, double HardLim);
  ReinforcingSteel(int tag);    
  ~ReinforcingSteel();

  const char *getClassType(void) const {return "ReinforcingSteel";};
  
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
  int recvSelf(int commitTag, Channel &theChannel, 
  	         FEM_ObjectBroker &theBroker);    
  
  void Print(OPS_Stream &s, int flag =0);
  
  virtual double getEnergy() { return Energy; } //by SAJalali

 protected:
  
 private:
	 double Energy; //by SAJalali

  #ifdef HelpDebugMat
    static int classCount;
    int thisClassNumber;
    int thisClassCommit;
    int thisClassStep;
  #endif

  double ZeroTol;
  double reduction;
  double fsu_fraction;
  double beta;
  int theBarFailed;

  // natural stress-strain variables
  double p;
  double Esp;   // natural Elastic Modulus
  double eshp;	// natural Hardening Strain
  double fshp;  // natural Hardening Stress
  double Eshp;  // natural Hardening Modulus
  double esup;  // natural Strain at Peak Stress
  double fsup;  // natural Peak Stress
  double Esup;  // natural peak stress Modulus
  double Eypp;  // natural Yield Plateau Modulus
  double fint;  // natural Stress yield plateau intersect
  double eyp;   // natural strain at yield
  double fyp;   // natural yield stress

  double fy;    // engineering yield strength (user input)
  double fsu;   // engineering ultimate strength (user input)
  double Es;    // engineering elastic modulus (user input)
  double Esh;   // engineering hardening slope (user input)
  double esh;   // engineering hardening strain (user input)
  double esu;   // engineering ultimate strain (user input)

  double eshpa;	// Curve smoothing Parameters (at SH transition)
  double Eshpb;	// These are used to eliminate a sudden discontinuity in stiffness

  double a1;    // Linear Hardening Constant (with relation to accumulated plastic strain)
  double hardLim;

  double THardFact;
  double CHardFact;

  // Strength degradation parameters
  double T_ePlastic[LastRule_RS/2+1];
  double C_ePlastic[LastRule_RS/2+1];
  //double Nbf;               // Cyclic Backbone factor used correct backbone proporsional to return strain
  double TFatDamage;
  double CFatDamage;
  double LDratio;
  double Fat1;
  double Fat2;
  double Deg1;
  int    BuckleModel;
  double BackStress;

  // Menegotto-Pinto Calibration Constants
  double RC1;
  double RC2;
  double RC3;

  // Menegotto-Pinto Equation parameters
  double TR;
  double Tfch;
  double TQ;
  double TEsec;
  double Tea;
  double Tfa;
  double TEa;
  double Teb;
  double Tfb;
  double TEb;

  double re;
  double rE1;
  double rE2;

  // Converged Menegotto-Pinto Equation parameters
  double CR[LastRule_RS/2+1];
  double Cfch[LastRule_RS/2+1];
  double CQ[LastRule_RS/2+1];
  double CEsec[LastRule_RS/2+1];
  double Cea[LastRule_RS/2+1];
  double Cfa[LastRule_RS/2+1];
  double CEa[LastRule_RS/2+1];
  double Ceb[LastRule_RS/2+1];
  double Cfb[LastRule_RS/2+1];
  double CEb[LastRule_RS/2+1];

  // Trial History Variables
  int    TBranchNum;
  int    TBranchMem;
  double Teo_p;
  double Teo_n;
  double Temax;
  double Temin;
  double TeAbsMax;
  double TeAbsMin;
  double TeCumPlastic;

  // Converged History Variables
  int    CBranchNum;
  double Ceo_p;
  double Ceo_n;
  double Cemax;
  double Cemin;
  double CeAbsMax;
  double CeAbsMin;
  double CeCumPlastic;

  // Trial State Variables
  double TStrain;           // Trial strain
  double TStress;           // Trial stress
  double TTangent;          // Trial tangent

  // Converged History Variables
  double CStrain;
  double CStress;
  double CTangent;

  // Private Functions
  int    inline Sign(double x);
  double Backbone_f(double ess);
  double Backbone_fNat(double essp);
  double Backbone_E(double ess);
  double Buckled_stress_Dhakal(double ess, double fss);
  double Buckled_stress_Gomes(double ess, double fss);
  double Buckled_mod_Gomes(double ess, double fss, double Ess);
  double Buckled_mod_Dhakal(double ess, double fss, double Ess);

  double inline MP_f(double e);
  double inline MP_E(double e);
  int    SetMP(void);
  double MPfunc(double a);
  void   inline SetTRp(void);
  void   inline SetTRn(void);
  void   inline SetTRp1(void);
  void   inline SetTRn1(void);
  void   SetPastCurve(int branch);

  int    BranchDriver(int res);
  int    Rule1(int res);
  int    Rule2(int res);
  int    Rule3(int res);
  int    Rule4(int res);
  int    Rule5(int res);
  int    Rule6(int res);
  int    Rule7(int res);
  int    Rule8(int res);
  int    Rule9(int res);
  int    Rule10(int res);
  int    Rule11(int res);
  int    Rule12(int res);

  double inline damage(double ehalfPlastic);
  double inline getPlasticStrain(double ehalf, double stressAmp);
  double scalefactor();
  double inline ReturnSlope(double dea);
  void updateHardeningLoaction(double PlasticStrain);
  void updateHardeningLoactionParams(void);
};

#endif
