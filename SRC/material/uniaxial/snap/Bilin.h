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

//Bilinear material with strength and stiffness deterioration
//Based on SNAP Stanford element (modified Ib.&Kraw. element)

//**********************************************************************                                                                     
// Written: Theodore Karavasilis 
// Lecturer
// Department of Engineering Science, University of Oxford, Oxford, U.K.
//**********************************************************************

#ifndef Bilin_h
#define Bilin_h

#include <UniaxialMaterial.h>

class Bilin : public UniaxialMaterial
{
  public:
  Bilin(int tag, 
	double Ke,         double As,         double AsNeg,   double My_pos,    double My_neg,
	double LamdaS,     double LamdaK,     double LamdaA,  double LamdaD,    double Cs,         
        double Ck,         double Ca,         double Cd,      double Thetap_pos,double Thetap_neg,
        double Thetapc_pos,double Thetapc_neg,double K,       double KNeg,      double Thetau_pos,
        double Thetau_neg, double PDPlus,     double PDNeg);
  Bilin(); 
  ~Bilin();
  const char *getClassType(void) const {return "Bilin";};
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
  
 protected:
  
 private:
  //my functions
  void interPoint(double& xInt, double& yInt, double x1,double y1,double m1,double x2,double y2,double m2);
  void snCalc(void);
  void spCalc(void);
  void envelPosCap2(double fy,double alphaPos,double alphaCap,double cpDsp,double& d,
		    double& f,double& ek,double elstk,double fyieldPos,double Resfac);
  void envelNegCap2(double fy,double alphaNeg,double alphaCap,double cpDsp,double& d,
		    double& f,double& ek,double elstk,double fyieldNeg,double Resfac);
  double boundPos(void);
  double boundNeg(void);
  void envHitsZero(double& f);


  //Fixed input material parameters 
  double Ke;
  double As;
  double AsNeg;
  double My_pos;
  double My_neg;
  double LamdaS;
  double LamdaK;
  double LamdaA;
  double LamdaD;
  double Cs;
  double Ck;
  double Ca;
  double Cd;
  double Thetap_pos;
  double Thetap_neg;
  double Thetapc_pos;
  double Thetapc_neg;
  double K;
  double KNeg;
  double Thetau_pos;
  double Thetau_neg;
  double PDPlus;
  double PDNeg;
  //State variables 
  double U, CU; //displacement
  double Force, CForce; //force
  double Tangent, CTangent; //tangent stiffness
  //History variables
  double dNewLoadPos,CdNewLoadPos;
  double dNewLoadNeg, CdNewLoadNeg;
  int flgstop,Cflgstop; 
  int flagdeg,Cflagdeg ; 
  int flagstopdeg,Cflagstopdeg; 
  double ekt, Cekt;
  int interup, Cinterup;  
  int kcode, Ckcode; 
  int kon, Ckon;
  int iCapNeg, CiCapNeg; 
  int iNoFneg, CiNoFneg; 
  int iNoFpos, CiNoFpos; 
  int iCapPos, CiCapPos; 
  int iDeg, CiDeg; 
  int LP, CLP;
  int LN, CLN; 
  double capSlope, CcapSlope; 
  double capDispPos, CcapDispPos; 
  double capDispNeg, CcapDispNeg; 
  double elstk, Celstk; 
  double fyieldPos, CfyieldPos; 
  double fyieldNeg, CfyieldNeg; 
  double alpha, Calpha;  
  double ecaps, Cecaps; 
  double ecapk, Cecapk; 
  double ecapd, Cecapd; 
  double cs, Ccs; 
  double ck, Cck; 
  double cd, Ccd; 
  double dmax, Cdmax; 
  double dmin, Cdmin; 
  double Enrgtot, CEnrgtot; 
  double Enrgc, CEnrgc; 
  double fyPos, CfyPos; 
  double fLimNeg, CfLimNeg; 
  double fyNeg, CfyNeg; 
  double ekP, CekP; 
  double ekunload, Cekunload; 
  double sp, Csp; 
  double sn, Csn; 
  double dP, CdP; 
  double fP, CfP; 
  double ek, Cek; 
  double stif, Cstif; 
  double dLimPos, CdLimPos; 
  double dLimNeg, CdLimNeg;  
  double vtot, Cvtot; 
  double ftot, Cftot; 
  double dtot, Cdtot; 
  double dn, Cdn; 
  double cpPos, CcpPos; 
  double cpNeg, CcpNeg; 
  double fLimPos, CfLimPos; 
  double dlstPos, CdlstPos; 
  double flstPos, CflstPos; 
  double dlstNeg, CdlstNeg; 
  double flstNeg, CflstNeg; 
  double ekexcurs, Cekexcurs; 
  double RSE, CRSE; 
  double fPeakPos, CfPeakPos; 
  double fPeakNeg, CfPeakNeg; 
  double dCap1Pos, CdCap1Pos; 
  double dCap2Pos, CdCap2Pos; 
  double dCap1Neg, CdCap1Neg; 
  double dCap2Neg, CdCap2Neg; 
  double alphaNeg, CalphaNeg; 
  double alphaPos, CalphaPos; 
  double ekhardNeg, CekhardNeg; 
  double ekhardPos, CekhardPos; 
  double fCapRefPos, CfCapRefPos; 
  double fCapRefNeg, CfCapRefNeg; 
  double Enrgts, CEnrgts; 
  double Enrgtk, CEnrgtk; 
  double Enrgtd, CEnrgtd; 
  double dyPos, CdyPos; 
  double dyNeg, CdyNeg; 
  double dyieldPos, CdyieldPos; 
  double dyieldNeg, CdyieldNeg; 
  double resSnHor, CresSnHor; 
  double fmax, Cfmax; 
  double fmin, Cfmin; 
  double resSp, CresSp; 
  double resSn, CresSn; 
  double fCapPos, CfCapPos; 
  double fCapNeg, CfCapNeg; 
  double snHor, CsnHor; 
  double spHor, CspHor; 
  double resSpHor, CresSpHor;
  double snEnv, CsnEnv;
  double resSnEnv, CresSnEnv; 
  double spEnv, CspEnv; 
  double resSpEnv, CresSpEnv; 
  double Resfac, CResfac; 
  double capSlopeOrig,CcapSlopeOrig;  
  double fracDispPos, CfracDispPos; 
  double fracDispNeg, CfracDispNeg; 
  double DPlus, CDPlus; 
  double DNeg, CDNeg; 
  double alphaN, CalphaN;
  double ecapa, Cecapa;
  double ca, Cca;
  double ResfacNeg, CResfacNeg;
  double myFcapping, CmyFcapping;
  double myFcappingNeg, CmyFcappingNeg;
  double capSlopeNeg, CcapSlopeNeg;
  int flagControlResponse, CflagControlResponse; 
  double capSlopeOrigNeg, CcapSlopeOrigNeg;
  double Uprev,CUprev;
};


#endif

