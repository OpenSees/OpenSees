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

//Modified Ibarra-Medina-Krawinkler with Bilinear hysteretic response

//**********************************************************************                                                                    
// Code Developed by: Dimitrios G. Lignos
// Assistant Professor, McGill University, Montreal Canada
// Originally Written by: Theodore Karavasilis
// Lecturer
// Department of Engineering Science, University of Oxford, Oxford, U.K.
// Re-Written by: D. G. Lignos, August 29th 2012
//**********************************************************************
// Adapted by: Filipe Ribeiro and Andre Barbosa, Sep 13th 2013
// Oregon State University, OR, USA
//**********************************************************************

#ifndef Bilin02_h
#define Bilin02_h

#include <UniaxialMaterial.h>

class Bilin02 : public UniaxialMaterial
{
  public:
  Bilin02(int tag,
        double Ke0,     double As,         double AsNeg,   double My_pos,    double My_neg,            
        double LamdaS,     double LamdaD, double LamdaA,  double LamdaK,    double Cs,        
        double Cd,         double Ca,         double Ck,      double Thetap_pos,double Thetap_neg,
        double Thetapc_pos,double Thetapc_neg,double K,       double KNeg,      double Thetau_pos,
        double Thetau_neg, double PDPlus,     double PDNeg, double nFactor);  // Updated: Filipe Ribeiro and Andre Barbosa    

     // Updated: Filipe Ribeiro and Andre Barbosa
  Bilin02(int tag,                                                                                      
        double Ke0,     double As,         double AsNeg,   double My_pos,    double My_neg,            
        double LamdaS,     double LamdaD, double LamdaA,  double LamdaK,    double Cs,                
        double Cd,         double Ca,         double Ck,      double Thetap_pos,double Thetap_neg,    
        double Thetapc_pos,double Thetapc_neg,double K,       double KNeg,      double Thetau_pos,    
        double Thetau_neg, double PDPlus,     double PDNeg);                                          
        //  Updated: Filipe Ribeiro and Andre Barbosa
  Bilin02();
  ~Bilin02();

  const char *getClassType(void) const {return "Bilin02";};
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
  double Ke0;               // Updated: Filipe Ribeiro and Andre Barbosa
  double nFactor;           // Updated: Filipe Ribeiro and Andre Barbosa
  double AsPos;
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
  double KPos;
  double KNeg;
  double Thetau_pos;
  double Thetau_neg;
  double PDPlus;
  double PDNeg;

  //State variables
  double U, CU; //displacement
  double Tangent, CTangent; //tangent stiffness

  //History variables
  double dNewLoadPos,CdNewLoadPos;
  double dNewLoadNeg, CdNewLoadNeg;
 
  int flagdeg,Cflagdeg ;
  int flagstopdeg,Cflagstopdeg;
  int interup, Cinterup;  
  int kon, Ckon;
  int iNoFneg, CiNoFneg;
  int iNoFpos, CiNoFpos;
  int LP, CLP;
  int LN, CLN;
  int flagControlResponse, CflagControlResponse;
  double capSlope, CcapSlope;
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
  double dLimPos, CdLimPos;
  double dLimNeg, CdLimNeg;  
  double cpPos, CcpPos;
  double cpNeg, CcpNeg;
  double fLimPos, CfLimPos;
  double ekexcurs, Cekexcurs;
  double RSE, CRSE;
  double fPeakPos, CfPeakPos;
  double fPeakNeg, CfPeakNeg;
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

  double capSlopeOrig,CcapSlopeOrig;  
  double capSlopeNeg, CcapSlopeNeg;
  double capSlopeOrigNeg, CcapSlopeOrigNeg;

  double Ke, CKe;                                      // Updated: Filipe Ribeiro and Andre Barbosa
  double capSlopeMember, CcapSlopeMember;              // Updated: Filipe Ribeiro and Andre Barbosa
  double capSlopeNegMember, CcapSlopeNegMember;        // Updated: Filipe Ribeiro and Andre Barbosa
  double prodBeta, CprodBeta;                          // Updated: Filipe Ribeiro and Andre Barbosa
  
  int commitCalledOnce;
};


#endif
  