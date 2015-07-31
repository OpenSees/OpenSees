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
                                                                     
// $Revision: 1 $
// $Date: 2011/02/01 12:35:01 $
// $Source: /usr/local/cvs/OpenSees/SRC/material/uniaxial/ModIMKPeakOriented02.h,v $
                                                                       
// Written: Dimitrios G. Lignos, PhD, Assistant Professor, McGill University
// Created: February 2011
// Revision: A
//
//**********************************************************************
// Adapted by: Filipe Ribeiro and Andre Barbosa, June 2015
// Oregon State University, OR, USA
//**********************************************************************
//  
// Description: This file contains the class interface for
// ModIMKPeakOriented02 Model.  ModIMKPeakOriented02 defines the modified IMK model with peak-oriented hysteretic response

#ifndef ModIMKPeakOriented02_h
#define ModIMKPeakOriented02_h

#include <UniaxialMaterial.h>
// Updated:Filipe Ribeiro and Andre Barbosa

class ModIMKPeakOriented02 : public UniaxialMaterial
{
  public:
  ModIMKPeakOriented02(int tag, double Ke0, double AlfanPos, double AlfanNeg, double My_pos, double My_neg,
		       double Ls, double Ld, double La, double Lk, double Cs, double Cd, double Ca, double Ck,
		       double ThetaPpos, double ThetaPneg, double ThetaPCpos, double ThetaPCneg,
		       double ResfacPos, double ResfacNeg, double FracDispPos, double FracDispNeg,
		       double DPos, double DNeg, double C_Fp, double C_Fn, double nFactor);               
  
  ModIMKPeakOriented02(int tag, double Ke0, double AlfanPos, double AlfanNeg, double My_pos, double My_neg,  
		       double Ls, double Ld, double La, double Lk, double Cs, double Cd, double Ca, double Ck,
		       double ThetaPpos, double ThetaPneg, double ThetaPCpos, double ThetaPCneg,
		       double ResfacPos, double ResfacNeg, double FracDispPos, double FracDispNeg,
		       double DPos, double DNeg, double C_Fp, double C_Fn);          
  
  ModIMKPeakOriented02(int tag, double Ke0, double AlfanPos, double AlfanNeg, double My_pos, double My_neg,  
		       double Ls, double Ld, double La, double Lk, double Cs, double Cd, double Ca, double Ck,
		       double ThetaPpos, double ThetaPneg, double ThetaPCpos, double ThetaPCneg,
		       double ResfacPos, double ResfacNeg, double FracDispPos, double FracDispNeg,
		       double DPos, double DNeg);    
  ModIMKPeakOriented02(int tag, double Ke0, double AlfanPos, double AlfanNeg, double My_pos, double My_neg,  
		       double Ls, double Ld, double La, double Lk, double Cs, double Cd, double Ca, double Ck,
		       double ThetaPpos, double ThetaPneg, double ThetaPCpos, double ThetaPCneg,
		       double ResfacPos, double ResfacNeg, double FracDispPos, double FracDispNeg,
		       double DPos, double DNeg, double nFactor);   
  
  ModIMKPeakOriented02();
  ~ModIMKPeakOriented02();
  
  const char *getClassType(void) const {return "ModIMKPeakOriented02";};
   
    int setTrialStrain(double strain, double strainRate = 0.0);
    double getStrain(void);
    double getStrainRate(void);
    double getStress(void);
   
    double getTangent(void);
    double getInitialTangent(void);
    double getDampTangent(void);


    int commitState(void);
    int revertToLastCommit(void);    
    int revertToStart(void);        

    UniaxialMaterial *getCopy(void);
   
    int sendSelf(int commitTag, Channel &theChannel);  
    int recvSelf(int commitTag, Channel &theChannel,
                 FEM_ObjectBroker &theBroker);    
   
    void Print(OPS_Stream &s, int flag =0);
    Response *setResponse (const char **argv, int argc, OPS_Stream &theOutputStream);
    int getResponse (int responseID, Information &matInformation);    

   
  protected:
   
  private:
       
    // Subroutines to be used inside the material
    void envelPosCap2(double fy, double alphaPos, double alphaCap, double cpDsp, double& d,
                      double& f, double& ek, double elstk, double fyieldPos, double Resfac,
                      double fracDisp, int& flagStop);
   
    void envelNegCap2(double fy, double alphaNeg, double alphaCap, double cpDsp, double& d,
                      double& f, double& ek, double elstk, double fyieldNeg, double Resfac,
                      double fracDisp, int& flagStop);
   
    // Fixed Input Material Variables
    double Ke0;                 // Updated:Filipe Ribeiro and Andre Barbosa    
    double nFactor;         // Updated:Filipe Ribeiro and Andre Barbosa
    double C_Fp;				// Updated:Filipe Ribeiro and Andre Barbosa
    double C_Fn;				// Updated:Filipe Ribeiro and Andre Barbosa
    double AlfanPos;    // Updated:Filipe Ribeiro and Andre Barbosa
    double AlfanNeg;    // Updated:Filipe Ribeiro and Andre Barbosa
    double My_pos;
    double My_neg;
   
    double Ls;
    double Ld;
    double La;
    double Lk;
   
    double Cs;
    double Cd;
    double Ca;
    double Ck;
   
    double ThetaPpos;
    double ThetaPneg;
    double ThetaPCpos;
    double ThetaPCneg;
   
    double ResfacPos;
    double ResfacNeg;
    double FracDispPos;
    double FracDispNeg;
   
    double DPos;
    double DNeg;
   
    // State Variables
    double Cstrain;     // Deformation
    double Cstress;     // Force
    double Tangent, Ctangent;   // Tangent Stiffness
   
    // Trial and Committeed State Variables
    double dP, CdP;         // Committed Strain
    double fP, CfP;         // Committed Stress
    double ek, Cek;         // Committed Tangent
   
    int Unl, CUnl;
    int kon, Ckon;
    int flagStop, CflagStop;
    int flagdeg, Cflagdeg;
   
    double dmax, Cdmax;
    double dmin, Cdmin;
    double fmin, Cfmin;
    double fmax, Cfmax;
   
    double fyPos, CfyPos;
    double fyNeg, CfyNeg;
   
    double dlstPos, CdlstPos;
    double dlstNeg, CdlstNeg;
   
    double flstPos, CflstPos;
    double flstNeg, CflstNeg;
   
    double sn, Csn;
    double sp, Csp;
   
    double Enrgc, CEnrgc;
    double Enrgtot, CEnrgtot;
   
    double Enrgts, CEnrgts;
    double Enrgtd, CEnrgtd;
    double Enrgtk, CEnrgtk;
    double Enrgta, CEnrgta;
   
    double fPeakPos, CfPeakPos;
    double fPeakNeg, CfPeakNeg;
   
    double capSlopePos, CcapSlopePos;
    double capSlopeNeg, CcapSlopeNeg;
   
    double fCapRefPos, CfCapRefPos;
    double fCapRefNeg, CfCapRefNeg;
   
    double ekunload, Cekunload;
   
    double cpNeg, CcpNeg;
    double cpPos, CcpPos;
   
    double ekhardPos, CekhardPos;
    double ekhardNeg, CekhardNeg;
    double ekexcurs, Cekexcurs;
    double ekP, CekP;
   
    double RSE, CRSE;
    double dres;
       
    double Ke, CKe;                         // Updated:Filipe Ribeiro and Andre Barbosa
    double AlfaPos, CAlfaPos;   // Updated:Filipe Ribeiro and Andre Barbosa
    double AlfaNeg, CAlfaNeg;   // Updated:Filipe Ribeiro and Andre Barbosa
    double prodBeta, CprodBeta; // Updated:Filipe Ribeiro and Andre Barbosa
    int commitCalledOnce;           // Updated:Filipe Ribeiro and Andre Barbosa
};


#endif
 

