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
// $Source: /usr/local/cvs/OpenSees/SRC/material/uniaxial/ModIMKPinching.h,v $
                                                                       
// Written: Dimitrios G. Lignos, PhD, Assistant Professor, McGill University
// Created: February 2011
// Revision: A
//
//**********************************************************************
// Adapted by: Filipe Ribeiro and Andre Barbosa, September 30th 2013
// Oregon State University, OR, USA
//**********************************************************************
// 
// Description: This file contains the class interface for
// ModIMKPinching Model.  ModIMKPinching defines the modified IMK model with pinched hysteretic response

#ifndef ModIMKPinching_h
#define ModIMKPinching_h

#include <UniaxialMaterial.h>

class ModIMKPinching : public UniaxialMaterial
{
  public:
    ModIMKPinching(int tag, double Ke0, double AlfanPos, double AlfanNeg, double My_pos, double My_neg, double FprPos, double FprNeg, double A_pinch,		// Updated: Filipe Ribeiro and Andre Barbosa
					   double Ls, double Ld, double La, double Lk, double Cs, double Cd, double Ca, double Ck,
					   double ThetaPpos, double ThetaPneg, double ThetaPCpos, double ThetaPCneg,
					   double ResfacPos, double ResfacNeg, double FracDispPos, double FracDispNeg,
					   double DPos, double DNeg, double nFactor);														// Updated: Filipe Ribeiro and Andre Barbosa
    
	ModIMKPinching(int tag, double Ke0, double AlfanPos, double AlfanNeg, double My_pos, double My_neg, double FprPos, double FprNeg, double A_pinch,		// Updated: Filipe Ribeiro and Andre Barbosa
					   double Ls, double Ld, double La, double Lk, double Cs, double Cd, double Ca, double Ck,
					   double ThetaPpos, double ThetaPneg, double ThetaPCpos, double ThetaPCneg,
					   double ResfacPos, double ResfacNeg, double FracDispPos, double FracDispNeg,
					   double DPos, double DNeg);	
    ModIMKPinching();
    ~ModIMKPinching();

    const char *getClassType(void) const {return "ModIMKPinching";};

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
    double Ke0;			// Updated: Filipe Ribeiro and Andre Barbosa
	double nFactor;		// Updated: Filipe Ribeiro and Andre Barbosa
    double AlfanPos;	// Updated: Filipe Ribeiro and Andre Barbosa
    double AlfanNeg;	// Updated: Filipe Ribeiro and Andre Barbosa
    double My_pos;
    double My_neg;
   
    double FprPos;
    double FprNeg;
    double A_pinch;
   
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
   
    // flags
    int kon, Ckon;
    int flagStop, CflagStop;
    int flagdeg, Cflagdeg;
   
    double dmax, Cdmax;
    double dmin, Cdmin;
    double fmin, Cfmin;
    double fmax, Cfmax;
   
    double fyPos, CfyPos;
    double fyNeg, CfyNeg;
   
    double sn, Csn;
    double sp, Csp;
   
    double Enrgc, CEnrgc;
    double Enrgtot, CEnrgtot;
   
    double Enrgts, CEnrgts;
    double Enrgtd, CEnrgtd;
    double Enrgtk, CEnrgtk;
    double Enrgta, CEnrgta;
   
    double capSlopePos, CcapSlopePos;
    double capSlopeNeg, CcapSlopeNeg;
   
    double fCapRefPos, CfCapRefPos;
    double fCapRefNeg, CfCapRefNeg;
   
    double fCapPos, CfCapPos;
    double fCapNeg, CfCapNeg;
   
    double fpDegPos, CfpDegPos;
    double fpDegNeg, CfpDegNeg;
   
    double ekunload, Cekunload;
   
    double cpNeg, CcpNeg;
    double cpPos, CcpPos;
   
    double ekhardPos, CekhardPos;
    double ekhardNeg, CekhardNeg;
    double ekexcurs, Cekexcurs;
    double ekP, CekP;
	
	double Ke, CKe;				// Updated: Filipe Ribeiro and Andre Barbosa
    double AlfaPos, CAlfaPos;	// Updated: Filipe Ribeiro and Andre Barbosa
    double AlfaNeg, CAlfaNeg;	// Updated: Filipe Ribeiro and Andre Barbosa
	double prodBeta, CprodBeta; // Updated: Filipe Ribeiro and Andre Barbosa
	int commitCalledOnce;		// Updated: Filipe Ribeiro and Andre Barbosa

       
};


#endif