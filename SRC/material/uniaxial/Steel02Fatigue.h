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
                                                                        
// $Revision: 1.3 $
// $Date: 2006-11-03 18:40:20 $
// $Source: /usr/local/cvs/OpenSees/SRC/material/uniaxial/Steel02.h,v $
                                                                      
// Written: fmk
// Created: 03/06
//
// Description: This file contains the class definition for 
// Steel02. Steel02 is based on an f2c of the FEDEAS material
// Steel02.f which is:
//-----------------------------------------------------------------------
// MENEGOTTO-PINTO STEEL MODEL WITH FILIPPOU ISOTROPIC HARDENING
//            written by MOHD YASSIN (1993)
//          adapted to FEDEAS material library
//    by E. Spacone, G. Monti and F.C. Filippou (1994)
//    Modified by Nasser A. Marafi (2018) to have strength and 
//   stiffness deterioration as per Kunnath et al. (2009)
//-----------------------------------------------------------------------


#ifndef Steel02Fatigue_h
#define Steel02Fatigue_h

#include <UniaxialMaterial.h>

class Steel02Fatigue : public UniaxialMaterial
{
  public:
	Steel02Fatigue(int tag,
	    double fy, double E0, double b,
	    double R0, double cR1, double cR2,
	    double a1, double a2, double a3, double a4, 
		double Cd , double Cf , double Alpha , double Beta, double minStrain , double maxStrain, double sigInit = 0.0 );

    
    // Constructor for no isotropic hardening
	Steel02Fatigue(int tag,
	    double fy, double E0, double b,
	    double R0, double cR1, double cR2,
		double Cd , double Cf , double Alpha , double Beta , double minStrain, double maxStrain);

    
    // Constructor for no isotropic hardening
    // Also provides default values for R0, cR1, and cR2
	Steel02Fatigue(int tag, double fy, double E0, double b,
		double Cd , double Cf , double Alpha , double Beta , double minStrain, double maxStrain);

	    
	Steel02Fatigue(void);
    virtual ~Steel02Fatigue();
    

    const char *getClassType(void) const {return "Steel02Fatigue";};

    double getInitialTangent(void);
    UniaxialMaterial *getCopy(void);

    int setTrialStrain(double strain, double strainRate = 0.0); 
    double getStrain(void);      
    double getStress(void);
    double getTangent(void);

	int commitSteel(void);
    
    int commitState(void);
    int revertToLastCommit(void);    
    int revertToStart(void);        
    
    int sendSelf(int commitTag, Channel &theChannel);  
    int recvSelf(int commitTag, Channel &theChannel, 
		 FEM_ObjectBroker &theBroker);    
    
    void Print(OPS_Stream &s, int flag =0);
    
 protected:
    
 private:
    // matpar : STEEL FIXED PROPERTIES
    double Fy;  //  = matpar(1)  : yield stress
    double E0;  //  = matpar(2)  : initial stiffness
    double b;   //  = matpar(3)  : hardening ratio (Esh/E0)
    double R0;  //  = matpar(4)  : exp transition elastic-plastic
    double cR1; //  = matpar(5)  : coefficient for changing R0 to R
    double cR2; //  = matpar(6)  : coefficient for changing R0 to R
    double a1;  //  = matpar(7)  : coefficient for isotropic hardening in compression
    double a2;  //  = matpar(8)  : coefficient for isotropic hardening in compression
    double a3;  //  = matpar(9)  : coefficient for isotropic hardening in tension
    double a4;  //  = matpar(10) : coefficient for isotropic hardening in tension
    double sigini; // initial 
    // hstvP : STEEL HISTORY VARIABLES
    double epsminP; //  = hstvP(1) : max eps in compression
    double epsmaxP; //  = hstvP(2) : max eps in tension
    double epsplP;  //  = hstvP(3) : plastic excursion
    double epss0P;  //  = hstvP(4) : eps at asymptotes intersection
    double sigs0P;  //  = hstvP(5) : sig at asymptotes intersection
    double epssrP;  //  = hstvP(6) : eps at last inversion point
    double sigsrP;  //  = hstvP(7) : sig at last inversion point
    int    konP;    //  = hstvP(8) : index for loading/unloading
    // hstv : STEEL HISTORY VARIABLES   
    double epsP;  //  = strain at previous converged step
    double sigP;  //  = stress at previous converged step
    double eP;    //   stiffness modulus at last converged step;

    double epsmin; 
    double epsmax; 
    double epspl;  
    double epss0;  
    double sigs0; 
    double epsr;  
    double sigr;  
    int    kon;    
    double sig;   
    double e;     
    double eps;   //  = strain at current step

	double Fatigue_DI; //Damage index
	double Fatigue_X; //Range in consideration
	double Fatigue_Y; //Previous Adjacent Range
	double Fatigue_A; //Peak or valley 1
	double Fatigue_B; //Peak or valley 2
	double Fatigue_C; //Peak or valley 2
	double Fatigue_D; //Peak or valley 4
	int   Fatigue_PCC; /*Previous Cycle counter flag if >1 then previous 'n'
			   cycles did not flag a complete cycle */
	int   Fatigue_R1F; //Flag for first  peak count
	int   Fatigue_R2F; //Flag for second peak count
	double Fatigue_cSlope; //Current Slope
	double Fatigue_PS; //Previous slope
	double Fatigue_EP; //Previous Strain
	int    Fatigue_SF; /*Start Flag = 0 if very first strain,
			   (i.e. when initializing)    = 1 otherwise */
	double Fatigue_DL; //Damage if current strain was last peak.

	double Fatigue_Dmax;
	double Fatigue_E0;
	double Fatigue_m;

	double minStrain;
	double maxStrain;

	bool Fatigue_Cfailed;
	double trialStrain;

	// added 6/9/2006
	// For recording strain ranges (SRXX) and Number of Cycles (NCXX)
	double Fatigue_SR1;  // Committed strain range at peak
	double Fatigue_NC1;  // Committed number of cycles at SR1 (i.e. 1.0 or 0.5)
	double Fatigue_SR2;  // Committed strain range 2 at PSUEDO peak - there are potentially two ranges
	double Fatigue_NC2;  // Committed number of cycles at SR2 2 (at PSUEDO peak) - there are potentially two ranges
	double Fatigue_SR3;  // Committed strain range 3 at PSUEDO peak - there are potentially two ranges
	double Fatigue_NC3;  // Committed number of cycles at SR2 3 (at PSUEDO peak) - there are potentially two ranges

	double Fatigue_epsmin;
	double Fatigue_epsmax;

	// Previous State
	double Fatigue_DIP; //Damage index
	double Fatigue_XP; //Range in consideration
	double Fatigue_YP; //Previous Adjacent Range
	double Fatigue_AP; //Peak or valley 1
	double Fatigue_BP; //Peak or valley 2
	double Fatigue_CP; //Peak or valley 2
	double Fatigue_DP; //Peak or valley 4
	int   Fatigue_PCCP; /*Previous Cycle counter flag if >1 then previous 'n'
					   cycles did not flag a complete cycle */
	int   Fatigue_R1FP; //Flag for first  peak count
	int   Fatigue_R2FP; //Flag for second peak count
	double Fatigue_cSlopeP; //Current Slope
	double Fatigue_PSP; //Previous slope
	double Fatigue_EPP; //Previous Strain
	int    Fatigue_SFP; /*Start Flag = 0 if very first strain,
					   (i.e. when initializing)    = 1 otherwise */
	double Fatigue_DLP; //Damage if current strain was last peak.

	double Fatigue_DmaxP;
	double Fatigue_E0P;
	double Fatigue_mP;

	double minStrainP;
	double maxStrainP;

	bool Fatigue_CfailedP;
	double Fatigue_trialStrainP;

	// added 6/9/2006
	// For recording strain ranges (SRXX) and Number of Cycles (NCXX)
	double Fatigue_SR1P;  // Committed strain range at peak
	double Fatigue_NC1P;  // Committed number of cycles at SR1 (i.e. 1.0 or 0.5)
	double Fatigue_SR2P;  // Committed strain range 2 at PSUEDO peak - there are potentially two ranges
	double Fatigue_NC2P;  // Committed number of cycles at SR2 2 (at PSUEDO peak) - there are potentially two ranges
	double Fatigue_SR3P;  // Committed strain range 3 at PSUEDO peak - there are potentially two ranges
	double Fatigue_NC3P;  // Committed number of cycles at SR2 3 (at PSUEDO peak) - there are potentially two ranges

	double Cf;
	double Cd;
	double Alpha;
	double Beta;
	
	double Zd;

	double Lambda_SR;
	double Lambda_SRP;

	double Fatigue_FyInitial;
	double Fatigue_Fy;
	double Fatigue_FyP;
};


#endif

