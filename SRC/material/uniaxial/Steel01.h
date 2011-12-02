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
                                                                        
// $Revision: 1.1.1.1 $
// $Date: 2000-09-15 08:23:22 $
// $Source: /usr/local/cvs/OpenSees/SRC/material/uniaxial/Steel01.h,v $
                                                                        
                                                                        
#ifndef Steel01_h
#define Steel01_h

// File: Steel01.h
//
// Written: MHS 
// Created: 06/99
// Revision: A
//
// Description: This file contains the class definition for 
// Steel01.h
// 
//
//
// What: "@(#) Steel01.h, revA"


#include <UniaxialMaterial.h>

// Default values for isotropic hardening parameters a1, a2, a3, and a4
#define STEEL_01_DEFAULT_A1        0.0
#define STEEL_01_DEFAULT_A2       55.0
#define STEEL_01_DEFAULT_A3        0.0
#define STEEL_01_DEFAULT_A4       55.0

class Steel01 : public UniaxialMaterial
{
  public:
    Steel01(int tag, double fy, double E0, double b,
       double a1 = STEEL_01_DEFAULT_A1, double a2 = STEEL_01_DEFAULT_A2,
       double a3 = STEEL_01_DEFAULT_A3, double a4 = STEEL_01_DEFAULT_A4,
       double min = NEG_INF_STRAIN, double max = POS_INF_STRAIN);
    Steel01();
    ~Steel01();

    int setTrialStrain(double strain, double strainRate = 0.0); 
    double getStrain(void);              
    double getStress(void);
    double getTangent(void);
    double getSecant (void);

    int commitState(void);
    int revertToLastCommit(void);    
    int revertToStart(void);        

    UniaxialMaterial *getCopy(void);
    
    int sendSelf(int commitTag, Channel &theChannel);  
    int recvSelf(int commitTag, Channel &theChannel, 
		 FEM_ObjectBroker &theBroker);    
    
    void Print(ostream &s, int flag =0);
    
  protected:
    
  private:
    /*** Material Properties ***/
    double fy;  // Yield stress
    double E0;  // Initial stiffness
    double b;   // Hardening ratio (b = Esh/E0)
    double a1;
    double a2;
    double a3;
    double a4;  // a1 through a4 are coefficients for isotropic hardening
    double epsmin;  // Minimum allowable strain
    double epsmax;  // Maximum allowable strain

    /*** Calculated Material Properties ***/
    double epsy;
    double Esh;

    /*** CONVERGED History Variables ***/
    double CminStrain;  // Minimum strain in compression
    double CmaxStrain;  // Maximum strain in tension
    double CshiftP;     // Shift in hysteresis loop for positive loading
    double CshiftN;     // Shift in hysteresis loop for negative loading
    int Cloading;       // Flag for loading/unloading
                        // 1 = loading (positive strain increment)
                        // -1 = unloading (negative strain increment)
                        // 0 initially
    int Cfailed;        // Flag indicating failure ...  a committed strain
                        // has fallen out of the range [epsmin, epsmax]


    /*** CONVERGED State Variables ***/    
    double Cstrain;
    double Cstress;
    double Ctangent;    

    /*** TRIAL History Variables ***/
    double TminStrain;
    double TmaxStrain;
    double TshiftP;
    double TshiftN;
    int Tloading;
    int Tfailed;
    
    /*** TRIAL State Variables ***/
    double Tstrain;
    double Tstress;
    double Ttangent; // Not really a state variable, but declared here
                     // for convenience

    // Sets trial history variables to converged history values
    void setHistoryVariables ();

    // Calculates the trial state variables based on the trial strain
    void determineTrialState (double dStrain);

    // Determines if a load reversal has occurred based on the trial strain
    void detectLoadReversal (double dStrain);
};

#endif
