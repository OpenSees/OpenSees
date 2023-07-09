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
// $Date: 2006-08-03 23:42:19 $
// $Source: /usr/local/cvs/OpenSees/SRC/material/uniaxial/Steel03.h,v $
                                                                        
                                                                        
#ifndef Steel03_h
#define Steel03_h

// File: Steel03.h
//
// Written: mackie
// Created: 06/2005
// Revision: A
//
// Description: This file contains the class definition for 
// Steel03.h. Steel03 is Steel01 verbatim but with added Giuffre-Menegotto-Pinto 
// transitions on the loading and unloading loops.  
// references:
// 1.) 	Menegotto, M., and Pinto, P.E. (1973). Method of analysis of cyclically loaded 
//	RC plane frames including changes in geometry and non-elastic behavior of 
//	elements under normal force and bending. Preliminary Report IABSE, vol 13. 
// 2.)	Dhakal, R.J., and Maekawa, K. (2002). Path-dependent cyclic stress-strain relationship
//	of reinforcing bar including buckling. Engineering Structures, 24(11): 1383-96.
// 3.)	Gomes, A., and Appleton, J. (1997). Nonlinear cyclic stress-strain relationship of 
//	reinforcing bars including buckling. Engineering Structures, 19(10): 822-6.

// 
//
//
// What: "@(#) Steel03.h, revA"


#include <UniaxialMaterial.h>

// Default values for isotropic hardening parameters a1, a2, a3, and a4
#define STEEL_03_DEFAULT_A1        0.0
#define STEEL_03_DEFAULT_A2       55.0
#define STEEL_03_DEFAULT_A3        0.0
#define STEEL_03_DEFAULT_A4       55.0

class Steel03 : public UniaxialMaterial
{
  public:
    Steel03(int tag, double fy, double E0, double b, double r, double cR1, double cR2,
       double a1 = STEEL_03_DEFAULT_A1, double a2 = STEEL_03_DEFAULT_A2,
       double a3 = STEEL_03_DEFAULT_A3, double a4 = STEEL_03_DEFAULT_A4);
    Steel03();
    ~Steel03();

    const char *getClassType(void) const {return "Steel03";};

    int setTrialStrain(double strain, double strainRate = 0.0); 
    int setTrial (double strain, double &stress, double &tangent, double strainRate = 0.0);
    double getStrain(void);              
    double getStress(void);
    double getTangent(void);
    double getInitialTangent(void) {return E0;};

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
    /*** Material Properties ***/
    double fy;  // Yield stress
    double E0;  // Initial stiffness
    double b;   // Hardening ratio (b = Esh/E0)
    double r;	// radius of rounded corners
    double cR1;
    double cR2;
    double a1;
    double a2;
    double a3;
    double a4;  // a1 through a4 are coefficients for isotropic hardening

    /*** CONVERGED History Variables ***/
    double CminStrain;  // Minimum strain in compression
    double CmaxStrain;  // Maximum strain in tension
    double CshiftP;     // Shift in hysteresis loop for positive loading
    double CshiftN;     // Shift in hysteresis loop for negative loading
    int Cloading;       // Flag for loading/unloading
                        // 1 = loading (positive strain increment)
                        // -1 = unloading (negative strain increment)
                        // 0 initially
    double CbStrain;
    double CbStress;
    double CrStrain;
    double CrStress;
    double Cplastic;

    /*** CONVERGED State Variables ***/    
    double Cstrain;
    double Cstress;
    double Ctangent;
    double CcurR;

    /*** TRIAL History Variables ***/
    double TminStrain;
    double TmaxStrain;
    double TshiftP;
    double TshiftN;
    int Tloading;
    double TbStrain;
    double TbStress;
    double TrStrain;
    double TrStress;
    double Tplastic;
    
    /*** TRIAL State Variables ***/
    double Tstrain;
    double Tstress;
    double Ttangent; // Not really a state variable, but declared here
                     // for convenience
    double TcurR;

    // Calculates the trial state variables based on the trial strain
    void determineTrialState (double dStrain);
    double getR (double x_in);

};

#endif
