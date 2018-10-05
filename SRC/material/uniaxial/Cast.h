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
// $Source: /usr/local/cvs/OpenSees/SRC/material/uniaxial/Cast.h,v $
                                                                      
// Written: Dimitrios Lignos
// Created: 04/2011
//
// Description: This file contains the class implementation of Cast. 
// This Cast is a material model that describes the Cast Fuse brace
// Developed by Gray et al. (2010).
// this model is based on:
// 1. Modified Menegotto-Pinto Steel Model with Filippou Isotropic Hardening
// 2. Closed Form solution for monotonic loading developed by Gray et al. (2010)
//-----------------------------------------------------------------------



#ifndef Cast_h
#define Cast_h

#include <UniaxialMaterial.h>

class Cast : public UniaxialMaterial
{
  public:
    Cast(int tag, double nLegs, double bo, 
		double h, double fy, double eo, double l, double b,
	    double R0, double cR1, double cR2,
	    double a1, double a2, double a3, double a4);
    
    // Constructor for no isotropic hardening
    Cast(int tag, double nLegs, double bo, double h, 
		 double fy, double eo, double l, double b,
	    double R0, double cR1, double cR2);
    
    // Constructor for no isotropic hardening
    // Also provides default values for R0, cR1, and cR2
    Cast(int tag, double nLegs, double bo, double h, 
		 double fy, double eo, double l, double b);
	    
    Cast(void);
    virtual ~Cast();
    

    const char *getClassType(void) const {return "Cast";};

    double getInitialTangent(void);
    UniaxialMaterial *getCopy(void);

    int setTrialStrain(double strain, double strainRate = 0.0); 
    double getStrain(void);      
    double getStress(void);
    double getTangent(void);
    
    int commitState(void);
    int revertToLastCommit(void);    
    int revertToStart(void);        
    
    int sendSelf(int commitTag, Channel &theChannel);  
    int recvSelf(int commitTag, Channel &theChannel, 
		 FEM_ObjectBroker &theBroker);    
    
    void Print(OPS_Stream &s, int flag =0);
    
 protected:
    
 private:
    // additional parameters addeded by Dimitrios Lignos
    double nLegs; // Number of fingers
	double bo;    // width of the individual yielding finger
	double h;     // thickness of individual finger
	double fy;    // Yield stress of the finger
	double eo;    // modulus of elasticity of a steel finger
    double l;     // Length of finger
	
	// standard parameters by original menegotto-pinto material
	double b;   // Hardening ratio (b = Esh/E0)
    double R0;  //  = matpar(4)  : exp transition elastic-plastic
    double cR1; //  = matpar(5)  : coefficient for changing R0 to R
    double cR2; //  = matpar(6)  : coefficient for changing R0 to R
    double a1;  //  = matpar(7)  : coefficient for isotropic hardening in compression
    double a2;  //  = matpar(8)  : coefficient for isotropic hardening in compression
    double a3;  //  = matpar(9)  : coefficient for isotropic hardening in tension
    double a4;  //  = matpar(10) : coefficient for isotropic hardening in tension
    
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
	
	// added by DL
	double kp;   // Initial stiffness of Cast Fuse
	double Pp;   // Monotonic yield strength of Cast Fuse
	
	double epsminr;   // added by MG
	double epsmaxr;   // added by MG
	double epsminrP;   // added by MG
	double epsmaxrP;   // added by MG
};


#endif

