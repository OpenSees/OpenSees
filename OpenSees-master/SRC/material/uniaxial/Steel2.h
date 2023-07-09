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
// $Source: /usr/local/cvs/OpenSees/SRC/material/uniaxial/Steel2.h,v $
                                                                      
// Written: fmk
// Created: 03/06
//
// Description: This file contains the class definition for 
// Steel2. Steel2 is based on an f2c of the FEDEAS material
// Steel02.f which is:
//-----------------------------------------------------------------------
// MENEGOTTO-PINTO STEEL MODEL WITH FILIPPOU ISOTROPIC HARDENING
//            written by MOHD YASSIN (1993)
//          adapted to FEDEAS material library
//    by E. Spacone, G. Monti and F.C. Filippou (1994)
//-----------------------------------------------------------------------


#ifndef Steel2_h
#define Steel2_h

#include <UniaxialMaterial.h>

class Steel2 : public UniaxialMaterial
{
  public:
    Steel2(int tag,
	    double fy, double E0, double b,
	    double R0, double cR1, double cR2,
	    double a1, double a2, double a3, double a4, double sigInit =0.0);
    
    // Constructor for no isotropic hardening
    Steel2(int tag,
	    double fy, double E0, double b,
	    double R0, double cR1, double cR2);
    
    // Constructor for no isotropic hardening
    // Also provides default values for R0, cR1, and cR2
    Steel2(int tag, double fy, double E0, double b);
	    
    Steel2(void);
    virtual ~Steel2();
    

    const char *getClassType(void) const {return "Steel2";};

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

    // 
    double epsminP;
    double epsmaxP;
    
    double epsplP;  //  = hstvP(3) : plastic excursion
    double epss0P;  //  = hstvP(4) : eps at asymptotes intersection
    double sigs0P;  //  = hstvP(5) : sig at asymptotes intersection
    double epssrP;  //  = hstvP(6) : eps at last inversion point
    double sigsrP;  //  = hstvP(7) : sig at last inversion point

	double TepsplP;
    double Tepss0P;
    double Tsigs0P;
    double TepssrP;
    double TsigsrP;

    double CepsplP;
    double Cepss0P;
    double Csigs0P;
    double CepssrP;
    double CsigsrP;

    int    konP;
    double epsP;
    double sigP;
    double eP;
	
    //
    double epsmin; 
    double epsmax; 

	double epspl;  
    double epss0;  
    double sigs0; 
    double epsr;  
    double sigr;  

    double Tepspl;  
    double Tepss0;  
    double Tsigs0; 
    double Tepsr;  
    double Tsigr;  

    double Cepspl;  
    double Cepss0;  
    double Csigs0; 
    double Cepsr;  
    double Csigr;  

    int    kon;    
    double sig;   
    double e;     
    double eps;   //  = strain at current step
};


#endif

