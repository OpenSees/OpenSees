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
                                                                       
// $Revision: 1.2                                                                              
// $Date: 03/2017                                                                              
// $Source:                                                                                    

// Written: R.Rahimi & R.Sepasdar & Dr. M. R. Banan                                            
// Created: 09/2012                                                                            
//                                                                                             
//                                                                                             
// Description: This file contains the class implementation of RambergOsgoodSteel.             
//-----------------------------------------------------------------------                      
//              RAMBERG-OSGOOD STEEL MODEL                                                     
//      Developed  by REZA RAHIMI, (Reza.Rahimi@Dal.Ca)   (2012)                               
//                                                                                             
//      Co-Developer: REZA SEPASDAR,                                                           
//      Supervisor: Dr. Mo. R. Banan,                                                          
//----------------------------------------------------------------------- 



#ifndef RambergOsgoodSteel_h
#define RambergOsgoodSteel_h

#include <UniaxialMaterial.h>

class RambergOsgoodSteel : public UniaxialMaterial
{
  public:

    RambergOsgoodSteel(int tag, double fy, double E0, double rezaA, double rezaN);
    RambergOsgoodSteel(int tag, double fy, double E0);
    RambergOsgoodSteel(void);

    virtual ~RambergOsgoodSteel();

    const char *getClassType(void) const {return "RambergOsgoodSteel";};

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

    int iii;
    int ki;
    // matpar : STEEL FIXED PROPERTIES
    double Fy;  //  = matpar(1)  : yield stress
    double E0;  //  = matpar(2)  : initial stiffness
    double rezaAA;   //  = matpar(3)  : hardening ratio (Esh/E0)
    double rezaNN;  //  = matpar(4)  : exp transition elastic-plastic
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
};


#endif
