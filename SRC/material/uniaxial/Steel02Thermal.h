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
                                                                        
//JZ. 05/2010


#ifndef Steel02Thermal_h
#define Steel02Thermal_h

#include <UniaxialMaterial.h>

class Steel02Thermal : public UniaxialMaterial
{
  public:
    Steel02Thermal(int tag,
		   double fy, double E0, double b,
		   double R0, double cR1, double cR2,
		   double a1, double a2, double a3, double a4, double sigInit =0.0);
    
    // Constructor for no isotropic hardening
    Steel02Thermal(int tag,
		   double fy, double E0, double b,
		   double R0, double cR1, double cR2);
    
    // Constructor for no isotropic hardening
    // Also provides default values for R0, cR1, and cR2
    Steel02Thermal(int tag, double fy, double E0, double b);
	    
    Steel02Thermal(void);
    virtual ~Steel02Thermal();
    

    const char *getClassType(void) const {return "Steel02Thermal";};

    double getInitialTangent(void);
    UniaxialMaterial *getCopy(void);

    //JZ this function is no use, just for the definiation of pure virtual function.
    int setTrialStrain(double strain, double strainRate);

    int setTrialStrain(double strain, double FiberTemperature, double strainRate); 

    double getStrain(void);      
    double getStress(void);
    double getTangent(void);
	
    double getThermalElongation(void); //***JZ
    double getElongTangent(double, double&, double&, double);//***JZ //PK add to include max temp
    
    int commitState(void);
    int revertToLastCommit(void);    
    int revertToStart(void);        
    
    int sendSelf(int commitTag, Channel &theChannel);  
    int recvSelf(int commitTag, Channel &theChannel, 
		 FEM_ObjectBroker &theBroker);    
    
    void Print(OPS_Stream &s, int flag =0);

    int getVariable(const char *variable, Information &);
    
 protected:
    
 private:

//JZ 11/10 /////////////////////////////////////////////////////////////start
   	double Temp;  // material temp
	//double steps;    //the amount of the steps. 
    double ThermalElongation; // eps(theata) = alpha * temperature 
	double FyT;
	double E0T;
//JZ 11/10 /////////////////////////////////////////////////////////////end 

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
    double E0P;  // Initial stiffness in last committed step;
	double FyP;   //Yield stress in last committed step;
	double FiberTP; //FiberTemperature in last committed step;

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

