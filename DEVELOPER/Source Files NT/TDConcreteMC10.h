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
   
 // Written: NT
 // Created: 2019
 //
 // Description: This file contains the class definition for TDConcreteMC10. 
 // TDConcreteMC10 is a time-dependent concrete material model that calculates
 // creep and shrinkage strains.
 /*-------------------------------
 ! Concrete Compression - Linear
 ! Concrete Tension - Tamai et al. (1988) "Average Stress-Strain Relationship in Post Yield Range of Steel Bar in Concrete"
 ! Concrete Creep - Linear superposition of creep coefficient, Model Code 2010 time function
 ! Concrete Shrinkage - Model Code 2010 time function
 -------------------------------*/
 //
 // The framework for this code was originally modified from Concrete02.

#ifndef TDConcreteMC10_h
#define TDConcreteMC10_h 

#include <UniaxialMaterial.h>
#include <Domain.h> //Added by AMK

class TDConcreteMC10 : public UniaxialMaterial //ntosic: changed name
{
  public:
    TDConcreteMC10(int tag, double _fc, double _ft, double _Ec, double _Ecm, double _beta, double _age, double _epsba, double _epsbb, double _epsda, double _epsdb, double _tcr, double _phiba, double _phibb, double _phida, double _phidb, double _phidc, double _tcast);

    TDConcreteMC10(void);

    virtual ~TDConcreteMC10();

    const char *getClassType(void) const {return "TDConcreteMC10";};    
    double getInitialTangent(void);
    UniaxialMaterial *getCopy(void);

    int setTrialStrain(double strain, double strainRate = 0.0); 
	double setCreepBasicStrain(double time, double stress); //Added by AMK //ntosic: split into basic and drying creep
	double setCreepDryingStrain(double time, double stress); //Added by AMK //ntosic: split into basic and drying creep
    double setStress(double strain, double &stiff); //Added by AMK
	double getCurrentTime(void); //Added by AMK
    double getStrain(void);
    double getPHIB_i(void); //Added by AMK //ntosic: split into basic and drying creep
	double getPHID_i(void); //Added by AMK //ntosic: split into basic and drying creep
    double getStress(void);
    double getTangent(void);
	double getCreepBasic(void); //Added by AMK //ntosic: split into basic and drying creep
	double getCreepDrying(void); //Added by AMK //ntosic: split into basic and drying creep
	double getMech(void); //Added by AMK 
	double setPhiBasic(double time, double tp); //Added by AMK //ntosic: split into basic and drying creep
	double setPhiDrying(double time, double tp); //Added by AMK //ntosic: split into basic and drying creep
	double setShrinkBasic(double time); //Added by AMK //ntosic: split into basic and drying shrinkage
	double setShrinkDrying(double time); //Added by AMK //ntosic: split into basic and drying shrinkage
	double getShrinkBasic(void); //Added by AMK //ntosic: split into basic and drying
	double getShrinkDrying(void); //Added by AMK //ntosic: split into basic and drying
    double getKappa(void);
    
    int commitState(void);
    int revertToLastCommit(void);    
    int revertToStart(void);        
    
    int sendSelf(int commitTag, Channel &theChannel);  
    int recvSelf(int commitTag, Channel &theChannel, 
		 FEM_ObjectBroker &theBroker);    
    
    void Print(OPS_Stream &s, int flag =0);

    int getVariable(const char *variable, Information &);
	
	
	//Added by AMK for recording Creep and Mechanical Strain:
	Response *setResponse(const char **argv, int argc,OPS_Stream &theOutput);
	int getResponse(int responseID, Information &matInfo);
    
 protected:
    
 private:
    void Tens_Envlp (double epsc, double &sigc, double &Ect);
    void Compr_Envlp (double epsc, double &sigc, double &Ect);

    // matpar : Concrete FIXED PROPERTIES
    double fc;    // concrete compression strength           : mp(1)
    //	double fcT;  //Time Dependent Strength
    double epsc0; // strain at compression strength          : mp(2)
    double fcu;   // stress at ultimate (crushing) strain    : mp(3)
    double epscu; // ultimate (crushing) strain              : mp(4)       
    double tcr;   // creep relationship age
    double ft;    // concrete tensile strength               : mp(6)
    //	double ftT;	//Time dependent strength (tension)
    double Ets;   // tension stiffening slope                : mp(7)
	double Ec;  //Concrete stiffness, Added by AMK
	double Ecm; //ntosic 28-day modulus, necessary for normalizing creep coefficient
	//	double EcT; //Time dependent stiffness
	double age;   // concrete age at first loading, Added by AMK
	double epsba; //ntosic
    double epsbb; //ntosic
	double epsda; //ntosic
	double epsdb; //ntosic
    double phiba; //ntosic
	double phibb; //ntosic
	double phida; //ntosic
	double phidb; //ntosic
	double phidc; //ntosic
    double sigCr; // stress that creep curve is based on //ntosic: CHANGE?
    double beta;
    double tcast;

    // hstvP : Concerete HISTORY VARIABLES last committed step
    double ecminP;  //  hstP(1)
    double ecmaxP;  // added by AMK
    double deptP;   //  hstP(2)
    double epsP;  //  = strain at previous converged step
    double sigP;  //  = stress at previous converged step
    double eP;    //   stiffness modulus at last converged step;

    // hstv : Concerete HISTORY VARIABLES  current step
    double ecmin;
    double ecmax; // added by AMK  
    double dept;   
    double sig;   
    double e;     
    double eps;
	
	//Added by AMK:
	int count;
	double epsInit;
	double sigInit;
	double eps_crb; //ntosic: split into basic and drying creep
	double eps_crd; //ntosic: split into basic and drying creep
	double eps_shb; //ntosic: split into basic and drying shrinkage
	double eps_shd; //ntosic: split into basic and drying shrinkage
	double eps_T;
	double eps_m;
	double epsP_m;
	double epsP_crb; //ntosic: split into basic and drying creep
	double epsP_crd; //ntosic: split into basic and drying creep
    double epsP_shb; //ntosic: split into basic and drying shrinkage
	double epsP_shd; //ntosic: split into basic and drying shrinkage
	double eps_total;
	double epsP_total;
	double e_total;
	double eP_total;
	double t; //Time
	double t_load; //loaded time
	double phib_i; //ntosic: split into basic and drying creep
	double phid_i; //ntosic: split into basic and drying creep
	double Et;
	int crack_flag;
	int crackP_flag;
    int iter; //Iteration number
    
    float PHIB_i[5000]; //ntosic: split into basic and drying creep
	float PHID_i[5000]; //ntosic: split into basic and drying creep
    float E_i[5000];
    float DSIG_i[5000];
    float dsig_i[5000];
    float TIME_i[5000]; //Time from the previous time step
    float DTIME_i[5000];
    float kappa[5000];
};


#endif

