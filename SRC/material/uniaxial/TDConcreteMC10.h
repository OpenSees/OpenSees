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
   
//----------------------------------------------------------------------------------------------------------------------------
// Developed by:
// Nikola D. Tosic (ntosic@imk.grf.bg.ac.rs)
// Department for Materials and Structure, Faculty of Civil Engineering, University of Belgrade, Serbia
// Adam M. Knaack (adam.knaack@schaefer-inc.com) 
// Schaefer-Inc, Cincinnati, Ohio, USA
// Yahya C. Kurama (ykurama@nd.edu)
// Department of Civil and Environmental Engineering and Earth Sciences, College of Engineering, University of Notre Dame, Notre Dame, Indiana, USA
//----------------------------------------------------------------------------------------------------------------------------

//----------------------------------------------------------------------------------------------------------------------------
// Created: 2019
// Last updated: 2019
//----------------------------------------------------------------------------------------------------------------------------

//----------------------------------------------------------------------------------------------------------------------------
// Description: This file contains the source code of TDConcreteMC10. 
// TDConcreteMC10 is a time-dependent concrete material model that calculates
// creep and shrinkage strains.
/*-------------------------------
! Concrete Compression - Linear
! Concrete Tension - Tamai, S., Shima, H., Izumo, J., Okamura, H. 1988. Average Stress-Strain Relationship in Post Yield Range of Steel Bar in Concrete, Concrete Library of JSCE, No. 11, 117-129.
! Concrete Creep - Linear superposition of creep coefficient, Model Code 2010 time function
! Concrete Shrinkage - Model Code 2010 time function
-------------------------------*/
// Detailed descriptions of the model and its implementation can be found in the following:
// (1) Knaack, A.M., Kurama, Y.C. 2018. Modeling Time-Dependent Deformations: Application for Reinforced Concrete Beams with 
//     Recycled Concrete Aggregates. ACI Structural J. 115, 175�190. doi:10.14359/51701153
// (2) Knaack, A.M., 2013. Sustainable concrete structures using recycled concrete aggregate: short-term and long-term behavior
//     considering material variability. PhD Dissertation, Civil and Environmental Engineering and Earth Sciences, University of Notre Dame, Notre Dame, Indiana, USA, 680 pp.
// A manual describing the use of the model and sample files can be found at:
// <https://data.mendeley.com/datasets/z4gxnhchky/3>
//----------------------------------------------------------------------------------------------------------------------------

//----------------------------------------------------------------------------------------------------------------------------
// Disclaimer: This software is provided �as is�, without any warranties, expressed or implied. In no event shall the developers be liable for any claim, damages, or liability arising from or in connection with this software.
//----------------------------------------------------------------------------------------------------------------------------

#ifndef TDConcreteMC10_h
#define TDConcreteMC10_h 

#include <UniaxialMaterial.h>
#include <Domain.h> //Added by AMK

class TDConcreteMC10 : public UniaxialMaterial //ntosic: changed name
{
  public:
    TDConcreteMC10(int tag, double _fc, double _ft, double _Ec, double _Ecm, double _beta, double _age, double _epsba, double _epsbb, double _epsda, double _epsdb, double _phiba, double _phibb, double _phida, double _phidb, double _tcast, double _cem);

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
    double sigCr; // stress that creep curve is based on //ntosic: CHANGE?
    double beta;
    double tcast;
	double cem; //ntosic

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

