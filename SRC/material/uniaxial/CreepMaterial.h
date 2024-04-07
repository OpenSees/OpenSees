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
 // Michael H. Scott
 //
 // Based on TDConcrete implementations by:
 // Adam M. Knaack (adam.knaack@schaefer-inc.com) 
 // Schaefer-Inc, Cincinnati, Ohio, USA
 // Nikola D. Tosic (ntosic@imk.grf.bg.ac.rs)
 // Department for Materials and Structure, Faculty of Civil Engineering, University of Belgrade, Serbia
 // Yahya C. Kurama (ykurama@nd.edu)
 // Department of Civil and Environmental Engineering and Earth Sciences, College of Engineering, University of Notre Dame, Notre Dame, Indiana, USA
 //----------------------------------------------------------------------------------------------------------------------------

 //----------------------------------------------------------------------------------------------------------------------------
 // Description: This file contains the source code of CreepMaterial. 
 // CreepMaterial is a wrapper that imposes creep and shrinkage evoluation equations
 // to any uniaxialMaterial.
 //----------------------------------------------------------------------------------------------------------------------------

#ifndef CreepMaterial_h
#define CreepMaterial_h

#include <UniaxialMaterial.h>
#include <Domain.h> //Added by AMK

class CreepMaterial : public UniaxialMaterial
{
public:
  CreepMaterial(int tag, double _fc, double _fcu, double _epscu, double _ft, double _Ec, double _beta, double _age, double _epsshu, double _epssha, double _tcr, double _epscru, double _epscra, double _epscrd, double _tcast);
  CreepMaterial(int tag, UniaxialMaterial &matl, double _age, double _epsshu, double _epssha, double _tcr, double _epscru, double _epscra, double _epscrd, double _tcast);
  
  CreepMaterial(void);
  
  virtual ~CreepMaterial();
  
  const char *getClassType(void) const {return "CreepMaterial";};    
  double getInitialTangent(void);
  UniaxialMaterial *getCopy(void);
  
  int setTrialStrain(double strain, double strainRate = 0.0); 
  double setCreepStrain(double time, double stress); //Added by AMK
  double getCurrentTime(void); //Added by AMK
  double getStrain(void);
  double getPHI_i(void); //Added by AMK      
  double getStress(void);
  double getTangent(void);
  double getCreep(void); //Added by AMK
  double getMech(void); //Added by AMK
  double setPhi(double time, double tp); //Added by AMK
  double setShrink(double time); //Added by AMK
  double getShrink(void); //Added by AMK
  
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
  UniaxialMaterial *wrappedMaterial;
  
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
  //	double EcT; //Time dependent stiffness
  double age;   // concrete age at first loading, Added by AMK
  double epsshu; // ultimate shrinkage
  double epssha; // shrinkage parameter
  double epscra; // concrete creep exponent parameter
  double epscru; // ultimate concrete creep
  //double sigCr; // stress that creep curve is based on
  double beta;
  double epscrd;
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
  double eps_cr;
  double eps_sh;
  double eps_T;
  double eps_m;
  double epsP_m;
  double epsP_cr;
  double epsP_sh;
  double eps_total;
  double epsP_total;
  double e_total;
  double eP_total;
  double t; //Time
  double t_load; //loaded time
  double phi_i;
  double Et;
  int crack_flag;
  int crackP_flag;
  int iter; //Iteration number

  enum{startSize=500, growSize=200};
  int maxSize;
  
  float *PHI_i;
  float *E_i;
  float *DSIG_i;
  float *TIME_i;
  float *DTIME_i;

  void expandArrays();
};


#endif

