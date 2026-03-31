// $Source: /usr/local/cvs/OpenSees/SRC/material/uniaxial/Concrete06.h,v $
// Created: 07/06
// Modified by: LMS 
// Description: This file contains the class implementation for Concrete06. Based on Concrete01

                                                                        
#ifndef Concrete06_h
#define Concrete06_h
#include <UniaxialMaterial.h>

class Concrete06 : public UniaxialMaterial
{
 public:
  Concrete06 (int tag, double fc, double eo, double r, double k, double alphaC, double fcr, double ecr, double b, double alphaT);
  Concrete06 ();
  ~Concrete06();

  const char *getClassType(void) const {return "Concrete06";};    

  int setTrialStrain(double strain, double strainRate = 0.0); 
  double getStrain(void);      
  double getStress(void);
  double getTangent(void);
  double getInitialTangent(void) { return fc/eo*r/(r-1.0);} //initial in compression
	  
  //  double getec(void) {return eo;};
  
  int commitState(void);
  int revertToLastCommit(void);    
  int revertToStart(void);        
  
  UniaxialMaterial *getCopy(void);
  
  int sendSelf(int commitTag, Channel &theChannel);  
  int recvSelf(int commitTag, Channel &theChannel, 
	       FEM_ObjectBroker &theBroker);    
  
  void Print(OPS_Stream &s, int flag =0);

    int getVariable(const char *variable, Information &);

   protected:

   private:
      /*** Material Properties ***/
	double ecr;
	double fcr;
	double b;
	double fc;
	double eo;
	double r;
	double k;
	double alphaC;
	double alphaT;

      /*** CONVERGED History Variables ***/
	double Cecmax;
	double Cet;
	double CetAccum;
	double Cscmax;
	double Cet1;
	double Cet2;
	double Cstmax;
	double Cetmax;
	double CEt;
	double CEr1;
	double CEr2;


      /*** CONVERGED State Variables ***/
    double Cstrain;
    double Cstress;   
    double Ctangent;	
					
      /*** TRIAL History Variables ***/
	double ecmax;
	double et;
	double etAccum;
	double scmax;
	double et1;
	double et2;
	double stmax;
	double etmax;
	double Et;
	double Er1;
	double Er2;

      /*** TRIAL State Variables ***/
    double Tstrain;
    double Tstress;
    double Ttangent; // 

	double Eci; 
    double Eti;
    double ecref;
    double scref ;
    double etref;
    double stref;

    void envelopeC(double e);
    void envelopeT(double e);
	void DefLoop(double Er);

// AddingSensitivity:BEGIN //////////////////////////////////////////
    int parameterID;
	Matrix *SHVs;
// AddingSensitivity:END ///////////////////////////////////////////
};


#endif


