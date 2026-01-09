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
//
// $Date: 2014/04/13 12:26
// 
// Written: Erkan Bicici 
//
// Description: The file contains the application of shear model 
//              developed by Sezen (2002). The material is defined 
//              with sezen point.  
//
                                                                      
#ifndef Shear01_h
#define Shear01_h

#include <UniaxialMaterial.h>

class Shear01 : public UniaxialMaterial
{
  public:
	Shear01(int tag, double ep, double crkp, double e2p, double eypp, double eyep, double eafp,  double en, double crkn, double e2n, double eypn, double eyen, double eafn);
   
    Shear01();    

    ~Shear01();

    int setTrialStrain(double strain, double strainRate = 0.0); 
    double getStrain(void);          
    double getStress(void);
    double getTangent(void);

    double getInitialTangent(void) {return Ep;};

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
   
    double ezero;	// initial strain
    double Ep;		// elastic modulus
    double ecrp;
    double E2p;
    double epp;		// plastic strain at last commit
    double eyp;
    double eup; 
    double Vcrp;
    double Vmaxp;
    double Edegp; 

    double En;		// elastic modulus
    double ecrn;
    double E2n;
    double epn;		// plastic strain at last commit
    double eyn;
    double eun; 
    double Vcrn;
    double Vmaxn;
    double Edegn;  
	double kpos;
	double kneg;
	double k2pos;
	double k2neg;
    double Edegp_2;
    double Edegn_2;

	int TIndicatorNcycle;
	int CIndicatorNcycle;
	//double EnvIndicator;
	int CEnvIndicator;
	int TEnvIndicator;


double ncycle;
double Cncycle;
double Tncycle;


    double TmaxStrain;   // maximum recorded strain
	double CmaxStrain;   // maximum recorded strain
	
	double TmaxStress;   // maximum recorded stress 
    double CmaxStress;   // maximum recorded stress 
    
	double TminStrain;   // minimum recorded strain
    double CminStrain;   // minimum recorded strain
    
	double TminStress;   // minimum recorded stress 
    double CminStress;   // minimum recorded stress 



	double TenvmaxStress;
	double CenvmaxStress;
	double TenvminStress;
	double CenvminStress;









    double trialStrain;	  // trial strain
    double trialStress;   // current trial stress
    double trialTangent;  // current trial tangent
    double commitStrain;  // last commited strain
    double commitStress;  // last commited stress
    double commitTangent; // last committed  tangent

    double Dif_Strain;
    double Cur_Dif_Strain;      // current difference between strain
    double Old_Dif_Strain;      // last committed difference between strains
	
	void setParameters(void);


	double PosEnvTangent(double trialStrain);
	double NegEnvTangent(double trialStrain);	
	double PosEnvStress(double trialStrain);
	double NegEnvStress(double trialStrain);	    
	void unloading(double Dif_Strain);
	void reloading(double Dif_Strain);

	void PositiveIncycle(double Dif_Strain);
	void NegativeIncycle(double Dif_Strain);

	void NumberCycle(void);
	double beta;

	void Envelope_Check(void);



};


#endif



