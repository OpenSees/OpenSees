                                                                       
#ifndef ConcreteD_h
#define ConcreteD_h

// Written: fmk 
//
// Description: This file contains the class definition for 
// ConcreteD. ConcreteD provides the abstraction
// of an elastic perfectly plastic uniaxial material, 
//
// What: "@(#) ConcreteD.h, revA"

#include <UniaxialMaterial.h>

class ConcreteD : public UniaxialMaterial
{
  public:
    ConcreteD(int tag, double fc0, double ec0,double ft0, double eptt0,double Ec0,
		double alphac0 ,double alphat0, double cesp0,double etap0);    
	ConcreteD(int tag, double fc0, double ec0,double ft0, double eptt0,double Ec0,
		double alphac0 ,double alphat0); 
    ConcreteD();    

    ~ConcreteD();

	const char *getClassType(void) const {return "ConcreteD";};

    int setTrialStrain(double strain, double strainRate = 0.0); 
	//int setTrial(double strain,double &stress,double &tangent,double strainRate=0.0);
	//int determineTrialState(double dStrain);

    double getStrain(void);          
    double getStress(void);
    double getTangent(void);
	double getSecant(void);
	void   envelope(void);
	void   unload(void);
    double getInitialTangent(void) {return Ec;};

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
//InPut Parameters
	double fcc;
	double epcc;
	double ft;
	double eptt;
	double Ec;
	double alphac;
	double alphat;
	double cesp;
	double etap;
//State Parameters
	int	   CLoadState;  //load state 0 loading,1 unloading,2 reloading
	double CStrain;     // last committed strain
    double CStress;     // last committed stress
    double CTangent;    // last committed tangent
	double CSecant;     // last committed secant
	double CDc;         // last committed damage parameter in compression
	double CDt;         // last committed damage parameter in tenssion
	double CDcp;
	double CDtp;
	double CEpp;        // last coomited plastic strain
	double CRc;         // last coomited compression threshold
	double CRt;         // last coomited tenssion threshold

//Trial Parameters
	int	   TLoadState;  // trial state 0 loading,1 unloading,2 reloading
	double TStrain;     // trial strain
    double TStress;     // trial stress
    double TTangent;    // trial  tangent
	double TDc;         // trial damage parameter in compression
	double TDt;         // trial damage parameter in tenssion
	double TDcp;
	double TDtp;
	double TEpp;        // trial plastic strain
	double TRc;         // trial compression threshold
	double TRt;         // trial tenssion threshold
	double TSecant;
};


#endif



