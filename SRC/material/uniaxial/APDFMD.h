
#ifndef APDFMD_h
#define APDFMD_h

#include <UniaxialMaterial.h>

class APDFMD : public UniaxialMaterial
{
  public:

    APDFMD(int tag, double fy1, double E1, double ad2, double fy2 , double E2, double rezaA2, double rezaN2);
    APDFMD(void);
 
   
    virtual ~APDFMD();

    const char *getClassType(void) const {return "APDFMD";};

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

     double findstress1(double targetStrain, double K, double mq, double yi, double rezaAA, double rezaNN);
     double findstressp(double targetStrain, double mz, double K, double mq, double rezaAA, double rezaNN);
     double findstressn(double targetStrain, double mz, double K, double mq, double rezaAA, double rezaNN);
   
   
    double ad2;
    double Fy1;  
    double Fy2;
    double E1;  
    double E2;
    double rezaAA2;
    double rezaNN2;
    double sigini; 


    double epsmin1;
    double epsmax1;
    double epsmin2;
    double epsmax2;
    double epspl;
    double epssA;
    double sigsA;
    double epssC;
    double sigsC;
    double stap;
    double stan;

    double epsr;
    double sigr;

    int    kon;
    double sig;
    double TrotMax;
    double eps;


    double epsminP1; 
    double epsmaxP1; 
    double epsminP2; 
    double epsmaxP2; 
    double epsplP;  
    double epssAP;  
    double sigsAP;  
    double epssCP;  
    double sigsCP;  
    double epssrP; 
    double sigsrP;  
    int    konP;    


    double epsP;  
    double sigP; 
    double CrotMax;

   
};


#endif
