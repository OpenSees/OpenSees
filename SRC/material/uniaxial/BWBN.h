                                                                        
#ifndef BWBN_h
#define BWBN_h


#include <UniaxialMaterial.h>
#include <Matrix.h>

//#define MAT_TAG_BWBN 5002

class BWBN : public UniaxialMaterial
{
  public:
    BWBN(int tag, 
	 double alpha,
	 double ko,
	 double n,
	 double gamma,
	 double beta,
	 double Ao,
	 double q,
	 double zetas,
	 double p,
	 double Shi,
	 double deltaShi,
	 double lamda,
	 double tolerance,
	 int maxNumIter);
    BWBN();	
    ~BWBN();

    const char *getClassType(void) const {return "BoucWenMaterial";};

    int setTrialStrain(double strain, double strainRate = 0.0); 
    double getStrain(void);          
    double getStress(void);
    double getTangent(void);
    double signum(double);
    int commitState(void);
    int revertToLastCommit(void);    
    int revertToStart(void);        
    UniaxialMaterial *getCopy(void);
    int sendSelf(int commitTag, Channel &theChannel);  
    int recvSelf(int commitTag, Channel &theChannel, 
		 FEM_ObjectBroker &theBroker);    
    void Print(OPS_Stream &s, int flag =0);
    
    double getInitialTangent(void);

  protected:
    
  private:

    // Material parameters
    double alpha;
    double ko;
    double n;
    double gamma;
    double beta;
    double Ao;
    double q;
    double zetas;
    double p;
    double Shi;
    double deltaShi;
    double lamda;
    
    // History variables (trial and committed)
    double Tstrain, Cstrain;
    double Tz, Cz;
    double Te, Ce;
    
    // Other variables
    double Tstress, Ttangent;
    
    double tolerance;
    int maxNumIter;
};


#endif



