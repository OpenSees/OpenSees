                                                                        
#ifndef DegradingPinchedBW_h
#define DegradingPinchedBW_h


#include <UniaxialMaterial.h>
#include <Matrix.h>



class DegradingPinchedBW : public UniaxialMaterial
{
  public:
    DegradingPinchedBW(int tag, 
	 double m,
	 double Fy,
	 double xu,
	 double alpha,
	 double ko,
	 double n,
	 double eta,
	 double beta,
	 double rhoeps,
	 double rhox,
	 double phi,
	 double deltak,
	 double deltaf,
	 double sigma,
	 double u,
	 double epsp,
	 double rhop,
	 double tolerance,
	 int maxNumIter);
	
    ~DegradingPinchedBW();

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
	double m;
	double Fy;
	double xu;
    double alpha;
    double ko;
    double n;
    double eta;
    double beta;
	double rhoeps;
	double rhox;
	double phi;
	double deltak;
	double deltaf;
	double sigma;
	double u;
	double epsp;
	double rhop;
    
    // History variables (trial and commited)
	double xmaxp;
    double xmax;
	double Tstrain, Cstrain;
    double Tz, Cz;
    double Te, Ce;
    
    // Ohter variables
    double Tstress, Ttangent;
    
    double tolerance;
    int maxNumIter;
};


#endif



