// Written: KK
// Revision: A                                                                 
                                                                        
#ifndef Concrete0_h
#define Concrete0_h

#include <UniaxialMaterial.h>

class Concrete0 : public UniaxialMaterial
{
  public:
    Concrete0(int tag, double E0C, double FPC, double E00); 
	Concrete0(int tag, double E0C, double FPC, double E00, double E0T, double FT); 
	Concrete0(int tag, double E0C, double FPC, double E00, double E0T, double FT, int MON); 
    Concrete0();    
    ~Concrete0();

	const char *getClassType(void) const {return "Concrete0";};

    int setTrialStrain(double strain, double strainRate); 
    int setTrial(double strain, double &stress, double &tangent, double strainRate); 
    double getStrain(void) {return Tstrain;};
    double getStress(void);
    double getTangent(void);
    double getInitialTangent(void);

	void determineTrialState(void);
	void compEnvelope(void);

    int commitState(void);
    int revertToLastCommit(void);    
    int revertToStart(void);        

    UniaxialMaterial *getCopy(void);
    
    int sendSelf(int commitTag, Channel &theChannel);  
    int recvSelf(int commitTag, Channel &theChannel, 
		 FEM_ObjectBroker &theBroker);    
    
    void Print(OPS_Stream &s, int flag =0);
    
    int setParameter(const char **argv, int argc, Parameter &param);
	int updateParameter(int parameterID, Information &info);

	Response *setResponse (const char **argv, int argc, OPS_Stream &theOutputStream); // KK
	int getResponse (int responseID, Information &matInformation);  // KK
	Vector getInputParameters(void); // KK

  protected:
    
  private:
	// Input parameters
	double e0c;
	double e0t;
	double fpc;
	double e00;
	double ft;

	int mon;

	// Compression envelope parameters
	double e0c033;
	double e0c066;
	double fpc033;
	double fpc066;
    double Ec0;
	double Ec033;
	double Ec066;
	double Ec00;

	// Tension envelope parameters
	double Et;

	// Tral state parameters
	double Tstrain;
	double Tstress;
	double Ttangent;

	// Committed state parameters
	double Cstrain;
	double Cstress;
	double Ctangent;
	double CompStiff;

	// History parameters
	int Cracking;
	int Crushing;

	double epsnmin;
	double signmin;

	double SMALL;

};


#endif

