#ifndef SPSW02_h
#define SPSW02_h
#include <UniaxialMaterial.h>
class SPSW02:public UniaxialMaterial
{
public:
	SPSW02(int tag, double fpy, double E0, double b, double t, double hs,
		double l, double R, double epsPCFac, double pstCapEFac, double gama, double c, double resFac);
	SPSW02(int tag, double E0, double b, double _FTS, double _FCS, double _cmpUnldngEFac,
		double _sigTEFac, double _sigTFfac, double _epsTFfac, double R, double epsPCFac,
		double pstCapEFac, double gama, double c, double resFac);
	
	SPSW02(void);
	virtual ~SPSW02();
	const char* getClassType(void) const {return "SPSW02";};
	double getInitialTangent();
	UniaxialMaterial* getCopy();
	virtual int setTrialStrain(double strain, double strainRate = 0.0);
	virtual double getStrain();
	virtual double getStress();
	virtual double getTangent();
	virtual int commitState();
	virtual int revertToLastCommit();
	virtual int revertToStart();
	virtual int sendSelf(int commitTag, Channel& theChannel);
	virtual int recvSelf(int commitTag, Channel& theChannel, FEM_ObjectBroker & theBroker);
	virtual void Print(OPS_Stream & s, int flag = 0);
    //virtual Response *setResponse (const char **argv, int argc, 
				//   OPS_Stream &theOutputStream);
    //virtual int getResponse (int responseID, Information &matInformation);
	virtual double getEnergy();
	virtual double getInitYieldStrain() {return FTS/E0;}

private:
	void updateDamage();
private:
	void Calc_sigcr (/*double& Fts, double & Fcs*/);
	void MenegottoPinto (double epsc, double bT, double R, double& sigc, double& ec);
	//void MenegottoPinto (double epsc, double bT, double R, double& sigc, double& ec, double er, double sr, double e0, double s0);
	//matpar : PLATE PROPERTIES
	double t;			//plate thickness
	double hs;			//plate height
	double l;			//plate length
	double fpy	;			//yield stress
	double E0	;			//initial stiffness
	double b	;			//hardening ratio
	double R;
	double Fts	;			//trial yield stress of a tension strip (damage not commited)
	double Fcs	;			//yield stress of a compression strip
	double FTS, FCS;		//yield stresses before damage
	double epsPCFac	;		//ratio between post cap strain and yield strain
	double pstcpEFac;		//ratio between post cap stifness and initial stiffness
	double gama, FailEnerg, c;			//damage parameters
	double resFac;
	//Commit HISTORY VARIABLES
	double epsmaxP	;		//max eps in tension
	double sigmaxP	;		//max sig in tension
	double epss0P	;		//eps at asymptotes intersection
	double sigs0P	;		//sig at ...
	double epssrP	;		//eps at last inversion point
	double sigsrP	;		//sig at ...
	double epsTFP	;		//eps at tension field redevelopment point
	double plstrP	;		//plastic strain
	int konP	;			//index for loading-unloading
	double epsP	;			//eps at previous converged step
	double sigP	;			//stress at ...
	double eP	;			//stiffness modulus at ...

	//trial history variables
	double epsmax	;		//
	double sigmax	;		//
	double epss0	;		//
	double sigs0	;		//
	double epsr		;		//
	double sigr		;		//
	double epsTF	;		//
	double plstr	;		//
	int    kon		;		//
	double sig		;		//
	double e		;		//
	double eps		;		//
	//bool capStrFlg;		//flag to indicate whether capping strain has been passed or not

	//Damage History Variables
	//trial version:
	double excurEnerg, totalEnerg, beta;
	//commit version:
	double excurEnergP, totalEnergP, betaP;

	bool givenParams;

	//behavior parameters
	double cmpUnldngEFac, sigTEFac, sigTFfac, epsTFfac;
};
#endif
