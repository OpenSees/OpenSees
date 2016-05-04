// Written/implemented by: 
// J.Y. LU,
// Southeast University, Nanjing, China
//	& 
// G. LIN (guan.lin@polyu.edu.hk),	
// The Hong Kong Polytechnic University, Hong Kong, China

/* *********************************************************************************************
** Compression Part:
** Cyclic behavior based on Ref: 
**		L. Lam and J.G. Teng (2009), 
**		"Stress-strain model for FRP-confined concrete under cyclic axial compression",
**		Engineering Structures, 31, 308-321.										
** Envelop curve based on Ref:
**		J.G. Teng, T. Jiang, L. Lam, Y.Z. Luo (2009),
**		"Refinement of a design-oriented stress-strain model for FRP-confined concrete",
**		Journal of Composites for Construction, ASCE, 13(4), 269-278
** Column simulation based on Ref:
**		J.G. Teng, L. Lam, G. Lin, J.Y. Lu and Q.G. Xiao (2015),
**		Numerical Simulation of FRP-Jacketed RC Columns Subjected to Cyclic and Seismic Loading,
**		Journal of Composites for Construction, ASCE, 20(1), 04015021.
** Tension Part:
** Based on Ref:
**		M.H.M. Yassin (1994),
**		Nonlinear analysis of prestressed concrete structures under monotonic and cyclic loads. 
**		Dissertation.University of California. Berkeley, California.
**		and made some improvements
** *********************************************************************************************/

#ifndef FRPConfinedConcrete02_h
#define FRPConfinedConcrete02_h

#include <UniaxialMaterial.h>
#define MAT_TAG_FRPConfinedConcrete02 5010

class FRPConfinedConcrete02 : public UniaxialMaterial
{
public:
	FRPConfinedConcrete02(int tag, double fc0, double Ec, double ec0, double t, double Efrp, double eps_h_rup, double R, double ft, double Ets, int Unit);
	FRPConfinedConcrete02(int tag, double fc0, double Ec, double ec0, double fcc, double ecu, double ft, double Ets, int Unit);  
	FRPConfinedConcrete02(int tag, double fc0, double Ec, double ec0, double ft, double Ets, int Unit);
	FRPConfinedConcrete02();    

	~FRPConfinedConcrete02();

	int setTrialStrain(double strain, double strainRate = 0.0); 
	double getStrain(void);          
	double getStress(void);
	double getTangent(void);

	double getInitialTangent(void) {return m_Ec;};

	int commitState(void);
	int revertToLastCommit(void);    
	int revertToStart(void);    

	UniaxialMaterial *getCopy(void);

	int sendSelf(int commitTag, Channel &theChannel);  
	int recvSelf(int commitTag, Channel &theChannel, 
		FEM_ObjectBroker &theBroker);    

	void Print(OPS_Stream &s, int flag =0);

	//////////////////////////////////////////////////////////////////////////
	void Tens_Envlp(double epsc, double &sigc, double &Ect);
	void Compr_Envlp(double epsc, double &sigc, double &Ect);
	void Compr_UnloadingPath(double epsc, double &sigc, double &Ect);
	void Compr_GetPlasticStrain();
	void Compr_GetStrainRecoveryRatio();
	void Compr_ReloadingPath(double epsc, double &sigc, double &Ect);
	void GetRefPoint();
	void GetDeterioratedStress();
	void GetStressDeteriorationRatio();

	// AddingSensitivity:BEGIN //////////////////////////////////////////
	int    setParameter             (const char **argv, int argc, Information &info);
	int    updateParameter          (int parameterID, Information &info);
	int    activateParameter        (int parameterID);
	double getStressSensitivity     (int gradNumber, bool conditional);
	int    commitSensitivity        (double strainGradient, int gradNumber, int numGrads);
	// AddingSensitivity:END ///////////////////////////////////////////
protected:

private:
	//////////////////////////////////////////////////////////////////////////
	//inputs
	double m_fc0;
	double m_Ec;
	double m_epsc0;
	double m_t;
	double m_Efrp;
	double m_eps_h_rup;
	double m_R;
	double m_Ets;
	double m_ft;

	int m_Unit;
	//////////////////////////////////////////////////////////////////////////
	//Assistant variable
	double m_epst;
	double m_epscu;
	double m_fl;
	double m_fcc;
	double m_E2;
	double m_Unitscale;
	//////////////////////////////////////////////////////////////////////////
	//Current state variable
	//Tension
	double m_epstn;
	double m_epstu;
	double m_Etr1;
	double m_Etr2;

	//Compression
	int m_n;
	int m_ne;
	int m_loadingflag;

	double m_Ere;
	double m_epsunenv;
	double m_sigunenv;
	double m_sigunn1;
	double m_epsretenv;
	double m_epsun;
	double m_sigun;
	double m_epsre;
	double m_sigre;
	double m_epsref;
	double m_sigref;
	double m_Eun0;

	double m_betaun;
	double m_fi;
	double m_fiful;
	double m_signew;

	double m_gamare;
	double m_omg;
	double m_omgful;
	double m_epspl;

	bool m_bSmallStress;
	bool m_bUltimate;

	double m_Tstrain;
	double m_Tstress;

	double m_trialTangent;

	//////////////////////////////////////////////////////////////////////////
	//History variable
	//Tension
	double m_epstnlast;
	double m_epstulast;
	double m_Etr1last;
	double m_Etr2last;

	//Compression
	int m_nlast;
	int m_nelast;
	int m_loadingflaglast;

	double m_Erelast;
	double m_epsunenvlast;
	double m_sigunenvlast;
	double m_sigunn1last;
	double m_epsretenvlast;
	double m_epsunlast;
	double m_sigunlast;
	double m_epsrelast;
	double m_sigrelast;
	double m_epsreflast;
	double m_sigreflast;
	double m_Eun0last;

	double m_betaunlast;
	double m_filast;
	double m_fifullast;
	double m_signewlast;

	double m_gamarelast;
	double m_omglast;
	double m_omgfullast;
	double m_epspllast;

	bool m_bSmallStresslast;
	bool m_bUltimatelast;

	double m_trialStrainlast;
	double m_trialStresslast;

	double m_trialTangentlast;

	// AddingSensitivity:BEGIN //////////////////////////////////////////
	int parameterID;
	Matrix *SHVs;
	// AddingSensitivity:END ///////////////////////////////////////////
};
#endif

