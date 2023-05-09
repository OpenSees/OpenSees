#ifndef SLModel_h
#define SLModel_h

#include <UniaxialMaterial.h>

class SLModel : public UniaxialMaterial
{
public:
    SLModel(int tag, double Dt, double E, double sgm_ini, double c, double gamma, double q, double beta, double sigmaC, double epsiC, double Ed1, double Ed2, double sigmaDM,
	double aSigma, double aE, double lambda1Degrad, double cDegrad);      
    SLModel();    
    ~SLModel();


	const char *getClassType(void) const {return "SLModel";};

    int setTrialStrain(double strain, double strainRate = 0.0); 
    double getStrain(void);          
    double getStress(void);
    double getTangent(void);

    double getInitialTangent(void) {return E;};

    int commitState(void);
    int revertToLastCommit(void);    
    int revertToStart(void);    

    UniaxialMaterial *getCopy(void);
    
    int sendSelf(int commitTag, Channel &theChannel);  
    int recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker);    
    
    void Print(OPS_Stream &s, int flag =0);

	// functions
	void StrainHardeningFunc(void);
	void YieldPointFunc(void);
	void BackBoneCompFunc(void);
	void BackBoneTenFunc(void);
	void BackBoneComp2Func(void);
	void BackBoneTen2Func(void);
	
	

protected:
    
private:
	
	
	////////////////////////////////////////////////////////
	double Dt;
	double E, sgm_ini, c, gamma, q, beta;
	double sigmaC, epsiC, Ed1, Ed2, sigmaDM, aSigma, aE, lambda1Degrad, cDegrad;
	double sigmaCDivSigmaY, epsiCDivEpsiY, Ed1DivE, Ed2DivE, sigmaDMDivSigmaC;
	double Dteq;
	/*double Dt, sgm_ini;
	double E, Dteq;
	double c, gamma, q, beta;*/
	double CapYieldStressM, CapYieldStrainM, Ed1EM, Ed2EM, DetCapStressM;
	
	double status;
	
	double p_teps, p_neps, p_teps_prev, p_neps_prev;
	double cum_p_teps;
	double sgm_0, alf_d, alf;
	double ytsgm_p, ytsgm_n, yteps_p, yteps_n;
	
	double teps, neps, tsgm, nsgm;
	double teps_prev, neps_prev, tsgm_prev, nsgm_prev;
	
	double cEu, cSgmy, cEpsy, cSgmc, cEpsc, cSgmd1, cEpsd1, cSgmd2, cEpsd2, cSgmb, cSgmd, cEsth, cEd1, cEd2;
	double cIniSgmy, cIniEpsy, cIniSgmc, cIniEpsc, cIniEsth, cIniEd1, cIniEd2, cIniSgmd1, cIniEpsd1, cIniSgmb, cIniSgmd, cIniSgmd2, cIniEpsd2;
	double tEu, tSgmy, tEpsy, tSgmp, tEpsp, tEpsp2, tEr, tEr2, refEps;
	
	double ay, au;
	double Lambda1, c1, Lambda2, c2, Lambda3, c3;
	double Et1, Et2, Et3;
	double Beta1, Beta2, Beta3;
	double Alpha1, Alpha2, Alpha3;
	double TotalE, DeltaE;
	
	double Tangent;
	double iInitial;
	
	
	
	////////////////////////////////////////////////////////
	double C_Dt, C_sgm_ini;
	double C_OP_Material;
	double C_E, C_Dteq;
	double C_c, C_gamma, C_q, C_beta;
	double C_CapYieldStressM, C_CapYieldStrainM, C_Ed1EM, C_Ed2EM, C_DetCapStressM;
	
	double C_status;
	
	double C_p_teps, C_p_neps, C_p_teps_prev, C_p_neps_prev;
	double C_cum_p_teps;
	double C_sgm_0, C_alf_d, C_alf;
	double C_ytsgm_p, C_ytsgm_n, C_yteps_p, C_yteps_n;
	
	double C_teps, C_neps, C_tsgm, C_nsgm;
	double C_teps_prev, C_neps_prev, C_tsgm_prev, C_nsgm_prev;
	
	double C_cEu, C_cSgmy, C_cEpsy, C_cSgmc, C_cEpsc, C_cSgmd1, C_cEpsd1, C_cSgmd2, C_cEpsd2, C_cSgmb, C_cSgmd, C_cEsth, C_cEd1, C_cEd2;
	double C_cIniSgmy, C_cIniEpsy, C_cIniSgmc, C_cIniEpsc, C_cIniEsth, C_cIniEd1, C_cIniEd2, C_cIniSgmd1, C_cIniEpsd1, C_cIniSgmb, C_cIniSgmd, C_cIniSgmd2, C_cIniEpsd2;
	double C_tEu, C_tSgmy, C_tEpsy, C_tSgmp, C_tEpsp, C_tEpsp2, C_tEr, C_tEr2, C_refEps;
	
	double C_ay, C_au;
	double C_Lambda1, C_c1, C_Lambda2, C_c2, C_Lambda3, C_c3;
	double C_Et1, C_Et2, C_Et3;
	double C_Beta1, C_Beta2, C_Beta3;
	double C_Alpha1, C_Alpha2, C_Alpha3;
	double C_TotalE, C_DeltaE;
	
	double C_Tangent;
	double C_iInitial;
	
};


#endif









































































