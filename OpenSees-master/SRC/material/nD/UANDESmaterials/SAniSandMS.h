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
                                                                       
// Written: José Abell, Haoyuan Liu & Federico Pisano
//
// Description: This file contains the implementation for the SAniSandMS class.
//
//  Alborz Ghofrani, Pedro Arduino
//			May 2013, University of Washington
//
// On revision 6602
//

#ifndef SAniSandMS_h
#define SAniSandMS_h

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <iostream>
#include <fstream>

#include <NDMaterial.h>
#include <Matrix.h>
#include <Vector.h>

#include <Information.h>
//#include <MaterialResponse.h>
#include <Parameter.h>

#include <Channel.h>
#include <FEM_ObjectBroker.h>
#include <string.h>


#include <elementAPI.h>

class SAniSandMS : public NDMaterial
{
  public:
	  
    // full constructor
  //   SAniSandMS(int tag, int classTag, double G0, double nu, double e_init, double Mc, 
  //   	double c, double lambda_c, double e0, double ksi,
  //   	double P_atm, double m, double h0, double ch, double nb, 
  //   	double A0, double nd, double zeta, double mu0, double beta, 
  //   	double w1, double w2, double mDen, 
		// int integrationScheme = 2, int tangentType = 2, int JacoType = 1, double TolF = 1.0e-7, double TolR = 1.0e-7);
   
    SAniSandMS(int tag, int classTag, double G0, double nu, double e_init, double Mc, 
    	double c, double lambda_c, double e0, double ksi,
    	double P_atm, double m, double h0, double ch, double nb, 
    	double A0, double nd, double zeta, double mu0, 
    		double beta,
			double mDen, 
		int integrationScheme = 2, int tangentType = 2, int JacoType = 1, double TolF = 1.0e-7, double TolR = 1.0e-7);

    // full constructor
    // SAniSandMS(int tag, double G0, double nu, double e_init, double Mc, double c, double lambda_c, double e0, double ksi,
				// 	double P_atm, double m, double h0, double ch, double nb, double A0, double nd, double zeta, double mu0, 
				// 	double beta, double w1, double w2, double mDen, 
				// 	int integrationScheme = 2, int tangentType = 2, int JacoType = 1, double TolF = 1.0e-7, double TolR = 1.0e-7);
    SAniSandMS(int tag, double G0, double nu, double e_init, double Mc, double c, double lambda_c, double e0, double ksi,
				double P_atm, double m, double h0, double ch, double nb, double A0, double nd, double zeta, double mu0, 
					double beta,
			 	 	double mDen, 
				int integrationScheme = 2, int tangentType = 2, int JacoType = 1, double TolF = 1.0e-7, double TolR = 1.0e-7);
    


    // null constructor
    SAniSandMS();
    // destructor
    ~SAniSandMS();
 
    NDMaterial *getCopy(const char *type);

    int commitState(void);
    int revertToLastCommit(void);
    int revertToStart(void);

    NDMaterial *getCopy(void);
    const char *getType(void) const;
    int        getOrder(void) const;

	// Recorder functions
	virtual const Vector& getStressToRecord() {return mSigma;};
	const Vector getState();
	const Vector getAlpha();
	const Vector getalpha_in();
	const Vector getAlphaM();
	double getMM();
	virtual const Vector& getEStrain();



    Response *setResponse (const char **argv, int argc, OPS_Stream &output);
    int getResponse (int responseID, Information &matInformation);

    int sendSelf(int commitTag, Channel &theChannel);  
    int recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker); 

    void Print(OPS_Stream &s, int flag =0);

	int setParameter(const char **argv, int argc, Parameter &param);
	int updateParameter(int responseID, Information &info);

	// send mass density to element in dynamic analysis
	double getRho(void) {return massDen;};
	double getVoidRatio(void) {return mVoidRatio;};
	int    getNumIterations(void) {return mIter;};

  protected:

	// Material constants
	double m_G0;
	double m_nu;
	double m_e_init;
	double m_Mc;
	double m_c;
	double m_lambda_c;
	double m_e0;
	double m_ksi;
	double m_P_atm;
	double m_m;
	double m_h0;
	double m_ch;
	double m_nb;
	double m_A0;
	double m_nd;


	double m_zeta;
	double m_mu0;
	// double m_beta;  old formulation
	double m_beta;   	//new formulation

	
	// internal variables
	Vector mEpsilon;    // strain tensor
	Vector mEpsilon_n;  // strain tensor (last committed)
	Vector mSigma;      // stress tensor
	Vector mSigma_n;    // stress tensor (last committed)
	Vector mEpsilonE;	// elastic strain tensor
	Vector mEpsilonE_n;	// elastic strain tensor (last committed)
	Vector mAlpha;		// back-stress ratio
	Vector mAlpha_n;	// back-stress ratio (last committed)
	Vector mAlphaM;		// Change!  back-stress ratio of the Memory Surface
	Vector mAlphaM_n;	// Change!  back-stress ratio of the Memory Surface (last committed)
	double mMM_plus;         // new ! SAniSandMS 2.0
	double mMM_plus_n;       // new ! SAniSandMS 2.0
	double mMM_minus;         // new ! SAniSandMS 2.0
	double mMM_minus_n;       // new ! SAniSandMS 2.0
	

	// Change!  Vector mAlpha_in;	// back-stress ratio at loading reversal
	Vector malpha_in;	// back-stress ratio at loading reversal
	//  Change! Vector mAlpha_in_n;	// back-stress ratio at loading reversal (last committed)
	Vector malpha_in_n;	// back-stress ratio at loading reversal (last committed)
	double mDGamma;		// plastic multiplier
	double mDGamma_n;	// plastic multiplier (last committed)
	// Vector mFabric;		// fabric tensor
	// Vector mFabric_n;	// fabric tensor (last committed)
	Matrix mCe;			// elastic tangent
	Matrix mCep;		// continuum elastoplastic tangent
	Matrix mCep_Consistent; // consistent elastoplastic tangent
	double mK;			// state dependent Bulk modulus
	double mG;			// state dependent Shear modulus
	double massDen;     // mass density for dynamic analysis
	double mVoidRatio;	// material void ratio

	int num_load_reversals;
	int num_load_reversals_n;

	bool m_firstLoading;

	double	mTolF;			// max drift from yield surface
	double	mTolR;			// tolerance for Newton iterations
	char unsigned mIter;	// number of iterations
	char unsigned mJacoType;// 0: Finite Difference Jacobian, 1: Analytical Jacobian
	char unsigned mScheme;	// 0: Forward Euler Explicit, 1: Backward Euler Implicit, 2: Backward Euler Implicit with considerations for stability, 
														// 3: FE Explicit with constrained strain increment
	char unsigned mTangType;// 0: Elastic Tangent, 1: Contiuum ElastoPlastic Tangent, 2: Consistent ElastoPlastic Tangent
	bool    mUseElasticTan;
	double	mEPS;			// machine epsilon (for FD jacobian)
	double	m_Pmin;			// Minimum allowable mean effective stress
	static char unsigned mElastFlag;	// 1: enforce elastic response

	static Vector mI1;			// 2nd Order Identity Tensor
	static Matrix mIIco;		// 4th-order identity tensor, covariant
	static Matrix mIIcon;		// 4th-order identity tensor, contravariant
	static Matrix mIImix;		// 4th-order identity tensor, mixed variant
	static Matrix mIIvol;		// 4th-order volumetric tensor, IIvol = I1 tensor I1 
	static Matrix mIIdevCon;	// 4th order deviatoric tensor, contravariant
	static Matrix mIIdevMix;	// 4th order deviatoric tensor, mixed variant
	static Matrix mIIdevCo;		// 4th order deviatoric tensor, covariant
	// initialize these Vector and Matrices:
	static class initTensors {
	public :
		initTensors() {
			// 2nd order identity tensor
			mI1.Zero();
			mI1(0) = 1.0;
			mI1(1) = 1.0;
			mI1(2) = 1.0;
			// 4th order mixed variant identity tensor
			mIImix.Zero();
			for (int i = 0; i<6; i++) {
				mIImix(i,i) = 1.0;
			}
			// 4th order covariant identity tensor
			mIIco = mIImix;
			mIIco(3,3) = 2.0;
			mIIco(4,4) = 2.0;
			mIIco(5,5) = 2.0;
			// 4th order contravariant identity tensor
			mIIcon = mIImix;
			mIIcon(3,3) = 0.5;
			mIIcon(4,4) = 0.5;
			mIIcon(5,5) = 0.5;
			// 4th order volumetric tensor, IIvol = I1 tensor I1
			mIIvol.Zero();
			for (int i = 0; i<3; i++) {
				mIIvol(i,0) = 1.0;
				mIIvol(i,1) = 1.0;
				mIIvol(i,2) = 1.0;
			}
			// 4th order contravariant deviatoric tensor
			mIIdevCon = mIIcon - one3*mIIvol;
			// 4th order covariant deviatoric tensor
			mIIdevCo  = mIIco  - one3*mIIvol;
			// 4th order mixed variant deviatoric tensor
			mIIdevMix = mIImix - one3*mIIvol;
		}
	} initTensorOps;

	// constant computation parameters
	static const double		one3;
	static const double		two3;
	static const double		root23;
	static const double		small;
	static const bool		debugFlag;
	static const double		maxStrainInc;
	static const char unsigned	mMaxSubStep; // Max number of substepping

	//Member Functions specific for SAniSandMS model
	void	initialize();

	void	integrate();
	void	elastic_integrator(const Vector& CurStress, const Vector& CurStrain, const Vector& CurElasticStrain,
					const Vector& NextStrain, Vector& NextElasticStrain, Vector& NextStress, Vector& NextAlpha,
					double& NextVoidRatio, double& G, double& K, Matrix& aC, Matrix& aCep, Matrix& aCep_Consistent, Vector& NextAlphaM);
	void	explicit_integrator(const Vector& CurStress, const Vector& CurStrain, const Vector& CurElasticStrain,
        			const Vector& CurAlpha, const Vector& CurAlphaM, const double CurMM_plus, double CurMM_minus, const Vector& Curalpha_in, const Vector& NextStrain,
        			Vector& NextElasticStrain, Vector& NextStress, Vector& NextAlpha, Vector& NextAlphaM, double& NextMM_plus, double& NextMM_minus,
        			double& NextDGamma, double& NextVoidRatio,  double& G, double& K, Matrix& aC, Matrix& aCep, Matrix& aCep_Consistent);
	void	MaxEnergyInc(const Vector& CurStress, const Vector& CurStrain, const Vector& CurElasticStrain,
        			const Vector& CurAlpha, const Vector& CurAlphaM, const double CurMM_plus, double CurMM_minus, const Vector& Curalpha_in, const Vector& NextStrain,
        			Vector& NextElasticStrain, Vector& NextStress, Vector& NextAlpha, Vector& NextAlphaM, double& NextMM_plus, double& NextMM_minus,
        			double& NextDGamma, double& NextVoidRatio,  double& G, double& K, Matrix& aC, Matrix& aCep, Matrix& aCep_Consistent);
	void	MaxStrainInc(const Vector& CurStress, const Vector& CurStrain, const Vector& CurElasticStrain,
        			const Vector& CurAlpha, const Vector& CurAlphaM, const double CurMM_plus, double CurMM_minus, const Vector& Curalpha_in, const Vector& NextStrain,
        			Vector& NextElasticStrain, Vector& NextStress, Vector& NextAlpha, Vector& NextAlphaM, double& NextMM_plus, double& NextMM_minus,
        			double& NextDGamma, double& NextVoidRatio,  double& G, double& K, Matrix& aC, Matrix& aCep, Matrix& aCep_Consistent);
	void	ForwardEuler(const Vector& CurStress, const Vector& CurStrain, const Vector& CurElasticStrain,
        			const Vector& CurAlpha, const Vector& CurAlphaM, const double CurMM_plus, double CurMM_minus, const Vector& Curalpha_in, const Vector& NextStrain,
        			Vector& NextElasticStrain, Vector& NextStress, Vector& NextAlpha, Vector& NextAlphaM, double& NextMM_plus, double& NextMM_minus,
        			double& NextDGamma, double& NextVoidRatio,  double& G, double& K, Matrix& aC, Matrix& aCep, Matrix& aCep_Consistent);
	void	ModifiedEuler(const Vector& CurStress, const Vector& CurStrain, const Vector& CurElasticStrain,
        			const Vector& CurAlpha, const Vector& CurAlphaM, const double CurMM_plus, double CurMM_minus, const Vector& Curalpha_in, const Vector& NextStrain,
        			Vector& NextElasticStrain, Vector& NextStress, Vector& NextAlpha, Vector& NextAlphaM, double& NextMM_plus, double& NextMM_minus,
        			double& NextDGamma, double& NextVoidRatio,  double& G, double& K, Matrix& aC, Matrix& aCep, Matrix& aCep_Consistent);
	void	RungeKutta4(const Vector& CurStress, const Vector& CurStrain, const Vector& CurElasticStrain,
        			const Vector& CurAlpha, const Vector& CurAlphaM, const double CurMM_plus, double CurMM_minus, const Vector& Curalpha_in, const Vector& NextStrain,
        			Vector& NextElasticStrain, Vector& NextStress, Vector& NextAlpha, Vector& NextAlphaM, double& NextMM_plus, double& NextMM_minus,
        			double& NextDGamma, double& NextVoidRatio,  double& G, double& K, Matrix& aC, Matrix& aCep, Matrix& aCep_Consistent);
	int		BackwardEuler_CPPM(const Vector& CurStress, const Vector& CurStrain, const Vector& CurElasticStrain,
        			const Vector& CurAlpha, const Vector& CurAlphaM, const double CurMM_plus, double CurMM_minus, const Vector& Curalpha_in, const Vector& NextStrain,
        			Vector& NextElasticStrain, Vector& NextStress, Vector& NextAlpha, Vector& NextAlphaM, double& NextMM_plus, double& NextMM_minus,
        			double& NextDGamma, double& NextVoidRatio,  double& G, double& K, Matrix& aC, Matrix& aCep, Matrix& aCep_Consistent);

	double	IntersectionFactor(const Vector& CurStress, const Vector& CurStrain, const Vector& NextStrain, const Vector& CurAlpha, 
				double a0, double a1);
	double	IntersectionFactor_Unloading(const Vector& CurStress, const Vector& CurStrain, const Vector& NextStrain, const Vector& CurAlpha);
	void	Stress_Correction( const Vector& CurStress, const Vector& CurStrain, const Vector& CurElasticStrain, 
		            const Vector& CurAlpha, const Vector& CurAlphaM,  const double CurMM_plus, const double CurMM_minus, const Vector& Curalpha_in, const Vector& NextStrain,
		            Vector& NextElasticStrain, Vector& NextStress, Vector& NextAlpha, Vector& NextAlphaM, double& NextMM_plus, double& NextMM_minus,
		            double& NextDGamma, double& NextVoidRatio, double& G, double& K, Matrix& aC, Matrix& aCep, Matrix& aCep_Consistent);
	
	int		NewtonIter(const Vector& xo, const Vector& inVar, Vector& x, Matrix& aCepPart);
	int		NewtonIter2(const Vector& xo, const Vector& inVar, Vector& sol, Matrix& aCepPart);
	int		NewtonSol(const Vector& x, const Vector &inVar, Vector& del, Matrix& Cep);
	int		NewtonIter3(const Vector& xo, const Vector& inVar, Vector& sol, Matrix& aCepPart);
	int		NewtonSol2(const Vector& x, const Vector &inVar, Vector& res, Vector& JRes, Vector& del, Matrix& Cep);
	int		NewtonIter2_negP(const Vector& xo, const Vector& inVar, Vector& sol, Matrix& aCepPart);
	int		NewtonSol_negP(const Vector &xo, const Vector &inVar, Vector& del, Matrix& Cep);
	Vector  NewtonRes(const Vector &xo, const Vector &inVar);
	Vector  NewtonRes_negP(const Vector &xo, const Vector &inVar);
	Vector	GetResidual(const Vector& x, const Vector& inVar);
	Matrix	GetJacobian(const Vector &x, const Vector &inVar);
	Matrix	GetFDMJacobian(const Vector &delta, const Vector &inVar);
	Vector	SetManzariComponent(const Vector& stress, const Vector& alpha,
				const Vector& fabric, const double& dGamma);
	Vector	SetManzariStateInVar(const Vector& nStrain, const Vector& cStrain, const Vector& cStress, 
				const Vector& cEStrain, const Vector& cAlpha, const Vector& cFabric,
				const double& cVoidRatio, const double& nVoidRatio, const Vector& Alpha_in);
	double	machineEPS();
	// Material Specific Methods
	double	Macauley(double x);
	double	MacauleyIndex(double x);
	double	g(const double cos3theta, const double c);
	double	GetF(const Vector& nStress, const Vector& nAlpha);
	double	GetPSI(const double& e, const double& p);
	double	GetLodeAngle(const Vector& n);
	void	GetElasticModuli(const Vector& sigma, const double& en, const double& en1,
				const Vector& nEStrain, const Vector& cEStrain, double &K, 
				double &G);
	void	GetElasticModuli(const Vector& sigma, const double& en, double &K, double &G);
	void	GetElasticModuli(const Vector& sigma, const double& en, double &K, double &G, const double& D);
	Matrix	GetStiffness(const double& K, const double& G);
	Matrix	GetCompliance(const double& K, const double& G);
	void	GetStateDependent(const Vector &stress, const Vector &alpha, 
    const Vector &alphaM, const double MM_plus, const double MM_minus, const Vector &alpha_in,  
     const double &e, Vector &n, Vector &d, Vector &b, Vector &bM, 
     double &cos3Theta, double &h, double& hM, double &psi, 
     double &rBtheta, double &rDtheta, double &b0, 
     double& A, double& D, double& B, double& C, Vector& R, double &Z);

	Matrix	GetElastoPlasticTangent(const Vector& NextStress, const double& NextDGamma, const Vector& CurStrain, const Vector& NextStrain,
				const double& G, const double& K, const double& B, const double& C,const double& D, const double& h, 
				const Vector& n, const Vector& d, const Vector& b) ;
	Vector	GetNormalToYield(const Vector &stress, const Vector &alpha);
	int	Check(const Vector& TrialStress, const Vector& stress, const Vector& CurAlpha, const Vector& NextAlpha);

	// Symmetric Tensor Operations
	double GetTrace(const Vector& v);
	Vector GetDevPart(const Vector& aV);
	Vector SingleDot(const Vector& v1, const Vector& v2);
	double DoubleDot2_2_Contr(const Vector& v1, const Vector& v2);
	double DoubleDot2_2_Cov(const Vector& v1, const Vector& v2);
	double DoubleDot2_2_Mixed(const Vector& v1, const Vector& v2);
	double GetNorm_Contr(const Vector& v);
	double GetNorm_Cov(const Vector& v);
	Matrix Dyadic2_2(const Vector& v1, const Vector& v2);
	Vector DoubleDot4_2(const Matrix& m1, const Vector& v1);
	Vector DoubleDot2_4(const Vector& v1, const Matrix& m1);
	Matrix DoubleDot4_4(const Matrix& m1, const Matrix& m2);
	Matrix SingleDot4_2(const Matrix& m1, const Vector& v1);
	Matrix SingleDot2_4(const Vector& v1, const Matrix& m1);
	Matrix Trans_SingleDot4T_2(const Matrix& m1, const Vector& v1);
	double Det(const Vector& aV);
	Vector Inv(const Vector& aV);
	Vector ToContraviant(const Vector& v1);
	Vector ToCovariant(const Vector& v1);
	Matrix ToContraviant(const Matrix& m1);
	Matrix ToCovariant(const Matrix& m1);

};

#endif
