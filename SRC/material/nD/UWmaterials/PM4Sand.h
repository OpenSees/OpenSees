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

// Computational Geomechanics Group
// Author: Long Chen, Pedro Arduino
// University of Washington
// Date:      Nov 2016
// Last Modified: Sep 2018

// Description: This file contains the implementation for the PM4Sand class.
// PM4Sand(Version 3.1): A Sand Plasticity Model For Earthquake Engineering Applications
// by R.W.Boulanger and K.Ziotopoulou
// Oct 2017

#ifndef PM4Sand_h
#define PM4Sand_h

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

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

class PM4Sand : public NDMaterial
{
public:
	// full constructor
	PM4Sand(int tag, int classTag, double Dr, double G0, double hp0, double mDen, double P_atm = 101.3, double h0 = -1, double emax = 0.8,
		double emin = 0.5, double nb = 0.5, double nd = 0.1, double Ado = -1, double z_max = -1, double cz = 250,
		double ce = -1, double phi_cv = 33.0, double nu = 0.3, double Cgd = 2.0, double Cdr = -1, double Ckaf = -1, double Q = 10,
		double R = 1.5, double m = 0.01, double Fsed_min = -1, double p_sdeo = -1, int integrationScheme = 1, int tangentType = 0,
		double TolF = 1.0e-7, double TolR = 1.0e-7);
	// full constructor
	PM4Sand(int tag, double Dr, double G0, double hp0, double mDen, double P_atm = 101.3, double h0 = -1, double emax = 0.8,
		double emin = 0.5, double nb = 0.5, double nd = 0.1, double Ado = -1, double z_max = -1, double cz = 250,
		double ce = -1, double phi_cv = 33.0, double nu = 0.3, double Cgd = 2.0, double Cdr = -1, double Ckaf = -1, double Q = 10,
		double R = 1.5, double m = 0.01, double Fsed_min = -1, double p_sdeo = -1, int integrationScheme = 1, int tangentType = 0,
		double TolF = 1.0e-7, double TolR = 1.0e-7);
	// null constructor
	PM4Sand();
	// destructor
	~PM4Sand();

	// send mass density to element in dynamic analysis
	double getRho(void) { return massDen; };
	double getVoidRatio(void) { return mVoidRatio; };
	int    getNumIterations(void) { return mIter; };

	int setTrialStrain(const Vector &v);
	int setTrialStrain(const Vector &v, const Vector &r);
	int initialize(Vector initStress);
	int initialize();
	NDMaterial *getCopy(const char *type);

	int commitState(void);
	int revertToLastCommit(void);
	int revertToStart(void);

	NDMaterial *getCopy(void);
	const char *getType(void) const;
	int        getOrder(void) const;

	// Recorder functions
	virtual const Vector& getStressToRecord() { return mSigma; };
	double getDGamma();
	const Vector getState();
	const Vector getAlpha();
	const Vector getFabric();
	const Vector getAlpha_in();
	const Vector getTracker();
	double getG();
	double getKp();
	const Vector getAlpha_in_p();
	const Matrix &getTangent();
	const Matrix &getInitialTangent();
	const Vector &getStress();
	const Vector &getStrain();
	const Vector &getElasticStrain();

	Response *setResponse(const char **argv, int argc, OPS_Stream &output);
	int getResponse(int responseID, Information &matInformation);

	int sendSelf(int commitTag, Channel &theChannel);
	int recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker);

	void Print(OPS_Stream &s, int flag = 0);

	int setParameter(const char **argv, int argc, Parameter &param);
	int updateParameter(int responseID, Information &info);


protected:

	// Material constants
	double m_Dr;
	double m_G0;
	double m_hpo;
	double massDen;     // mass density for dynamic analysis
	double m_P_atm;
	double m_h0;
	double m_emax;
	double m_emin;
	double m_e_init;
	double m_nb;
	double m_nd;
	double m_Ado;
	double m_cz;
	double m_ce;
	double m_Mc;
	double m_nu;
	double m_Cgd;
	double m_Cdr;
	double m_Ckaf;
	double m_Q;
	double m_R;
	double m_m;
	double m_z_max;
	double m_Fsed_min;
	double m_p_sedo;
	int m_FirstCall;
	int m_PostShake;

	// internal variables
	Vector mEpsilon;    // strain tensor
	Vector mEpsilon_n;  // strain tensor (last committed)
	Vector mEpsilon_r;  // negative strain tensor for returning
	Vector mSigma;      // stress tensor
	Vector mSigma_n;    // stress tensor (last committed)
	Vector mSigma_r;    // negative stress tensor for returning
	Vector mSigma_b;    // stress tensor offset from initial stress state outside bounding surface correction
	Vector mEpsilonE;	// elastic strain tensor
	Vector mEpsilonE_n;	// elastic strain tensor (last committed)
	Vector mEpsilonE_r; // negative elastic strain tensor for returning
	Vector mAlpha;		// back-stress ratio
	Vector mAlpha_n;	// back-stress ratio (last committed)
	Vector mAlpha_in;	// back-stress ratio at loading reversal
	Vector mAlpha_in_n;	// back-stress ratio at loading reversal (last committed)
	Vector mAlpha_in_p; // previous back-stress ratio at loading reversal
	Vector mAlpha_in_p_n; // previous back-stress ratio at loading reversal (last committed)
	Vector mAlpha_in_true;  // true initial back stress ratio tensor
	Vector mAlpha_in_true_n;  // true initial back stress ratio tensor (last committed)
	Vector mAlpha_in_max; // Maximum value of initial back stress ratio
	Vector mAlpha_in_max_n; // Maximum value of initial back stress ratio (last committed)
	Vector mAlpha_in_min; // Minimum value of initial back stress ratio
	Vector mAlpha_in_min_n; // Minimum value of initial back stress ratio (last committed)
	double mDGamma;		// plastic multiplier
	double mDGamma_n;	// plastic multiplier (last committed)
	Vector mFabric;		// fabric tensor
	Vector mFabric_n;	// fabric tensor (last committed)
	Vector mFabric_in;  // fabric tensor at loading reversal
	Vector mFabric_in_n;  // fabric tensor at loading reversal (last committed)
	Matrix mCe;			// elastic tangent
	Matrix mCep;		// continuum elastoplastic tangent
	Matrix mCep_Consistent; // consistent elastoplastic tangent
	double mK;			// state dependent Bulk modulus
	double mG;			// state dependent Shear modulus
	double mVoidRatio;	// material void ratio
	double mKp;         // plastic mudulus
	double mzcum;       // current cumulated fabric
	double mzpeak;      // current peak fabric
	double mpzp;
	double mzxp;        // product of z and p
	double mMb;
	double mMd;
	double mMcur;       // current stress ratio
	Vector mTracker;      // internal paramter tracker

	double	mTolF;			// max drift from yield surface
	double	mTolR;			// tolerance for Newton iterations
	char unsigned mIter;	// number of iterations
	char unsigned mScheme;	// 1: Forward Euler Explicit, 2: Modified Euler Explicit
	char unsigned mTangType;// 0: Elastic Tangent, 1: Contiuum ElastoPlastic Tangent, 2: Consistent ElastoPlastic Tangent
	double	m_Pmin;			// Minimum allowable mean effective stress
	double  m_Pmin2;        // Minimum p for Cpzp2 and Cpmin
	bool    m_pzpFlag;          // flag for updating pzp
	static char unsigned   me2p;	// 0: enforce elastic response

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
	public:
		initTensors() {
			// 2nd order identity tensor
			mI1.Zero();
			mI1(0) = 1.0;
			mI1(1) = 1.0;
			// 4th order mixed variant identity tensor
			mIImix.Zero();
			for (int i = 0; i < 3; i++) {
				mIImix(i, i) = 1.0;
			}
			// 4th order covariant identity tensor
			mIIco = mIImix;
			mIIco(2, 2) = 2.0;
			// 4th order contravariant identity tensor
			mIIcon = mIImix;
			mIIcon(2, 2) = 0.5;
			// 4th order volumetric tensor, IIvol = I1 tensor I1
			mIIvol.Zero();
			for (int i = 0; i < 2; i++) {
				mIIvol(i, 0) = 1.0;
				mIIvol(i, 1) = 1.0;
			}
			// 4th order contravariant deviatoric tensor
			mIIdevCon = mIIcon - 0.5*mIIvol;
			// 4th order covariant deviatoric tensor
			mIIdevCo = mIIco - 0.5*mIIvol;
			// 4th order mixed variant deviatoric tensor
			mIIdevMix = mIImix - 0.5*mIIvol;
		}
	} initTensorOps;

	// constant computation parameters
	static const double		one3;
	static const double		two3;
	static const double		root23;
	static const double     root12;
	static const double		small;
	static const bool		debugFlag;
	static const double		maxStrainInc;
	static const char unsigned	mMaxSubStep; // Max number of substepping

											 //Member Functions specific for PM4Sand model
											 //void	initialize();
	void	integrate();
	void	elastic_integrator(const Vector& CurStress, const Vector& CurStrain, const Vector& CurElasticStrain,
		const Vector& NextStrain, Vector& NextElasticStrain, Vector& NextStress, Vector& NextAlpha,
		double& NextVoidRatio, double& G, double& K, Matrix& aC, Matrix& aCep, Matrix& aCep_Consistent);
	void	explicit_integrator(const Vector& CurStress, const Vector& CurStrain, const Vector& CurElasticStrain,
		const Vector& CurAlpha, const Vector& CurFabric, const Vector& alpha_in, const Vector& alpha_in_p, const Vector& NextStrain,
		Vector& NextElasticStrain, Vector& NextStress, Vector& NextAlpha, Vector& NextFabric,
		double& NextDGamma, double& NextVoidRatio, double& G, double& K, Matrix& aC, Matrix& aCep, Matrix& aCep_Consistent);
	void	ForwardEuler(const Vector& CurStress, const Vector& CurStrain, const Vector& CurElasticStrain,
		const Vector& CurAlpha, const Vector& CurFabric, const Vector& alpha_in, const Vector& alpha_in_p, const Vector& NextStrain,
		Vector& NextElasticStrain, Vector& NextStress, Vector& NextAlpha, Vector& NextFabric,
		double& NextDGamma, double& NextVoidRatio, double& G, double& K, Matrix& aC, Matrix& aCep, Matrix& aCep_Consistent);
	void	ModifiedEuler(const Vector& CurStress, const Vector& CurStrain, const Vector& CurElasticStrain,
		const Vector& CurAlpha, const Vector& CurFabric, const Vector& alpha_in, const Vector& alpha_in_p, const Vector& NextStrain,
		Vector& NextElasticStrain, Vector& NextStress, Vector& NextAlpha, Vector& NextFabric,
		double& NextDGamma, double& NextVoidRatio, double& G, double& K, Matrix& aC, Matrix& aCep, Matrix& aCep_Consistent);
	void	RungeKutta4(const Vector& CurStress, const Vector& CurStrain, const Vector& CurElasticStrain,
		const Vector& CurAlpha, const Vector& CurFabric, const Vector& alpha_in, const Vector& alpha_in_p, const Vector& NextStrain,
		Vector& NextElasticStrain, Vector& NextStress, Vector& NextAlpha, Vector& NextFabric,
		double& NextDGamma, double& NextVoidRatio, double& G, double& K, Matrix& aC, Matrix& aCep, Matrix& aCep_Consistent);
	void	MaxStrainInc(const Vector& CurStress, const Vector& CurStrain, const Vector& CurElasticStrain,
		const Vector& CurAlpha, const Vector& CurFabric, const Vector& alpha_in, const Vector& alpha_in_p, const Vector& NextStrain,
		Vector& NextElasticStrain, Vector& NextStress, Vector& NextAlpha, Vector& NextFabric,
		double& NextDGamma, double& NextVoidRatio, double& G, double& K, Matrix& aC, Matrix& aCep, Matrix& aCep_Consistent);

	double	IntersectionFactor(const Vector& CurStress, const Vector& CurStrain, const Vector& NextStrain, const Vector& CurAlpha,
		double a0, double a1);
	double	IntersectionFactor_Unloading(const Vector& CurStress, const Vector& CurStrain, const Vector& NextStrain, const Vector& CurAlpha);
	void Stress_Correction(Vector& NextStress, Vector& NextAlpha, const Vector& alpha_in, const Vector& alpha_in_p, const Vector& CurFabric, double& NextVoidRatio);
	void Stress_Correction(Vector& NextStress, Vector& NextAlpha, const Vector& dAlpha, const double m, const Vector& R, const Vector& n, const Vector& r);
	// Material Specific Methods
	double	Macauley(double x);
	double	MacauleyIndex(double x);
	double	GetF(const Vector& nStress, const Vector& nAlpha);
	double	GetKsi(const double& e, const double& p);
	void	GetElasticModuli(const Vector& sigma, double &K, double &G);
	void	GetElasticModuli(const Vector& sigma, double &K, double &G, double &Mcur, const double& zcum);
	Matrix	GetStiffness(const double& K, const double& G);
	Matrix	GetCompliance(const double& K, const double& G);
	void	GetStateDependent(const Vector &stress, const Vector &alpha, const Vector &alpha_in, const Vector& alpha_in_p
		, const Vector &fabric, const Vector &fabric_in, const double &G, const double &zcum, const double &zpeak
		, const double &pzp, const double &Mcur, const double &dr, Vector &n, double &D, Vector &R, double &K_p
		, Vector &alphaD, double &Cka, double &h, Vector &b, double &AlphaAlphaBDotN);
	Matrix	GetElastoPlasticTangent(const Vector& NextStress, const Matrix& aCe, const Vector& R, const Vector& n, const double K_p);
	Vector	GetNormalToYield(const Vector &stress, const Vector &alpha);
	int	Check(const Vector& TrialStress, const Vector& stress, const Vector& CurAlpha, const Vector& NextAlpha);

	// Symmetric Tensor Operations
	double GetTrace(const Vector& v);
	Vector GetDevPart(const Vector& aV);
	double DoubleDot2_2_Contr(const Vector& v1, const Vector& v2);
	double DoubleDot2_2_Cov(const Vector& v1, const Vector& v2);
	double DoubleDot2_2_Mixed(const Vector& v1, const Vector& v2);
	double GetNorm_Contr(const Vector& v);
	double GetNorm_Cov(const Vector& v);
	Matrix Dyadic2_2(const Vector& v1, const Vector& v2);
	Vector DoubleDot4_2(const Matrix& m1, const Vector& v1);
	Vector DoubleDot2_4(const Vector& v1, const Matrix& m1);
	Matrix DoubleDot4_4(const Matrix& m1, const Matrix& m2);
	Vector ToContraviant(const Vector& v1);
	Vector ToCovariant(const Vector& v1);
};
#endif
