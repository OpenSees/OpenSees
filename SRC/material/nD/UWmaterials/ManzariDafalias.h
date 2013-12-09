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
                                                                       
// Written: Alborz Ghofrani, Pedro Arduino
//			May 2013, University of Washington
                                                                      
// Description: This file contains the implementation for the ManzariDafalias class.

#ifndef ManzariDafalias_h
#define ManzariDafalias_h

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include <NDMaterial.h>
#include <Matrix.h>
#include <Vector.h>

#include <Information.h>
#include <MaterialResponse.h>
#include <Parameter.h>

#include <Channel.h>
#include <FEM_ObjectBroker.h>
#include <string.h>


#include <elementAPI.h>

class ManzariDafalias : public NDMaterial
{
  public:

    // full constructor
    ManzariDafalias(int tag, int classTag, double G0, double nu, double e_init, double Mc, double c, double lambda_c, double e0, double ksi,
	double P_atm, double m, double h0, double ch, double nb, double A0, double nd, double z_max, double cz, double mDen, int integrationScheme = 2,
	int tangentType = 2, int JacoType = 1, double TolF = 1.0e-7, double TolR = 1.0e-7);
    // null constructor
    ManzariDafalias();
    // destructor
    ~ManzariDafalias();
 
    NDMaterial *getCopy(const char *type);

    int commitState(void);
    int revertToLastCommit(void);
    int revertToStart(void);

    NDMaterial *getCopy(void);
    const char *getType(void) const;
    int        getOrder(void) const;

	virtual const Vector& getStressToRecord() {return mSigma;};
	const Vector getState();
	const Vector getAlpha();
	const Vector getFabric();
	const Vector getAlpha_in();
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
	double m_z_max;
	double m_cz;
	
	// internal variables
    Vector mEpsilonE;
	Vector mEpsilonE_n;
	Vector mAlpha;
	Vector mAlpha_n;
	Vector mAlpha_in;
	Vector mAlpha_in_n;
	double mDGamma_n;
	double mDGamma;
	Vector mFabric;
	Vector mFabric_n;
	Matrix mCe;
	Matrix mCep;
	Matrix mCep_Consistent;
	double mK;  // state dependent Bulk modulus
	double mG;  // state dependent Shear modulus
	double massDen;         // mass density for dynamic analysis
	double mVoidRatio;

	Vector mEpsilon;        // strain tensor, step n+1
	Vector mEpsilon_n;      // strain tensor, step n
	Vector mSigma;          // stress tensor, step n+1
	Vector mSigma_n;        // stress tensor, step n

	double	mTolF;
	double	mTolR;
	double  mEPS;
	int  	mIter;
	int     mJacoType;          // 0: FDM Jacobian, 1: Analytical Jacobian
	int		mScheme;            // 0: Forward Euler Explicit, 1: Backward Euler Implicit, 2: Backward Euler Implicit with considerations for stability
	int     mTangType;          // 0: Elastic Tangent, 1: Contiuum ElastoPlastic Tangent, 2: Consistent ElastoPlastic Tangent
	int     mOrgTangType;
	int     mLoadUnloadFlag;
	static int  mElastFlag;     // 1: enforce elastic response
	bool initializeState;

	Vector mI1;				// 2nd Order Identity Tensor
	Matrix mIIco;			// 4th-order identity tensor, covariant
	Matrix mIIcon;			// 4th-order identity tensor, contravariant
	Matrix mIImix;          // 4th-order identity tensor, mixed variant
	Matrix mIIvol;			// 4th-order volumetric tensor, IIvol = I1 tensor I1 
	Matrix mIIdevCon;       // 4th order deviatoric tensor, contravariant
	Matrix mIIdevMix;       // 4th order deviatoric tensor, mixed variant
	Matrix mIIdevCo;		// 4th order deviatoric tensor, covariant

	// constant computation parameters
	static const double one3;
	static const double two3;
	static const double root23;
	static const double small;
	static const bool   debugFlag;
	static const int    mMaxSubStep; // Max number of implicit substepping

	//Member Functions specific for ManzariDafalias model
	void initialize();
	void plastic_integrator();
	void elastic_integrator(const Vector& CurStress, const Vector& CurStrain, const Vector& CurElasticStrain,
		const Vector& NextStrain, Vector& NextElasticStrain, Vector& NextStress, Vector& NextAlpha,
		double& NextVoidRatio, double& G, double& K, Matrix& aC, Matrix& aCep, Matrix& aCep_Consistent);
	void explicit_integrator(const Vector& CurStress, const Vector& CurStrain, const Vector& CurElasticStrain,
		const Vector& CurAlpha, const Vector& CurFabric, const Vector& alpha_in, const Vector& NextStrain,
		Vector& NextElasticStrain, Vector& NextStress, Vector& NextAlpha, Vector& NextFabric, Vector& NextAlpha_in, 
		double& NextDGamma, double& NextVoidRatio, double& G, double& K, Matrix& aC, Matrix& aCep, Matrix& aCep_Consistent) ;
	int  implicit_integrator(const Vector& CurStress, const Vector& CurStrain, const Vector& CurElasticStrain,
		const Vector& CurAlpha, const Vector& CurFabric, const Vector& alpha_in, const Vector& NextStrain,
		Vector& NextElasticStrain, Vector& NextStress, Vector& NextAlpha, Vector& NextFabric, Vector& NextAlpha_in,
		double& NextDGamma,	double& NextVoidRatio, double& G, double& K, Matrix& aC, Matrix& aCep, Matrix& aCep_Consistent, int implicitLevel = 1) ;
	int    NewtonSolve(const Vector& xo, const Vector& inVar, Vector& x, Matrix& aCepPart);
	Vector GetResidual(const Vector& x, const Vector& inVar);
	Matrix GetJacobian(const Vector &x, const Vector &inVar);
	Matrix GetFDMJacobian(const Vector delta, const Vector inVar);
	Vector SetManzariComponent(const Vector& stress, const Vector& alpha,
							 const Vector& fabric, const double& dGamma);
	Vector SetManzariStateInVar(const Vector& nStrain, const Vector& cStrain, const Vector& cStress, 
									  const Vector& cEStrain, const Vector& cAlpha, const Vector& cFabric,
									  const double& cVoidRatio, const double& nVoidRatio, const Vector& Alpha_in);
	Vector NormalizeJacobian(Matrix& Jaco);
	void   DenormalizeJacobian(Matrix& JInv, const Vector& norms);
	double machineEPS();
	// Material Specific Methods
	double Macauley(double x);
	double MacauleyIndex(double x);
	double g(const double cos3theta, const double c);
	double GetF(const Vector& nStress, const Vector& nAlpha);
	double GetPSI(const double& e, const double& p);
	double GetLodeAngle(const Vector& n);
	void   GetElasticModuli(const Vector& sigma, const double& en, const double& en1,
								const Vector& nEStrain, const Vector& cEStrain, double &K, double &G);
	Matrix GetStiffness(const double& K, const double& G);
	Matrix GetCompliance(const double& K, const double& G);
	void   GetStateDependent( const Vector &stress, const Vector &alpha, const double &e
							, const Vector &alpha_in, Vector &n, Vector &d, Vector &b
							, double &cos3Theta, double &h, double &psi, double &alphaBtheta
							, double &alphaDtheta, double &b0);
	Matrix GetElastoPlasticTangent(const Vector& NextStress, const double& NextDGamma, const double& G
		                    , const double& K, const double& B, const double& C, const double& D
							, const double& h, const Vector& n, const Vector& d, const Vector& b, const Vector& InVar);
	Vector GetNormalToYield(const Vector &stress, const Vector &alpha);
	int    Check(const Vector& TrialStress, const Vector& stress, const Vector& CurAlpha, const Vector& NextAlpha);

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

};
#endif
