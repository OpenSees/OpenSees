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
                                                                       
// Created: Pedro Arduino, UW, 11.2011
//
// Description: This file contains the class definition for ManzariDafalias.

#ifndef ManzariDafalias_h
#define ManzariDafalias_h

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include <NDMaterial.h>
#include <Matrix.h>
#include <Vector.h>

class ManzariDafalias : public NDMaterial
{
  public:

    // full constructor
    ManzariDafalias(int tag, int classTag, double Ko, double Go, double v, double b, double Patm,
	                                       double Ao, double ho, double Cm, double Me, double Mc,
										   double kBE, double kBC, double kDE, double kDC, double ecRef,
										   double lambda, double Pref, double m, double Fmax, double Cf,
										   double eo, double mDen);
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
    int getOrder(void) const;

	Vector GetState();
    Response *setResponse (const char **argv, int argc, OPS_Stream &output);
    int getResponse (int responseID, Information &matInformation);

    int sendSelf(int commitTag, Channel &theChannel);  
    int recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker); 

    void Print(OPS_Stream &s, int flag =0);

	int setParameter(const char **argv, int argc, Parameter &param);
  	int updateParameter(int responseID, Information &eleInformation);

	// send mass density to element in dynamic analysis
	double getRho(void) {return massDen;};

  protected:

	// memeber variables and functions
    Vector mParam;
	// internal variables
    Vector mEpsilonE;
	Vector mAlpha;
	double mDGamma;
	Vector mFabric;
	Matrix mCe;
	Matrix mCep;
	double mM;
	double mK;  // state dependent Bulk modulus
	double mG;  // state dependent Shear modulus
	double massDen;         // mass density for dynamic analysis
	Vector mState;          // vector of state param for recorders

	Vector mEpsilon;        // strain tensor, step n+1
	Vector mEpsilon_n;      // strain tensor, step n
	Vector mSigma;          // stress tensor, step n+1
	Vector mSigma_n;        // stress tensor, step n

	double	mTolF;
	double	mTolR;
	double	mIter;
	double  mEPS;
	double  mElastFlag;
	bool initializeState;

	Vector mI1;				// 2nd Order Identity Tensor
	Matrix mIIco;			// 4th-order identity tensor, covariant
	Matrix mIIcon;			// 4th-order identity tensor, contravariant
	Matrix mIImix;          // 4th-order identity tensor, mixed variant
	Matrix mIIvol;			// 4th-order volumetric tensor, IIvol = I1 tensor I1 
	Matrix mIIdevCon;       // 4th order deviatoric tensor, contravariant
	Matrix mIIdevMix;       // 4th order deviatoric tensor, mixed variant

	// constant computation parameters
	static const double one3;
	static const double two3;
	static const double root23;
	static const double PI;

	//Member Functions specific for ManzariDafalias model
	void initialize();
	void plastic_integrator();
	
	Vector SetManzariComponent(const Vector& eStrain, const Vector& alpha, const double& m,
							   const Vector& fabric, const double& dGamma);
	Vector SetManzariStateInVar(const Vector& nStrain, const Vector& cStrain, const Vector& cStress, 
							    const Vector& cEStrain, const Vector& cAlpha, const double& cM,
							    const Vector& cFabric);
	Vector HypoElastic(const Vector& cStress, const Vector& nEStrain, 
					   const Vector& cEStrain, Matrix& C);
	double F(Vector& nextStress, Vector& nAlpha, double m);
	Vector NewtonSolve(const Vector& xo, const Vector& inVar, Matrix& aCepPart);
	Vector GetResidual(const Vector&  x, const Vector& inVar);

	void StateDepend( const Vector &stress, const Vector &alpha, const double &m,
					  const Vector &strain, Vector &n, Vector &d, Vector &b,
					  double &bref, double &psi, double &Theta, double &alphaBtheta);

	double PSI(const Vector& strain, const double& p);
	Matrix GetJacobian(const Vector &x, const Vector &inVar);
	Vector getdCosThetaOverdEE(const Vector& r, const double &p, const Vector& s, double K, double G, double theta);
	Vector getdCosThetaOverdAlpha(const Vector& rBar, double theta);
	Vector getdAlphaDOverdEE(double c, double theta, const Vector& dCosThetaOverdEE,
							 double psi, double cd, double K, double p);
	Vector getdAlphaDOverdAlpha(double c, double theta, const Vector& dCosThetaOverdAlpha,
							    double psi, double cd);
	Vector getdAlphaBOverdEE(double c, double theta, const Vector& dCosThetaOverdEE,
							 double psi, double cb, double K, double p);
	Vector getdAlphaBOverdAlpha(double c, double theta, const Vector dCosThetaOverdAlpha,
							    double psi, double cb, double K, double p);
	double getdgOverdCosTheta(const double theta, const double c);
	int sign(double x);

	// tensor operations
	double GetContraNorm(const Vector& v);
	double GetCovariantNorm(const Vector& v); 
	double GetTrace(const Vector& v);
	Vector GetDevPart(const Vector& v);
	double GetJ2(const Vector& v);
	double GetJ3(const Vector& v);
	double GetLodeAngle(const Vector& v);
	double Det(const Vector& v);
	double g(const double Theta, const double c);
	double Macauley(double x);
	int    MacauleyIndex(double x);
	double machineEPS();

	Vector Dot(const Vector& v1);
	double DoubleDot2_2(const Vector& v1, const Vector& v2);       // doubledot product vector-vector
	Vector DoubleDot2_4(const Vector& v1, const Matrix& m1);       // doubledot product vector-matrix
	Vector DoubleDot4_2(const Matrix& m1, const Vector& v2);       // doubledot product matrix-vector
	Matrix DoubleDot4_4(const Matrix& m1, const Matrix& m2);       // doubledot product matrix-matrix
	Matrix Dyadic2_2(const Vector& v1, const Vector& v2);          // dyadic product vector-vector

};
#endif
