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

#ifndef BoundingCamClay_h
#define BoundingCamClay_h

// Written: Kathryn Petek	
//          December 2004
// Modified: Chris McGann
//           January 2011

// Description: This file contains the class definition for BoundingCamClay. 


#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include <NDMaterial.h>
#include <Matrix.h>
#include <Vector.h>

class BoundingCamClay : public NDMaterial
{
  public:

    // full constructor
    BoundingCamClay(int tag, int classTag, double massDen, double C, double bulk, double OCR,  
							 double mu_o, double Alpha, double lambda, double h, double m);
    // null constructor
    BoundingCamClay();
    
    // destructor
    ~BoundingCamClay();
 
    NDMaterial *getCopy(const char *type);

    int commitState(void);
    int revertToLastCommit(void);
    int revertToStart(void);

    NDMaterial *getCopy(void);
    const char *getType(void) const;
    int getOrder(void) const;

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

	// input material parameters
	double iC;			    // ellipsoildal axis ratio
	double mBulk;           // initial bulk modulus
	double iOCR;            // overconsolidation ratio
	double ikappa;			// elastic compressibility index
	double imu_o;			// elastic shear modulus
	double ialpha;			// pressure-dependent parameter
	double ilambda;			// compressibility soil index for virgin loading
	double ih;				// hardening parameter
	double im;				// hardening parameter
	double iepsE_vo;		// initial volumetric strain invariant
	double massDen;         // mass density for dynamic analysis

	// internal variables
	double mr;              // hardening response variable for load funct, step n+1
	double mr_n;            // hardening response variable for load funct, step n
	double mR;              // hardening response variable for bounding funct, step n+1
	double mR_n;            // hardening response variable for bounding funct, step n
	double mKappa;          // loading/bounding function relation, step n+1
	double mKappa_n;        // loading/bounding function relation, step n
	double mp_o;			// initial mean stress
	double mTHETA;          // compressiblity constant, theta = 1/(ilam - ikap)
	double mStressRatio;    // computational variable

	Vector mEpsilon;        // strain tensor
	Vector mEpsilon_P;      // plastic strain tensor, step n+1
	Vector mEpsilon_n_P;    // plastic strain tensor, step n
	Vector mSigma;          // stress tensor, step n+1
	Vector mSigma_n;        // stress tensor, step n
	Vector mSIGMAo;         // normalized stress tensor, step n+1
	Vector mSIGMAo_n;       // normalized stress tensor, step n

	bool flagReversal;      // flag indicating switch to unloading response
	static double mElastFlag;
	bool initializeState;

	Matrix mCe;				// elastic tangent stiffness matrix
	Matrix mCep;			// elastoplastic tangent stiffness matrix
	Vector mI1;				// 2nd Order Identity Tensor
	Matrix mIIco;			// 4th-order identity tensor, covariant
	Matrix mIIcon;			// 4th-order identity tensor, contravariant
	Matrix mIImix;          // 4th-order identity tensor, mixed variant
	Matrix mIIvol;			// 4th-order volumetric tensor, IIvol = I1 tensor I1 
	Matrix mIIdevCon;       // 4th order deviatoric tensor, contravariant
	Matrix mIIdevMix;       // 4th order deviatoric tensor, mixed variant
	Matrix mM;              // 4th Order yield function tensor, covariant
	Vector mState;          // vector of state param for recorders
	  

	// member functions
	void initialize();
	void plastic_integrator();
	Matrix GetElasticOperator(double p, double ev, double es, Vector n);
	Matrix GetCep(double kappa, double r, double R, double dgamma, double rho, double eta,
                        double nu, Vector SIGMAo, Vector xi, Vector df_dSigma, Matrix Dep);
	Matrix GetComplianceOperator(double p, double ev, double es, Vector n);
	double GetContraNorm(Vector v);
	double GetCovariantNorm(Vector v); 
	double GetTrace(Vector v);
	Vector GetState(); 
	Vector GetCenter();

	// tensor functions
	double DoubleDot2_2(Vector v1, Vector v2);       // doubledot product vector-vector
	Vector DoubleDot2_4(Vector v1, Matrix m1);       // doubledot product vector-matrix
	Vector DoubleDot4_2(Matrix m1, Vector v2);       // doubledot product matrix-vector
	Matrix DoubleDot4_4(Matrix m1, Matrix m2);       // doubledot product matrix-matrix
	Matrix Dyadic2_2(Vector v1, Vector v2);          // dyadic product vector-vector

	// constant computation parameters
	static const double one3 ;
	static const double two3 ;
	static const double root23 ;
};
#endif
