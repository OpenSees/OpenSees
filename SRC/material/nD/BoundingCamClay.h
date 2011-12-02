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
**                                                                    **
**                                                                    **
**                                                                    **
** ****************************************************************** */

#ifndef BoundingCamClay_h
#define BoundingCamClay_h

// $Revision: 1.1 $
// $Date: 2010-02-17 20:52:18 $
// $Source: /usr/local/cvs/OpenSees/SRC/material/nD/BoundingCamClay.h,v $

// Written: kap	
// Created: 12/04
//
// Description: This file contains the class definition for BoundingCamClay. 
//

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include <NDMaterial.h>
#include <Matrix.h>
#include <Vector.h>
#include <Tensor.h>

class BoundingCamClay : public NDMaterial
{
  public:
      // Full Constructor
	  BoundingCamClay(int tag, int classTag, double C, double r, double R, double p_o, double kappa,
							 double mu_o, double alpha, double lambda, double h, 
							 double m, double epsE_vo);

	  //Null Constructor
	  BoundingCamClay();
    
	  //Destructor
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
    int recvSelf(int commitTag, Channel &theChannel, 
		  FEM_ObjectBroker &theBroker); 

	void Print(OPS_Stream &s, int flag =0);

  protected:

	  //input material parameters
	  double iC;			// ellipsoildal axis ratio
	  double ip_o;			// preconsolidation pressure
	  double ikappa;		// hyperelastic compressibility parameter
	  double imu_o;			// elastic shear modulus
	  double ialpha;		// hyperelastic parameter
	  double ilambda;		// compressibility parameter
	  double ih;			// hardening parameter
	  double im;			// hardening parameter
	  double iepsE_vo;		// initial volumetric strain invariant

	  //internal variables
	  double mr;
	  double mr_n;
	  double mR;
	  double mR_n;
	  double mKappa;
	  double mKappa_n;
	  double mTHETA;

	  Vector mEpsilon;
	  //Vector mEpsilon_n;
	  Vector mEpsilon_P;
	  Vector mEpsilon_n_P;
	  Vector mSigma;
	  Vector mSigma_n;
	  Vector mSIGMAo;
	  Vector mSIGMAo_n;
	  //Vector mBeta;
	  //Vector mBeta_n;

	  bool flagReversal;

	  Matrix mM;
	  Matrix mphi;

	  Matrix mCe;			// elastic tangent stiffness matrix
	  Matrix mCep;			// elastoplastic tangent stiffness matrix
	  Vector mI1;			// 2nd Order Identity Tensor
	  Matrix mII;			// 4th Order Identity Tensor
	  Matrix mIIvol;		// IIvol = I1 tensor I1  
	    

	  //functions
	  void initialize();	// initializes variables

	  //plasticity integration routine
	  void plastic_integrator( ) ;

	  Matrix GetElasticModulus(double p, double q, double ev, double es, Vector n);

	  Matrix GetConsistentModulus(double kappa, double dgamma, double rho, double eta, double nu, double r, double R, double p, double ev, double es, Vector SIGMAo, Vector sigmaMinusAlpha, Vector dfdSigma, Vector n);

	  Matrix GetCompliantModulus(double p, double ev, double es, Vector n);
	  //
	  double Norm_EngStrain(Vector v); 
	  double GetVolInv(Vector v);
	  double GetDevInv(Vector v);

	  // Tensor functions
	  double DoubleDot2_2(Vector v1, Vector v2);
	  Vector DoubleDot2_4(Vector v1, Matrix m1);
	  Vector DoubleDotStrain2_4(Vector v1, Matrix m1);
	  Vector DoubleDot4_2(Matrix m1, Vector v2);
	  Matrix DoubleDot4_4(Matrix m1, Matrix m2);
	  Matrix DoubleDotStrain4_4(Matrix m1, Matrix m2);
	  Matrix Dyadic2_2(Vector v1, Vector v2);


	  //parameters
	  static const double one3 ;
	  static const double two3 ;
	  static const double root23 ;


};


#endif
