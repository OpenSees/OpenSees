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

#ifndef DruckerPrager_h
#define DruckerPrager_h

// $Revision: 1.1 $
// $Date: 2010-02-04 00:44:04 $
// $Source: /usr/local/cvs/OpenSees/SRC/material/nD/DruckerPrager.h,v $
                                                                      
// Written: Kathryn Petek, Peter Mackenzie-Helnwein, and Pedro Arduino
// Created: 12/04
//
// Description: This file contains the class definition for DruckerPrager. 
//

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include <NDMaterial.h>
#include <Matrix.h>
#include <Vector.h>
#include <ID.h>

class DruckerPrager : public NDMaterial
{
  public:
    // Full Constructor
    DruckerPrager(int tag, int classTag, double bulk, double shear,
		  double s_y, double r, double r_bar, double Kinfinity, double Kinit, 
		  double d1, double d2, double H, double t, double massDen = 0.0, double atm = 101.0);

  // Elastic Constructor
  //	  DruckerPrager(int tag, double bulk, double shear);
  
  //Null Constructor
  DruckerPrager();
  
  //Destructor
  ~DruckerPrager();
  
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

  double getRho(void) {return massDen;};
  
 protected:
  
  //material parameters
  double mKref;			// reference Bulk Modulus 
  double mGref;			// reference Shear Modulus
  double mPatm;		    // reference stress first invariant (pressure)
  double mK;			// bulk modulus 
  double mG;			// shear modulus
  double msigma_y;		// yield strength 
  double mrho;			// volumetric term
  double mrho_bar;		// nonassociative flow term
  double mKinf;			// nonlinear isotropic hardening term
  double mKo;			// nonlinear isotropic hardening term
  double mdelta1; 		// exponential hardening term for drucker prager surface
  double mdelta2;       // exponential hardening term for tension cutoff surface
  double mHard;			// hardening constant
  double mtheta;		// hardening constant
  double mTo;           // initial tension cutoff strength

  double massDen;
  
  //internal variables
  Vector mEpsilon;
  Vector mEpsilon_n_p;	// plastic strain vector at step n, trail e_p
  Vector mEpsilon_n1_p;	// plastic strain vector at step n+1 
  Vector mSigma;
  
  Vector mBeta_n;		// backstress at step n, beta_np1_trial = beta_n
  Vector mBeta_n1;		// backstress at step n+1
  
  double mHprime;		// derivative of linear kinematic hardening term 
  
  double mAlpha1_n;		// alpha1_n
  double mAlpha1_n1;	// alpha1_n+1
  double mAlpha2_n;		// alpha2_n
  double mAlpha2_n1;	// alpha2_n+1
  
  int mElastFlag;    // Flag to determine elastic behavior
  int mFlag;
  
  Matrix mCe;			// elastic tangent stiffness matrix
  Matrix mCep;			// elastoplastic tangent stiffness matrix
  Vector mI1;			// 2nd Order Identity Tensor	
  Matrix mIIvol;		// IIvol = I1 tensor I1  
  Matrix mIIdev;		// 4th Order Deviatoric Tensor
  
  Vector mState;		// state vector for output
  
  //functions
  void initialize();	// initializes variables
  int  updateElasticParam(void); //updated Elastic Parameters based on mean stress 
  
  //plasticity integration routine
  void plastic_integrator(void);
  
  double Kiso(double alpha1);		// isotropic hardening function
  double Kisoprime(double alpha1);	//
  double T(double alpha2);
  double deltaH(double dGamma);
  
  Vector getState();		// fills vector of state variables for output
  
  //parameters
  static const double one3 ;
  static const double two3 ;
  static const double root23 ;
};


#endif
