/* ****************************************************************** **
**    OpenSees - Open System for Earthquake Engineering Simulation    **
**          Pacific Earthquake Engineering Research Center            **
**                                                                    **
**                                                                    **
** (C) Copyright 2001, The Regents of the University of California    **
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
** Reliability module developed by:                                   **
**   Terje Haukaas (haukaas@ce.berkeley.edu)                          **
**   Armen Der Kiureghian (adk@ce.berkeley.edu)                       **
**                                                                    **
** ****************************************************************** */
                                                                        
// $Revision: 1.9 $
// $Date: 2007-11-06 19:32:36 $
// $Source: /usr/local/cvs/OpenSees/SRC/reliability/analysis/transformation/NatafProbabilityTransformation.h,v $


//
// Written by Terje Haukaas (haukaas@ce.berkeley.edu)
//

#ifndef NatafProbabilityTransformation_h
#define NatafProbabilityTransformation_h

#include <Vector.h>
#include <Matrix.h>
#include <ReliabilityDomain.h>
#include <MatrixOperations.h>

class NatafProbabilityTransformation : public ProbabilityTransformation
{

public:
	NatafProbabilityTransformation(ReliabilityDomain *passedReliabilityDomain,
						  int printFlag);
	~NatafProbabilityTransformation();

	int transform_x_to_u(Vector &u);
	int transform_u_to_x(const Vector &u, Vector &x);
	int getJacobian_x_to_u(Matrix &Jxu);
	int getJacobian_u_to_x(const Vector &u, Matrix &Jux);

	Vector meanSensitivityOf_x_to_u(const Vector &x, int rvTag);
	Vector stdvSensitivityOf_x_to_u(const Vector &x, int rvTag);

protected:

private:

	// Private data members
	ReliabilityDomain *theReliabilityDomain;
	int nrv;
	Matrix *correlationMatrix;
	
	// LAPACK representation from DPOTRF
	// Use DPOTRS to get "inverse"
	double *lapackA;
	double *lapackB;
	int lapackCholesky(void);
	bool isFactorized;

	int printFlag;

	// Private member functions
	void setCorrelationMatrix(int pertMeanOfThisRV, int pertStdvOfThisRV, double h);
	// Representing diag matrix as vector
	int getJacobian_z_x(const Vector &z, Vector &Jzx);
	int z_to_x(const Vector &z, Vector &x);
	int x_to_z(Vector &z);

	// Auxiliary member functions for manual evaluation of 
	// the integral equation to find Nataf correlation
	double phi2(double z_i, 
				double z_j, 
				double rho);
	double integrand(int rv_i,
				     double z_i, 
				     double mean_i,
				     double stdv_i, 
				     int rv_j,
				     double z_j,
				     double mean_j, 
				     double stdv_j,
				     double rho);
	double doubleIntegral(int rv_i,
						  double mean_i,
						  double stdv_i, 
						  int rv_j,
						  double mean_j, 
						  double stdv_j,
						  double rho);
	double residualFunction(double rho_original, 
						    double rho,
						    int rv_i, 
						    double mean_i, 
						    double stdv_i, 
						    int rv_j, 
						    double mean_j, 
						    double stdv_j);
	double solveForCorrelation(int rv_i, int rv_j, double rho_original);
};

#endif
