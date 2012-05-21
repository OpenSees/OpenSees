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
                                                                        
// $Revision: 1.6 $
// $Date: 2008-02-29 19:47:19 $
// $Source: /usr/local/cvs/OpenSees/SRC/reliability/analysis/direction/SQPsearchDirectionMeritFunctionAndHessian.cpp,v $


//
// Written by Terje Haukaas (haukaas@ce.berkeley.edu) 
//

#include <SQPsearchDirectionMeritFunctionAndHessian.h>
#include <SearchDirection.h>
#include <MeritFunctionCheck.h>
//#include <HessianApproximation.h>
#include <Vector.h>
#include <Matrix.h>
#include <math.h>


SQPsearchDirectionMeritFunctionAndHessian::SQPsearchDirectionMeritFunctionAndHessian(
															double pc_bar, 
															double pe_bar)
:SearchDirection(), MeritFunctionCheck()
{
	//theHessianApproximation = 0;

	// Parameters
	alpha = 0.0;
	c_bar = pc_bar;
	e_bar = pe_bar;


	// History variables
	B = 0;
	delta = 1.0;
	c = c_bar;
	lambda = 1.0;


	// Check parameters
	if (c_bar < 1.0) {
		opserr << "ERROR: Parameter c_bar in SQP algorithm is invalid." << endln;
	}
	if (e_bar > 1.0) {
		opserr << "ERROR: Parameter e_bar in SQP algorithm is invalid." << endln;
	}
}

SQPsearchDirectionMeritFunctionAndHessian::~SQPsearchDirectionMeritFunctionAndHessian()
{
	if (B != 0) 
		delete B;
}


//int
//SQPsearchDirectionMeritFunctionAndHessian::setHessianApproximation(HessianApproximation *passedHessianApproximation)
//{
	//theHessianApproximation = passedHessianApproximation;
	//return 0;
//}




// SEARCH DIRECTION METHODS
const Vector&
SQPsearchDirectionMeritFunctionAndHessian::getSearchDirection()
{
	return searchDirection;
}

int
SQPsearchDirectionMeritFunctionAndHessian::computeSearchDirection(
							int stepNumber,
							const Vector &u, 
							double g, 
							const Vector &gradG )
{
	// Initial declarations
	int i,j;


	// Number of random variables in the problem
	int nrv = u.Size();


	// If it's the first step; set the Hessian approximation to unity
	// and initialize history parameters. 
	if (stepNumber == 1) {
		//if (theHessianApproximation != 0) {
		//	theHessianApproximation->setHessianToIdentity(nrv);
		//}
		//else {
		//	opserr << "WARNING: SQPsearchDirectionMeritFunctionAndHessian::computeSearchDirection() -- " << endln
		//		<< "theHessianApproximation has not been set! " << endln;
		//}
		delta = 1.0;
		c = c_bar;
		lambda = 1.0;
	}


	// Get the Hessian approximation
	Matrix Hessian;
	//if (theHessianApproximation != 0) {
	//	Hessian = theHessianApproximation->getHessianApproximation();
	//}
	//else {
	//	opserr << "WARNING: SQPsearchDirectionMeritFunctionAndHessian::computeSearchDirection() -- " << endln
	//		<< "theHessianApproximation has not been set! " << endln;
	//}

	
	// Establish coefficient matrix
	Matrix A((nrv+1),(nrv+1));
	for (i=0; i<(nrv+1); i++) {
		for (j=0; j<(nrv+1); j++) {
			if (i==nrv && j==nrv) {
				A(i,j) = 0.0;
			}
			else {
				if (i<nrv && j<nrv) {
					A(i,j) = Hessian(i,j);
				}
				else if (i==nrv) {
					A(i,j) = gradG(j);
				}
				else {
					A(i,j) = gradG(i);
				}
			}
		}
	}


	// Establish right-hand side vector
	Vector b(nrv+1);
	for (i=0; i<nrv; i++) {
		b(i) = -u(i);
	}
	b(nrv) = -g;


	// Solve the system of equations by computing inverse (this should be done more efficiently!)
	//MatrixOperations theMatrixOperations(A);
	//theMatrixOperations.computeInverse();
	//Matrix invA = theMatrixOperations.getInverse();
	//Vector dAndKappa = invA ^ b; // why transpose A?  According to PEER report, this must be a bug -- MHS

	Vector dAndKappa(nrv+1);
	A.Solve(b, dAndKappa);

	// Separate out direction vector and kappa value
	Vector d(nrv);
	for (i=0; i<nrv; i++) {
		d(i) = dAndKappa(i);
	}
	kappa = dAndKappa(nrv);
	searchDirection = d;
	





	// Compute product d*B*d and dd
	Vector temp2 = Hessian ^ searchDirection;
	double dBd = temp2 ^ searchDirection;
	double dd = searchDirection^searchDirection;


	// Update the history parameter 'delta'
	double term1 = dBd/dd;
	if ( term1 < delta ) {
		delta = term1;
	}


	// Compute 'e' (needed to compute 'i')
	double e;
	if ( fabs(kappa-lambda)<1.0e-9 ) {
		e = e_bar;
	}
	else {
		e = dd / ((kappa-lambda)*(kappa-lambda));
	}


	// Compute 'i' (needed to compute 'c')
	double iii = -ceil(log(0.25*e*delta*(1-0.25*delta))/log(c_bar));
	if (iii < 0) {
		iii = 0;
	}


	// Update history parameter 'c'
	term1 = pow(c_bar,iii);
	if ( c < term1 ) {
		c = term1;
	}

	return 0;
}


double 
SQPsearchDirectionMeritFunctionAndHessian::getMeritFunctionValue(const Vector &u,
								 double g,
								 const Vector &grad_G)
{
	opserr << "WARNING: SQPsearchDirectionMeritFunctionAndHessian::getMeritFunctionValue() --" << endln
		<< " no explicit merit function value is computed." << endln;
	return 0.0;
}



int
SQPsearchDirectionMeritFunctionAndHessian::updateMeritParameters(const Vector &u,
								 double g,
								 const Vector &grad_G, int reschk)
{
	opserr << "WARNING: SQPsearchDirectionMeritFunctionAndHessian::updateMeritParameters() --" << endln
		<< " no explicit merit function value is computed." << endln;

	return 0;
}



int
SQPsearchDirectionMeritFunctionAndHessian::check(const Vector &u_old, 
						 double g_old, 
						 const Vector &grad_G_old, 
						 double stepSize,
						 const Vector &stepDirection,
						 double g_new, int reschk)
{
	// Have 'c' and 'delta' and 'lambda' as history parameters
	// and 'kappa' stored to be used in this method


	// Initial declarations
	int i;

	// Problem size
	int nrv = u_old.Size();


	// New point in standard normal space
	Vector u_new = u_old + stepSize*stepDirection;


	// Check that the factor alpha is set (since this is a dual purpose class...)
	if (alpha == 0.0) {
		opserr << "ERROR: SQPsearchDirectionMeritFunctionAndHessian::check()" << endln
			<< " the alpha factor is not set! " << endln;
	}


	// Left side of merit function check
	double tempNewLambda = lambda+stepSize*(kappa-lambda);
	double LHS = 0.5*(u_new^u_new) + tempNewLambda*g_new + 0.5*c*g_new*g_new;


	// Right side of merit function check
	double term1 = 0.5*(u_old^u_old) + lambda*g_old + 0.5*c*g_old*g_old;
	Vector gradient = u_old + lambda*grad_G_old + (c*g_old)*grad_G_old;
	Vector extendedGradient(nrv+1);
	for (i=0; i<nrv; i++) {
		extendedGradient(i) = gradient(i);
	}
	extendedGradient(nrv) = g_old;
	Vector extendedDirection(nrv+1);
	for (i=0; i<nrv; i++) {
		extendedDirection(i) = stepDirection(i);
	}
	extendedDirection(nrv) = (kappa-lambda);
//	double RHS = term1 + alpha*stepSize*(extendedGradient^extendedDirection);
//	(Maybe use absolute value of 'g' and 'kappa-lambda'?)
	double RHS = term1 + alpha*stepSize*(gradient^stepDirection);

	// Do the check
	if (  LHS  <=  RHS ) {
		return 0;  // ok
	}
	else {
		return -1; // not ok
	}

}





int
SQPsearchDirectionMeritFunctionAndHessian::setAlpha(double palpha)
{
	alpha = palpha;

	if (alpha > 0.5) {
		opserr << "ERROR: Parameter alpha in SQP algorithm is invalid." << endln;
	}

	return 0;
}



const Matrix&
SQPsearchDirectionMeritFunctionAndHessian::getHessianApproximation()
{
	return (*B);
}

int
SQPsearchDirectionMeritFunctionAndHessian::updateHessianApproximation(const Vector &u_old,
								      double g_old,
								      const Vector &gradG_old,
								      double stepSize,
								      const Vector &searchDirection,
								      double g_new,
								      const Vector &gradG_new)
{
	if (B == 0) {
		this->setHessianToIdentity(u_old.Size());
	}

	// Initial declarations
	int i,j;

	// Re-compute the new point in the standard normal space
	Vector u_new = u_old + stepSize * searchDirection;


	// The gradient of the Lagrangian (the old one)
	Vector gradL_old = u_old + lambda*gradG_old;


	// At this point we should update the lambda (the Lagrangian multiplier)
	lambda = lambda + stepSize * (kappa-lambda);


	// The gradient of the Lagrangian (the new one)
	Vector gradL_new = u_new + lambda*gradG_new;


	// Computations to obtain the vector 'q'
	Vector qPrime = gradL_new - gradL_old;
	
	Vector temp2 = (*B) ^ searchDirection;

	double dBd = temp2 ^ searchDirection ;

	double dqPrime = searchDirection ^ qPrime;

	double theta;
	if (dqPrime >= 0.2*stepSize*dBd) {
		theta = 1.0;
	}
	else {
		theta = (0.8*stepSize*dBd) / (stepSize*dBd-dqPrime);
	}
	Vector q = theta*qPrime + (1-theta)*stepSize*((*B)^searchDirection);


	// Update the Hessian approximation
	Vector Bd = (*B)^searchDirection;
	int nrv = Bd.Size();
	Matrix Bdd(nrv,nrv);
	for (i=0; i<nrv; i++) {
		for (j=0; j<nrv; j++) {
			Bdd(i,j) = Bd(i)*searchDirection(j);
		}
	}
	Matrix BddB = Bdd*(*B);
	Matrix qq(nrv,nrv);
	for (i=0; i<nrv; i++) {
		for (j=0; j<nrv; j++) {
			qq(i,j) = q(i)*q(j);
		}
	}
	Matrix newB = (*B) + 1.0/(stepSize*(q^searchDirection))*(qq) - (1.0/(dBd))*(BddB);
	(*B) = newB;


	return 0;
}




int
SQPsearchDirectionMeritFunctionAndHessian::setHessianToIdentity(int size)
{

	// Check that is exists
	if (B == 0) {
		B = new Matrix(size,size);
	}

	// Set it to unity
	for (int i=0; i<size; i++) {
		for (int j=0; j<size; j++) {
			if (i==j) {
				(*B)(i,j) = 1.0;
			}
			else {
				(*B)(i,j) = 0.0;
			}
		}
	}

	return 0;
}
