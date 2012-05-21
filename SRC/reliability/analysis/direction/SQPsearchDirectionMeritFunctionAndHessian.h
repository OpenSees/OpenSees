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
                                                                        
// $Revision: 1.4 $
// $Date: 2008-02-29 19:47:19 $
// $Source: /usr/local/cvs/OpenSees/SRC/reliability/analysis/direction/SQPsearchDirectionMeritFunctionAndHessian.h,v $


//
// Written by Terje Haukaas (haukaas@ce.berkeley.edu) 
//

#ifndef SQPsearchDirectionMeritFunctionAndHessian_h
#define SQPsearchDirectionMeritFunctionAndHessian_h

#include <SearchDirection.h>
#include <MeritFunctionCheck.h>
//#include <HessianApproximation.h>
#include <Vector.h>
#include <Matrix.h>

class SQPsearchDirectionMeritFunctionAndHessian : public SearchDirection, public MeritFunctionCheck
{

public:
	SQPsearchDirectionMeritFunctionAndHessian(double c_bar, double e_bar);
	~SQPsearchDirectionMeritFunctionAndHessian();

	// METHODS FOR SEARCH DIRECTION
	int computeSearchDirection(int stepNumber, 
				   const Vector &passed_u, 
				   double passed_gFunctionValue, 
				   const Vector &passedGradientInStandardNormalSpace);
	const Vector &getSearchDirection();

	// METHODS FOR MERIT FUNCTION CHECK
	int check(const Vector &u_old, 
		  double g_old, 
		  const Vector &grad_G_old, 
		  double stepSize,
		  const Vector &stepDirection,
		  double g_new,
			  /////S added by K Fujimura ////
			  int reschk=0);
	//double getMeritFunctionValue(Vector u, double g, Vector grad_G);
	//int updateMeritParameters(Vector u, double g, Vector grad_G, int reschk=0);

	int setAlpha(double alpha);

	double getMeritFunctionValue(const Vector &u, double g,
				     const Vector &grad_G);
	int updateMeritParameters(const Vector &u, double g,
				  const Vector &grad_G, int reschk=0);
			  /////E added by K Fujimura ////


	// METHODS FOR HESSIAN APPROXIMATION
    // KRM 5-19-2012
    // this needs to be updated so classes are not being mixed, use proper OO
	const Matrix &getHessianApproximation();
	int     setHessianToIdentity(int size);
	//int     setHessianApproximation(HessianApproximation *theHessianApproximation);
	int     updateHessianApproximation(const Vector &u_old,
					   double g_old,
					   const Vector &gradG_old,
					   double stepSize,
					   const Vector &searchDirection,
					   double g_new,
					   const Vector &grad_G_new);

protected:

private:
	//HessianApproximation *theHessianApproximation;

	// Parameters
	double c_bar, e_bar, alpha;

	// To be returned...
	Vector searchDirection;
	double stepSize;

	// History data
	Matrix *B;
	double delta, c, lambda; 
	double kappa;
};

#endif
