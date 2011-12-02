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
                                                                        
// $Revision: 1.2 $
// $Date: 2003-10-27 23:45:42 $
// $Source: /usr/local/cvs/OpenSees/SRC/reliability/analysis/direction/SQPsearchDirectionMeritFunctionAndHessian.h,v $


//
// Written by Terje Haukaas (haukaas@ce.berkeley.edu) 
//

#ifndef SQPsearchDirectionMeritFunctionAndHessian_h
#define SQPsearchDirectionMeritFunctionAndHessian_h

#include <SearchDirection.h>
#include <MeritFunctionCheck.h>
#include <HessianApproximation.h>
#include <Vector.h>

class SQPsearchDirectionMeritFunctionAndHessian : public SearchDirection, public MeritFunctionCheck, public HessianApproximation
{

public:
	SQPsearchDirectionMeritFunctionAndHessian(double c_bar, double e_bar);
	~SQPsearchDirectionMeritFunctionAndHessian();

	// METHODS FOR SEARCH DIRECTION
	int computeSearchDirection(	int stepNumber, 
								Vector passed_u, 
								double passed_gFunctionValue, 
								Vector passedGradientInStandardNormalSpace);
	Vector getSearchDirection();

	// METHODS FOR MERIT FUNCTION CHECK
	int	check(Vector u_old, 
			  double g_old, 
			  Vector grad_G_old, 
			  double stepSize,
			  Vector stepDirection,
			  double g_new);
	double getMeritFunctionValue(Vector u, double g, Vector grad_G);
	int updateMeritParameters(Vector u, double g, Vector grad_G);

	int setAlpha(double alpha);

	// METHODS FOR HESSIAN APPROXIMATION
	Matrix  getHessianApproximation();
	int     setHessianToIdentity(int size);
	int     setHessianApproximation(HessianApproximation *theHessianApproximation);
	int     updateHessianApproximation(Vector u_old,
									   double g_old,
									   Vector gradG_old,
									   double stepSize,
									   Vector searchDirection,
									   double g_new,
									   Vector grad_G_new);

protected:

private:
	HessianApproximation *theHessianApproximation;

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
