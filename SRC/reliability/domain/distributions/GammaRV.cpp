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
                                                                        
// $Revision: 1.11 $
// $Date: 2008-05-08 15:32:54 $
// $Source: /usr/local/cvs/OpenSees/SRC/reliability/domain/distributions/GammaRV.cpp,v $


//
// Written by Terje Haukaas (haukaas@ce.berkeley.edu) 
//

#include <GammaRV.h>
#include <NormalRV.h>
#include <cmath>

GammaRV::GammaRV(int passedTag, double passedMean, double passedStdv)
:RandomVariable(passedTag, RANDOM_VARIABLE_gamma)
{
	int setp = setParameters(passedMean,passedStdv);
	if (setp < 0)
		opserr << "Error setting parameters in Gamma RV with tag " << this->getTag() << endln;
}


GammaRV::GammaRV(int passedTag, const Vector &passedParameters)
:RandomVariable(passedTag, RANDOM_VARIABLE_gamma)
{
	if (passedParameters.Size() != 2) {
		opserr << "Gamma RV requires 2 parameters, k and lambda, for RV with tag " <<
		this->getTag() << endln;
		
		// this will create terminal errors
		k = 0;
		lambda = 0;
		
	} else {
		
		k = passedParameters(0);
		lambda = passedParameters(1);
	}
}


GammaRV::~GammaRV()
{
}


const char *
GammaRV::getType()
{
	return "GAMMA";
}


double 
GammaRV::getMean()
{
	return k/lambda;
}


double 
GammaRV::getStdv()
{
	return sqrt(k)/lambda;
}


const Vector &
GammaRV::getParameters(void) {
	static Vector temp(2);
	temp(0) = k;
	temp(1) = lambda;
	return temp;
}


int 
GammaRV::setParameters(double mean, double stdv)
{
	k = (mean/stdv) * (mean/stdv);
	lambda = mean / (stdv*stdv);
	
	return 0;
}


double
GammaRV::getPDFvalue(double rvValue)
{
	double result;
	if ( 0.0 < rvValue ) {
		result = lambda*pow((lambda*rvValue),(k-1.0))*exp(-lambda*rvValue) / gammaFunction(k);
	}
	else {
		result = 0.0;
	}

	return result;
}


double
GammaRV::getCDFvalue(double rvValue)
{
	double result;
	if ( 0.0 < rvValue ) {
		result = incompleteGammaFunction(k,(lambda*rvValue));
	}
	else {
		result = 0.0;
	}
	
	return result;
}


double
GammaRV::getInverseCDFvalue(double probValue)
{
	// Here we want to solve the nonlinear equation:
	//         probValue = getCDFvalue(x)
	// with respect to x. 
	// A Newton scheme to find roots - f(x)=0 - looks something like:
	//         x(i+1) = x(i) - f(xi)/f'(xi)
	// In our case the function f(x) is: f(x) = probValue - getCDFvalue(x)
	// The derivative of the function can be found approximately by a
	// finite difference scheme where e.g. stdv/200 is used as perturbation.
	//
	// Note to writer of above: Newton does not converge for all cases, particularly those
	// that look like CDFs when the start point is not close to the zero point.  
	
	// bootstrap with estimate from equivalent normal distribution
	static NormalRV normRV( 1, getMean(), getStdv());
    double x_old = normRV.getInverseCDFvalue(probValue);
	
	double tol = 1.0e-8;
	double result = x_old;

	double x_new;
	double f;
	double df;
	for (int i = 1; i <= 50;  i++ )  {
		// Evaluate function
		f = getCDFvalue(x_old) - probValue;
		
		// evaluate derivative
		df = getPDFvalue(x_old);
		
		if ( fabs(df) > 1e-8 ) {
			// Take a Newton step
			x_new = x_old - f/df;
			
			// Check convergence; quit or continue
			if (fabs(f/df) < tol)
				return x_new;
		} else {
			// gradient zero
			return x_old;
		}
		
		x_old = x_new;
	}
	
	opserr << "WARNING: Did not converge to find inverse CDF!" << endln;
	return result;
}


void
GammaRV::Print(OPS_Stream &s, int flag)
{
	s << "Gamma RV #" << this->getTag() << endln;
	s << "\tk = " << k << endln;
	s << "\tlambda = " << lambda << endln;
}
