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
                                                                        
// $Revision: 1.7 $
// $Date: 2007-02-17 21:27:23 $
// $Source: /usr/local/cvs/OpenSees/SRC/reliability/domain/distributions/BetaRV.cpp,v $


//
// Written by Terje Haukaas (haukaas@ce.berkeley.edu) 
//

#include <BetaRV.h>
#include <Vector.h>
#include <cmath>

BetaRV::BetaRV(int passedTag, 
			   const Vector &passedParameters)
:RandomVariable(passedTag, RANDOM_VARIABLE_beta), startValue(0)
{
	if (passedParameters.Size() != 4) {
		opserr << "Beta RV requires 4 parameters, a b q and r, for RV with tag " <<
		this->getTag() << endln;
		
		// this will create terminal errors
		a = 0;
		b = 0;
		q = 0;
		r = 0;
		
	} else {
		
		a = passedParameters(0);
		b = passedParameters(1);
		q = passedParameters(2);
		r = passedParameters(3);
	}
	
}


BetaRV::~BetaRV()
{
}


const char *
BetaRV::getType()
{
	return "BETA";
}


double 
BetaRV::getMean()
{
	return (a*r+b*q)/(q+r);
}


double 
BetaRV::getStdv()
{
	return ((b-a)/(q+r)) * sqrt(q*r/(q+r+1));
}


const Vector &
BetaRV::getParameters(void) {
	static Vector temp(4);
	temp(0) = a;
	temp(1) = b;
	temp(2) = q;
	temp(3) = r;
	return temp;
}


double
BetaRV::getPDFvalue(double rvValue)
{
	double result;
	if ( a <= rvValue && rvValue <= b ) {
		double par1 = pow(rvValue-a,q-1.0);
		double par2 = pow(b-rvValue,r-1.0);
		double par3 = betaFunction(q,r);
		double par4 = pow(b-a,q+r-1.0);
		result = par1*par2/(par3*par4);
	}
	else {
		result = 0.0;
	}
	return result;
}


double
BetaRV::getCDFvalue(double rvValue)
{

	double result = 0.0;

	if ( a < rvValue && rvValue < b ) {
		// There exists no closed form expression for the Beta CDF.
		// In this preliminary implementation of the Beta random variable,
		// numerical integration - using Simpsons rule - is employed.
		// The aim is to integrate the PDF from 'a' to 'rvValue'.
		int n_2 = 100; // Half the number of intervals
		double h = rvValue-a;
		double fa = getPDFvalue(a);
		double fb = getPDFvalue(rvValue);
		double sum_fx2j = 0.0;
		double sum_fx2j_1 = 0.0;
		for (int j=1;  j<=n_2;  j++) {
			sum_fx2j = sum_fx2j + getPDFvalue(   (double) (a+(j*2)*h/(2*n_2))   );
			sum_fx2j_1 = sum_fx2j_1 + getPDFvalue(   (double)(a+(j*2-1)*h/(2*n_2))   );
		}
		sum_fx2j = sum_fx2j - getPDFvalue((double)(rvValue));
		result = h/(2*n_2)/3.0*(fa + 2.0*sum_fx2j + 4.0*sum_fx2j_1 + fb);
	}
	else if (rvValue<=a) {
		result = 0.0;
	}
	else {
		result = 1.0;
	}

	return result;
}


double
BetaRV::getInverseCDFvalue(double probValue)
{
	double result = 0.0;
	// Here we want to solve the nonlinear equation:
	//         probValue = getCDFvalue(x)
	// with respect to x. 
	// A Newton scheme to find roots - f(x)=0 - looks something like:
	//         x(i+1) = x(i) - f(xi)/f'(xi)
	// In our case the function f(x) is: f(x) = probValue - getCDFvalue(x)
	// The derivative of the function can be found approximately by a
	// finite difference scheme where e.g. stdv/200 is used as perturbation.
	double tol = 0.000001;
	double x_old = getMean();   // Start at the mean of the random variable
	double x_new;
	double f;
	double df;
	double h;
	double perturbed_f;
	for (int i=1;  i<=100;  i++ )  {

		// Evaluate function
		f = probValue - getCDFvalue(x_old);
		// Evaluate perturbed function
		h = getStdv()/200.0;
		perturbed_f = probValue - getCDFvalue(x_old+h);

		// Evaluate derivative of function
		df = ( perturbed_f - f ) / h;

		if ( fabs(df) < 1.0e-15) {
			opserr << "WARNING: BetaRV::getInverseCDFvalue() -- zero derivative " << endln
				<< " in Newton algorithm. " << endln;
		}
		else {

			// Take a Newton step
			x_new = x_old - f/df;
			
			// Check convergence; quit or continue
			if (fabs(1.0-fabs(x_old/x_new)) < tol) {
				return x_new;
			}
			else {
				if (i==100) {
					opserr << "WARNING: Did not converge to find inverse CDF!" << endln;
					return 0.0;
				}
				else {
					x_old = x_new;
				}
			
			}
		}
	}

	return result;
}


void
BetaRV::Print(OPS_Stream &s, int flag)
{
	s << "Beta RV #" << this->getTag() << endln;
	s << "\ta = " << a << ", b = " << b << endln;
	s << "\tq = " << q << ", r = " << r << endln;
}
