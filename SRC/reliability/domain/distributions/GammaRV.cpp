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
// $Date: 2003-03-04 00:44:33 $
// $Source: /usr/local/cvs/OpenSees/SRC/reliability/domain/distributions/GammaRV.cpp,v $


//
// Written by Terje Haukaas (haukaas@ce.berkeley.edu) 
//

#include <GammaRV.h>
#include <math.h>
#include <string.h>
#include <Vector.h>
#include <classTags.h>

GammaRV::GammaRV(int passedTag, 
		 double passedMean,
		 double passedStdv,
		 double passedStartValue)
:RandomVariable(passedTag, RANDOM_VARIABLE_gamma)
{
	tag = passedTag ;
	k = (passedMean/passedStdv) * (passedMean/passedStdv);
	lambda = passedMean / (passedStdv*passedStdv);
	startValue = passedStartValue;
}
GammaRV::GammaRV(int passedTag, 
		 double passedParameter1,
		 double passedParameter2,
		 double passedParameter3,
		 double passedParameter4,
		 double passedStartValue)
:RandomVariable(passedTag, RANDOM_VARIABLE_gamma)
{
	tag = passedTag ;
	k = passedParameter1;
	lambda = passedParameter2;
	startValue = passedStartValue;
}
GammaRV::GammaRV(int passedTag, 
		 double passedMean,
		 double passedStdv)
:RandomVariable(passedTag, RANDOM_VARIABLE_gamma)
{
	tag = passedTag ;
	k = (passedMean/passedStdv) * (passedMean/passedStdv);
	lambda = passedMean / (passedStdv*passedStdv);
	startValue = getMean();
}
GammaRV::GammaRV(int passedTag, 
		 double passedParameter1,
		 double passedParameter2,
		 double passedParameter3,
		 double passedParameter4)
:RandomVariable(passedTag, RANDOM_VARIABLE_gamma)
{
	tag = passedTag ;
	k = passedParameter1;
	lambda = passedParameter2;
	startValue = getMean();
}


GammaRV::~GammaRV()
{
}


void
GammaRV::Print(OPS_Stream &s, int flag)
{
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
		// Take a Newton step
		x_new = x_old - f/df;
		// Check convergence; quit or continue
		if (fabs(1.0-fabs(x_old/x_new)) < tol) {
			result = x_new;
		}
		else {
			if (i==100) {
				opserr << "WARNING: Did not converge to find inverse CDF!" << endln;
				result = 0.0;
			}
			else {
			x_old = x_new;
			}
		}
	}
	return result;
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


double 
GammaRV::getStartValue()
{
	return startValue;
}



double GammaRV::getParameter1()  {return k;}
double GammaRV::getParameter2()  {return lambda;}
double GammaRV::getParameter3()  {opserr<<"No such parameter in r.v. #"<<tag<<endln; return 0.0;}
double GammaRV::getParameter4()  {opserr<<"No such parameter in r.v. #"<<tag<<endln; return 0.0;}





double
GammaRV::gammaFunction(double x)
{
	double res;

	if (x==0 || ( x < 0.0 && floor(x)==x ) )  {
		opserr << "Invalid input to the gamma function" << endln;
	}
	else {
		Vector p(9);
		Vector q(9);
		Vector c(8);

		p(0) = 0.0;
		p(1) = -1.71618513886549492533811e+0; 
		p(2) = 2.47656508055759199108314e+1;
		p(3) = -3.79804256470945635097577e+2; 
		p(4) = 6.29331155312818442661052e+2;
		p(5) = 8.66966202790413211295064e+2; 
		p(6) = -3.14512729688483675254357e+4;
		p(7) = -3.61444134186911729807069e+4; 
		p(8) = 6.64561438202405440627855e+4;
		
		q(0) = 0.0;
		q(1) = -3.08402300119738975254353e+1; 
		q(2) = 3.15350626979604161529144e+2;
		q(3) = -1.01515636749021914166146e+3; 
		q(4) = -3.10777167157231109440444e+3;
		q(5) = 2.25381184209801510330112e+4; 
		q(6) = 4.75584627752788110767815e+3;
		q(7) = -1.34659959864969306392456e+5; 
		q(8) = -1.15132259675553483497211e+5;
		
		c(0) = 0.0;
		c(1) = -1.910444077728e-03; 
		c(2) = 8.4171387781295e-04;
		c(3) = -5.952379913043012e-04; 
		c(4) = 7.93650793500350248e-04;
		c(5) = -2.777777777777681622553e-03; 
		c(6) = 8.333333333333333331554247e-02;
		c(7) = 5.7083835261e-03;

		double pi = 3.14159265358979;
		double y;
		double y1;
		double fact;
		double x1;
		double xn;
		double ysq;
		double sum;
		double spi;
		bool flag01 = false;
		bool flag1_12 = false;
		bool flagNegative = false;

		// If x is negative
		if (x<0.0) {
			y = -x;
			y1 = floor(y);
			res = y - y1;
			fact = -pi / sin(pi*res) * (1 - 2*fmod(y1,2));
			x = y + 1;
			flagNegative = true;
		}
		// Now x is positive

		// Map x in interval [0,1] to [1,2]
		if (x<1.0) {
			x1 = x;
			x = x1 + 1.0;
			flag01 = true;
		}

		// Map x in interval [1,12] to [1,2]
		if (x<12.0) {
			xn = floor(x) - 1;
			x = x - xn;
			// Evaluate approximation for 1 < x < 2
			double z = x - 1.0;
			double xnum = 0.0;
			double xden = xnum + 1.0;

			
			for (int i = 1 ; i<=8; i++ ) {
				xnum = (xnum + p(i)) * z;
				xden = xden * z + q(i);
			}

			res = xnum / xden + 1.0;
			flag1_12 = true;
		}

		// Adjust result for case  0.0 < x < 1.0
		if (flag01)  {
			res = res / x1;
		}
		else if (flag1_12){   // Adjust result for case  2.0 < x < 12.0
			double max_xn = xn;
			for (int m=1;  m<=max_xn; m++) {
				res = res * x;
				x = x + 1;
				xn = xn - 1;
			}
		}

		// Evaluate approximation for x >= 12
		if (x>=12.0) {
			y = x;
			ysq = y * y;
			sum = c(7);
			for (int i = 1; i<=6; i++ ) {
				sum = sum / ysq + c(i);
			}

			spi = 0.9189385332046727417803297;
			sum = sum / y - y + spi;
			sum = sum + (y-0.5)*log(y);
			res = exp(sum);
		}

		if (flagNegative) {
			res = fact / res;
		}
	}
	
	return res;
}



double 
GammaRV::incompleteGammaFunction(double x, double a)
{
	double gam = log(gammaFunction(a));
	double b = x;
	if (x==0.0) {
		b = 0.0;
	}
	if (a==0.0) {
		b = 1.0;
	}
	// Series expansion for x < a+1
	if (a!=0.0 && x!=0.0 && x<a+1.0) {
		double ap = a;
		double sum = 1.0/ap;
		double del = sum;
		while (fabs(del) >= 1.0e-10*fabs(sum)) {
			ap = ap + 1.0;
			del = x * del / ap;
			sum = sum + del;
		}
		b = sum * exp(-x + a*log(x) - gam);
	}

	// Continued fraction for x >= a+1
	if ((a != 0) && (x != 0) && (x >= a+1.0))  {
		double a0 = 1.0;
		double a1 = x;
		double b0 = 0.0;
		double b1 = a0;
		double fac = 1.0;
		double n = 1;
		double g = b1;
		double gold = b0;
		while (fabs(g-gold) >= 1.0e-10*fabs(g))  {
			gold = g;
			double ana = n - a;
			a0 = (a1 + a0 *ana) * fac;
			b0 = (b1 + b0 *ana) * fac;
			double anf = n*fac;
			a1 = x * a0 + anf * a1;
			b1 = x * b0 + anf * b1;
			fac = 1.0 / a1;
			g = b1 * fac;
			n = n + 1.0;
		}
		b = 1 - exp(-x + a*log(x) - gam) * g;
	}

	return b;
}
