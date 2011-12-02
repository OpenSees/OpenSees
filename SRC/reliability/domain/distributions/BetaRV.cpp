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
// $Source: /usr/local/cvs/OpenSees/SRC/reliability/domain/distributions/BetaRV.cpp,v $


//
// Written by Terje Haukaas (haukaas@ce.berkeley.edu) 
//

#include <BetaRV.h>
#include <GammaRV.h>
#include <math.h>
#include <string.h>
#include <classTags.h>
#include <OPS_Globals.h>

BetaRV::BetaRV(int passedTag, 
		 double passedParameter1,
		 double passedParameter2,
		 double passedParameter3,
		 double passedParameter4,
		 double passedStartValue)
:RandomVariable(passedTag, RANDOM_VARIABLE_beta)
{
	tag = passedTag ;
	a = passedParameter1;
	b = passedParameter2;
	q = passedParameter3;
	r = passedParameter4;
	startValue = passedStartValue;
}
BetaRV::BetaRV(int passedTag, 
		 double passedParameter1,
		 double passedParameter2,
		 double passedParameter3,
		 double passedParameter4)
:RandomVariable(passedTag, RANDOM_VARIABLE_beta)
{
	tag = passedTag ;
	a = passedParameter1;
	b = passedParameter2;
	q = passedParameter3;
	r = passedParameter4;
	startValue = getMean();
}


BetaRV::~BetaRV()
{
}


void
BetaRV::Print(OPS_Stream &s, int flag)
{
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


double 
BetaRV::getStartValue()
{
	return startValue;
}

double BetaRV::getParameter1()  {return a;}
double BetaRV::getParameter2()  {return b;}
double BetaRV::getParameter3()  {return q;}
double BetaRV::getParameter4()  {return r;}



double 
BetaRV::betaFunction(double q, double r)
{
/*	OLD CODE: 
	GammaRV *aGammaRV = new GammaRV(1, 0.0, 1.0, 0.0);
	double par1,par2,par3;
	par1 = aGammaRV->gammaFunction(q);
	par2 = aGammaRV->gammaFunction(q);
	par3 = aGammaRV->gammaFunction(q+r);
	delete aGammaRV;
	return par1*par2/par3;
*/

	// Matlab definition of the beta function:
	//    y = exp(gammaln(q)+gammaln(r)-gammaln(q+r));
	//    ... where gammaln(.) = ln(gamma(.))
	// So, try this instead:
	GammaRV *aGammaRV = new GammaRV(1, 0.0, 1.0, 0.0);
	double gammaq,gammar,gammaqpr;
	gammaq = aGammaRV->gammaFunction(q);
	gammar = aGammaRV->gammaFunction(r);
	gammaqpr = aGammaRV->gammaFunction(q+r);
	delete aGammaRV;
	double loggammaq,loggammar,loggammaqpr;
	loggammaq = log(gammaq);
	loggammar = log(gammar);
	loggammaqpr = log(gammaqpr);
	return exp(loggammaq+loggammar-loggammaqpr);
}
