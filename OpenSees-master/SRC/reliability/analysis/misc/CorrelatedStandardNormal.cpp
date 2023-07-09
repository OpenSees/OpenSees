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
                                                                        
// $Revision: 1.1 $
// $Date: 2008-05-11 19:52:54 $
// $Source: /usr/local/cvs/OpenSees/SRC/reliability/analysis/misc/CorrelatedStandardNormal.cpp,v $

#include <CorrelatedStandardNormal.h>
#include <ReliabilityDomain.h>
#include <NormalRV.h>
#include <math.h>
#include <string.h>
#include <float.h>
#include <stdlib.h>


CorrelatedStandardNormal::CorrelatedStandardNormal(double corrcoeff)
{
	setCorrelation(corrcoeff);
	
	// test getCDF using different integration schemes
	//testCDF();
}


CorrelatedStandardNormal::~CorrelatedStandardNormal()
{

}


int
CorrelatedStandardNormal::setCorrelation(double corrcoeff)
{
	rho = corrcoeff;
	
	if (rho > 1 || rho < -1) {
		opserr << "Warning: correlation coefficient in CorrelatedStandardNormal must be [-1,1], setting to zero"
			<< endln;
		rho = 0;
	}
	
	return 0;
}


double
CorrelatedStandardNormal::bivariatePDF(double b1, double b2, double r)
{
	double pdf = 0.0;
	static const double pi = acos(-1.0);
	
	double oneMinusr2 = 1.0 - r*r;
	if ( oneMinusr2 <= DBL_EPSILON )
		return 0.0;
	
	pdf = 1.0/(2.0*pi*sqrt(oneMinusr2)) * exp( -1.0/2.0/(oneMinusr2)*( b1*b1 + b2*b2 - 2.0*r*b1*b2 ) );
	if (pdf < 0)
		pdf = 0.0;

	return pdf;
}


double 
CorrelatedStandardNormal::exponentialForm(double b1, double b2, double th)
{
	if ( th <= DBL_EPSILON )
		return 0.0;
	
	return exp(-1.0/2.0/sin(th)/sin(th) * (b1*b1 + b2*b2 - 2.0*b1*b2*cos(th)));
}


double
CorrelatedStandardNormal::getPDF(double b1, double b2)
{
	return bivariatePDF(b1,b2,rho);
}


double
CorrelatedStandardNormal::SimpsonOwen(double a, double b, double beta1, double beta2)
{
	// Composite Simpson's numerical integration
	// n must be even, endpoints a and b
	int n = 1600;
	double h, fa, fb;
	double integral = 0.0;
	
	h = (b-a)/n;
	fa = bivariatePDF(beta1,beta2,a);
	fb = bivariatePDF(beta1,beta2,b);
	
	double sum_fx2j_1 = 0.0;
	double sum_fx2j_2 = 0.0;
	for (int j=2; j <= n/2;  j++)
		sum_fx2j_1 = sum_fx2j_1 + bivariatePDF(beta1,beta2,a+h*(2*j-2));
	for (int j=1; j <= n/2;  j++)
		sum_fx2j_2 = sum_fx2j_2 + bivariatePDF(beta1,beta2,a+h*(2*j-1));

	integral = h/3.0*(fa + 2.0*sum_fx2j_1 + 4.0*sum_fx2j_2 + fb);

	return integral;
}


double
CorrelatedStandardNormal::SimpsonSheppard(double r, double beta1, double beta2)
{
	// Composite Simpson's numerical integration
	// n must be even, endpoints a and b
	int n = 1600;
	double h, fa, fb;
	double integral = 0.0;
	static const double pi = acos(-1.0);
	double b = pi;
	double a = acos(r);
	
	h = (b-a)/n;
	fa = exponentialForm(beta1,beta2,a);
	fb = exponentialForm(beta1,beta2,b);
	
	double sum_fx2j_1 = 0.0;
	double sum_fx2j_2 = 0.0;
	for (int j=2; j <= n/2;  j++)
		sum_fx2j_1 = sum_fx2j_1 + exponentialForm(beta1,beta2,a+h*(2*j-2));
	for (int j=1; j <= n/2;  j++)
		sum_fx2j_2 = sum_fx2j_2 + exponentialForm(beta1,beta2,a+h*(2*j-1));

	integral = h/3.0*(fa + 2.0*sum_fx2j_1 + 4.0*sum_fx2j_2 + fb);

	return integral/2.0/pi;
}


double
CorrelatedStandardNormal::getCDFowen(double b1, double b2, int popt)
{	
	// single integral form by Owen 1956
    static NormalRV uRV(1, 0.0, 1.0);
    double cdf = 0.0;
	double thresh = 0.99;
	double integral = 0.0;
	
	if ( rho > thresh ) {
		// create a refined domain near +1
		integral = SimpsonOwen(0,thresh,b1,b2);
		integral += SimpsonOwen(thresh,rho,b1,b2);
	} else if ( rho < -thresh ) {
		// create a refined domain near -1
		integral = SimpsonOwen(0,-thresh,b1,b2);
		integral += SimpsonOwen(-thresh,rho,b1,b2);
	} else {
		integral = SimpsonOwen(0,rho,b1,b2);
	}
	
	cdf = uRV.getCDFvalue(b1) * uRV.getCDFvalue(b2) + integral;
	if (popt == 1)
		opserr << " Owen = " << cdf;
		
	return cdf;
}


double
CorrelatedStandardNormal::getCDFsheppard(double b1, double b2, int popt)
{
	// single integral form by Sheppard 1900
    static NormalRV uRV(1, 0.0, 1.0);
    double cdf = 0.0;
    
	if (b1 >= 0 && b2 >= 0)
		cdf = 1.0 + SimpsonSheppard(rho,b1,b2) - uRV.getCDFvalue(-b1) - uRV.getCDFvalue(-b2);
	else if (b1 < 0 && b2 < 0)
		cdf = SimpsonSheppard(rho,-b1,-b2);
	else if (b1 < 0 && b2 >= 0)
		cdf = 1.0 - uRV.getCDFvalue(-b1) - SimpsonSheppard(-rho,-b1,b2);
	else
		cdf = 1.0 - uRV.getCDFvalue(-b2) - SimpsonSheppard(-rho,b1,-b2);
	
	if (popt == 1)
		opserr << " Sheppard = " << cdf;

	return cdf;
}


double
CorrelatedStandardNormal::getCDFadaptive(double b1, double b2, int popt)
{
	// adaptive quadrature from Quan
    static NormalRV uRV(1, 0.0, 1.0);
	double cdf = 0.0;
	double integral = 0.0;
	
	double fa = bivariatePDF(b1,b2,0);
	double fb = bivariatePDF(b1,b2,rho/2.0);
	double fc = bivariatePDF(b1,b2,rho);
	integral = getAdaptiveIntegralValue(1.0e-12, 0.0,rho, fa,fb,fc, b1,b2);
	cdf = uRV.getCDFvalue(b1) * uRV.getCDFvalue(b2) + integral;
	
	if (popt == 1)
		opserr << " Adaptive = " << cdf;
	
	return cdf;
}


double
CorrelatedStandardNormal::getCDF(double b1, double b2)
{
    static NormalRV uRV(1, 0.0, 1.0);
	double cdf = 0.0;
	int popt = 0;
	
	// treat several special cases before attempting to integrate numerically
	if ( 1.0 - rho <= DBL_EPSILON ) {
		// rho = 1
		if (b1 >= b2 && b1 >= 0)
			cdf = 2.0 - uRV.getCDFvalue(b1) - uRV.getCDFvalue(-b1) - uRV.getCDFvalue(-b2);
		else if (b2 > b1 && b2 >= 0)
			cdf = 2.0 - uRV.getCDFvalue(b2) - uRV.getCDFvalue(-b1) - uRV.getCDFvalue(-b2);
		else if (-b1 >= -b2)
			cdf = 1.0 - uRV.getCDFvalue(-b1);
		else
			cdf = 1.0 - uRV.getCDFvalue(-b2);
		
	} else if ( 1.0 + rho <= DBL_EPSILON ) {
		// rho = -1
		if (b1 + b2 >= 0)
			cdf = 1.0 - uRV.getCDFvalue(-b1) - uRV.getCDFvalue(-b2);
		else
			cdf = 0.0;
		
	} else {
		// single integral form by Owen 1956
		//cdf = getCDFowen(b1,b2,popt);
		
		// single integral form by Sheppard 1900
		cdf = getCDFsheppard(b1,b2,popt);
		
		// adaptive quadrature from Quan
		//cdf = getCDFadaptive(b1,b2,popt);
	}
	
	if (cdf > 1)
		cdf = 1.0;
	if (cdf < 0)
		cdf = 0.0;
	
	return cdf;
}


double 
CorrelatedStandardNormal::getAdaptiveIntegralValue(double tol, double lowerBound, double upperBound, 
	double fa, double fb, double fc, double beta1, double beta2)
{
	double b = (lowerBound+upperBound)/2.0;
	double fd = bivariatePDF(beta1,beta2, (lowerBound+b)/2.0);
	double fe = bivariatePDF(beta1,beta2, (upperBound+b)/2.0);

	double S1 = (upperBound-lowerBound)*(fa+4*fb+fc)/6.0;
	double S2 = (upperBound-lowerBound)*(fa+4*fd+2*fb+4*fe+fc)/12.0;
	if (fabs(S2-S1)/15.0 <= tol) 
		return (16.0*S2-S1)/15.0;
		
	else {
		double int1 = getAdaptiveIntegralValue( tol/2.0, lowerBound, b, fa, fd, fb, beta1, beta2);
		double int2 = getAdaptiveIntegralValue( tol/2.0, b, upperBound, fb, fe, fc, beta1, beta2);
		return int1+int2;
	}
}


void
CorrelatedStandardNormal::testCDF()
{
	setCorrelation(0.35);
	opserr << "rho = " << 0.35 << endln;
	opserr << getCDF(-1.56,1.55) << endln;
	opserr << getCDF(-2.05,-1.07) << endln;
	opserr << getCDF(3.05,1.69) << endln << endln;
	
	setCorrelation(-0.35);
	opserr << "rho = " << -0.35 << endln;
	opserr << getCDF(-1.56,1.55) << endln;
	opserr << getCDF(-2.05,-1.07) << endln;
	opserr << getCDF(3.05,1.69) << endln << endln;
	
	setCorrelation(-0.99);
	opserr << "rho = " << -0.99 << endln;
	opserr << getCDF(-1.56,1.55) << endln;
	opserr << getCDF(-2.05,1.07) << endln;
	opserr << getCDF(3.05,-1.69) << endln << endln;
	
	setCorrelation(0.99);
	opserr << "rho = " << 0.99 << endln;
	opserr << getCDF(-1.56,1.55) << endln;
	opserr << getCDF(-2.05,1.07) << endln;
	opserr << getCDF(3.05,-1.69) << endln << endln;
	
	setCorrelation(-0.999);
	opserr << "rho = " << -0.999 << endln;
	opserr << getCDF(-1.56,1.55) << endln;
	opserr << getCDF(-2.05,1.07) << endln;
	opserr << getCDF(3.05,-1.69) << endln << endln;
	
	setCorrelation(0.9999);
	opserr << "rho = " << 0.9999 << endln;
	opserr << getCDF(1.56,1.55) << endln;
	opserr << getCDF(-2.05,-1.07) << endln;
	opserr << getCDF(3.05,1.69) << endln;
	opserr << getCDF(0.65,-0.09) << endln << endln;
	
	setCorrelation(-0.999999);
	opserr << "rho = " << -0.999999 << endln;
	opserr << getCDF(-1.56,1.55) << endln;
	opserr << getCDF(-2.05,1.07) << endln;
	opserr << getCDF(3.05,-1.69) << endln << endln;
	
	setCorrelation(-1);
	opserr << "rho = " << -1 << endln;
	opserr << getCDF(-1.56,1.55) << endln;
	opserr << getCDF(-2.05,1.07) << endln;
	opserr << getCDF(3.05,-1.69) << endln << endln;
	
	setCorrelation(1);
	opserr << "rho = " << 1 << endln;
	opserr << getCDF(1.56,1.55) << endln;
	opserr << getCDF(-2.05,-1.07) << endln;
	opserr << getCDF(3.05,1.69) << endln;
	opserr << getCDF(0.65,-0.09) << endln << endln;

}
