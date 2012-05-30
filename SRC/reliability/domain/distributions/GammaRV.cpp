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

	// KRM -- Note to writer of above: Newton does not converge for all cases, particularly those
	// that look like CDFs when the start point is not close to the zero point.  
	// Bootstrap with estimate from equivalent normal distribution
	static NormalRV normRV( 1, getMean(), getStdv());
    double x_old = normRV.getInverseCDFvalue(probValue);
	
	double tol = 1.0e-6;
	double x_new = x_old;
	double dx = x_old;
    int step = 1;
    int nmax = 50;
    
    if (probValue > 1.0e-3) {
        while ( (step < nmax) && (fabs(dx/x_old) > tol) ) {
            // compute dx
            dx = (getCDFvalue(x_new) - probValue) / getPDFvalue(x_new);

            // Take a Newton step
            x_new = x_new - dx;
            step++;
        }
    }
    
    // use bisection for small probabilities (x > 0 always for gamma). Here we (arbitrarily)
    // limit the left extreme x to be 0.1
    else {
        double left = 0.1;
        double right = 2.5*x_old;
        int change = 0;
        double fa = 0;
        double fc = 0;
        
        while (change == 0 && step <= nmax) {
            x_new = (left+right)/2.0;
            fc = getCDFvalue(x_new) - probValue;
            fa = getCDFvalue(left) - probValue;
            if (fabs(fc) < tol || (right-left)/2.0 < tol)
                change = 1;
            
            if ( fa>0 && fc>0 )
                left = x_new;
            else if ( fa<0 && fc<0 )
                left = x_new;
            else
                right = x_new;
            
            step++;
        }
	}
	
    if ( step >= 50 ) {
        opserr << "WARNING: GammaRV did not converge to find inverse CDF!" << endln;
        // decide later what to do with this terminal case
    }
    
	return x_new;
}


int 
GammaRV::getCDFparameterSensitivity(Vector &dFdP)
{
    // returns gradient of F(x) with respect to distribution parameters
    double rvValue = this->getCurrentValue();
    double cdfValue = getCDFvalue(rvValue);
    
    // dFdk
    // no closed form expression for this one, use finite differences
    double k_old = k;
    double dh = k/1000.0;
    k += dh;
    dFdP(0) = ( getCDFvalue(rvValue) - cdfValue )/dh;
    k = k_old;
    
    // dFdlambda
    dFdP(1) = rvValue/lambda * getPDFvalue(rvValue);
    
    return 0;
}


int
GammaRV::getParameterMeanSensitivity(Vector &dPdmu)
{
    // returns gradient of distribution parameters with respect to the mean
    double mu = getMean();
    double sig = getStdv();
    
    // dkdmu
    dPdmu(0) = 2*mu/sig/sig;
    
    // dlambdadmu
    dPdmu(1) = 1/sig/sig;
    
    return 0;
}


int
GammaRV::getParameterStdvSensitivity(Vector &dPdstdv)
{
    // returns gradient of distribution parameters with respect to the stdv
    double mu = getMean();
    double sig = getStdv();
    
    // dkdsig
    dPdstdv(0) = -2*mu*mu/sig/sig/sig;
    
    // dlambdadsig
    dPdstdv(1) = -2*mu/sig/sig/sig;
    
    return 0;
}


void
GammaRV::Print(OPS_Stream &s, int flag)
{
	s << "Gamma RV #" << this->getTag() << endln;
	s << "\tk = " << k << endln;
	s << "\tlambda = " << lambda << endln;
}
