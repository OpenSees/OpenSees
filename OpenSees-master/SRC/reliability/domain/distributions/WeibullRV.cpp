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
// $Source: /usr/local/cvs/OpenSees/SRC/reliability/domain/distributions/WeibullRV.cpp,v $


//
// Written by Terje Haukaas (haukaas@ce.berkeley.edu)
//

#include <WeibullRV.h>
#include <Vector.h>
#include <cmath>

WeibullRV::WeibullRV(int passedTag, 
		     double passedMean, double passedStdv)
:RandomVariable(passedTag, RANDOM_VARIABLE_weibull)
{
	int setp = setParameters(passedMean,passedStdv);
	if (setp < 0)
		opserr << "Error setting parameters in Weibull RV with tag " << this->getTag() << endln;
}


WeibullRV::WeibullRV(int passedTag,
		     const Vector &passedParameters)
:RandomVariable(passedTag, RANDOM_VARIABLE_weibull)
{
	
	if (passedParameters.Size() != 2) {
		opserr << "Weibull RV requires 2 parameters, u and k, for RV with tag " <<
		this->getTag() << endln;
		
		// this will create terminal errors
		u = 0;
		k = 0;
		
	} else {
		
		u = passedParameters(0);
		k = passedParameters(1);
		
	}
}


WeibullRV::~WeibullRV()
{
}


const char *
WeibullRV::getType()
{
	return "WEIBULL";
}


double 
WeibullRV::getMean()
{
	double result = u * gammaFunction(1.0+1.0/k);
	return result;
}


double 
WeibullRV::getStdv()
{
	double a = gammaFunction(1.0+2.0/k);
	double b = gammaFunction(1.0+1.0/k);
	double result = u*sqrt(a-b*b);
	return result;
}


const Vector &
WeibullRV::getParameters(void) {
	static Vector temp(2);
	temp(0) = u;
	temp(1) = k;
	return temp;
}


int
WeibullRV::setParameters(double mean, double stdv)
{
	double cov = stdv/mean;
	double c = 1+cov*cov;
	double k_prev = log(1/(c-1));
	double del = 1.0;
    double del_prev = 10.0;
	
	// now use Newtons method with nice f/f' function
	int ncount = 1;
	int nmax = 100;
    double tolk = 1.0e-8;
	
	while (del > tolk && ncount <= nmax) {
		k = k_prev - (c*pow(gammaFunction(1/k_prev),2)-2*k_prev*gammaFunction(2/k_prev)) / 
            ( 2*gammaFunction(2/k_prev+1)*(-harmonicNumber(1/k_prev)+harmonicNumber(2/k_prev)) );
		del = fabs(k-k_prev);
        
        if (del > del_prev) {
            // be careful with non-convergent cases, exit to bisection method
            break;
        }
        
        k_prev = k;
        del_prev = del;
		ncount++;
	}
    
    // check whether we need to use bisection method to supplement
    if (del > tolk) {
        double left = 0.05;
        double right = 1000;
        int change = 0;
        double fa = 0;
        double fc = 0;
        ncount = 1;
        
        while (change == 0 && ncount <= nmax) {
            k = (left+right)/2.0;
            fc = -c + gammaFunction(1+2/k)/gammaFunction(1+1/k)/gammaFunction(1+1/k);
            fa = -c + gammaFunction(1+2/left)/gammaFunction(1+1/left)/gammaFunction(1+1/left);
            if (fabs(fc) < tolk || (right-left)/2.0 < tolk) {
                change = 1;
                del = tolk;
            }
            
            if ( fa>0 && fc>0 )
                left = k;
            else if ( fa<0 && fc<0 )
                left = k;
            else
                right = k;
            
            ncount++;
        }
    }
    
    // unlikely but possible that k was not found using either method above
    if (del > tolk) {
		opserr << "Warning: Weibull distribution did not converge during setParameters()" << endln;
        return -1;
    }
    
	u = mean/gammaFunction(1+1/k);
	
	return 0;
}


double
WeibullRV::getPDFvalue(double rvValue)
{
	double result;
	if ( 0.0 < rvValue) {
		result = k/u * pow(rvValue/u,k-1) * exp(-pow(rvValue/u,k));
	}
	else {
		result = 0.0;
	}
	return result;
}


double
WeibullRV::getCDFvalue(double rvValue)
{
	double result;
	if ( 0.0 < rvValue) {
		result = 1 - exp(-pow(rvValue/u,k));
	}
	else {
		result = 0.0;
	}
	return result;
}


double
WeibullRV::getInverseCDFvalue(double probValue)
{
	return u * pow(-log(1-probValue),1/k);
}


int 
WeibullRV::getCDFparameterSensitivity(Vector &dFdP)
{
    // returns gradient of F(x) with respect to distribution parameters
    double rvValue = this->getCurrentValue();
    
    // dFdu
    dFdP(0) = -rvValue/u * getPDFvalue(rvValue);
    
    // dFdk
    dFdP(1) = rvValue/k*log(rvValue/u) * getPDFvalue(rvValue);
    
    return 0;
}


int
WeibullRV::getParameterMeanSensitivity(Vector &dPdmu)
{
    // returns gradient of distribution parameters with respect to the mean
    double mu = getMean();
    double sig = getStdv();
    
    // dudmu
    // no closed form expression for this one, use finite differences
    double u_old = u;
    double k_old = k;
    double dh = mu/1000.0;
    setParameters(mu+dh,sig);
    dPdmu(0) = (u-u_old)/dh;

    // dkdmu
    // use finite differences from above, or check with closed form
    dPdmu(1) = (k-k_old)/dh;
    
    // reset original parameters
    k = k_old;
    u = u_old;
   
    return 0;
}


int
WeibullRV::getParameterStdvSensitivity(Vector &dPdstdv)
{
    // returns gradient of distribution parameters with respect to the stdv
    double mu = getMean();
    double sig = getStdv();
    
    // dudsig
    // no closed form expression for this one, use finite differences
    double u_old = u;
    double k_old = k;
    double dh = sig/1000.0;
    setParameters(mu,sig+dh);
    dPdstdv(0) = (u-u_old)/dh;
    
    // dkdsig
    // use finite differences from above, or check with closed form
    dPdstdv(1) = (k-k_old)/dh;
    
    // reset original parameters
    k = k_old;
    u = u_old;
    
    return 0;
}


void
WeibullRV::Print(OPS_Stream &s, int flag)
{
	s << "Weibull RV #" << this->getTag() << endln;
	s << "\tu = " << u << endln;
	s << "\tk = " << k << endln;
}

