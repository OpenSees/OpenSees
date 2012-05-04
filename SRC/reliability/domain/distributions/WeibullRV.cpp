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
	
	// now use Newtons method with nice f/f' function
	int ncount = 1;
	int nmax = 100;
	
	while (del > 1.0e-8 && ncount <= nmax) {
		k = k_prev - (c*pow(gammaFunction(1/k_prev),2)-2*k_prev*gammaFunction(2/k_prev)) / 
            ( 2*gammaFunction(2/k_prev+1)*(-harmonicNumber(1/k_prev)+harmonicNumber(2/k_prev)) );
		del = fabs(k-k_prev);
		k_prev = k;
		ncount++;
	}
	if (ncount >= nmax) {
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
	return u * pow((-log(1-probValue)),(1/k));
}


void
WeibullRV::Print(OPS_Stream &s, int flag)
{
	s << "Weibull RV #" << this->getTag() << endln;
	s << "\tu = " << u << endln;
	s << "\tk = " << k << endln;
}

