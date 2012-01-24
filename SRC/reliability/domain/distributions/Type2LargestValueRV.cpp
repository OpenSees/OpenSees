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
// $Source: /usr/local/cvs/OpenSees/SRC/reliability/domain/distributions/Type2LargestValueRV.cpp,v $


//
// Written by Terje Haukaas (haukaas@ce.berkeley.edu)
//

#include <Type2LargestValueRV.h>
#include <Vector.h>
#include <cmath>

Type2LargestValueRV::Type2LargestValueRV(int passedTag, 
					 double passedMean, double passedStdv)
:RandomVariable(passedTag, RANDOM_VARIABLE_type2largestvalue)
{
	int setp = setParameters(passedMean,passedStdv);
	if (setp < 0)
		opserr << "Error setting parameters in Type2LargestValue RV with tag " << this->getTag() << endln;
}


Type2LargestValueRV::Type2LargestValueRV(int passedTag,
					 const Vector &passedParameters)
:RandomVariable(passedTag, RANDOM_VARIABLE_type2largestvalue)
{
	
	if (passedParameters.Size() != 2) {
		opserr << "Type2LargestValue RV requires 2 parameters, u and k, for RV with tag " <<
		this->getTag() << endln;
		
		// this will create terminal errors
		u = 0;
		k = 0;
		
	} else {
		
		u = passedParameters(0);
		k = passedParameters(1);
	}
}


Type2LargestValueRV::~Type2LargestValueRV()
{
}


const char *
Type2LargestValueRV::getType()
{
	return "TYPE2LARGESTVALUE";
}


double 
Type2LargestValueRV::getMean()
{
	double result;
	result = u * gammaFunction(1.0-1.0/k);
	return result;
}


double 
Type2LargestValueRV::getStdv()
{
	double a = gammaFunction(1.0-2.0/k);
	double b = gammaFunction(1.0-1.0/k);
	double result = u*sqrt(a-b*b);
	return result;
}


const Vector &
Type2LargestValueRV::getParameters(void) {
	static Vector temp(2);
	temp(0) = u;
	temp(1) = k;
	return temp;
}


int
Type2LargestValueRV::setParameters(double mean, double stdv)
{
	double cov = stdv/mean;
	double xk = 1.5;
	function131(xk,cov,mean);
	
	return 0;
}


double
Type2LargestValueRV::getPDFvalue(double rvValue)
{
	double result;
	if ( 0.0 < rvValue ) {
		result = k/u * pow(u/rvValue,k+1) * exp(-pow(u/rvValue,k));
	}
	else {
		result = 0.0;
	}
	return result;
}


double
Type2LargestValueRV::getCDFvalue(double rvValue)
{
	double result;
	if ( 0.0 < rvValue ) {
		result = exp(-pow(u/rvValue,k));
	}
	else {
		result = 0.0;
	}
	return result;
}


double
Type2LargestValueRV::getInverseCDFvalue(double probValue)
{
	return u * pow((-log(probValue)),(-1.0/k));
}


void
Type2LargestValueRV::Print(OPS_Stream &s, int flag)
{
	s << "Type2LargestValue RV #" << this->getTag() << endln;
	s << "\tu = " << u << endln;
	s << "\tk = " << k << endln;
}


void
Type2LargestValueRV::function131(double xk, double cov, double mean)
{
	xk = xk + 1.0;
	double x1 = 1.0 - 1.0/xk;
	double x2 = 1.0 - 2.0/xk;
	double gm1 = gammaFunction(x1);
	double gm2 = gammaFunction(x2);
	double vy = sqrt(gm2/gm1/gm1 - 1.0);
	if (cov-vy < 0.0) 
		function131(xk,cov,mean);
	else if (cov-vy == 0.0) 
		function134(xk,cov,mean);
	else 
		function132(xk,cov,mean);
}


void
Type2LargestValueRV::function132(double xk, double cov, double mean)
{
	xk = xk - 0.1;
	double x1 = 1.0 - 1.0/xk;
	double x2 = 1.0 - 2.0/xk;
	double gm1 = gammaFunction(x1);
	double gm2 = gammaFunction(x2);
	double vy = sqrt(gm2/gm1/gm1 - 1.0);
	if (cov-vy < 0.0) 
		function133(xk,cov,mean);
	else if (cov-vy == 0.0) 
		function134(xk,cov,mean);
	else 
		function132(xk,cov,mean);
}


void
Type2LargestValueRV::function133(double xk, double cov, double mean)
{
	xk = xk + 0.01;
	double x1 = 1.0 - 1.0/xk;
	double x2 = 1.0 - 2.0/xk;
	double gm1 = gammaFunction(x1);
	double gm2 = gammaFunction(x2);
	double vy = sqrt(gm2/gm1/gm1 - 1.0);
	if (cov-vy < 0.0) 
		function133(xk,cov,mean);
	else if (cov-vy == 0.0) 
		function134(xk,cov,mean);
	else 
		function134(xk,cov,mean);
}


void
Type2LargestValueRV::function134(double xk, double gm1, double mean)
{
	u = mean/gm1;
	k = xk;
}
