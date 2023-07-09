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
// $Source: /usr/local/cvs/OpenSees/SRC/reliability/domain/distributions/ExponentialRV.cpp,v $


//
// Written by Terje Haukaas (haukaas@ce.berkeley.edu) 
//

#include <ExponentialRV.h>
#include <Vector.h>
#include <cmath>

ExponentialRV::ExponentialRV(int passedTag, double passedMean, double passedStdv)
:RandomVariable(passedTag, RANDOM_VARIABLE_exponential)
{
	int setp = setParameters(passedMean,passedStdv);
	if (setp < 0)
		opserr << "Error setting parameters in Exponential RV with tag " << this->getTag() << endln;
}


ExponentialRV::ExponentialRV(int passedTag, const Vector &passedParameters)
:RandomVariable(passedTag, RANDOM_VARIABLE_exponential)
{
	
	if (passedParameters.Size() != 1) {
		opserr << "Exponential RV requires 1 parameter, lambda, for RV with tag " <<
		this->getTag() << endln;
		
		// this will create terminal errors
		lambda = 0;
		
	} else {
		
		lambda = passedParameters(0);
		
	}
}


ExponentialRV::~ExponentialRV()
{
}


const char *
ExponentialRV::getType()
{
	return "EXPONENTIAL";
}


double 
ExponentialRV::getMean()
{
	return 1.0/lambda;
}


double 
ExponentialRV::getStdv()
{
	return 1.0/lambda;
}


const Vector &
ExponentialRV::getParameters(void) {
	static Vector temp(1);
	temp(0) = lambda;
	return temp;
}


int
ExponentialRV::setParameters(double mean, double stdv)
{
	lambda = 1.0/mean;
	return 0;
}


double
ExponentialRV::getPDFvalue(double rvValue)
{
	double result;
	if ( 0.0 <= rvValue ) {
		result = lambda * exp(-lambda * rvValue);
	}
	else {
		result = 0.0;
	}
	return result;
}


double
ExponentialRV::getCDFvalue(double rvValue)
{
	double result;
	if ( 0.0 <= rvValue ) {
		result = 1 - exp(-lambda*rvValue);
	}
	else {
		result = 0.0;
	}
	return result;
}


double
ExponentialRV::getInverseCDFvalue(double probValue)
{
	return -log(1.0-probValue)/lambda;
}


int 
ExponentialRV::getCDFparameterSensitivity(Vector &dFdP)
{
    // returns gradient of F(x) with respect to distribution parameters
    double rvValue = this->getCurrentValue();
    
    // dFdlambda
    dFdP(0) = rvValue * exp(-lambda*rvValue);
    
    return 0;
}


int
ExponentialRV::getParameterMeanSensitivity(Vector &dPdmu)
{
    // returns gradient of distribution parameters with respect to the mean
    double mu = getMean();
    
    // dlambdadmu
    dPdmu(0) = -1/mu/mu;
    
    return 0;
}


int
ExponentialRV::getParameterStdvSensitivity(Vector &dPdstdv)
{
    // returns gradient of distribution parameters with respect to the stdv
    double sig = getStdv();
    
    // dlambdadsig
    dPdstdv(0) = -1/sig/sig;
    
    return 0;
}


void
ExponentialRV::Print(OPS_Stream &s, int flag)
{
	s << "Exponential RV #" << this->getTag() << endln;
	s << "\tlambda = " << lambda << endln;
}

