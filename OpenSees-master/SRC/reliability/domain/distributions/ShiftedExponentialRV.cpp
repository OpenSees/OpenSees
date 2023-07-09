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
// $Source: /usr/local/cvs/OpenSees/SRC/reliability/domain/distributions/ShiftedExponentialRV.cpp,v $


//
// Written by Terje Haukaas (haukaas@ce.berkeley.edu) 
//

#include <ShiftedExponentialRV.h>
#include <Vector.h>
#include <cmath>

ShiftedExponentialRV::ShiftedExponentialRV(int passedTag, 
					   double passedMean, double passedStdv)
:RandomVariable(passedTag, RANDOM_VARIABLE_shiftedexponential)
{
	int setp = setParameters(passedMean,passedStdv);
	if (setp < 0)
		opserr << "Error setting parameters in ShiftedExponential RV with tag " << this->getTag() << endln;
}


ShiftedExponentialRV::ShiftedExponentialRV(int passedTag,
					   const Vector &passedParameters)
:RandomVariable(passedTag, RANDOM_VARIABLE_shiftedexponential)
{
    
	if (passedParameters.Size() != 2) {
		opserr << "ShiftedExponential RV requires 2 parameters, lambda and x0, for RV with tag " <<
		this->getTag() << endln;
		
		// this will create terminal errors
		lambda = 0;
		x0 = 0;
		
	} else {
		
		lambda = passedParameters(0);
		x0 = passedParameters(1);
		
	}
}


ShiftedExponentialRV::~ShiftedExponentialRV()
{
}


const char *
ShiftedExponentialRV::getType()
{
	return "SHIFTEDEXPONENTIAL";
}


double 
ShiftedExponentialRV::getMean()
{
	return x0 + 1/lambda;
}


double 
ShiftedExponentialRV::getStdv()
{
	return 1/lambda;
}


const Vector &
ShiftedExponentialRV::getParameters(void) {
	static Vector temp(2);
	temp(0) = lambda;
	temp(1) = x0;
	return temp;
}


int
ShiftedExponentialRV::setParameters(double mean, double stdv)
{
	lambda = 1/stdv;
	x0 = mean - stdv;
	
	return 0;
}


double
ShiftedExponentialRV::getPDFvalue(double rvValue)
{
	double result;
	if ( x0 <= rvValue ) {
		result = lambda * exp(-lambda * (rvValue-x0));
	}
	else {
		result = 0.0;
	}
	return result;
}


double
ShiftedExponentialRV::getCDFvalue(double rvValue)
{
	double result;
	if ( x0 <= rvValue ) {
		result = 1 - exp(-lambda*(rvValue-x0));
	}
	else {
		result = 0.0;
	}
	return result;
}


double
ShiftedExponentialRV::getInverseCDFvalue(double probValue)
{
	return x0 - (log(1.0-probValue)) / lambda;
}


int 
ShiftedExponentialRV::getCDFparameterSensitivity(Vector &dFdP)
{
    // returns gradient of F(x) with respect to distribution parameters
    double rvValue = this->getCurrentValue();
    
    // dFdlambda
    dFdP(0) = (rvValue-x0)/lambda * getPDFvalue(rvValue);
    
    // dFdx0
    dFdP(1) = -1 * getPDFvalue(rvValue);
    
    return 0;
}


int
ShiftedExponentialRV::getParameterMeanSensitivity(Vector &dPdmu)
{
    // returns gradient of distribution parameters with respect to the mean
    
    // dlambdadmu
    dPdmu(0) = 0;
    
    // dx0dmu
    dPdmu(1) = 1;
    
    return 0;
}


int
ShiftedExponentialRV::getParameterStdvSensitivity(Vector &dPdstdv)
{
    // returns gradient of distribution parameters with respect to the stdv
    double sig = getStdv();
    
    // dlambdadsig
    dPdstdv(0) = -1/sig/sig;
    
    // dx0dsig
    dPdstdv(1) = -1;
    
    return 0;
}


void
ShiftedExponentialRV::Print(OPS_Stream &s, int flag)
{
	s << "ShiftedExponential RV #" << this->getTag() << endln;
	s << "\tlambda = " << lambda << endln;
	s << "\tx0 = " << x0 << endln;
}
