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
#include <math.h>
#include <string.h>
#include <classTags.h>
#include <OPS_Globals.h>

ExponentialRV::ExponentialRV(int passedTag, 
		 double passedMean,
		 double passedStdv,
		 double passedStartValue)
:RandomVariable(passedTag, RANDOM_VARIABLE_exponential, passedStartValue)
{
	// Note: this constructor is void.
	opserr << "WARNING: This type of random variable is not uniquely defined by mean and stdv." << endln;
}
ExponentialRV::ExponentialRV(int passedTag, 
		 double passedParameter1,
		 double passedParameter2,
		 double passedParameter3,
		 double passedParameter4,
		 double passedStartValue)
:RandomVariable(passedTag, RANDOM_VARIABLE_exponential, passedStartValue)
{
	lambda = passedParameter1;
}
ExponentialRV::ExponentialRV(int passedTag, 
		 double passedMean,
		 double passedStdv)
:RandomVariable(passedTag, RANDOM_VARIABLE_exponential, passedMean)
{
	// Note: this constructor is void.
	opserr << "WARNING: This type of random variable is not uniquely defined by mean and stdv." << endln;
}
ExponentialRV::ExponentialRV(int passedTag, 
		 double passedParameter1,
		 double passedParameter2,
		 double passedParameter3,
		 double passedParameter4)
:RandomVariable(passedTag, RANDOM_VARIABLE_exponential)
{
	lambda = passedParameter1;
	this->setStartValue(getMean());
}


ExponentialRV::~ExponentialRV()
{
}


void
ExponentialRV::Print(OPS_Stream &s, int flag)
{
}


double
ExponentialRV::getPDFvalue(double rvValue)
{
	double result;
	if ( 0.0 < rvValue ) {
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
	if ( 0.0 < rvValue ) {
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


const char *
ExponentialRV::getType()
{
	return "EXPONENTIAL";
}


double 
ExponentialRV::getMean()
{
	return 1/lambda;
}



double 
ExponentialRV::getStdv()
{
	return 1/lambda;
}


double
ExponentialRV::getParameter1()
{
  return lambda;
}
