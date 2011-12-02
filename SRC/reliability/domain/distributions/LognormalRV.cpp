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
// $Date: 2003-10-27 23:04:39 $
// $Source: /usr/local/cvs/OpenSees/SRC/reliability/domain/distributions/LognormalRV.cpp,v $


//
// Written by Terje Haukaas (haukaas@ce.berkeley.edu)
//

#include <LognormalRV.h>
#include <NormalRV.h>
#include <math.h>
#include <string.h>
#include <classTags.h>
#include <OPS_Globals.h>

LognormalRV::LognormalRV(int passedTag, 
		 double passedMean,
		 double passedStdv,
		 double passedStartValue)
:RandomVariable(passedTag, RANDOM_VARIABLE_lognormal)
{
	if (passedMean<0.0) {
		isPositive = false;
		passedMean = -passedMean;
	}
	else {
		isPositive = true;
	}

	tag = passedTag ;
	zeta = sqrt(   log(   1+(passedStdv/passedMean)*(passedStdv/passedMean)   )   );
	lambda = log(passedMean) - 0.5*zeta*zeta;
	startValue = passedStartValue;
}



LognormalRV::LognormalRV(int passedTag, 
		 double passedParameter1,
		 double passedParameter2,
		 double passedParameter3,
		 double passedParameter4,
		 double passedStartValue)
:RandomVariable(passedTag, RANDOM_VARIABLE_lognormal)
{
	if (passedParameter2<0.0) {
		isPositive = false;
		passedParameter2 = -passedParameter2;
	}
	else {
		isPositive = true;
	}

	tag = passedTag ;
	lambda = passedParameter1;
	zeta = passedParameter2;
	startValue = passedStartValue;
}



LognormalRV::LognormalRV(int passedTag, 
		 double passedMean,
		 double passedStdv)
:RandomVariable(passedTag, RANDOM_VARIABLE_lognormal)
{
	if (passedMean<0.0) {
		isPositive = false;
		passedMean = -passedMean;
	}
	else {
		isPositive = true;
	}

	tag = passedTag ;
	zeta = sqrt(   log(   1+(passedStdv/passedMean)*(passedStdv/passedMean)   )   );
	lambda = log(passedMean) - 0.5*zeta*zeta;
	startValue = getMean();
}
LognormalRV::LognormalRV(int passedTag, 
		 double passedParameter1,
		 double passedParameter2,
		 double passedParameter3,
		 double passedParameter4)
:RandomVariable(passedTag, RANDOM_VARIABLE_lognormal)
{
	if (passedParameter2<0.0) {
		isPositive = false;
		passedParameter2 = -passedParameter2;
	}
	else {
		isPositive = true;
	}

	tag = passedTag ;
	lambda = passedParameter1;
	zeta = passedParameter2;
	startValue = getMean();
}


LognormalRV::~LognormalRV()
{
}


void
LognormalRV::Print(OPS_Stream &s, int flag)
{
}


double
LognormalRV::getPDFvalue(double rvValue)
{
	if (!isPositive) {
		// The formal answer is: f(x) = f_pos(x+2|x|), but let's do it simple here
		rvValue = -rvValue;
	}

	double result;
	if ( 0.0 < rvValue ) {
		double pi = 3.14159265358979;
		result = 1/(sqrt(2*pi)*zeta*rvValue) * exp(-0.5* pow ( (log(rvValue)-lambda) / zeta, 2 )  );
	}
	else {
		result = 0.0;
	}
	return result;
}


double
LognormalRV::getCDFvalue(double rvValue)
{
	double result;

	static NormalRV aStandardNormalRV( 1, 0.0, 1.0, 0.0);
	
	if (isPositive) {
		if ( 0.0 < rvValue ) {
			result = aStandardNormalRV.getCDFvalue((log(rvValue)-lambda)/zeta);
		}
		else {
			result = 0.0;
		}
	}
	else {
		if ( rvValue < 0.0 ) {
			result = aStandardNormalRV.getCDFvalue((log(fabs(rvValue))-lambda)/zeta);
			result = 1.0-result;
		}
		else {
			result = 1.0;
		}
	}

	// Return result depending on type of random variable
	if (isPositive) {
		return result;
	}
	else {
		return 1-result;
	}



/*
	// First, flip it around if it's a negative lognormal
	if (!isPositive) {
		rvValue = -rvValue;
	}

	// Compute the ordinary CDF
	double result;
	if ( 0.0 < rvValue ) {
		RandomVariable *aStandardNormalRV;
		aStandardNormalRV= new NormalRV( 1, 0.0, 1.0, 0.0);
		result = aStandardNormalRV->getCDFvalue((log(rvValue)-lambda)/zeta);
		delete aStandardNormalRV;	
	}
	else {
		result = 0.0;
	}

	// Return result depending on type of random variable
	if (isPositive) {
		return result;
	}
	else {
		return 1-result;
	}
*/
}


double
LognormalRV::getInverseCDFvalue(double probValue)
{
	if ( probValue > 1.0 || probValue < 0.0) {
		opserr << "Invalid probability value input to inverse CDF function" << endln;
		return 0.0;
	}
	else {
		static NormalRV aStandardNormalRV( 1, 0.0, 1.0, 0.0);

		if (isPositive) {
			double inverseNormal = aStandardNormalRV.getInverseCDFvalue(probValue);
			return exp(inverseNormal*zeta + lambda);
		}
		else {
			double inverseNormal = aStandardNormalRV.getInverseCDFvalue(1.0-probValue);
			return (-exp(inverseNormal*zeta + lambda));
		}
	}
}


const char *
LognormalRV::getType()
{
	return "LOGNORMAL";
}


double 
LognormalRV::getMean()
{
	if (isPositive) {
		return exp(lambda+0.5*zeta*zeta);
	}
	else {
		return -exp(lambda+0.5*zeta*zeta);
	}
}



double 
LognormalRV::getStdv()
{
	return exp(lambda+0.5*zeta*zeta)*sqrt(exp(zeta*zeta)-1);
}


double 
LognormalRV::getStartValue()
{
	return startValue;
}


double LognormalRV::getParameter1()  {return lambda;}

double LognormalRV::getParameter2()  
{
	if (isPositive) {
		return zeta;
	}
	else {
		return -zeta;
	}
}

double LognormalRV::getParameter3()  {opserr<<"No such parameter in r.v. #"<<tag<<endln; return 0.0;}

double LognormalRV::getParameter4()  {opserr<<"No such parameter in r.v. #"<<tag<<endln; return 0.0;}
