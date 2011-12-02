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
// $Date: 2003-03-04 00:44:35 $
// $Source: /usr/local/cvs/OpenSees/SRC/reliability/domain/distributions/ShiftedExponentialRV.cpp,v $


//
// Written by Terje Haukaas (haukaas@ce.berkeley.edu) 
//

#include <ShiftedExponentialRV.h>
#include <math.h>
#include <string.h>
#include <classTags.h>
#include <OPS_Globals.h>

ShiftedExponentialRV::ShiftedExponentialRV(int passedTag, 
		 double passedMean,
		 double passedStdv,
		 double passedStartValue)
:RandomVariable(passedTag, RANDOM_VARIABLE_shiftedexponential)
{
	tag = passedTag ;
	lambda = 1/passedStdv;
	x0 = passedMean - passedStdv;
	startValue = passedStartValue;
}
ShiftedExponentialRV::ShiftedExponentialRV(int passedTag, 
		 double passedParameter1,
		 double passedParameter2,
		 double passedParameter3,
		 double passedParameter4,
		 double passedStartValue)
:RandomVariable(passedTag, RANDOM_VARIABLE_shiftedexponential)
{
	tag = passedTag ;
	lambda = passedParameter1;
	x0 = passedParameter2;
	startValue = passedStartValue;
}
ShiftedExponentialRV::ShiftedExponentialRV(int passedTag, 
		 double passedMean,
		 double passedStdv)
:RandomVariable(passedTag, RANDOM_VARIABLE_shiftedexponential)
{
	tag = passedTag ;
	lambda = 1/passedStdv;
	x0 = passedMean - passedStdv;
	startValue = getMean();
}
ShiftedExponentialRV::ShiftedExponentialRV(int passedTag, 
		 double passedParameter1,
		 double passedParameter2,
		 double passedParameter3,
		 double passedParameter4)
:RandomVariable(passedTag, RANDOM_VARIABLE_shiftedexponential)
{
	tag = passedTag ;
	lambda = passedParameter1;
	x0 = passedParameter2;
	startValue = getMean();
}


ShiftedExponentialRV::~ShiftedExponentialRV()
{
}


void
ShiftedExponentialRV::Print(OPS_Stream &s, int flag)
{
}


double
ShiftedExponentialRV::getPDFvalue(double rvValue)
{
	double result;
	if ( x0 < rvValue ) {
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
	if ( x0 < rvValue ) {
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


double 
ShiftedExponentialRV::getStartValue()
{
	return startValue;
}

double ShiftedExponentialRV::getParameter1()  {return lambda;}
double ShiftedExponentialRV::getParameter2()  {return x0;}
double ShiftedExponentialRV::getParameter3()  {opserr<<"No such parameter in r.v. #"<<tag<<endln; return 0.0;}
double ShiftedExponentialRV::getParameter4()  {opserr<<"No such parameter in r.v. #"<<tag<<endln; return 0.0;}
