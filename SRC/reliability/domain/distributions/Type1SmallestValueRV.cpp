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
                                                                        
// $Revision: 1.8 $
// $Date: 2007-02-17 21:27:23 $
// $Source: /usr/local/cvs/OpenSees/SRC/reliability/domain/distributions/Type1SmallestValueRV.cpp,v $


//
// Written by Terje Haukaas (haukaas@ce.berkeley.edu)
//

#include <Type1SmallestValueRV.h>
#include <math.h>
#include <string.h>
#include <classTags.h>
#include <OPS_Globals.h>

Type1SmallestValueRV::Type1SmallestValueRV(int passedTag, 
		 double passedMean,
		 double passedStdv,
		 double passedStartValue)
:RandomVariable(passedTag, RANDOM_VARIABLE_type1smallestvalue, passedStartValue)
{
	double gamma = 0.5772156649;
	double pi = 3.14159265358979;
	u = passedMean + gamma * passedStdv * sqrt(6.0) / pi;
	alpha = pi / (passedStdv*sqrt(6.0));
}
Type1SmallestValueRV::Type1SmallestValueRV(int passedTag, 
		 double passedParameter1,
		 double passedParameter2,
		 double passedParameter3,
		 double passedParameter4,
		 double passedStartValue)
:RandomVariable(passedTag, RANDOM_VARIABLE_type1smallestvalue, passedStartValue)
{
	u = passedParameter1;
	alpha = passedParameter2;
}
Type1SmallestValueRV::Type1SmallestValueRV(int passedTag, 
		 double passedMean,
		 double passedStdv)
:RandomVariable(passedTag, RANDOM_VARIABLE_type1smallestvalue, passedMean)
{
	double gamma = 0.5772156649;
	double pi = 3.14159265358979;
	u = passedMean + gamma * passedStdv * sqrt(6.0) / pi;
	alpha = pi / (passedStdv*sqrt(6.0));
}
Type1SmallestValueRV::Type1SmallestValueRV(int passedTag, 
		 double passedParameter1,
		 double passedParameter2,
		 double passedParameter3,
		 double passedParameter4)
:RandomVariable(passedTag, RANDOM_VARIABLE_type1smallestvalue)
{
	u = passedParameter1;
	alpha = passedParameter2;
	this->setStartValue(getMean());
}


Type1SmallestValueRV::~Type1SmallestValueRV()
{
}


void
Type1SmallestValueRV::Print(OPS_Stream &s, int flag)
{
}


double
Type1SmallestValueRV::getPDFvalue(double rvValue)
{
	return alpha*exp(alpha*(rvValue-u)-exp(alpha*(rvValue-u)));
}


double
Type1SmallestValueRV::getCDFvalue(double rvValue)
{
	return 1-exp(-exp(alpha*(rvValue-u)));
}


double
Type1SmallestValueRV::getInverseCDFvalue(double probValue)
{
	return (alpha*u + log(-log(1.0-probValue))) / alpha;
}


const char *
Type1SmallestValueRV::getType()
{
	return "TYPE1SMALLESTVALUE";
}


double 
Type1SmallestValueRV::getMean()
{
	double gamma = 0.5772156649;
	return u - gamma/alpha;
}



double 
Type1SmallestValueRV::getStdv()
{
	double pi = 3.14159265358979;
	return pi/(sqrt(6.0)*alpha);
}




double
Type1SmallestValueRV::getParameter1()  
{
  return u;
}

double
Type1SmallestValueRV::getParameter2()  
{
  return alpha;
}
