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
// $Date: 2004-08-27 17:51:50 $
// $Source: /usr/local/cvs/OpenSees/SRC/reliability/domain/distributions/GumbelRV.cpp,v $


//
// Written by Terje Haukaas (haukaas@ce.berkeley.edu)
//

#include <GumbelRV.h>
#include <math.h>
#include <string.h>
#include <classTags.h>
#include <OPS_Globals.h>

GumbelRV::GumbelRV(int passedTag, 
		 double passedMean,
		 double passedStdv,
		 double passedStartValue)
:RandomVariable(passedTag, RANDOM_VARIABLE_gumbel)
{
	tag = passedTag ;
	double pi = 3.14159265358979;
	double gamma = 0.5772156649;
	u = passedMean - (gamma*sqrt(6.0)*passedStdv)/pi;
	alpha = pi / ( sqrt(6.0) * passedStdv );
	startValue = passedStartValue;
}
GumbelRV::GumbelRV(int passedTag, 
		 double passedParameter1,
		 double passedParameter2,
		 double passedParameter3,
		 double passedParameter4,
		 double passedStartValue)
:RandomVariable(passedTag, RANDOM_VARIABLE_gumbel)
{
	tag = passedTag ;
	u = passedParameter1;
	alpha = passedParameter2;
	startValue = passedStartValue;
}
GumbelRV::GumbelRV(int passedTag, 
		 double passedMean,
		 double passedStdv)
:RandomVariable(passedTag, RANDOM_VARIABLE_gumbel)
{
	tag = passedTag ;
	double pi = 3.14159265358979;
	double gamma = 0.5772156649;
	u = passedMean - (gamma*sqrt(6.0)*passedStdv)/pi;
	alpha = pi / ( sqrt(6.0) * passedStdv );
	startValue = getMean();
}
GumbelRV::GumbelRV(int passedTag, 
		 double passedParameter1,
		 double passedParameter2,
		 double passedParameter3,
		 double passedParameter4)
:RandomVariable(passedTag, RANDOM_VARIABLE_gumbel)
{
	tag = passedTag ;
	u = passedParameter1;
	alpha = passedParameter2;
	startValue = getMean();
}


GumbelRV::~GumbelRV()
{
}


void
GumbelRV::Print(OPS_Stream &s, int flag)
{
}


double
GumbelRV::getPDFvalue(double rvValue)
{
	return alpha * exp(-alpha*(rvValue-u)-exp(-alpha*(rvValue-u)));
}


double
GumbelRV::getCDFvalue(double rvValue)
{
	return exp(-exp(-alpha*(rvValue-u)));
}


double
GumbelRV::getInverseCDFvalue(double probValue)
{
	return (alpha*u - log(-log(probValue))) / alpha;
}


const char *
GumbelRV::getType()
{
	return "GUMBEL";
}


double 
GumbelRV::getMean()
{
	double gamma = 0.5772156649;
	return u+gamma/alpha;
}



double 
GumbelRV::getStdv()
{
	double pi = 3.14159265358979;
	return pi/(sqrt(6.0)*alpha);
}


double 
GumbelRV::getStartValue()
{
	return startValue;
}

double GumbelRV::getParameter1()  {return u;}
double GumbelRV::getParameter2()  {return alpha;}
double GumbelRV::getParameter3()  {opserr<<"No such parameter in r.v. #"<<tag<<endln; return 0.0;}
double GumbelRV::getParameter4()  {opserr<<"No such parameter in r.v. #"<<tag<<endln; return 0.0;}
