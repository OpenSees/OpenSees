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
// $Date: 2003-03-04 00:44:34 $
// $Source: /usr/local/cvs/OpenSees/SRC/reliability/domain/distributions/LaplaceRV.cpp,v $


//
// Written by Terje Haukaas (haukaas@ce.berkeley.edu) 
//

#include <LaplaceRV.h>
#include <math.h>
#include <string.h>
#include <classTags.h>
#include <OPS_Globals.h>

LaplaceRV::LaplaceRV(int passedTag, 
		 double passedMean,
		 double passedStdv,
		 double passedStartValue)
:RandomVariable(passedTag, RANDOM_VARIABLE_laplace)
{
	tag = passedTag ;
	alpha = passedMean;
	beta = sqrt(2.0)/passedStdv;
	startValue = passedStartValue;
}
LaplaceRV::LaplaceRV(int passedTag, 
		 double passedParameter1,
		 double passedParameter2,
		 double passedParameter3,
		 double passedParameter4,
		 double passedStartValue)
:RandomVariable(passedTag, RANDOM_VARIABLE_laplace)
{
	tag = passedTag ;
	alpha = passedParameter1;
	beta = passedParameter2;
	startValue = passedStartValue;
}
LaplaceRV::LaplaceRV(int passedTag, 
		 double passedMean,
		 double passedStdv)
:RandomVariable(passedTag, RANDOM_VARIABLE_laplace)
{
	tag = passedTag ;
	alpha = passedMean;
	beta = sqrt(2.0)/passedStdv;
	startValue = getMean();
}
LaplaceRV::LaplaceRV(int passedTag, 
		 double passedParameter1,
		 double passedParameter2,
		 double passedParameter3,
		 double passedParameter4)
:RandomVariable(passedTag, RANDOM_VARIABLE_laplace)
{
	tag = passedTag ;
	alpha = passedParameter1;
	beta = passedParameter2;
	startValue = getMean();
}


LaplaceRV::~LaplaceRV()
{
}


void
LaplaceRV::Print(OPS_Stream &s, int flag)
{
}


double
LaplaceRV::getPDFvalue(double rvValue)
{
	return 0.5*beta*exp(-beta*fabs(rvValue-alpha));
}


double
LaplaceRV::getCDFvalue(double rvValue)
{
	double result;
	if (rvValue < alpha)  {
		result = 0.5*exp(-beta*fabs(rvValue-alpha));
	}
	else  {
		result = 1 - 0.5*exp(-beta*fabs(rvValue-alpha));
	}
	return result;
}


double
LaplaceRV::getInverseCDFvalue(double rvValue)
{
	return 0.0;
}


const char *
LaplaceRV::getType()
{
	return "LAPLACE";
}


double 
LaplaceRV::getMean()
{
	return alpha;
}



double 
LaplaceRV::getStdv()
{
	return sqrt(2.0)/beta;
}


double 
LaplaceRV::getStartValue()
{
	return startValue;
}

double LaplaceRV::getParameter1()  {return alpha;}
double LaplaceRV::getParameter2()  {return beta;}
double LaplaceRV::getParameter3()  {opserr<<"No such parameter in r.v. #"<<tag<<endln; return 0.0;}
double LaplaceRV::getParameter4()  {opserr<<"No such parameter in r.v. #"<<tag<<endln; return 0.0;}
