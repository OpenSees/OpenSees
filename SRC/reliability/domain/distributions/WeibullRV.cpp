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
// $Date: 2003-03-04 00:44:37 $
// $Source: /usr/local/cvs/OpenSees/SRC/reliability/domain/distributions/WeibullRV.cpp,v $


//
// Written by Terje Haukaas (haukaas@ce.berkeley.edu)
//

#include <WeibullRV.h>
#include <math.h>
#include <string.h>
#include <GammaRV.h>
#include <classTags.h>
#include <OPS_Globals.h>

WeibullRV::WeibullRV(int passedTag, 
		 double passedMean,
		 double passedStdv,
		 double passedStartValue)
:RandomVariable(passedTag, RANDOM_VARIABLE_weibull)
{
	tag = passedTag;
	setParameters(passedMean,passedStdv);
	startValue = passedStartValue;
}
WeibullRV::WeibullRV(int passedTag, 
		 double passedParameter1,
		 double passedParameter2,
		 double passedParameter3,
		 double passedParameter4,
		 double passedStartValue)
:RandomVariable(passedTag, RANDOM_VARIABLE_weibull)
{
	tag = passedTag ;
	u = passedParameter1;
	k = passedParameter2;
	startValue = passedStartValue;
}
WeibullRV::WeibullRV(int passedTag, 
		 double passedMean,
		 double passedStdv)
:RandomVariable(passedTag, RANDOM_VARIABLE_weibull)
{
	tag = passedTag;
	setParameters(passedMean,passedStdv);
	startValue = getMean();
}
WeibullRV::WeibullRV(int passedTag, 
		 double passedParameter1,
		 double passedParameter2,
		 double passedParameter3,
		 double passedParameter4)
:RandomVariable(passedTag, RANDOM_VARIABLE_weibull)
{
	tag = passedTag ;
	u = passedParameter1;
	k = passedParameter2;
	startValue = getMean();
}


WeibullRV::~WeibullRV()
{
}


void
WeibullRV::Print(OPS_Stream &s, int flag)
{
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


const char *
WeibullRV::getType()
{
	return "WEIBULL";
}


double 
WeibullRV::getMean()
{
	GammaRV *aGammaRV = new GammaRV(1, 0.0, 1.0, 0.0);
	double result = u * aGammaRV->gammaFunction(1.0+1.0/k);
	delete aGammaRV;
	return result;
}



double 
WeibullRV::getStdv()
{
	GammaRV *aGammaRV = new GammaRV(1, 0.0, 1.0, 0.0);
	double a = aGammaRV->gammaFunction(1.0+2.0/k);
	double b = aGammaRV->gammaFunction(1.0+1.0/k);
	delete aGammaRV;
	double result = u*sqrt(a-b*b);
	return result;
}


double 
WeibullRV::getStartValue()
{
	return startValue;
}

double WeibullRV::getParameter1()  {return u;}
double WeibullRV::getParameter2()  {return k;}
double WeibullRV::getParameter3()  {opserr<<"No such parameter in r.v. #"<<tag<<endln; return 0.0;}
double WeibullRV::getParameter4()  {opserr<<"No such parameter in r.v. #"<<tag<<endln; return 0.0;}

void
WeibullRV::setParameters(double mean, double stdv)
{
	double cov = stdv/mean;
	double xk = 0.0;
	function141(xk, cov,mean);
}


void
WeibullRV::function141(double xk, double cov, double mean)
{
	GammaRV *aGammaRV = new GammaRV(1, 0.0, 1.0, 0.0);
	xk = xk + 1.0;
	double x1 = 1.0 + 1.0/xk;
	double x2 = 1.0 + 2.0/xk;
	double gm1 = aGammaRV->gammaFunction(x1);
	double gm2 = aGammaRV->gammaFunction(x2);
	delete aGammaRV;
	double vy = sqrt(gm2/gm1/gm1 - 1.0);
	if (cov-vy < 0.0) 
		function141(xk,cov,mean);
	else if (cov-vy == 0.0) 
		function144(xk,cov,mean);
	else 
		function142(xk,cov,mean);
}
void
WeibullRV::function142(double xk, double cov, double mean)
{
	xk = xk - 0.1;
	double x1 = 1.0 + 1.0/xk;
	double x2 = 1.0 + 2.0/xk;
	GammaRV *aGammaRV = new GammaRV(1, 0.0, 1.0, 0.0);
	double gm1 = aGammaRV->gammaFunction(x1);
	double gm2 = aGammaRV->gammaFunction(x2);
	delete aGammaRV;
	double vy = sqrt(gm2/gm1/gm1 - 1.0);
	if (cov-vy < 0.0) 
		function143(xk,cov,mean);
	else if (cov-vy == 0.0) 
		function144(xk,cov,mean);
	else 
		function142(xk,cov,mean);
}
void
WeibullRV::function143(double xk, double cov, double mean)
{
	xk = xk + 0.01;
	double x1 = 1.0 + 1.0/xk;
	double x2 = 1.0 + 2.0/xk;
	GammaRV *aGammaRV = new GammaRV(1, 0.0, 1.0, 0.0);
	double gm1 = aGammaRV->gammaFunction(x1);
	double gm2 = aGammaRV->gammaFunction(x2);
	delete aGammaRV;
	double vy = sqrt(gm2/gm1/gm1 - 1.0);
	if (cov-vy < 0.0) 
		function143(xk,cov,mean);
	else if (cov-vy == 0.0) 
		function144(xk,cov,mean);
	else 
		function144(xk,cov,mean);


}
void
WeibullRV::function144(double xk, double gm1, double mean)
{
	u = mean/gm1;
	k = xk;
}
