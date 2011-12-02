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
// $Source: /usr/local/cvs/OpenSees/SRC/reliability/domain/distributions/ParetoRV.cpp,v $


//
// Written by Terje Haukaas (haukaas@ce.berkeley.edu) 
//

#include <ParetoRV.h>
#include <math.h>
#include <string.h>
#include <classTags.h>
#include <OPS_Globals.h>

ParetoRV::ParetoRV(int passedTag, 
		 double passedParameter1,
		 double passedParameter2,
		 double passedParameter3,
		 double passedParameter4,
		 double passedStartValue)
:RandomVariable(passedTag, RANDOM_VARIABLE_pareto)
{
	tag = passedTag ;
	k = passedParameter1;
	u = passedParameter2;
	startValue = passedStartValue;
}
ParetoRV::ParetoRV(int passedTag, 
		 double passedParameter1,
		 double passedParameter2,
		 double passedParameter3,
		 double passedParameter4)
:RandomVariable(passedTag, RANDOM_VARIABLE_pareto)
{
	tag = passedTag ;
	k = passedParameter1;
	u = passedParameter2;
	startValue = getMean();
}


ParetoRV::~ParetoRV()
{
}


void
ParetoRV::Print(OPS_Stream &s, int flag)
{
}


double
ParetoRV::getPDFvalue(double rvValue)
{
	double result;
	if ( u <= rvValue ) {
		result = k/u * pow(u/rvValue,k+1);
	}
	else {
		result = 0.0;
	}
	return result;
}


double
ParetoRV::getCDFvalue(double rvValue)
{
	double result;
	if ( u <= rvValue ) {
		result = 1-pow(u/rvValue,k);
	}
	else {
		result = 0.0;
	}
	return result;
}


double
ParetoRV::getInverseCDFvalue(double rvValue)
{
	return 0.0;
}


const char *
ParetoRV::getType()
{
	return "PARETO";
}


double 
ParetoRV::getMean()
{
	return k*u/(k-1);
}



double 
ParetoRV::getStdv()
{
	return sqrt(k/(k-2))*(u/(k-1));
}


double 
ParetoRV::getStartValue()
{
	return startValue;
}

double ParetoRV::getParameter1()  {return k;}
double ParetoRV::getParameter2()  {return u;}
double ParetoRV::getParameter3()  {opserr<<"No such parameter in r.v. #"<<tag<<endln; return 0.0;}
double ParetoRV::getParameter4()  {opserr<<"No such parameter in r.v. #"<<tag<<endln; return 0.0;}

