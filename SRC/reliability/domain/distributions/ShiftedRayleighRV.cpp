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
// $Source: /usr/local/cvs/OpenSees/SRC/reliability/domain/distributions/ShiftedRayleighRV.cpp,v $


//
// Written by Terje Haukaas (haukaas@ce.berkeley.edu)
//

#include <ShiftedRayleighRV.h>
#include <math.h>
#include <string.h>
#include <classTags.h>
#include <OPS_Globals.h>

ShiftedRayleighRV::ShiftedRayleighRV(int passedTag, 
		 double passedMean,
		 double passedStdv,
		 double passedStartValue)
:RandomVariable(passedTag, RANDOM_VARIABLE_shiftedrayleigh, passedStartValue)
{
	double pi = 3.14159265358979;
	u = 2.0 * passedStdv / sqrt(3.0*pi+4.0);
	x0 = passedMean + passedStdv*sqrt(pi) / sqrt(3.0*pi+4.0);
}
ShiftedRayleighRV::ShiftedRayleighRV(int passedTag, 
		 double passedParameter1,
		 double passedParameter2,
		 double passedParameter3,
		 double passedParameter4,
		 double passedStartValue)
:RandomVariable(passedTag, RANDOM_VARIABLE_shiftedrayleigh, passedStartValue)
{
	u = passedParameter1;
	x0 = passedParameter2;
}
ShiftedRayleighRV::ShiftedRayleighRV(int passedTag, 
		 double passedMean,
		 double passedStdv)
:RandomVariable(passedTag, RANDOM_VARIABLE_shiftedrayleigh, passedMean)
{
	double pi = 3.14159265358979;
	u = 2.0 * passedStdv / sqrt(3.0*pi+4.0);
	x0 = passedMean + passedStdv*sqrt(pi) / sqrt(3.0*pi+4.0);
}
ShiftedRayleighRV::ShiftedRayleighRV(int passedTag, 
		 double passedParameter1,
		 double passedParameter2,
		 double passedParameter3,
		 double passedParameter4)
:RandomVariable(passedTag, RANDOM_VARIABLE_shiftedrayleigh)
{
	u = passedParameter1;
	x0 = passedParameter2;
	this->setStartValue(getMean());
}


ShiftedRayleighRV::~ShiftedRayleighRV()
{
}


void
ShiftedRayleighRV::Print(OPS_Stream &s, int flag)
{
}


double
ShiftedRayleighRV::getPDFvalue(double rvValue)
{
	double result;
	if ( x0 < rvValue ) {
		result = 2.0*(rvValue-x0)/(u*u) * exp(-pow(((rvValue-x0)/u),2));
	}
	else {
		result = 0.0;
	}
	return result;
}


double
ShiftedRayleighRV::getCDFvalue(double rvValue)
{
	double result;
	if ( x0 < rvValue ) {
		result = 1 - exp(-pow(((rvValue-x0)/u),2));
	}
	else {
		result = 0.0;
	}
	return result;
}


double
ShiftedRayleighRV::getInverseCDFvalue(double probValue)
{
	return x0 + u * sqrt(-log(1-probValue));
}


const char *
ShiftedRayleighRV::getType()
{
	return "SHIFTEDRAYLEIGH";
}


double 
ShiftedRayleighRV::getMean()
{
	double pi = 3.14159265358979;
	return x0 + 0.5 * u * sqrt(pi);
}



double 
ShiftedRayleighRV::getStdv()
{
	double pi = 3.14159265358979;
	return 0.5 * u * sqrt(4.0-pi);
}


double
ShiftedRayleighRV::getParameter1()
{
  return u;
}

double
ShiftedRayleighRV::getParameter2()
{
  return x0;
}
