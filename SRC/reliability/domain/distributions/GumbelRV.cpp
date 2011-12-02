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
                                                                        
// $Revision: 1.9 $
// $Date: 2008-05-08 15:32:54 $
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
:RandomVariable(passedTag, RANDOM_VARIABLE_gumbel, passedStartValue)
{
	setParameters(passedMean,passedStdv);
}
GumbelRV::GumbelRV(int passedTag, 
		 double passedParameter1,
		 double passedParameter2,
		 double passedParameter3,
		 double passedParameter4,
		 double passedStartValue)
:RandomVariable(passedTag, RANDOM_VARIABLE_gumbel, passedStartValue)
{
	u = passedParameter1;
	alpha = passedParameter2;
}
GumbelRV::GumbelRV(int passedTag, 
		 double passedMean,
		 double passedStdv)
:RandomVariable(passedTag, RANDOM_VARIABLE_gumbel, passedMean)
{
	setParameters(passedMean,passedStdv);
}
GumbelRV::GumbelRV(int passedTag, 
		 double passedParameter1,
		 double passedParameter2,
		 double passedParameter3,
		 double passedParameter4)
:RandomVariable(passedTag, RANDOM_VARIABLE_gumbel)
{
	u = passedParameter1;
	alpha = passedParameter2;
	this->setStartValue(getMean());
}


GumbelRV::~GumbelRV()
{
}

void
GumbelRV::setParameters(double mean, double stdv)
{
	double pi = acos(-1.0);
	double eulergamma = 0.57721566490153286061;
	u = mean - (eulergamma*sqrt(6.0)*stdv)/pi;
	alpha = pi / ( sqrt(6.0) * stdv );
}

void
GumbelRV::Print(OPS_Stream &s, int flag)
{
  s << "Gumbel RV #" << this->getTag() << endln;
  s << "\tu = " << u << endln;
  s << "\talpha = " << alpha << endln;
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
	double eulergamma = 0.57721566490153286061;
	return u+eulergamma/alpha;
}



double 
GumbelRV::getStdv()
{
	double pi = acos(-1.0);
	return pi/(sqrt(6.0)*alpha);
}


double
GumbelRV::getParameter1()
{
  return u;
}

double
GumbelRV::getParameter2()
{
  return alpha;
}
