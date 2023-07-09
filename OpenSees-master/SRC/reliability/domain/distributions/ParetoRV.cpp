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
// $Source: /usr/local/cvs/OpenSees/SRC/reliability/domain/distributions/ParetoRV.cpp,v $


//
// Written by Terje Haukaas (haukaas@ce.berkeley.edu) 
//

#include <ParetoRV.h>
#include <cmath>
#include <Vector.h>

ParetoRV::ParetoRV(int passedTag, const Vector &passedParameters)
:RandomVariable(passedTag, RANDOM_VARIABLE_pareto)
{
	if (passedParameters.Size() != 2) {
		opserr << "Pareto RV requires 2 parameters, k and u, for RV with tag " <<
		this->getTag() << endln;
		
		// this will create terminal errors
		k = 0;
		u = 0;
		
	} else {
		
		k = passedParameters(0);
		u = passedParameters(1);
	}
	
}


ParetoRV::~ParetoRV()
{
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


const Vector &
ParetoRV::getParameters(void) {
	static Vector temp(2);
	temp(0) = k;
	temp(1) = u;
	return temp;
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
ParetoRV::getInverseCDFvalue(double probValue)
{
	if (k <= 0) {
		// shape should be greater than 0
		return 0.0;
	}
	else {
		return pow((1 - probValue) / pow(u, k), -1 / k);
	}
}


void
ParetoRV::Print(OPS_Stream &s, int flag)
{
	s << "Pareto RV #" << this->getTag() << endln;
	s << "\tk = " << k << endln;
	s << "\tu = " << u << endln;
}
