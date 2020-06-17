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
#include <Vector.h>
#include <cmath>

ShiftedRayleighRV::ShiftedRayleighRV(int passedTag, 
				     double passedMean, double passedStdv)
:RandomVariable(passedTag, RANDOM_VARIABLE_shiftedrayleigh)
{
	int setp = setParameters(passedMean,passedStdv);
	if (setp < 0)
		opserr << "Error setting parameters in ShiftedRayleigh RV with tag " << this->getTag() << endln;
}


ShiftedRayleighRV::ShiftedRayleighRV(int passedTag,
				     const Vector &passedParameters)
:RandomVariable(passedTag, RANDOM_VARIABLE_shiftedrayleigh)
{
	
	if (passedParameters.Size() != 2) {
		opserr << "ShiftedRayleigh RV requires 2 parameters, u and x0, for RV with tag " <<
		this->getTag() << endln;
		
		// this will create terminal errors
		u = 0;
		x0 = 0;
		
	} else {
		
		u = passedParameters(0);
		x0 = passedParameters(1);
		
	}
}


ShiftedRayleighRV::~ShiftedRayleighRV()
{
}


const char *
ShiftedRayleighRV::getType()
{
	return "SHIFTEDRAYLEIGH";
}


double 
ShiftedRayleighRV::getMean()
{
	return x0 + 0.5 * u * sqrt(pi);
}


double 
ShiftedRayleighRV::getStdv()
{
	return 0.5 * u * sqrt(4.0-pi);
}


const Vector &
ShiftedRayleighRV::getParameters(void) {
	static Vector temp(2);
	temp(0) = u;
	temp(1) = x0;
	return temp;
}


int
ShiftedRayleighRV::setParameters(double mean, double stdv)
{
	u = 2.0 * stdv / sqrt(4.0 - pi);
	x0 = mean - stdv*sqrt(pi) / sqrt(4.0 - pi);
	
	return 0;
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


int 
ShiftedRayleighRV::getCDFparameterSensitivity(Vector &dFdP)
{
    // returns gradient of F(x) with respect to distribution parameters
    double rvValue = this->getCurrentValue();
    
    // dFdu
    dFdP(0) = -(rvValue-x0)/u * getPDFvalue(rvValue);
    
    // dFdx0
    dFdP(1) = -1 * getPDFvalue(rvValue);
    
    return 0;
}


int
ShiftedRayleighRV::getParameterMeanSensitivity(Vector &dPdmu)
{
    // returns gradient of distribution parameters with respect to the mean
    
    // dudmu
    dPdmu(0) = 0;
    
    // dx0dmu
    dPdmu(1) = 1;
    
    return 0;
}


int
ShiftedRayleighRV::getParameterStdvSensitivity(Vector &dPdstdv)
{
    // returns gradient of distribution parameters with respect to the stdv
    
    // dudsig
    dPdstdv(0) = 2/sqrt(4.0+3.0*pi);
    
    // dx0dsig
    dPdstdv(1) = sqrt(pi)/sqrt(4.0+3.0*pi);
    
    return 0;
}


void
ShiftedRayleighRV::Print(OPS_Stream &s, int flag)
{
	s << "ShiftedRayleigh RV #" << this->getTag() << endln;
	s << "\tu = " << u << endln;
	s << "\tx0 = " << x0 << endln;
}
