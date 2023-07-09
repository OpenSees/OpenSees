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
// $Date: 2008-05-08 15:32:54 $
// $Source: /usr/local/cvs/OpenSees/SRC/reliability/domain/distributions/RayleighRV.cpp,v $


//
// Written by Terje Haukaas (haukaas@ce.berkeley.edu)
//

#include <RayleighRV.h>
#include <Vector.h>
#include <cmath>

RayleighRV::RayleighRV(int passedTag, const Vector &passedParameters)
:RandomVariable(passedTag, RANDOM_VARIABLE_rayleigh)
{
	
	if (passedParameters.Size() != 1) {
		opserr << "Rayleigh RV requires 1 parameter, u, for RV with tag " <<
		this->getTag() << endln;
		
		// this will create terminal errors
		u = 0;
		
	} else {
		
		u = passedParameters(0);
		
	}
}


RayleighRV::~RayleighRV()
{
}


const char *
RayleighRV::getType()
{
	return "RAYLEIGH";
}


double 
RayleighRV::getMean()
{
	//double pi = acos(-1.0);
	return 0.5*sqrt(pi) * u;
}


double 
RayleighRV::getStdv()
{
	//double pi = acos(-1.0);
	return 0.5 * sqrt(4.0-pi) * u;
}


const Vector &
RayleighRV::getParameters(void) {
	static Vector temp(1);
	temp(0) = u;
	return temp;
}


double
RayleighRV::getPDFvalue(double rvValue)
{
	double result;
	if ( 0.0 < rvValue ) {
		result = 2.0*rvValue/(u*u) * exp(-pow(rvValue/u,2));
	}
	else {
		result = 0.0;
	}
	return result;
}


double
RayleighRV::getCDFvalue(double rvValue)
{
	double result;
	if ( 0.0 < rvValue ) {
		result = 1-exp(-pow(rvValue/u,2));
	}
	else {
		result = 0.0;
	}
	return result;
}


double
RayleighRV::getInverseCDFvalue(double probValue)
{
	return u*sqrt(-log(1-probValue));
}


int 
RayleighRV::getCDFparameterSensitivity(Vector &dFdP)
{
    // returns gradient of F(x) with respect to distribution parameters
    double rvValue = this->getCurrentValue();
    
    // dFdu
    dFdP(0) = -rvValue/u * getPDFvalue(rvValue);
    
    return 0;
}


int
RayleighRV::getParameterMeanSensitivity(Vector &dPdmu)
{
    // returns gradient of distribution parameters with respect to the mean
    
    // dudmu
    dPdmu(0) = 2/sqrt(pi);
    
    return 0;
}


int
RayleighRV::getParameterStdvSensitivity(Vector &dPdstdv)
{
    // returns gradient of distribution parameters with respect to the stdv
    
    // dudsig
    dPdstdv(0) = 2/sqrt(4-pi);
    
    return 0;
}


void
RayleighRV::Print(OPS_Stream &s, int flag)
{
	s << "Rayleigh RV #" << this->getTag() << endln;
	s << "\tu = " << u << endln;
}
