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
#include <Vector.h>
#include <cmath>

GumbelRV::GumbelRV(int passedTag, double passedMean, double passedStdv)
:RandomVariable(passedTag, RANDOM_VARIABLE_gumbel)
{
	int setp = setParameters(passedMean,passedStdv);
	if (setp < 0)
		opserr << "Error setting parameters in Gumbel RV with tag " << this->getTag() << endln;
}


GumbelRV::GumbelRV(int passedTag, const Vector &passedParameters)
:RandomVariable(passedTag, RANDOM_VARIABLE_gumbel)
{
	
	if (passedParameters.Size() != 2) {
		opserr << "Gumbel RV requires 2 parameters, u and alpha, for RV with tag " <<
		this->getTag() << endln;
		
		// this will create terminal errors
		u = 0;
		alpha = 0;
		
	} else {
		
		u = passedParameters(0);
		alpha = passedParameters(1);
		
	}
}


GumbelRV::~GumbelRV()
{
}


const char *
GumbelRV::getType()
{
	return "GUMBEL";
}


double 
GumbelRV::getMean()
{
	//double eulergamma = 0.57721566490153286061;
	return u+euler/alpha;
}


double 
GumbelRV::getStdv()
{
	//static const double pi = std::acos(-1.0);
	return pi/(sqrt(6.0)*alpha);
}


const Vector &
GumbelRV::getParameters(void) {
	static Vector temp(2);
	temp(0) = u;
	temp(1) = alpha;
	return temp;
}


int
GumbelRV::setParameters(double mean, double stdv)
{
	//static const double pi = std::acos(-1.0);
	//double eulergamma = 0.57721566490153286061;
	
	u = mean - (euler*sqrt(6.0)*stdv)/pi;
	alpha = pi / ( sqrt(6.0) * stdv );
	
	return 0;
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


int 
GumbelRV::getCDFparameterSensitivity(Vector &dFdP)
{
    // returns gradient of F(x) with respect to distribution parameters
    double rvValue = this->getCurrentValue();
    
    // dFdu
    dFdP(0) = -1 * getPDFvalue(rvValue);
    
    // dFdalpha
    dFdP(1) = -(u-rvValue)/alpha * getPDFvalue(rvValue);
    
    return 0;
}


int
GumbelRV::getParameterMeanSensitivity(Vector &dPdmu)
{
    // returns gradient of distribution parameters with respect to the mean
    
    // dudmu
    dPdmu(0) = 1;
    
    // dalphadmu
    dPdmu(1) = 0;
    
    return 0;
}


int
GumbelRV::getParameterStdvSensitivity(Vector &dPdstdv)
{
    // returns gradient of distribution parameters with respect to the stdv
    double sig = getStdv();
    
    // dudsig
    dPdstdv(0) = -sqrt(6.0)/pi*euler;
    
    // dalphadsig
    dPdstdv(1) = -pi/sqrt(6.0)/sig/sig;
    
    return 0;
}


void
GumbelRV::Print(OPS_Stream &s, int flag)
{
	s << "Gumbel RV #" << this->getTag() << endln;
	s << "\tu = " << u << endln;
	s << "\talpha = " << alpha << endln;
}
