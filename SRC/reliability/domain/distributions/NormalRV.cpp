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
                                                                        
// $Revision: 1.15 $
// $Date: 2008-05-08 15:32:54 $
// $Source: /usr/local/cvs/OpenSees/SRC/reliability/domain/distributions/NormalRV.cpp,v $


//
// Written by Terje Haukaas (haukaas@ce.berkeley.edu) 
//

#include <NormalRV.h>
#include <Vector.h>
#include <cmath>
#include <float.h>

NormalRV::NormalRV(int passedTag, double passedMean, double passedStdv)
:RandomVariable(passedTag, RANDOM_VARIABLE_normal)
{
	int setp = setParameters(passedMean,passedStdv);
	if (setp < 0)
		opserr << "Error setting parameters in Normal RV with tag " << this->getTag() << endln;

}


NormalRV::NormalRV(int passedTag, const Vector &passedParameters)
:RandomVariable(passedTag, RANDOM_VARIABLE_normal)
{
	if (passedParameters.Size() != 2) {
		opserr << "Normal RV requires 2 parameters, mu and sigma, for RV with tag " <<
		this->getTag() << endln;
		
		// this will create terminal errors
		mu = 0;
		sigma = 0;
		
	} else {
		
		mu = passedParameters(0);
		sigma = passedParameters(1);
	}
}


NormalRV::~NormalRV()
{
}


const char *
NormalRV::getType()
{
	return "NORMAL";
}


double 
NormalRV::getMean()
{
	return mu;
}


double 
NormalRV::getStdv()
{
	return sigma;
}


const Vector &
NormalRV::getParameters(void) {
	static Vector temp(2);
	temp(0) = mu;
	temp(1) = sigma;
	return temp;
}


int 
NormalRV::setParameters(double mean, double stdv)
{
	mu = mean;
	sigma = stdv;
	
	return 0;
}


double
NormalRV::getPDFvalue(double rvValue)
{
	//static const double pi = std::acos(-1.0);
	static const double oneOverSigmaRootTwoPi = 1.0/(sigma * sqrt(2.0*pi));
	
	//return 1 / sqrt ( 2.0 * pi ) * exp ( - 0.5 * pow ( ( ( rvValue - mu ) / sigma ), 2.0 ) );
	return oneOverSigmaRootTwoPi * exp ( - 0.5 * pow ( ( ( rvValue - mu ) / sigma ), 2.0 ) );
}


double
NormalRV::getCDFvalue(double rvValue)
{
    // should add error checking for when sigma is zero
    
    // now use base class functions
	double result = standardNormalPhi( (rvValue-mu)/sigma );

	return result;
}


double
NormalRV::getInverseCDFvalue(double probValue)
{
	// now use base class functions
	double result = getMean() + getStdv() * standardNormalInversePhi(probValue);
    
	return result;
}


double 
NormalRV::transform_x_to_u(void)
{
    return (this->getCurrentValue()-mu)/sigma;
}


double
NormalRV::transform_u_to_x(double uVal)
{
    return uVal*sigma + mu;
}


double
NormalRV::gradient_x_to_u(double uVal)
{
    return sigma;
    
}


int 
NormalRV::getCDFparameterSensitivity(Vector &dFdP)
{
    // returns gradient of F(x) with respect to distribution parameters
    double rvValue = this->getCurrentValue();
    
    // dFdmu
    dFdP(0) = -1 * getPDFvalue(rvValue);
    
    // dFdsig
    dFdP(1) = -(rvValue-mu)/sigma * getPDFvalue(rvValue);
    
    return 0;
}


int
NormalRV::getParameterMeanSensitivity(Vector &dPdmu)
{
    // returns gradient of distribution parameters with respect to the mean
    dPdmu(0) = 1;
    dPdmu(1) = 0;
    
    return 0;
}


int
NormalRV::getParameterStdvSensitivity(Vector &dPdstdv)
{
    // returns gradient of distribution parameters with respect to the stdv
    dPdstdv(0) = 0;
    dPdstdv(1) = 1;
    
    return 0;
}


void
NormalRV::Print(OPS_Stream &s, int flag)
{
	s << "Normal RV # " << this->getTag() << endln;
	s << "\tmean = " << mu << endln;
	s << "\tst.dev. = " << sigma << endln;
}
