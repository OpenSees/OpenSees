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
                                                                        
// $Revision: 1.10 $
// $Date: 2008-04-10 16:23:32 $
// $Source: /usr/local/cvs/OpenSees/SRC/reliability/domain/distributions/LognormalRV.cpp,v $


//
// Written by Terje Haukaas (haukaas@ce.berkeley.edu)
//

#include <LognormalRV.h>
#include <NormalRV.h>
#include <cmath>
#include <Vector.h>

LognormalRV::LognormalRV(int passedTag, double passedMean, double passedStdv)
:RandomVariable(passedTag, RANDOM_VARIABLE_lognormal)
{
	if (passedMean<0.0) {
		isPositive = false;
		passedMean = -passedMean;
	}
	else {
		isPositive = true;
	}

	int setp = setParameters(passedMean,passedStdv);
	if (setp < 0)
		opserr << "Error setting parameters in Lognormal RV with tag " << this->getTag() << endln;
	
}


LognormalRV::LognormalRV(int passedTag, const Vector &passedParameters)
:RandomVariable(passedTag, RANDOM_VARIABLE_lognormal)
{
	if (passedParameters.Size() != 2) {
		opserr << "Lognormal RV requires 2 parameters, lambda and zeta, for RV with tag " <<
			this->getTag() << endln;
		
		// this will create terminal errors
		lambda = 0;
		zeta = 0;
		isPositive = true;
		
	} else {
		
		lambda = passedParameters(0);
		zeta = passedParameters(1);

		if (lambda<0.0) {
			isPositive = false;
			lambda = -lambda;
		}
		else {
			isPositive = true;
		}
	}
}


LognormalRV::~LognormalRV()
{
}


const char *
LognormalRV::getType()
{
	return "LOGNORMAL";
}


double 
LognormalRV::getMean()
{
	if (isPositive)
		return exp(lambda+0.5*zeta*zeta);
	else
		return -exp(lambda+0.5*zeta*zeta);
}


double 
LognormalRV::getStdv()
{
	return exp(lambda+0.5*zeta*zeta)*sqrt(exp(zeta*zeta)-1);
}


const Vector &
LognormalRV::getParameters(void) {
	static Vector temp(2);
	temp(0) = lambda;
	temp(1) = zeta;
	return temp;
}


int 
LognormalRV::setParameters(double mean, double stdv)
{
	zeta = sqrt(   log(   1+(stdv/mean)*(stdv/mean)   )   );
	lambda = log(mean) - 0.5*zeta*zeta;
	
	return 0;
}


double
LognormalRV::getPDFvalue(double rvValue)
{
	//static const double pi = std::acos(-1.0);
	
	if (!isPositive) {
		// The formal answer is: f(x) = f_pos(x+2|x|), but let's do it simple here
		rvValue = -rvValue;
	}

	double result;
	if ( 0.0 < rvValue ) {
		result = 1/(sqrt(2*pi)*zeta*rvValue) * exp(-0.5* pow ( (log(rvValue)-lambda) / zeta, 2 )  );
	}
	else {
		result = 0.0;
	}
	return result;
}


double
LognormalRV::getCDFvalue(double rvValue)
{
	double result;

    // now use base class functions
	if (isPositive) {
		if ( 0.0 < rvValue ) {
			result = standardNormalPhi( (log(rvValue)-lambda)/zeta );
		}
		else {
			result = 0.0;
		}
	}
	else {
		if ( rvValue < 0.0 ) {
			result = standardNormalPhi( (log(fabs(rvValue))-lambda)/zeta );
			result = 1.0-result;
		}
		else {
			result = 1.0;
		}
	}

	// Return result depending on type of random variable
	if (isPositive) {
		return result;
	}
	else {
		return 1-result;
	}



/*
	// First, flip it around if it's a negative lognormal
	if (!isPositive) {
		rvValue = -rvValue;
	}

	// Compute the ordinary CDF
	double result;
	if ( 0.0 < rvValue ) {
		RandomVariable *aStandardNormalRV;
		aStandardNormalRV= new NormalRV( 1, 0.0, 1.0, 0.0);
		result = aStandardNormalRV->getCDFvalue((log(rvValue)-lambda)/zeta);
		delete aStandardNormalRV;	
	}
	else {
		result = 0.0;
	}

	// Return result depending on type of random variable
	if (isPositive) {
		return result;
	}
	else {
		return 1-result;
	}
*/
}


double
LognormalRV::getInverseCDFvalue(double probValue)
{
    double inverseNormal = 0;
    
    // now use base class functions
    if (isPositive) {
        inverseNormal = standardNormalInversePhi(probValue);
        return exp(inverseNormal*zeta + lambda);
    }
    else {
        inverseNormal = standardNormalInversePhi(1.0-probValue);
        return (-exp(inverseNormal*zeta + lambda));
    }
}


double 
LognormalRV::transform_x_to_u(void)
{
    if (isPositive)
        return (log(this->getCurrentValue())-lambda)/zeta;
    else
        return -(log(fabs(this->getCurrentValue())) - lambda)/zeta;
}


double
LognormalRV::transform_u_to_x(double uVal)
{
    if (isPositive)
        return exp(uVal*zeta + lambda);
    else
        return -exp(-uVal*zeta + lambda);
}


double
LognormalRV::gradient_x_to_u(double uVal)
{
    // check this works with negative LN
    return this->getCurrentValue()*zeta;
    
}


int 
LognormalRV::getCDFparameterSensitivity(Vector &dFdP)
{
    // returns gradient of F(x) with respect to distribution parameters
    double rvValue = this->getCurrentValue();
    
    // dFdlambda
    dFdP(0) = -rvValue * getPDFvalue(rvValue);
    
    // dFdzeta
    dFdP(1) = -1/zeta*rvValue*(log(rvValue)-lambda) * getPDFvalue(rvValue);
    
    return 0;
}


int
LognormalRV::getParameterMeanSensitivity(Vector &dPdmu)
{
    // returns gradient of distribution parameters with respect to the mean
    double mu = getMean();
    double sig = getStdv();
    
    // dlambdadmu
    dPdmu(0) = (mu*mu+2*sig*sig)/mu/(mu*mu+sig*sig);
    
    // dzetadmu
    dPdmu(1) = -sig*sig/mu/(mu*mu+sig*sig)/zeta;
    
    return 0;
}


int
LognormalRV::getParameterStdvSensitivity(Vector &dPdstdv)
{
    // returns gradient of distribution parameters with respect to the stdv
    double mu = getMean();
    double sig = getStdv();
    
    // dlambdadsig
    dPdstdv(0) = -sig/(mu*mu+sig*sig);
    
    // dzetadsig
    dPdstdv(1) = sig/(mu*mu+sig*sig)/zeta;
    
    return 0;
}


void
LognormalRV::Print(OPS_Stream &s, int flag)
{
  s << "Lognormal RV #" << this->getTag() << endln;
  s << "\tlambda = " << lambda << endln;
  s << "\tzeta = " << zeta << endln;
}
