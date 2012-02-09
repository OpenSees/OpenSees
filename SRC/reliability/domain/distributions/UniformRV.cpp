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
// $Source: /usr/local/cvs/OpenSees/SRC/reliability/domain/distributions/UniformRV.cpp,v $


//
// Written by Terje Haukaas (haukaas@ce.berkeley.edu)
//

#include <UniformRV.h>
#include <Vector.h>
#include <cmath>

UniformRV::UniformRV(int passedTag, double passedMean, double passedStdv)
:RandomVariable(passedTag, RANDOM_VARIABLE_uniform)
{
	int setp = setParameters(passedMean,passedStdv);
	if (setp < 0)
		opserr << "Error setting parameters in Uniform RV with tag " << this->getTag() << endln;
}


UniformRV::UniformRV(int passedTag, const Vector &passedParameters)
:RandomVariable(passedTag, RANDOM_VARIABLE_uniform)
{
	
	if (passedParameters.Size() != 2) {
		opserr << "Uniform RV requires 2 parameters, a and b, for RV with tag " <<
		this->getTag() << endln;
		
		// this will create terminal errors
		a = 0;
		b = 0;
		
	} else {
		
		a = passedParameters(0);
		b = passedParameters(1);
		
	}
}


UniformRV::~UniformRV()
{
}


const char *
UniformRV::getType()
{
	return "UNIFORM";
}


double 
UniformRV::getMean()
{
	//return (a+b)/2.0;
	return 0.5*(a+b);
}


double 
UniformRV::getStdv()
{
	static const double oneOverTwoRootThree = 1.0/(2.0*sqrt(3.0));
	
	//return (b-a)/(2.0*sqrt(3.0));
	return (b-a)*oneOverTwoRootThree;
}


const Vector &
UniformRV::getParameters(void) {
	static Vector temp(2);
	temp(0) = a;
	temp(1) = b;
	return temp;
}


int
UniformRV::setParameters(double mean, double stdv)
{
	a = mean - sqrt(3.0)*stdv;
	b = mean + sqrt(3.0)*stdv;
	
	return 0;
}


double
UniformRV::getPDFvalue(double rvValue)
{
	double result;
	if ( rvValue >= a && rvValue <= b ) {
		result = 1.0/(b-a);
	}
	else {
		result = 0.0;
	}
	return result;
}


double
UniformRV::getCDFvalue(double rvValue)
{
	double result;
	if ( rvValue >= a && rvValue <= b ) {
		result = (rvValue-a)/(b-a);
	}
	else if (rvValue > b) {
		result = 1.0;
	}
	else {
		result = 0.0;
	}
	return result;
}

double
UniformRV::getInverseCDFvalue(double probValue)
{
  //return probValue * b - probValue * a + a;
  return probValue * (b-a) + a;
}


double 
UniformRV::transform_x_to_u(void)
{
    double cdf = getCDFvalue(this->getCurrentValue());
    return standardNormalInversePhi(cdf);
}


double
UniformRV::transform_u_to_x(double uVal)
{
    return getInverseCDFvalue(standardNormalPhi(uVal));
}


double
UniformRV::gradient_x_to_u(double uVal)
{
    static const double oneOverRootTwoPi = 1.0/sqrt(2.0*pi);
	double pdf = oneOverRootTwoPi * exp ( -0.5 * pow(uVal, 2.0) );
    return pdf/getPDFvalue(this->getCurrentValue());
    
}


void
UniformRV::Print(OPS_Stream &s, int flag)
{
	s << "Uniform RV #" << this->getTag() << endln;
	s << "\ta = " << a << endln;
	s << "\tb = " << b << endln;
}
