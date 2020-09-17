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
// $Source: /usr/local/cvs/OpenSees/SRC/reliability/domain/distributions/ChiSquareRV.cpp,v $


//
// Written by Terje Haukaas (haukaas@ce.berkeley.edu) 
//

#include <ChiSquareRV.h>
#include <Vector.h>
#include <cmath>

ChiSquareRV::ChiSquareRV(int passedTag, double passedMean, double passedStdv)
:RandomVariable(passedTag, RANDOM_VARIABLE_chisquare)
{
	int setp = setParameters(passedMean,passedStdv);
	if (setp < 0)
		opserr << "Error setting parameters in ChiSquare RV with tag " << this->getTag() << endln;
}


ChiSquareRV::ChiSquareRV(int passedTag, const Vector &passedParameters)
:RandomVariable(passedTag, RANDOM_VARIABLE_chisquare)
{
	
	if (passedParameters.Size() != 1) {
		opserr << "ChiSquare RV requires 1 parameter, nu, for RV with tag " <<
		this->getTag() << endln;
		
		// this will create terminal errors
		nu = 0;
		
	} else {
		
		nu = passedParameters(0);
		
	}
}


ChiSquareRV::~ChiSquareRV()
{
}


const char *
ChiSquareRV::getType()
{
	return "CHISQUARE";
}


double 
ChiSquareRV::getMean()
{
	return nu;
}


double 
ChiSquareRV::getStdv()
{
	return sqrt(2*nu);
}


const Vector &
ChiSquareRV::getParameters(void) {
	static Vector temp(1);
	temp(0) = nu;
	return temp;
}


int
ChiSquareRV::setParameters(double mean, double stdv)
{
	nu = mean;
	
	return 0;
}


double
ChiSquareRV::getPDFvalue(double rvValue)
{
	double result;
	if ( 0.0 < rvValue ) {
		double a = gammaFunction(0.5*nu);
		result = pow(0.5*rvValue,0.5*nu-1.0)*exp(-0.5*rvValue)/(2.0*a);
	}
	else {
		result = 0.0;
	}
	return result;
}


double
ChiSquareRV::getCDFvalue(double rvValue)
{
	double result;
	if ( 0.0 < rvValue ) {
		result = incompleteGammaFunction(0.5*nu,0.5*rvValue);
	}
	else {
		result = 0.0;
	}
	return result;
}


double
ChiSquareRV::getInverseCDFvalue(double probValue)
{
	// Here we want to solve the nonlinear equation:
	//         probValue = getCDFvalue(x)
	// with respect to x. 
	// A Newton scheme to find roots - f(x)=0 - looks something like:
	//         x(i+1) = x(i) - f(xi)/f'(xi)
	// In our case the function f(x) is: f(x) = probValue - getCDFvalue(x)

	double x_old = nu;

	double tol = 1.0e-7;
	double x_new = x_old;
	double dx = x_old;
	int step = 1;
	int nmax = 50;

	for (int i = 1; i <= nmax; i++) {

		dx = (getCDFvalue(x_new) - probValue) / getPDFvalue(x_new);
		
		// Take a Newton step
		x_new = x_old - dx;

		if (fabs(getCDFvalue(x_new) - probValue) < tol) {
			return x_new;
		}
		else {
			if (i == nmax) {
				opserr << "WARNING: ChiSquaredRV did not converge to find inverse CDF!" << endln;
				return 0.0;
			}
			else {
				x_old = x_new;
			}
		}
	}

	return x_new;
}


int 
ChiSquareRV::getCDFparameterSensitivity(Vector &dFdP)
{
    // returns gradient of F(x) with respect to distribution parameters
    double rvValue = this->getCurrentValue();
    double cdfValue = getCDFvalue(rvValue);
    
    // dFdnu
    // no closed form expression for this one, use finite differences
    double nu_old = nu;
    double dh = nu/1000.0;
    nu += dh;
    dFdP(0) = ( getCDFvalue(rvValue) - cdfValue )/dh;
    nu = nu_old;
    
    return 0;
}


int
ChiSquareRV::getParameterMeanSensitivity(Vector &dPdmu)
{
    // returns gradient of distribution parameters with respect to the mean
    
    // dnudmu
    dPdmu(0) = 0.5;
    
    return 0;
}


int
ChiSquareRV::getParameterStdvSensitivity(Vector &dPdstdv)
{
    // returns gradient of distribution parameters with respect to the stdv
    double sig = getStdv();
    
    // dnudsig
    dPdstdv(0) = sig;
    
    return 0;
}


void
ChiSquareRV::Print(OPS_Stream &s, int flag)
{
	s << "ChiSquare RV #" << this->getTag() << endln;
	s << "\tnu = " << nu << endln;
}
