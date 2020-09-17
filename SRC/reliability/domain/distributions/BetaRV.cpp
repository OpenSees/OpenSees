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
// $Source: /usr/local/cvs/OpenSees/SRC/reliability/domain/distributions/BetaRV.cpp,v $


//
// Written by Terje Haukaas (haukaas@ce.berkeley.edu) 
//

#include <BetaRV.h>
#include <Vector.h>
#include <cmath>

BetaRV::BetaRV(int passedTag, const Vector &passedParameters)
:RandomVariable(passedTag, RANDOM_VARIABLE_beta)
{
	if (passedParameters.Size() != 4) {
		opserr << "Beta RV requires 4 parameters, a b q and r, for RV with tag " <<
		this->getTag() << endln;
		
		// this will create terminal errors
		a = 0;
		b = 1;
		q = 0;
		r = 0;
		
	} else {
		
		a = passedParameters(0);
		b = passedParameters(1);
		q = passedParameters(2);
		r = passedParameters(3);
	}
	
}


BetaRV::~BetaRV()
{
}


const char *
BetaRV::getType()
{
	return "BETA";
}


double 
BetaRV::getMean()
{
	return (a*r+b*q)/(q+r);
}


double 
BetaRV::getStdv()
{
	return ((b-a)/(q+r)) * sqrt(q*r/(q+r+1));
}


const Vector &
BetaRV::getParameters(void) {
	static Vector temp(4);
	temp(0) = a;
	temp(1) = b;
	temp(2) = q;
	temp(3) = r;
	return temp;
}


double
BetaRV::getPDFvalue(double rvValue)
{
	double result;
	if ( ((q<1.0 && rvValue > a) || (q >= 1.0 && rvValue >= a)) && ((r < 1.0 && rvValue < b) || (r >= 1.0 && rvValue <= b))) {
		double par1 = pow(rvValue-a,q-1.0);
		double par2 = pow(b-rvValue,r-1.0);
		double par3 = betaFunction(q,r);
		double par4 = pow(b-a,q+r-1.0);
		result = par1*par2/(par3*par4);
	}
	else {
		result = 0.0;
	}
	return result;
}

double
BetaRV::getCDFvalue(double rvValue)
{
	if (rvValue <= a) {
		return 0.0;
	}
	if (rvValue >= b) {
		return 1.0;
	}

	double z = (rvValue - a) / (b - a);
	double Ix = incompleteBetaFunction(q, r, z);
	//opserr << "x = " << rvValue << " z = " << z << " Ix = " << Ix;

	return Ix;
}


double
BetaRV::getInverseCDFvalue(double probValue)
{
	double result = 0.0;
	// Here we want to solve the nonlinear equation:
	//         probValue = getCDFvalue(x)
	// with respect to x. 
	// A Newton scheme to find roots - f(x)=0 - looks something like:
	//         x(i+1) = x(i) - f(xi)/f'(xi)
	// In our case the function f(x) is: f(x) = probValue - getCDFvalue(x)
	// The derivative of the function can be found approximately by a
	// finite difference scheme where e.g. stdv/200 is used as perturbation.

	// Newton method does not converge in general starting from mean particularly when q and/or r are < 1.0
	double tol = 1e-6;
	double x_old = getMean();   // Start at the mean of the random variable
	double x_new;
	double f;
	double df;
	double h;
	double perturbed_f;
	for (int i=1;  i<=100;  i++ )  {

		// Evaluate function
		f = probValue - getCDFvalue(x_old);
		// Evaluate perturbed function
		h = getStdv()/200.0;
		perturbed_f = probValue - getCDFvalue(x_old+h);

		// Evaluate derivative of function
		df = ( perturbed_f - f ) / h;

		if ( fabs(df) < 1.0e-15) {
			/*opserr << "WARNING: BetaRV::getInverseCDFvalue() -- zero derivative " << " in Newton algorithm. " << endln;
			opserr << "Switch bisection method" << endln;*/
			// In case of zero derivative in Newton MEthod (This often happen because we went outside the boundaries)
			// We use slower but robust bisection method
			double x_A, x_B;
			double f_A, f_B, f_M;
			x_A = a;
			x_B = b;
			f_A = probValue - getCDFvalue(x_A);
			f_B = probValue - getCDFvalue(x_B);
			for (int j = 1; j <= 200; j++) {
				x_new = (x_A + x_B) / 2.0;
				f_M = probValue - getCDFvalue(x_new);
				if (fabs(f_M) < tol) {
					return x_new;
				}
				else {
					if (i == 200) {
						opserr << "WARNING: Did not converge to find inverse CDF!" << endln;
						return 0.0;
					}
					if (f_M * f_A < 0) {
						// The zero is between A and M
						x_B = x_new;
						f_B = f_M;
					}
					else {
						// The zero is between M and B
						x_A = x_new;
						f_A = f_M;
					}
				}
			}

		}
		else {

			// Take a Newton step
			x_new = x_old - f/df;
			
			// Check convergence; quit or continue
			// Changed the tolerance check because if x_new is zero it doesn't work
			if (fabs(probValue - getCDFvalue(x_new)) < tol) {
				return x_new;
			}
			else {
				if (i==100) {
					opserr << "WARNING: Did not converge to find inverse CDF!" << endln;
					return 0.0;
				}
				else {
					x_old = x_new;
				}
			
			}
		}
	}

	return result;
}


int 
BetaRV::getCDFparameterSensitivity(Vector &dFdP)
{
    // returns gradient of F(x) with respect to distribution parameters
    double rvValue = this->getCurrentValue();
    double cdfValue = getCDFvalue(rvValue);
    
    // dFda
    // no closed form expression for this one, use finite differences
    double a_old = a;
    double dh = a/1000.0;
    a += dh;
    dFdP(0) = ( getCDFvalue(rvValue) - cdfValue )/dh;
    a = a_old;
    
    // dFdb
    // no closed form expression for this one, use finite differences
    double b_old = b;
    dh = b/1000.0;
    b += dh;
    dFdP(1) = ( getCDFvalue(rvValue) - cdfValue )/dh;
    b = b_old;
    
    // dFdq
    // no closed form expression for this one, use finite differences
    double q_old = q;
    dh = q/1000.0;
    q += dh;
    dFdP(2) = ( getCDFvalue(rvValue) - cdfValue )/dh;
    q = q_old;
    
    // dFdr
    // no closed form expression for this one, use finite differences
    double r_old = r;
    dh = r/1000.0;
    r += dh;
    dFdP(3) = ( getCDFvalue(rvValue) - cdfValue )/dh;
    r = r_old;    
    
    return 0;
}


int
BetaRV::getParameterMeanSensitivity(Vector &dPdmu)
{
    // returns gradient of distribution parameters with respect to the mean
    
    // dadmu
    dPdmu(0) = 0;
    
    // dbdmu
    dPdmu(1) = 0;
    
    // dqdmu
    dPdmu(2) = (q+q*q-r-q*r-2*r*r) / (a*r-b*r);
    
    // drdmu
    dPdmu(3) = (2*q*q+q+q*r-r-r*r) /q/(a-b);
    
    return 0;
}


int
BetaRV::getParameterStdvSensitivity(Vector &dPdstdv)
{
    // returns gradient of distribution parameters with respect to the stdv
    
    // dadsig
    dPdstdv(0) = 0;
    
    // dbdsig
    dPdstdv(1) = 0;
    
    // dqdsig
    dPdstdv(2) = 2*q*q*r/(a-b)/pow(q*r/(1+q+r),1.5);
    
    // drdsig
    dPdstdv(3) = 2*q*r*r/(a-b)/pow(q*r/(1+q+r),1.5);
    
    return 0;
}


void
BetaRV::Print(OPS_Stream &s, int flag)
{
	s << "Beta RV #" << this->getTag() << endln;
	s << "\ta = " << a << ", b = " << b << endln;
	s << "\tq = " << q << ", r = " << r << endln;
}
