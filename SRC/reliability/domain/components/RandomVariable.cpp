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
// $Date: 2008-04-10 16:22:57 $
// $Source: /usr/local/cvs/OpenSees/SRC/reliability/domain/components/RandomVariable.cpp,v $


//
// Written by Terje Haukaas (haukaas@ce.berkeley.edu) 
//

#include <RandomVariable.h>
#include <cmath>
#include <classTags.h>
#include <OPS_Globals.h>
#include <float.h>

// define common math constants needed by the random variables
const double RandomVariable::pi = std::acos(-1.0);
const double RandomVariable::euler = 0.57721566490153286061;
const double RandomVariable::zeta3 = 1.2020569031595942854;
const double RandomVariable::zeta5 = 1.0369277551433699263;
#define isnan(x) ((x)!=(x))


RandomVariable::RandomVariable(int tag, int classTag)
  :ReliabilityDomainComponent(tag,classTag),
   startValue(0.0), rvValue(0.0)
{

}


RandomVariable::~RandomVariable()
{
}


double 
RandomVariable::transform_x_to_u(void)
{
    double cdf = getCDFvalue(this->getCurrentValue());
    return standardNormalInversePhi(cdf);
}


double
RandomVariable::transform_u_to_x(double uVal)
{
    return getInverseCDFvalue(standardNormalPhi(uVal));
}


double
RandomVariable::gradient_x_to_u(double uVal)
{
    static const double oneOverRootTwoPi = 1.0/sqrt(2.0*pi);
	double pdf = oneOverRootTwoPi * exp ( -0.5 * pow(uVal, 2.0) );
    return pdf/getPDFvalue(this->getCurrentValue());
    
}


double 
RandomVariable::getCDFMeanSensitivity(void)
{
    // size vector based on number of parameters in random variable
    Vector temp = this->getParameters();
    Vector dFdP(temp.Size());
    Vector dPdmu(temp.Size());
    
    // returns dF/dmu
    this->getCDFparameterSensitivity(dFdP);
    this->getParameterMeanSensitivity(dPdmu);
    
    return dFdP^dPdmu;
}


double 
RandomVariable::getCDFStdvSensitivity(void)
{
    // size vector based on number of parameters in random variable
    Vector temp = this->getParameters();
    Vector dFdP(temp.Size());
    Vector dPdsig(temp.Size());
    
    // returns dF/dsigma
    this->getCDFparameterSensitivity(dFdP);
    this->getParameterStdvSensitivity(dPdsig);
    
    return dFdP^dPdsig;
}


void
RandomVariable::Print(OPS_Stream &s, int flag)
{
  s << "RV #" << this->getTag() << ' ' << this->getType() << endln;
  s << "\tParameters: " << this->getParameters() << endln;

}


int 
RandomVariable::setNewTag(int newTag)
{
  TaggedObject::setTag(newTag);
  
  return 0;
}


double
RandomVariable::gammaFunction(double x)
{
#ifdef _LINUX
	return tgamma(x);
#endif
    
	double res;
	//static const double pi = std::acos(-1.0);
	
	if (x==0 || ( x < 0.0 && floor(x)==x ) )  {
		opserr << "Invalid input to the gamma function" << endln;
	}
	else {
		Vector p(9);
		Vector q(9);
		Vector c(8);
		
		p(0) = 0.0;
		p(1) = -1.71618513886549492533811e+0; 
		p(2) = 2.47656508055759199108314e+1;
		p(3) = -3.79804256470945635097577e+2; 
		p(4) = 6.29331155312818442661052e+2;
		p(5) = 8.66966202790413211295064e+2; 
		p(6) = -3.14512729688483675254357e+4;
		p(7) = -3.61444134186911729807069e+4; 
		p(8) = 6.64561438202405440627855e+4;
		
		q(0) = 0.0;
		q(1) = -3.08402300119738975254353e+1; 
		q(2) = 3.15350626979604161529144e+2;
		q(3) = -1.01515636749021914166146e+3; 
		q(4) = -3.10777167157231109440444e+3;
		q(5) = 2.25381184209801510330112e+4; 
		q(6) = 4.75584627752788110767815e+3;
		q(7) = -1.34659959864969306392456e+5; 
		q(8) = -1.15132259675553483497211e+5;
		
		c(0) = 0.0;
		c(1) = -1.910444077728e-03; 
		c(2) = 8.4171387781295e-04;
		c(3) = -5.952379913043012e-04; 
		c(4) = 7.93650793500350248e-04;
		c(5) = -2.777777777777681622553e-03; 
		c(6) = 8.333333333333333331554247e-02;
		c(7) = 5.7083835261e-03;
		
		double y, y1;
		double fact;
		double x1, xn;
		double ysq;
		double sum;
		double spi;
		bool flag01 = false;
		bool flag1_12 = false;
		bool flagNegative = false;
		
		// If x is negative
		if (x<0.0) {
			y = -x;
			y1 = floor(y);
			res = y - y1;
			fact = -pi / sin(pi*res) * (1 - 2*fmod(y1,2));
			x = y + 1;
			flagNegative = true;
		}
		// Now x is positive
		
		// Map x in interval [0,1] to [1,2]
		if (x<1.0) {
			x1 = x;
			x = x1 + 1.0;
			flag01 = true;
		}
		
		// Map x in interval [1,12] to [1,2]
		if (x<12.0) {
			xn = floor(x) - 1;
			x = x - xn;
			// Evaluate approximation for 1 < x < 2
			double z = x - 1.0;
			double xnum = 0.0;
			double xden = xnum + 1.0;
			
			
			for (int i = 1 ; i<=8; i++ ) {
				xnum = (xnum + p(i)) * z;
				xden = xden * z + q(i);
			}
			
			res = xnum / xden + 1.0;
			flag1_12 = true;
		}
		
		// Adjust result for case  0.0 < x < 1.0
		if (flag01)  {
			res = res / x1;
		}
		else if (flag1_12){   // Adjust result for case  2.0 < x < 12.0
			double max_xn = xn;
			for (int m=1;  m<=max_xn; m++) {
				res = res * x;
				x = x + 1;
				xn = xn - 1;
			}
		}
		
		// Evaluate approximation for x >= 12
		if (x>=12.0) {
			y = x;
			ysq = y * y;
			sum = c(7);
			for (int i = 1; i<=6; i++ ) {
				sum = sum / ysq + c(i);
			}
			
			spi = 0.9189385332046727417803297;
			sum = sum / y - y + spi;
			sum = sum + (y-0.5)*log(y);
			res = exp(sum);
		}
		
		if (flagNegative) {
			res = fact / res;
		}
	}
	
	return res;
}


double 
RandomVariable::incompleteGammaFunction(double a, double x)
{
	double b = x;
	if (x == 0.0) {
		b = 0.0;
	}
	if (a == 0.0) {
		b = 1.0;
	}
	// Series expansion for x < a+1
	if ( (a != 0.0) && (x != 0.0) && (x < a+1.0) ) {
		double ap = a;
		double del = 1.0;
		double sum = del;
		while (fabs(del) >= 1.0e-12*fabs(sum)) {
			ap = ap + 1.0;
			del = x * del / ap;
			sum = sum + del;
		}
		b = sum * exp(-x + a*log(x) - log(gammaFunction(a+1)) );
		if ( x > 0 && b > 1 )
			b = 1.0;
		
	} else {
		// Continued fraction for x >= a+1
		double a0 = 1.0;
		double a1 = x;
		double b0 = 0.0;
		double b1 = a0;
		double fac = 1.0/a1;
		double n = 1;
		double g = b1*fac;
		double gold = b0;
		while (fabs(g-gold) >= 1.0e-12*fabs(g))  {
			gold = g;
			double ana = n - a;
			a0 = (a1 + a0 *ana) * fac;
			b0 = (b1 + b0 *ana) * fac;
			double anf = n*fac;
			a1 = x * a0 + anf * a1;
			b1 = x * b0 + anf * b1;
			fac = 1.0 / a1;
			g = b1 * fac;
			n = n + 1.0;
		}
		b = 1.0 - exp(-x + a*log(x) - log(gammaFunction(a)) ) * g;
	}
	
	return b;
}


double 
RandomVariable::betaFunction(double q, double r)
{
	/*	OLD CODE: 
	 GammaRV *aGammaRV = new GammaRV(1, 0.0, 1.0, 0.0);
	 double par1,par2,par3;
	 par1 = aGammaRV->gammaFunction(q);
	 par2 = aGammaRV->gammaFunction(q);
	 par3 = aGammaRV->gammaFunction(q+r);
	 delete aGammaRV;
	 return par1*par2/par3;
	 */
	
	// Matlab definition of the beta function:
	//    y = exp(gammaln(q)+gammaln(r)-gammaln(q+r));
	//    ... where gammaln(.) = ln(gamma(.))
	// So, try this instead:
	double gammaq,gammar,gammaqpr;
	gammaq = gammaFunction(q);
	gammar = gammaFunction(r);
	gammaqpr = gammaFunction(q+r);
	
	double loggammaq,loggammar,loggammaqpr;
	loggammaq = log(gammaq);
	loggammar = log(gammar);
	loggammaqpr = log(gammaqpr);
	return exp(loggammaq+loggammar-loggammaqpr);
}

double 
RandomVariable::incompleteBetaFunction(double q, double r, double x) {
	// Regularized Incomplete Beta Function
	// Original implementation 
	// Lewis Van Winkle
	// http://CodePlea.com

	const double stop = 1e-8;
	const double tiny = 1e-30;

	if (x < 0.0 || x > 1.0) {
		return 0.0;
	}
	if (x > ((q + 1) / (q + r + 2.0))) {
		return 1.0 - incompleteBetaFunction(r, q, 1 - x);
	}

	const double lbeta_qr = log(gammaFunction(q)) + log(gammaFunction(r)) - log(gammaFunction(q + r));

	const double front = exp(log(x) * q + log(1.0 - x) * r - lbeta_qr) / q;

	// Lentz's algorithm to evaluate the continued fraction
	double f = 1.0, c = 1.0, d = 0.0;

	int i, m;
	for (i = 0; i <= 200; ++i) {
		m = i / 2;

		double numerator;
		if (i == 0) {
			numerator = 1.0;
		}
		else if (i % 2 == 0) {
			// Even terms
			numerator = (m * (r - m) * x) / ((q + 2.0 * m - 1.0) * (q + 2.0 * m));
		}
		else {
			// Odd terms
			numerator = -((q+m)*(q+r+m)*x) / ((q+2.0*m)*(q+2.0*m+1.0));
		}

		// Perform an iteration of Lentz algorithm
		d = 1.0 + numerator * d;
		if (fabs(d) < tiny) {
			d = tiny;
		}
		d = 1.0 / d;

		c = 1.0 + numerator / c;
		if (fabs(c) < tiny) {
			c = tiny;
		}

		const double cd = c * d;
		f *= cd;

		// check for stop
		if (fabs(1.0 - cd) < stop) {
			return front * (f - 1.0);
		}
	}

	return 0.0;
}


double 
RandomVariable::errorFunction(double x)
{
#ifdef _LINUX
	return erf(x);
#endif

	// ErrorFunction(x) = 2/sqrt(pi) * integral from 0 to x of exp(-t^2) dt.
	double a1,a2,a3,a4,a5;
	double b1,b2,b3,b4;
	double c1,c2,c3,c4,c5,c6,c7,c8,c9;
	double d1,d2,d3,d4,d5,d6,d7,d8;
	double p1,p2,p3,p4,p5,p6;
	double q1,q2,q3,q4,q5;
	double y,z,xnum,xden,del,result;
	//double pi = acos(-1.0);
	
	// Evaluate  errorFunction  for  |x| <= 0.46875
    if ( fabs(x) <= 0.46875 )
	{
		a1 = 3.16112374387056560e00;
		a2 = 1.13864154151050156e02;
		a3 = 3.77485237685302021e02; 
		a4 = 3.20937758913846947e03;
		a5 = 1.85777706184603153e-1;
		b1 = 2.36012909523441209e01; 
		b2 = 2.44024637934444173e02;
		b3 = 1.28261652607737228e03; 
		b4 = 2.84423683343917062e03;
		y = fabs(x);
		z = y * y;
		xnum = a5 * z;
		xden = z;
		xnum = (xnum + a1) * z;
		xden = (xden + b1) * z;
		xnum = (xnum + a2) * z;
		xden = (xden + b2) * z;
		xnum = (xnum + a3) * z;
		xden = (xden + b3) * z;
		result = x * (xnum + a4) / (xden + b4);
	}
	
	
	// Evaluate  errorFunction  for 0.46875 <= |x| <= 4.0
	else if ( fabs(x) > 0.46875 && fabs(x) <= 4.0 )
	{
		c1 = 5.64188496988670089e-1; 
		c2 = 8.88314979438837594e00;
		c3 = 6.61191906371416295e01; 
		c4 = 2.98635138197400131e02;
		c5 = 8.81952221241769090e02; 
		c6 = 1.71204761263407058e03;
		c7 = 2.05107837782607147e03; 
		c8 = 1.23033935479799725e03;
		c9 = 2.15311535474403846e-8;
		d1 = 1.57449261107098347e01; 
		d2 = 1.17693950891312499e02;
		d3 = 5.37181101862009858e02; 
		d4 = 1.62138957456669019e03;
		d5 = 3.29079923573345963e03; 
		d6 = 4.36261909014324716e03;
		d7 = 3.43936767414372164e03; 
		d8 = 1.23033935480374942e03;
		y = fabs(x);
		xnum = c9 * y;
		xden = y;
		xnum = (xnum + c1) * y;
		xden = (xden + d1) * y;
		xnum = (xnum + c2) * y;
		xden = (xden + d2) * y;
		xnum = (xnum + c3) * y;
		xden = (xden + d3) * y;
		xnum = (xnum + c4) * y;
		xden = (xden + d4) * y;
		xnum = (xnum + c5) * y;
		xden = (xden + d5) * y;
		xnum = (xnum + c6) * y;
		xden = (xden + d6) * y;
		xnum = (xnum + c7) * y;
		xden = (xden + d7) * y;
		result = (xnum + c8) / (xden + d8);
		z = floor(y*16)/16;
		del = (y-z) * (y+z);
		result = exp(-z * z) * exp(-del) * result;
	}
	
	
	// Evaluate  erfc  for |x| > 4.0
	else if ( fabs(x) > 4.0 )
	{
		p1 = 3.05326634961232344e-1; 
		p2 = 3.60344899949804439e-1;
		p3 = 1.25781726111229246e-1; 
		p4 = 1.60837851487422766e-2;
		p5 = 6.58749161529837803e-4; 
		p6 = 1.63153871373020978e-2;
		q1 = 2.56852019228982242e00; 
		q2 = 1.87295284992346047e00;
		q3 = 5.27905102951428412e-1; 
		q4 = 6.05183413124413191e-2;
		q5 = 2.33520497626869185e-3;
		y = fabs(x);
		z = 1 / (y * y);
		xnum = p6 * z;
		xden = z;
		xnum = (xnum + p1) * z;
		xden = (xden + q1) * z;
		xnum = (xnum + p2) * z;
		xden = (xden + q2) * z;
		xnum = (xnum + p3) * z;
		xden = (xden + q3) * z;
		xnum = (xnum + p4) * z;
		xden = (xden + q4) * z;
		result = z * (xnum + p5) / (xden + q5);
		result = (1/sqrt(pi) -  result) / y;
		z = floor(y*16)/16;
		del = (y-z) * (y+z);
		result = exp(-z * z) * exp(-del) * result;
		if ( isnan(result) )
			result = 0;
	}
	
	else {
		// should never reach here unless the input argument was nan
		opserr << "RandomVariable::errorFunction WARNING: invalid input (" << x << ")" << endln;
		result = 0;
	}
	
	
	// Final computations
    if ( x > 0.46875 )
		result = (0.5 - result) + 0.5;
	if ( x < -0.46875 )
		result = (-0.5 + result) - 0.5;
	
	return result;
}


double 
RandomVariable::inverseErrorFunction(double y)
{
	double a1,a2,a3,a4;
	double b1,b2,b3,b4;
	double c1,c2,c3,c4;
	double d1,d2;
	double x,z;
	
	// Coefficients in rational approximations.
	a1 = 0.886226899; 
	a2 = -1.645349621;  
	a3 = 0.914624893;
	a4 = 0.140543331;
	b1 = -2.118377725;
	b2 = 1.442710462;
	b3 = -0.329097515;
	b4 = 0.012229801;
	c1 = -1.970840454;
	c2 = -1.624906493;
	c3 = 3.429567803;
	c4 = 1.641345311;
	d1 = 3.543889200;  
	d2 = 1.637067800;
	
	// Central range
	if ( fabs(y) <= 0.7 )
	{
		z = y * y;
		x = y * (((a4*z+a3)*z+a2)*z+a1) / ((((b4*z+b3)*z+b2)*z+b1)*z+1);
	}
	
	// Near end-points of range
	else if ( y > 0.7  &&  y < 1 )
	{
		z = sqrt(-log((1-y)/2));
		x = (((c4*z+c3)*z+c2)*z+c1) / ((d2*z+d1)*z+1);
	}
	
	else if ( y < -0.7  &&  y > -1 )
	{
		z = sqrt(-log((1+y)/2));
		x = -(((c4*z+c3)*z+c2)*z+c1) / ((d2*z+d1)*z+1);
	}
	
	// domain input errors: either abs(y) >= 1, or nan was input
	else 
	{
		// this is the -1 limit
		if ( y < 0 )
			x = -DBL_MAX;
		
		// this is the +1 limit
		else if ( y > 0 )
			x = DBL_MAX;
		
		else {
			opserr << "RandomVariable::inverseErrorFunction WARNING: input (" << y << ") outside of [-1,1] domain." << endln;
			x = 0;
		}
		
		// do not perform Newton on these boundary values
		return x;
	}
	
	//static const double pi = acos(-1.0);
	static const double twoOverRootPi = 2.0/sqrt(pi);
	// Two steps of Newton-Raphson correction to full accuracy.
	// Without these steps, erfinv(y) would be about 3 times
	// faster to compute, but accurate to only about 6 digits.
	//x = x - (errorFunction(x) - y) / (2.0/sqrt(pi) * exp(-pow(x,2.0)));
	//x = x - (errorFunction(x) - y) / (2.0/sqrt(pi) * exp(-pow(x,2.0)));
	x = x - (errorFunction(x) - y) / (twoOverRootPi * exp(-pow(x,2.0)));
	x = x - (errorFunction(x) - y) / (twoOverRootPi * exp(-pow(x,2.0)));
	
	return x;
	
}


double
RandomVariable::standardNormalPhi(double uValue)
{
	static const double oneOverRootTwo = 1.0/sqrt(2.0);
    
	//Phi(x) = 0.5 * erfc(-x/sqrt(2))
	double result = 0.5 * (1.0 + errorFunction( uValue*oneOverRootTwo ));
    
	return result;
}


double
RandomVariable::standardNormalInversePhi(double probValue)
{
	double trval = probValue;
	if (trval <= 0.0) {
        // technically this should return negative infinity
		//opserr << "WARNING: Invalid probability value (" << trval << ") input <= 0 to RandomVariable::standardNormalInversePhi()" << endln;
		trval = 2.0*DBL_EPSILON;
	} else if (trval >= 1.0) {
        // technically this should return positive infinity
		//opserr << "WARNING: Invalid probability value (" << trval << ") input >= 1 to RandomVariable::standardNormalInversePhi()" << endln;
		trval = 1.0-2.0*DBL_EPSILON;
	}
	static const double rootTwo = sqrt(2.0);
	double result = rootTwo * inverseErrorFunction(2.0*trval-1.0);
	return result;
}


double 
RandomVariable::harmonicNumber(double n)
{
	double Hn;
	
	if (n > 1.25) {
		// asymptotic harmonic number series approximation
		Hn = log(n) + euler + 1/(2*n) - 1/(12*n*n) + 1/120/pow(n,4) - 1/252/pow(n,6);
	}
    else if (n > 0.8) {
        // Taylor series expansion about n = 1
        Hn = 1+1/6*(pi*pi-6)*(n-1) + (1-zeta3)*pow(n-1,2) + (pow(pi,4)/90-1)*pow(n-1,3) +
        (1-zeta5)*pow(n-1,4);
    }
	else if (n > 0.2) {
		// Taylor series expansion about n = 1/2
		Hn = 2+2*pow(1-2*n,2) + (-4+pi*pi/2)*(n-0.5) + (-16+pow(pi,4)/6)*pow(n-0.5,3) + 
        32*pow(n-0.5,4) - log(4.0) - 7/4*pow(1-2*n,2)*zeta3 - 31*pow(n-0.5,4)*zeta5;
	}
	else {
		// Taylor series expansion about the origin
		Hn = pow(pi,2)*n/6 + pow(pi,4)*pow(n,3)/90 - zeta3*pow(n,2) - zeta5*pow(n,4);
	}
	
	return Hn;
}
