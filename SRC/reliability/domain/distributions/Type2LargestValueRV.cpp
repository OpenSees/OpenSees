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
// $Source: /usr/local/cvs/OpenSees/SRC/reliability/domain/distributions/Type2LargestValueRV.cpp,v $


//
// Written by Terje Haukaas (haukaas@ce.berkeley.edu)
//

#include <Type2LargestValueRV.h>
#include <Vector.h>
#include <cmath>

Type2LargestValueRV::Type2LargestValueRV(int passedTag, 
					 double passedMean, double passedStdv)
:RandomVariable(passedTag, RANDOM_VARIABLE_type2largestvalue)
{
	int setp = setParameters(passedMean,passedStdv);
	if (setp < 0)
		opserr << "Error setting parameters in Type2LargestValue RV with tag " << this->getTag() << endln;
}


Type2LargestValueRV::Type2LargestValueRV(int passedTag,
					 const Vector &passedParameters)
:RandomVariable(passedTag, RANDOM_VARIABLE_type2largestvalue)
{
	
	if (passedParameters.Size() != 2) {
		opserr << "Type2LargestValue RV requires 2 parameters, u and k, for RV with tag " <<
		this->getTag() << endln;
		
		// this will create terminal errors
		u = 0;
		k = 0;
		
	} else {
		
		u = passedParameters(0);
		k = passedParameters(1);
	}
}


Type2LargestValueRV::~Type2LargestValueRV()
{
}


const char *
Type2LargestValueRV::getType()
{
	return "TYPE2LARGESTVALUE";
}


double 
Type2LargestValueRV::getMean()
{
	double result;
	result = u * gammaFunction(1.0-1.0/k);
	return result;
}


double 
Type2LargestValueRV::getStdv()
{
	double a = gammaFunction(1.0-2.0/k);
	double b = gammaFunction(1.0-1.0/k);
	double result = u*sqrt(a-b*b);
	return result;
}


const Vector &
Type2LargestValueRV::getParameters(void) {
	static Vector temp(2);
	temp(0) = u;
	temp(1) = k;
	return temp;
}


int
Type2LargestValueRV::setParameters(double mean, double stdv)
{
	double cov = stdv/mean;
	double c = 1+cov*cov;
	double k_prev = log(1/(c-1));
	double del = 1.0;
    double del_prev = 10.0;
	
	// now use Newtons method with nice f/f' function
	int ncount = 1;
	int nmax = 100;
    double tolk = 1.0e-8;
	
	while (del > tolk && ncount <= nmax) {
		k = k_prev - (-c*pow(gammaFunction(-1/k_prev),2)-2*k_prev*gammaFunction(-2/k_prev)) / 
            ( 2*gammaFunction(-2/k_prev+1)*(-harmonicNumber(-1/k_prev)+harmonicNumber(-2/k_prev)) );
		del = fabs(k-k_prev);
        
        if (del > del_prev) {
            // be careful with non-convergent cases, exit to bisection method
            break;
        }
        
		k_prev = k;
        del_prev = del;
		ncount++;
	}
    
    // check whether we need to use bisection method to supplement
    if (del > tolk) {
        double left = 2.05;
        double right = 1000;
        int change = 0;
        double fa = 0;
        double fc = 0;
        ncount = 1;
        
        while (change == 0 && ncount <= nmax) {
            k = (left+right)/2.0;
            fc = -c + gammaFunction(1-2/k)/gammaFunction(1-1/k)/gammaFunction(1-1/k);
            fa = -c + gammaFunction(1-2/left)/gammaFunction(1-1/left)/gammaFunction(1-1/left);
            if (fabs(fc) < tolk || (right-left)/2.0 < tolk) {
                change = 1;
                del = tolk;
            }
            
            if ( fa>0 && fc>0 )
                left = k;
            else if ( fa<0 && fc<0 )
                left = k;
            else
                right = k;
            
            ncount++;
        }
    }
    
    // unlikely but possible that k was not found using either method above
    if (del > tolk) {
		opserr << "Warning: Type2Largest distribution did not converge during setParameters()" << endln;
        return -1;
    }
	
	u = mean/gammaFunction(1-1/k);
	
	return 0;
}


double
Type2LargestValueRV::getPDFvalue(double rvValue)
{
	double result;
	if ( 0.0 < rvValue ) {
		result = k/u * pow(u/rvValue,k+1) * exp(-pow(u/rvValue,k));
	}
	else {
		result = 0.0;
	}
	return result;
}


double
Type2LargestValueRV::getCDFvalue(double rvValue)
{
	double result;
	if ( 0.0 < rvValue ) {
		result = exp(-pow(u/rvValue,k));
	}
	else {
		result = 0.0;
	}
	return result;
}


double
Type2LargestValueRV::getInverseCDFvalue(double probValue)
{
	return u * pow((-log(probValue)),(-1.0/k));
}


void
Type2LargestValueRV::Print(OPS_Stream &s, int flag)
{
	s << "Type2LargestValue RV #" << this->getTag() << endln;
	s << "\tu = " << u << endln;
	s << "\tk = " << k << endln;
}
