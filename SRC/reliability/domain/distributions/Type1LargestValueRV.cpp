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
// $Date: 2007-02-17 21:27:23 $
// $Source: /usr/local/cvs/OpenSees/SRC/reliability/domain/distributions/Type1LargestValueRV.cpp,v $


//
// Written by Terje Haukaas (haukaas@ce.berkeley.edu) 
//

#include <Type1LargestValueRV.h>
#include <Vector.h>
#include <cmath>

Type1LargestValueRV::Type1LargestValueRV(int passedTag, 
					 double passedMean, double passedStdv)
:RandomVariable(passedTag, RANDOM_VARIABLE_type1largestvalue)
{
	int setp = setParameters(passedMean,passedStdv);
	if (setp < 0)
		opserr << "Error setting parameters in Type1LargestValue RV with tag " << this->getTag() << endln;
}


Type1LargestValueRV::Type1LargestValueRV(int passedTag,
					 const Vector &passedParameters)
:RandomVariable(passedTag, RANDOM_VARIABLE_type1largestvalue)
{
	
	if (passedParameters.Size() != 2) {
		opserr << "Type1LargestValue RV requires 2 parameters, u and alpha, for RV with tag " <<
		this->getTag() << endln;
		
		// this will create terminal errors
		u = 0;
		alpha = 0;
		
	} else {
		
		u = passedParameters(0);
		alpha = passedParameters(1);
		
	}
}


Type1LargestValueRV::~Type1LargestValueRV()
{
}


const char *
Type1LargestValueRV::getType()
{
	return "TYPE1LARGESTVALUE";
}


double 
Type1LargestValueRV::getMean()
{
	//double gamma = 0.5772156649;
	return u+euler/alpha;
}


double 
Type1LargestValueRV::getStdv()
{
	//double pi = 3.14159265358979;
	return pi/(sqrt(6.0)*alpha);
}


const Vector &
Type1LargestValueRV::getParameters(void) {
	static Vector temp(2);
	temp(0) = u;
	temp(1) = alpha;
	return temp;
}


int
Type1LargestValueRV::setParameters(double mean, double stdv)
{
	//static const double pi = std::acos(-1.0);
	//double eulergamma = 0.57721566490153286061;
	
	u = mean - (euler*sqrt(6.0)*stdv)/pi;
	alpha = pi / ( sqrt(6.0) * stdv );
	
	return 0;
}


double
Type1LargestValueRV::getPDFvalue(double rvValue)
{
	return alpha * exp(-alpha*(rvValue-u)-exp(-alpha*(rvValue-u)));
}


double
Type1LargestValueRV::getCDFvalue(double rvValue)
{
	return exp(-exp(-alpha*(rvValue-u)));
}


double
Type1LargestValueRV::getInverseCDFvalue(double probValue)
{
	return (alpha*u - log(-log(probValue))) / alpha;
}


void
Type1LargestValueRV::Print(OPS_Stream &s, int flag)
{
	s << "Type1LargestValue RV #" << this->getTag() << endln;
	s << "\tu = " << u << endln;
	s << "\talpha = " << alpha << endln;
}
