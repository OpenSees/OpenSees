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
// $Source: /usr/local/cvs/OpenSees/SRC/reliability/domain/distributions/Type3SmallestValueRV.cpp,v $


//
// Written by Terje Haukaas (haukaas@ce.berkeley.edu)
//

#include <Type3SmallestValueRV.h>
#include <Vector.h>
#include <cmath>

Type3SmallestValueRV::Type3SmallestValueRV(int passedTag, 
					   const Vector &passedParameters)
:RandomVariable(passedTag, RANDOM_VARIABLE_type3smallestvalue)
{
	
	if (passedParameters.Size() != 3) {
		opserr << "Type3SmallestValue RV requires 3 parameters, epsilon and u and k, for RV with tag " <<
		this->getTag() << endln;
		
		// this will create terminal errors
		epsilon = 0;
		u = 0;
		k = 0;
		
	} else {
		
		epsilon	= passedParameters(0);
		u = passedParameters(1);
		k = passedParameters(2);
		
	}
}


Type3SmallestValueRV::~Type3SmallestValueRV()
{
}


const char *
Type3SmallestValueRV::getType()
{
	return "TYPE3SMALLESTVALUE";
}


double 
Type3SmallestValueRV::getMean()
{
	double result;
	result = epsilon + (u-epsilon) * gammaFunction(1.0+1.0/k);
	return result;
}



double 
Type3SmallestValueRV::getStdv()
{
	double a = gammaFunction(1.0+2.0/k);
	double b = gammaFunction(1.0+1.0/k);
	return (u-epsilon) * sqrt(a-b*b);
}


const Vector &
Type3SmallestValueRV::getParameters(void) {
	static Vector temp(3);
	temp(0) = epsilon;
	temp(1) = u;
	temp(2) = k;
	return temp;
}


double
Type3SmallestValueRV::getPDFvalue(double rvValue)
{
	double result;
	if ( epsilon < rvValue ) {
		result = k/(u-epsilon) * pow(((rvValue-epsilon)/(u-epsilon)),k-1) * exp(-pow(((rvValue-epsilon)/(u-epsilon)),k));
	}
	else {
		result = 0.0;
	}
	return result;
}


double
Type3SmallestValueRV::getCDFvalue(double rvValue)
{
	double result;
	if ( epsilon < rvValue ) {
		result = 1 - exp(-pow(((rvValue-epsilon)/(u-epsilon)),k));
	}
	else {
		result = 0.0;
	}
	return result;
}


double
Type3SmallestValueRV::getInverseCDFvalue(double probValue)
{
	return (u-epsilon) * ( epsilon/(u-epsilon) + pow((-log(1-probValue)),(1/k)) );
}


void
Type3SmallestValueRV::Print(OPS_Stream &s, int flag)
{
	s << "Type3SmallestValue RV #" << this->getTag() << endln;
	s << "\tepsilon = " << epsilon << endln;
	s << "\tu = " << u << endln;
	s << "\tk = " << k << endln;
}


