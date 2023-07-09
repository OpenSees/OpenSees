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
                                                                        
// $Revision: 1.2 $
// $Date: 2007-02-17 21:27:23 $
// $Source: /usr/local/cvs/OpenSees/SRC/reliability/domain/distributions/UserDefinedRV.cpp,v $


//
// Written by Terje Haukaas (haukaas@ce.berkeley.edu) 
//

#include <UserDefinedRV.h>
#include <Vector.h>
#include <cmath>

UserDefinedRV::UserDefinedRV(int passedTag, const Vector &pxPoints,
			     const Vector &pPDFpoints)
  :RandomVariable(passedTag, RANDOM_VARIABLE_userdefined),
   xPoints(pxPoints), PDFpoints(pPDFpoints)
{
	// This used to be in TclReliabilityBuilder
	int numPoints = xPoints.Size();
	
	// Normalize the PDF
	double sum = 0.0;
	for ( int i=1; i<numPoints; i++ ) {
		sum += 0.5 * (PDFpoints(i)+PDFpoints(i-1)) * (xPoints(i)-xPoints(i-1));
	}
	double diff;
	if ( fabs(sum-1.0) < 1.0e-6) {
		// It's normalized enough...!
	}
	else if (sum < 1.0) {
		diff = (1.0-sum)/(xPoints(numPoints-1)-xPoints(0));
		opserr << "WARNING: The PDF of random variable " << this->getTag() << " is normalized by " << endln
		<< "         uniformly increasing the PDF values by " << diff << endln;
		for (int i=0; i<numPoints; i++) {
			PDFpoints(i) = PDFpoints(i) + diff;
		}
	}
	else {
		diff = (sum-1.0)/(xPoints(numPoints-1)-xPoints(0));
		opserr << "WARNING: The PDF of random variable " << this->getTag() << " is normalized by " << endln
		<< "         uniformly decreasing the PDF values by " << diff << endln;
		for (int i=0; i<numPoints; i++) {
			PDFpoints(i) = PDFpoints(i) - diff;
		}
	}
	
	// Check if the PDF became negative somewhere
	for (int i=0; i<numPoints; i++) {
		if ( PDFpoints(i) < 0.0 ) {
			opserr << "ERROR: Some PDF values of random variable " << this->getTag() << endln
			<< "became negative after normalization. " << endln;
		}
	}
	
	// Check that it has been normalized
	sum = 0.0;
	for ( int i=1; i<numPoints; i++ ) {
		sum += 0.5 * (PDFpoints(i)+PDFpoints(i-1)) * (xPoints(i)-xPoints(i-1));
	}
	if ( fabs(1.0-sum) > 1.0e-6 ) {
		opserr << "ERROR: Did not succeed in normalizing user-defined distribution." << endln;
	}
	
}


UserDefinedRV::~UserDefinedRV()
{

}


const char *
UserDefinedRV::getType()
{
	return "USERDEFINED";
}


double 
UserDefinedRV::getMean()
{
	double sum = 0.0;
	double a, b;
	for ( int i=1; i<xPoints.Size(); i++ ) {
		a = (PDFpoints(i)-PDFpoints(i-1))/(xPoints(i)-xPoints(i-1));
		b = PDFpoints(i-1) - a * xPoints(i-1);
		sum += a*(xPoints(i))*(xPoints(i))*(xPoints(i))/3.0
			+  0.5*b*(xPoints(i))*(xPoints(i))
			-  a*(xPoints(i-1))*(xPoints(i-1))*(xPoints(i-1))/3.0
			-  0.5*b*(xPoints(i-1))*(xPoints(i-1));
	}
	return sum;
}


double 
UserDefinedRV::getStdv()
{
	double sum = 0.0;
	double a, b;
	double mu = getMean();
	for ( int i=1; i<xPoints.Size(); i++ ) {
		a = (PDFpoints(i)-PDFpoints(i-1))/(xPoints(i)-xPoints(i-1));
		b = PDFpoints(i-1) - a * xPoints(i-1);
		double x1 = xPoints(i-1);
		double x2 = xPoints(i);
		sum += 0.25*a*x2*x2*x2*x2
			-  2.0/3.0*mu*a*x2*x2*x2
			+ mu*mu*b*x2
			- mu*b*x2*x2
			+ b*x2*x2*x2/3.0
			+ 0.5*mu*mu*a*x2*x2
			- 0.25*a*x1*x1*x1*x1
			+ 2.0/3.0*mu*a*x1*x1*x1
			- mu*mu*b*x1
			+ mu*b*x1*x1
			- b*x1*x1*x1/3.0
			- 0.5*mu*mu*a*x1*x1;
	}
	return sqrt(sum);
}


void
UserDefinedRV::Print(OPS_Stream &s, int flag)
{
}


double
UserDefinedRV::getPDFvalue(double x)
{

	double ok = -1.0;
	if ( x < xPoints(0) ) {
		return 0.0;
	}
	else if ( x > xPoints(xPoints.Size()-1) ) {
		return 0.0;
	}
	else {
		for ( int i=1; i<xPoints.Size(); i++ ) {
			if ( x <= xPoints(i) ) {
				double a = (PDFpoints(i)-PDFpoints(i-1))/(xPoints(i)-xPoints(i-1));
				ok = ( a * (x-xPoints(i-1)) + PDFpoints(i-1) );
				break;
			}
		}
		if (ok<0.0) {
			opserr << "ERROR: UserDefinedRV::getPDFvalue() -- distribution function seems to be invalid." << endln;
		}
		return ok;
	}
}


double
UserDefinedRV::getCDFvalue(double x)
{
	double ok = -1.0;
	if ( x < xPoints(0) ) {
		return 0.0;
	}
	else if ( x > xPoints(xPoints.Size()-1) ) {
		return 1.0;
	}
	else {
		double sum = 0.0;
		for ( int i=1; i<xPoints.Size(); i++ ) {
			if ( x <= xPoints(i) ) {
				
				ok = ( sum + 0.5 * (getPDFvalue(x)+PDFpoints(i-1)) * (x-xPoints(i-1)) );
				break;
			}
			sum += 0.5 * (PDFpoints(i)+PDFpoints(i-1)) * (xPoints(i)-xPoints(i-1));
		}
		if (ok<0.0) {
			opserr << "ERROR: UserDefinedRV::getPDFvalue() -- distribution function seems to be invalid." << endln;
		}
		return ok;
	}
}


double
UserDefinedRV::getInverseCDFvalue(double p)
{
	if (p < 0.0 || p > 1.0) {
		opserr << "ERROR: Illegal probability value input to UserDefinedRV::getInverseCDFvalue()" << endln;
		return 0.0;
	}
	else if (p == 0.0) {
		return (xPoints(0));
	}


	double ok;
	for ( int i=1; i<xPoints.Size(); i++ ) {
		if ( getCDFvalue(xPoints(i)) > p) {
			double a = (getPDFvalue(xPoints(i))-getPDFvalue(xPoints(i-1)))/(xPoints(i)-xPoints(i-1));
			if ( a==0.0 && getPDFvalue(xPoints(i-1))==0.0 ) {
				opserr << "ERROR: An inside region of PDF is constant zero..." << endln;
			}
			else if (a==0.0) {
				ok = (p-getCDFvalue(xPoints(i-1))+getPDFvalue(xPoints(i-1))*xPoints(i-1))/getPDFvalue(xPoints(i-1));
			}
			else {

				double A = 0.5*a;
				double B = getPDFvalue(xPoints(i-1));
				double C = getCDFvalue(xPoints(i-1)) - p;
				double x_minus_x_i_1 = (-B+sqrt(B*B-4.0*A*C))/(2.0*A);
				ok = x_minus_x_i_1 + (xPoints(i-1));
			}
			break;
		}
	}
	return ok;
}


