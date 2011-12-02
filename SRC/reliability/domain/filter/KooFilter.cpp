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
// $Date: 2008-03-13 22:36:37 $
// $Source: /usr/local/cvs/OpenSees/SRC/reliability/domain/filter/KooFilter.cpp,v $


//
// Written by Terje Haukaas (haukaas@ce.berkeley.edu)
//

#include <math.h>
#include <KooFilter.h>
#include <Filter.h>
#include <classTags.h>


KooFilter::KooFilter(int tag, double period, double dampingRatio)
:Filter(tag,FILTER_standardLinearOscillator)
{
	double pi = 3.14159265358979;
	wn = 2*pi/period;
	xi = dampingRatio;
}

KooFilter::~KooFilter()
{
}

double
KooFilter::getAmplitude(double time, double dT)
{
	if (time<0.0) {

		return 0.0;

	}
	else {

		double mysqrt = sqrt(1.0-pow(xi,2.0));

		double wd = wn * mysqrt;

		double term1 = wn/mysqrt * sin(wd*time);

		double term2 = 2.0*xi*wn*cos(wd*time);

		double theWholeThing = -(term1+term2)*exp(-xi*wn*time);

		return theWholeThing;
	}
}

double
KooFilter::getMaxAmplitude()
{
	double wd = wn * sqrt(1.0-pow(xi,2.0));

	double result = wd/(xi*wn*sqrt((xi*xi*wn*wn+wd*wd)/(xi*xi*wn*wn)))
		*exp(-xi*wn*(atan(wd/(xi*wn))/wd));

	return result;
}

double
KooFilter::getTimeOfMaxAmplitude()
{

	double wd = wn * sqrt(1.0-pow(xi,2.0));

	return (atan(wd/(xi*wn))/wd);
}

void
KooFilter::Print(OPS_Stream &s, int flag)  
{
}
