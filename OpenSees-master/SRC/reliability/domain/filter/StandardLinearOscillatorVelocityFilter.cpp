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
                                                                        
// $Revision: 1.4 $
// $Date: 2008-03-13 22:36:38 $
// $Source: /usr/local/cvs/OpenSees/SRC/reliability/domain/filter/StandardLinearOscillatorVelocityFilter.cpp,v $


//
// Written by Terje Haukaas (haukaas@ce.berkeley.edu)
//

#include <StandardLinearOscillatorVelocityFilter.h>
#include <Filter.h>
#include <classTags.h>
#include <math.h>

StandardLinearOscillatorVelocityFilter::StandardLinearOscillatorVelocityFilter(int tag, double period, double dampingRatio)
:Filter(tag,FILTER_standardLinearOscillator)
{
	double pi = 3.14159265358979;
	wn = 2*pi/period;
	xi = dampingRatio;
}

StandardLinearOscillatorVelocityFilter::~StandardLinearOscillatorVelocityFilter()
{
}

double
StandardLinearOscillatorVelocityFilter::getAmplitude(double time, double dT)
{
	if (time<0.0) {
		return 0.0;
	}
	else {
		double wd = wn * sqrt(1.0-pow(xi,2.0));
		return (  ( wd*cos(wd*time) - xi*wn*sin(wd*time) ) * exp(-xi*wn*time)  );
		// Should maybe include the 1/mwd factor too
	}
}

double
StandardLinearOscillatorVelocityFilter::getMaxAmplitude()
{
	double wd = wn * sqrt(1.0-pow(xi,2.0));

	opserr << "ERROR: The getMaxAmplitude() method is not implemented for velocity filter." << endln;

	double result = 0.0;

	return result;
}

double
StandardLinearOscillatorVelocityFilter::getTimeOfMaxAmplitude()
{
	double wd = wn * sqrt(1.0-pow(xi,2.0));

	opserr << "ERROR: The getTimeOfMaxAmplitude() method is not implemented for velocity filter." << endln;

	return 0.0;
}

void
StandardLinearOscillatorVelocityFilter::Print(OPS_Stream &s, int flag)  
{
}
