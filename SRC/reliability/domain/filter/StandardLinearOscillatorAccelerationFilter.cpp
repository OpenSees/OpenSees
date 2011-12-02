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
// $Date: 2008-03-13 22:36:37 $
// $Source: /usr/local/cvs/OpenSees/SRC/reliability/domain/filter/StandardLinearOscillatorAccelerationFilter.cpp,v $


//
// Written by Terje Haukaas (haukaas@ce.berkeley.edu)
//

#include <StandardLinearOscillatorAccelerationFilter.h>
#include <Filter.h>
#include <classTags.h>
#include <math.h>

StandardLinearOscillatorAccelerationFilter::StandardLinearOscillatorAccelerationFilter(int tag, double period, double dampingRatio)
:Filter(tag,FILTER_standardLinearOscillator)
{
	double pi = 3.14159265358979;
	wn = 2*pi/period;
	xi = dampingRatio;
}

StandardLinearOscillatorAccelerationFilter::~StandardLinearOscillatorAccelerationFilter()
{
}

double
StandardLinearOscillatorAccelerationFilter::getAmplitude(double time, double dT)
{
	if (time<0.0) {
		return 0.0;
	}
	else {
		double wd = wn * sqrt(1.0-pow(xi,2.0));
		return (  ( xi*xi*wn*wn*sin(wd*time) - 2.0*xi*wn*wd*cos(wd*time) - wd*wd*sin(wd*time) ) * exp(-xi*wn*time)  );
	}
}

double
StandardLinearOscillatorAccelerationFilter::getMaxAmplitude()
{
	double wd = wn * sqrt(1.0-pow(xi,2.0));

	opserr << "ERROR: The getMaxAmplitude() method is not implemented for acceleration filter." << endln;

	double result = 0.0;

	return result;
}

double
StandardLinearOscillatorAccelerationFilter::getTimeOfMaxAmplitude()
{
	double wd = wn * sqrt(1.0-pow(xi,2.0));

	opserr << "ERROR: The getTimeOfMaxAmplitude() method is not implemented for acceleration filter." << endln;

	return 0.0;
}

void
StandardLinearOscillatorAccelerationFilter::Print(OPS_Stream &s, int flag)  
{
}
