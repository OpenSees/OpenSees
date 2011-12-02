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
** ******  ************************************************************ */

// $Revision: 1.2 $
// $Date: 2008-03-13 22:25:26 $
// $Source: /usr/local/cvs/OpenSees/SRC/reliability/analysis/telm/NewStandardLinearOscillatorAccelerationFilter.cpp,v $

#include <NewStandardLinearOscillatorAccelerationFilter.h>
#include <Filter.h>
#include <classTags.h>
#include <math.h>

NewStandardLinearOscillatorAccelerationFilter::NewStandardLinearOscillatorAccelerationFilter(int tag, double period, double dampingRatio, double pdtpulse)
:Filter(tag,FILTER_standardLinearOscillator)
{
	double pi = 4.0*atan(1.0);
//3.14159265358979;
	if(period!=0.0){
		wn = 2*pi/period;
		xi = dampingRatio;
        wd = wn * sqrt(1.0-pow(xi,2.0));
		coefc = -2.0*xi*wn;
		coefs = (xi*xi*wn*wn-wd*wd)/wd;
	}else{
		wn = 0.0;
		xi = 0.0;
		wd = 0.0;
		coefc = 1.0;
		coefs = 1.0;
	}
    dtpulse=pdtpulse;
}

NewStandardLinearOscillatorAccelerationFilter::~NewStandardLinearOscillatorAccelerationFilter()
{
}

double
NewStandardLinearOscillatorAccelerationFilter::getAmplitude(double time, double dt)
{
/*	if (time<0.0) {
		return 0.0;
	}
	else {
		if(wd!=0.0){
			return (  ( coefs*sin(wd*time) + coefc*cos(wd*time) ) * exp(-xi*wn*time)  );
		}else{
			if(time < 1.0e-8) {
				return 1.0;
			}else{
				return 0.0;
			}
		}
	}
*/
	if (time<-1.0e-7) {
		return 0.0;
	}else{
		double t1=time;
	    if(t1<0.0) t1=0.0;
	    double t2=t1+dtpulse;
        double sss1=exp(-xi*wn*t1)*((xi*wn/wd)*sin(wd*t1)-cos(wd*t1));
        double sss2=exp(-xi*wn*t2)*((xi*wn/wd)*sin(wd*t2)-cos(wd*t2));
        double ddd=sss2-sss1;
		return ddd/dtpulse;
	}
}

double
NewStandardLinearOscillatorAccelerationFilter::getMaxAmplitude()
{
	double wd = wn * sqrt(1.0-pow(xi,2.0));

	opserr << "ERROR: The getMaxAmplitude() method is not implemented for acceleration filter." << endln;

	double result = 0.0;

	return result;
}

double
NewStandardLinearOscillatorAccelerationFilter::getTimeOfMaxAmplitude()
{
	double wd = wn * sqrt(1.0-pow(xi,2.0));

	opserr << "ERROR: The getTimeOfMaxAmplitude() method is not implemented for acceleration filter." << endln;

	return 0.0;
}

void
NewStandardLinearOscillatorAccelerationFilter::Print(OPS_Stream &s, int flag)  
{
}
