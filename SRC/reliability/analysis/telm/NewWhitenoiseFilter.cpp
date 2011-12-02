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
// $Date: 2008-03-13 22:25:26 $
// $Source: /usr/local/cvs/OpenSees/SRC/reliability/analysis/telm/NewWhitenoiseFilter.cpp,v $

#include <NewWhitenoiseFilter.h>
#include <Filter.h>
#include <classTags.h>


NewWhitenoiseFilter::NewWhitenoiseFilter(int tag, double pdelta)
:Filter(tag,FILTER_standardLinearOscillator)
{
	dt_analysis=pdelta;
	cdt_analysis=1.0/dt_analysis;

}

NewWhitenoiseFilter::~NewWhitenoiseFilter()
{
}

double
NewWhitenoiseFilter::getAmplitude(double time, double dt)
{
	if (time<-1.0e-7) {
		return 0.0;
//	}else if( time < dt_analysis ){
	}else if( time < 1.0e-5){
//		return 1.0;
		return cdt_analysis;

	}else {
		return 0.0;
	}
}

double
NewWhitenoiseFilter::getMaxAmplitude()
{
//	double wd = wn * sqrt(1.0-pow(xi,2.0));

	opserr << "ERROR: The getMaxAmplitude() method is not implemented for acceleration filter." << endln;

	double result = 0.0;

	return result;
}

double
NewWhitenoiseFilter::getTimeOfMaxAmplitude()
{
//	double wd = wn * sqrt(1.0-pow(xi,2.0));

	opserr << "ERROR: The getTimeOfMaxAmplitude() method is not implemented for acceleration filter." << endln;

	return 0.0;
}

void
NewWhitenoiseFilter::Print(OPS_Stream &s, int flag)  
{
}
