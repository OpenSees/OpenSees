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
// $Date: 2003-10-27 23:04:40 $
// $Source: /usr/local/cvs/OpenSees/SRC/reliability/domain/modulatingFunction/TrapezoidalModulatingFunction.cpp,v $


//
// Written by Terje Haukaas (haukaas@ce.berkeley.edu)
//

#include <TrapezoidalModulatingFunction.h>
#include <ModulatingFunction.h>
#include <classTags.h>


TrapezoidalModulatingFunction::TrapezoidalModulatingFunction(int tag,
												 Filter *theFilt, 
												 double pt1,
												 double pt2,
												 double pt3,
												 double pt4,
												 double pamplitude)
:ModulatingFunction(tag,MODULATING_FUNCTION_trapezoidal)
{
	t1 = pt1;
	t2 = pt2;
	t3 = pt3;
	t4 = pt4;
	amplitude = pamplitude;

	if (t1>t2 || t2>t3 || t3>t4) {
		opserr << "WARNING: Inconsistent input to Trapezoidal Modulating Function" << endln;
	}

	theFilter = theFilt;
}

TrapezoidalModulatingFunction::~TrapezoidalModulatingFunction()
{
}

double
TrapezoidalModulatingFunction::getAmplitude(double time)
{
	if (time < t1) {
		return 0.0;
	}
	else if (time < t2) {
		double a=amplitude/(t2-t1); 
		return (a*(time-t1));
	}
	else if (time < t3) {
		return amplitude;
	}
	else if (time < t4) {
		double a=-amplitude/(t4-t3);
		return (amplitude+a*(time-t3));
	}
	else {
		return 0.0;
	}
}

Filter *
TrapezoidalModulatingFunction::getFilter()
{
	return theFilter;
}

double
TrapezoidalModulatingFunction::getMaxAmplitude()
{
	return amplitude;
}

void
TrapezoidalModulatingFunction::Print(OPS_Stream &s, int flag)  
{
}

