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
                                                                        
// $Revision: 1.1 $
// $Date: 2003-10-27 23:04:40 $
// $Source: /usr/local/cvs/OpenSees/SRC/reliability/domain/modulatingFunction/KooModulatingFunction.cpp,v $


//
// Written by Terje Haukaas (haukaas@ce.berkeley.edu)
//
#include <math.h>
#include <KooModulatingFunction.h>
#include <ModulatingFunction.h>
#include <classTags.h>


KooModulatingFunction::KooModulatingFunction(int tag,
												 Filter *theFilt, 
												 double pt1,
												 double pt2)
:ModulatingFunction(tag,MODULATING_FUNCTION_trapezoidal)
{
	t1 = pt1;
	t2 = pt2;

	if (t1>t2 ) {
		opserr << "WARNING: Inconsistent input to Koo Modulating Function" << endln;
	}

	theFilter = theFilt;
}

KooModulatingFunction::~KooModulatingFunction()
{
}

double
KooModulatingFunction::getAmplitude(double time)
{
	double amplitude;
	if (time < 0.0) {
		amplitude = 0.0;
	}
	else if (time < t1) {
		amplitude = time*time/25.0; 
	}
	else if (time < t2) {
		amplitude = 1.0;
	}
	else {
		amplitude = exp(-0.5*(time-10.0));
	}

	return amplitude;
}

Filter *
KooModulatingFunction::getFilter()
{
	return theFilter;
}

double
KooModulatingFunction::getMaxAmplitude()
{
	return 1.0;
}

void
KooModulatingFunction::Print(OPS_Stream &s, int flag)  
{
}

