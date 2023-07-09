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
                                                                        
// $Revision: 1.3 $
// $Date: 2003-10-27 23:04:40 $
// $Source: /usr/local/cvs/OpenSees/SRC/reliability/domain/modulatingFunction/GammaModulatingFunction.cpp,v $


//
// Written by Terje Haukaas (haukaas@ce.berkeley.edu)
//

#include <GammaModulatingFunction.h>
#include <ModulatingFunction.h>
#include <classTags.h>
#include <math.h>

GammaModulatingFunction::GammaModulatingFunction(int tag,
												 Filter *theFilt,
												 double pa,
												 double pb,
												 double pc)
:ModulatingFunction(tag,MODULATING_FUNCTION_gamma)
{
	a = pa;
	b = pb;
	c = pc;
	theFilter = theFilt;
}

GammaModulatingFunction::~GammaModulatingFunction()
{
}

double
GammaModulatingFunction::getAmplitude(double time)
{
	return ( a * pow(time,b) * exp(-c*time)  );
}

Filter *
GammaModulatingFunction::getFilter()
{
	return theFilter;
}

double
GammaModulatingFunction::getMaxAmplitude()
{
	return ( a * pow( (b/c), b) * exp(-b) );
}

void
GammaModulatingFunction::Print(OPS_Stream &s, int flag)  
{
}

