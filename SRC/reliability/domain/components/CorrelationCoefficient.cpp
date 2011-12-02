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
// $Date: 2003-03-04 00:44:23 $
// $Source: /usr/local/cvs/OpenSees/SRC/reliability/domain/components/CorrelationCoefficient.cpp,v $


//
// Written by Terje Haukaas (haukaas@ce.berkeley.edu)
//

#include <CorrelationCoefficient.h>
#include <classTags.h>

CorrelationCoefficient::CorrelationCoefficient(int passedTag,
							int passedRv1,
							int passedRv2,
							double passedCorrelation)
:ReliabilityDomainComponent(passedTag, CORRELATION_COEFFICIENT)
{
	tag = passedTag;
	rv1 = passedRv1;
	rv2 = passedRv2;
	correlation = passedCorrelation;
}


CorrelationCoefficient::~CorrelationCoefficient()
{
}



void
CorrelationCoefficient::Print(OPS_Stream &s, int flag)  
{
}



int
CorrelationCoefficient::getRv1()
{
	return rv1;
}

int
CorrelationCoefficient::getRv2()
{
	return rv2;
}

double
CorrelationCoefficient::getCorrelation()
{
	return correlation;
}
