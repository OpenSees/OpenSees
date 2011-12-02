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
// $Date: 2003-03-04 00:44:25 $
// $Source: /usr/local/cvs/OpenSees/SRC/reliability/domain/components/ReliabilityDomainComponent.cpp,v $


//
// Written by Terje Haukaas (haukaas@ce.berkeley.edu)
//

#include <ReliabilityDomain.h>
#include <ReliabilityDomainComponent.h>

ReliabilityDomainComponent::ReliabilityDomainComponent(int passedTag, int passedClassTag)
:TaggedObject(passedTag)
{
}


ReliabilityDomainComponent::~ReliabilityDomainComponent()
{
    // does nothing
}


void
ReliabilityDomainComponent::setReliabilityDomain(ReliabilityDomain *model)
{
    // sets the pointer 
    theReliabilityDomain = model;
}


ReliabilityDomain *
ReliabilityDomainComponent::getReliabilityDomain(void) const
{
    // returns the current pointer
    return theReliabilityDomain;
}


