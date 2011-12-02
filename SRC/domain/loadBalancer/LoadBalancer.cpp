/* ****************************************************************** **
**    OpenSees - Open System for Earthquake Engineering Simulation    **
**          Pacific Earthquake Engineering Research Center            **
**                                                                    **
**                                                                    **
** (C) Copyright 1999, The Regents of the University of California    **
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
** ****************************************************************** */
                                                                        
// $Revision: 1.2 $
// $Date: 2003-02-14 23:00:58 $
// $Source: /usr/local/cvs/OpenSees/SRC/domain/loadBalancer/LoadBalancer.cpp,v $
                                                                        
                                                                        
 // File: ~/domain/loadBalancer/LoadBalancer.C
// 
// Written: fmk 
// Created: Fri Aug 29 17:43:25 1997
// Revision: A
//
// Description: This file contains the class definition for LoadBalancer.
// A LoadBalancer is an object used to partition a PartitionedDomain.
//
// What: "@(#) LoadBalancer.C, revA"

#include <LoadBalancer.h>
 
LoadBalancer::LoadBalancer()
:theDomainPartitioner(0)
{
    
}

LoadBalancer::~LoadBalancer()
{
    
}

void
LoadBalancer::setLinks(DomainPartitioner &thePartitioner)
{
    theDomainPartitioner = &thePartitioner;
}

DomainPartitioner *
LoadBalancer::getDomainPartitioner(void)
{
    if (theDomainPartitioner == 0) {
	opserr << "WARNING LoadBalancer::getDomainPartitioner() ";
	opserr << " no DomainPartitioner is set - has setLinks() been called?\n";
    }
    
    return theDomainPartitioner;
}

