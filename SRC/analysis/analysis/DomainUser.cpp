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
                                                                        
// $Revision: 1.1.1.1 $
// $Date: 2000-09-15 08:23:16 $
// $Source: /usr/local/cvs/OpenSees/SRC/analysis/analysis/DomainUser.cpp,v $
                                                                        
                                                                        
// File: ~/analysis/analysis/DomainUser.C
// 
// Written: fmk 
// Created: Sun Sept 15 11:47:47: 1996
// Revision: A
//
// Description: This file contains the implementation of DomainUser.
// DomainUser is an abstract base class, i.e. no objects of it's
// type can be created. 
//
// What: "@(#) DomainUser.C, revA"

#include <DomainUser.h>
#include <Analysis.h>
#include <Domain.h>

// DomainUser(theDomain &theDomain);
// All DomainUser are associated with a single domain. 


DomainUser::DomainUser(Domain &theDomain)
:myDomain(theDomain)
{

}

// ~DomainUser();
// All DomainUser are associated with a single domain, the destructor
// removes the link in the domain by invoking {\em removeDomainUser(*this)}
// on the domain. 

DomainUser::~DomainUser()
{
//    myDomain.removeDomainUser(this);
}

// Domain *getDomainPtr(void) const;
// A const method which returns a pointer to the Domain object on which
// the DomainUser performs its DomainUser.

Domain &
DomainUser::getDomain(void) const
{
    return myDomain;
}


