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
// $Date: 2000-09-15 08:23:18 $
// $Source: /usr/local/cvs/OpenSees/SRC/domain/domain/partitioned/PartitionedDomainSubIter.cpp,v $
                                                                        
                                                                        
// File: ~/OOP/domain/domain/partitioned/PartitionedDomainSubIter.C
//
// Written: fmk 
// Created: Fri Sep 20 15:27:47: 1996
// Revision: A
//
// Description: This file contains the method definitions for class 
// PartitionedDomainSubIter. PartitionedDomainSubIter is a 
// class for iterating through the subdomains of an object of type
// PartitionedDomain.

#include <PartitionedDomainSubIter.h>

#include <Subdomain.h>
#include <TaggedObjectIter.h>
#include <TaggedObjectStorage.h>


// PartitionedDomainSubIter(PartitionedDomainain &theDomain):
//	constructor that takes the model, just the basic iter

PartitionedDomainSubIter::PartitionedDomainSubIter(TaggedObjectStorage *theStorage)
  :myIter(theStorage->getComponents())
{

}


PartitionedDomainSubIter::~PartitionedDomainSubIter()
{
}    

void
PartitionedDomainSubIter::reset(void)
{
    myIter.reset();
}    


Subdomain *
PartitionedDomainSubIter::operator()(void)
{
    // check if we still have Subdomains in the model
    // if not return 0, indicating we are done
    TaggedObject *theComponent = myIter();
    if (theComponent == 0)
	return 0;
    else {
	Subdomain *result = (Subdomain *)theComponent;
	return result;
    }
}

    
    
