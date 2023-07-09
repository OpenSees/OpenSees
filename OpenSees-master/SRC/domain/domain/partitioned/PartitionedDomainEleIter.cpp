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
// $Source: /usr/local/cvs/OpenSees/SRC/domain/domain/partitioned/PartitionedDomainEleIter.cpp,v $
                                                                        
                                                                        
// File: ~/OOP/domain/domain/partitioned/PartitionedDomainEleIter.C
//
// Written: fmk 
// Created: Fri Sep 20 15:27:47: 1996
// Revision: A
//
// Description: This file contains the method definitions for class 
// PartitionedDomainEleIter. PartitionedDomainEleIter is a 
// class for iterating through the elements of a ParitionedDomain
// domain. 
//
// What: "@(#) PartitionedDomainEleIter.C, revA"

#include <PartitionedDomainEleIter.h>

#include <Subdomain.h>
#include <Element.h>
#include <SingleDomEleIter.h>
#include <PartitionedDomain.h>
#include <ArrayOfTaggedObjectsIter.h>
#include <ArrayOfTaggedObjects.h>



// PartitionedDomainEleIter(SingleDomain &theDomain):
//	constructor that takes the model, just the basic iter

PartitionedDomainEleIter::
PartitionedDomainEleIter(PartitionedDomain *partitionedDomain)
  :subdomainIter(0), currentIter(0), currentSubdomain(0),
   thePartitionedDomain(partitionedDomain)
{
    mainEleIter = new SingleDomEleIter(thePartitionedDomain->elements); 
    subdomainIter = new ArrayOfTaggedObjectsIter(
		     *(thePartitionedDomain->theSubdomains));
}


PartitionedDomainEleIter::~PartitionedDomainEleIter()
{
    delete subdomainIter;
}    

void
PartitionedDomainEleIter::reset(void)
{
    mainDomain = true;
    mainEleIter->reset();
    subdomainIter->reset();
    currentIter = mainEleIter;

    TaggedObject *currentObject = (*subdomainIter)();
    if (currentObject != 0)
	currentSubdomain = (Subdomain *)currentObject;
    else
	currentSubdomain = 0;
}

Element *
PartitionedDomainEleIter::operator()(void)
{
    Element *theEle;

    while ((currentSubdomain != 0 || mainDomain == true)) {
	if (mainDomain == true) {
	    theEle = (*currentIter)();

	    if (theEle != 0) {
		return theEle;
            }
	    else {
		mainDomain = false;
		Element *res = currentSubdomain;
		TaggedObject *currentObject = (*subdomainIter)();
		currentSubdomain = (Subdomain *)currentObject;
		return res;
	    }
	} else {
	    Element *res = currentSubdomain;
	    TaggedObject *currentObject = (*subdomainIter)();
	    currentSubdomain = (Subdomain *)currentObject;
	    return res;
	}
    }
    
    // we will only get here if we are done 
    return 0;
}
