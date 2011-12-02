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
// $Date: 2000-09-15 08:23:19 $
// $Source: /usr/local/cvs/OpenSees/SRC/domain/subdomain/SubdomainNodIter.cpp,v $
                                                                        
                                                                        
// File: ~/OOP/domain/subdomain/SubdomainNodIter.C
//
// Written: fmk 
// Created: Fri Sep 20 15:27:47: 1996
// Revision: A
//
// Description: This file contains the method definitions for class 
// SubdomainNodIter. SubdomainNodIter is a 
// class for iterating through the Nodes of a ParitionedDomain
// domain. 
//
// What: "@(#) SubdomainNodIter.C, revA"

#include <SubdomainNodIter.h>

#include <Subdomain.h>
#include <Node.h>
#include <ArrayOfTaggedObjectsIter.h>
#include <ArrayOfTaggedObjects.h>


// SubdomainNodIter(SingleDomain &theDomain):
//	constructor that takes the model, just the basic iter

SubdomainNodIter::SubdomainNodIter(Subdomain &theSub)
  :currentIter(0), theSubdomain(&theSub), external(true)
{
}


SubdomainNodIter::~SubdomainNodIter()
{
}    

void
SubdomainNodIter::reset(void)
{
    currentIter = &(theSubdomain->getExternalNodeIter());
    external = true;
}

Node *
SubdomainNodIter::operator()(void)
{
    Node *theNod;
    
    theNod = (*currentIter)();
    if (theNod != 0)
	return theNod;
    else  
	if (external == true) {
	    currentIter = &(theSubdomain->getInternalNodeIter()); 
	    external = false;
	    return (*currentIter)();
	}     
	else 
	    return 0;
}	


    
    



