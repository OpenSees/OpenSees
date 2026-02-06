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
                                                                        
// $Revision: 1.0 $
// $Date: 2025-05-29$
// $Source: /usr/local/cvs/OpenSees/SRC/domain/domain/single/SingleDomEQ_Iter.cpp,v $
                                                                        
                                                                        
// File: ~/OOP/domain/domain/SingleDomEQ_Iter.C
//
// Written: Yuli Huang (yulee@berkeley.edu)
// Created: 05/2020
// Revision: A
//
// Description: This file contains the method definitions for class 
// SingleDomEQ_Iter. SingleDomEQ_Iter is a class for iterating through the 
// elements of a domain. 

#include "SingleDomEQ_Iter.h"

#include <EQ_Constraint.h>
#include <TaggedObjectIter.h>
#include <TaggedObjectStorage.h>


// SingleDomEQ_Iter(SingleDomain &theDomain):
//	constructor that takes the model, just the basic iter

SingleDomEQ_Iter::SingleDomEQ_Iter(TaggedObjectStorage *theStorage)
  :myIter(theStorage->getComponents())
{
}

SingleDomEQ_Iter::~SingleDomEQ_Iter()
{
}    


void
SingleDomEQ_Iter::reset(void)
{
    myIter.reset();
}



EQ_Constraint *
SingleDomEQ_Iter::operator()(void)
{
    // check if we still have EQ_Constraints in the model
    // if not return 0, indicating we are done
    TaggedObject *theComponent = myIter();
    if (theComponent == 0)
	return 0;
    else {
	EQ_Constraint *result = (EQ_Constraint *)theComponent;
	return result;
    }
}

    
    


    
    
