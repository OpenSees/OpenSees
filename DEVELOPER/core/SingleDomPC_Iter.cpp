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
// $Date: 2012-08-22 12:09:28 $
// $Source: /usr/local/cvs/OpenSees/SRC/domain/domain/single/SingleDomPC_Iter.cpp,v $
                                                                        
                                                                        
// File: ~/OOP/domain/domain/SingleDomPC_Iter.C
//
// Written: Minjie
// Created: Aug 22 2012
// Revision: A
//
// Description: This file contains the method definitions for class 
// SingleDomPC_Iter. SingleDomPC_Iter is a class for iterating through the 
// elements of a domain. 

#include "SingleDomPC_Iter.h"

#include <Pressure_Constraint.h>
#include <TaggedObjectIter.h>
#include <TaggedObjectStorage.h>


// SingleDomPC_Iter(SingleDomain &theDomain):
//	constructor that takes the model, just the basic iter

SingleDomPC_Iter::SingleDomPC_Iter(TaggedObjectStorage *theStorage)
  :myIter(theStorage->getComponents())
{
}

SingleDomPC_Iter::~SingleDomPC_Iter()
{
}    


void
SingleDomPC_Iter::reset(void)
{
    myIter.reset();
}


Pressure_Constraint *
SingleDomPC_Iter::operator()(void)
{
    // check if we still have Pressure_Constraints in the model
    // if not return 0, indicating we are done
    TaggedObject *theComponent = myIter();
    if (theComponent == 0)
	return 0;
    else {
	Pressure_Constraint *result = (Pressure_Constraint *)theComponent;
	return result;
    }
}

    
    


    
    
