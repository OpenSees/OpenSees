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
// $Source: /usr/local/cvs/OpenSees/SRC/domain/domain/single/SingleDomLC_Iter.cpp,v $
                                                                        
                                                                        
// File: ~/OOP/domain/domain/SingleDomLC_Iter.C
//
// Written: fmk 
// Created: Fri Sep 20 15:27:47: 1996
// Revision: A
//
// Description: This file contains the method definitions for class 
// SingleDomLC_Iter. SingleDomLC_Iter is a class for iterating through the 
// elements of a domain. 

#include "SingleDomLC_Iter.h"

#include <LoadCase.h>
#include <TaggedObjectIter.h>
#include <TaggedObjectStorage.h>

// SingleDomLC_Iter(SingleDomain &theDomain):
//	constructor that takes the model, just the basic iter

SingleDomLC_Iter::SingleDomLC_Iter(TaggedObjectStorage *theStorage)
  :myIter(theStorage->getComponents())
{
}

SingleDomLC_Iter::~SingleDomLC_Iter()
{
}    

void
SingleDomLC_Iter::reset(void)
{
    myIter.reset();
}


LoadCase *
SingleDomLC_Iter::operator()(void)
{
    // check if we still have LoadCases in the model
    // if not return 0, indicating we are done
    TaggedObject *theComponent = myIter();
    if (theComponent == 0)
	return 0;
    else {
	LoadCase *result = (LoadCase *)theComponent;
	return result;
    }
}

    
    
