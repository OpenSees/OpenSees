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
// $Date: 2006-12-07 00:43:27 $
// $Source: /usr/local/cvs/OpenSees/SRC/reliability/domain/components/LimitStateFunctionIter.cpp,v $

// Description: This file contains the method definitions for class 
// LimitStateFunctionIter. LimitStateFunctionIter is a class for iterating through the 
// elements of a domain. 

#include <LimitStateFunctionIter.h>

#include <Element.h>
#include <TaggedObjectIter.h>
#include <TaggedObjectStorage.h>


// LimitStateFunctionIter(SingleDomain &theDomain):
//	constructor that takes the model, just the basic iter

LimitStateFunctionIter::LimitStateFunctionIter(TaggedObjectStorage *theStorage)
  :myIter(theStorage->getComponents())
{
}


LimitStateFunctionIter::~LimitStateFunctionIter()
{
}    

void
LimitStateFunctionIter::reset(void)
{
    myIter.reset();
}    


LimitStateFunction *
LimitStateFunctionIter::operator()(void)
{
    // check if we still have elements in the model
    // if not return 0, indicating we are done
    TaggedObject *theComponent = myIter();
    if (theComponent == 0)
	return 0;
    else {
		LimitStateFunction *result = (LimitStateFunction *)theComponent;
		return result;
    }
}
