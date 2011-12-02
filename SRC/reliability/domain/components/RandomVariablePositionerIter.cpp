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
                                                                        
// $Revision: 1.1 $
// $Date: 2006-12-06 23:03:50 $
// $Source: /usr/local/cvs/OpenSees/SRC/reliability/domain/components/RandomVariablePositionerIter.cpp,v $

// Description: This file contains the method definitions for class 
// RandomVariablePositionerIter. RandomVariablePositionerIter is a class for iterating through the 
// elements of a domain. 

#include <RandomVariablePositionerIter.h>

#include <Element.h>
#include <TaggedObjectIter.h>
#include <TaggedObjectStorage.h>


// RandomVariablePositionerIter(SingleDomain &theDomain):
//	constructor that takes the model, just the basic iter

RandomVariablePositionerIter::RandomVariablePositionerIter(TaggedObjectStorage *theStorage)
  :myIter(theStorage->getComponents())
{
}


RandomVariablePositionerIter::~RandomVariablePositionerIter()
{
}    

void
RandomVariablePositionerIter::reset(void)
{
    myIter.reset();
}    


RandomVariablePositioner *
RandomVariablePositionerIter::operator()(void)
{
    // check if we still have elements in the model
    // if not return 0, indicating we are done
    TaggedObject *theComponent = myIter();
    if (theComponent == 0)
	return 0;
    else {
	RandomVariablePositioner *result = (RandomVariablePositioner *)theComponent;
	return result;
    }
}
