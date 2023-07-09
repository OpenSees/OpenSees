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
// $Date: 2005-11-28 22:07:24 $
// $Source: /usr/local/cvs/OpenSees/SRC/analysis/model/FE_EleIter.cpp,v $
                                                                        
// Written: fmk 
// Created: 10/05
// Revision: A
//

#include "FE_EleIter.h"

#include <FE_Element.h>
#include <TaggedObjectIter.h>
#include <TaggedObjectStorage.h>


// FE_EleIter(SingleDomain &theDomain):
//	constructor that takes the model, just the basic iter

FE_EleIter::FE_EleIter(TaggedObjectStorage *theStorage)
  :myIter(&(theStorage->getComponents()))
{
}


FE_EleIter::~FE_EleIter()
{
}    

void
FE_EleIter::reset(void)
{
    myIter->reset();
}    


FE_Element *
FE_EleIter::operator()(void)
{
    // check if we still have elements in the model
    // if not return 0, indicating we are done
    TaggedObject *theComponent = (*myIter)();
    if (theComponent == 0)
	return 0;
    else {
	FE_Element *result = (FE_Element *)theComponent;
	return result;
    }
}

    
    
