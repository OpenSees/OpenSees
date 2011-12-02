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
// $Date: 2006-12-13 18:17:37 $
// $Source: /usr/local/cvs/OpenSees/SRC/domain/domain/single/SingleDomParamIter.cpp,v $

// Description: This file contains the method definitions for class 
// SingleDomParamIter. SingleDomParamIter is a class for iterating through the 
// parameters of a domain. 

#include <SingleDomParamIter.h>

#include <Parameter.h>
#include <TaggedObjectIter.h>
#include <TaggedObjectStorage.h>


// SingleDomParamIter(SingleDomain &theDomain):
//	constructor that takes the model, just the basic iter

SingleDomParamIter::SingleDomParamIter(TaggedObjectStorage *theStorage)
  :myIter(theStorage->getComponents())
{
}


SingleDomParamIter::~SingleDomParamIter()
{
}    

void
SingleDomParamIter::reset(void)
{
    myIter.reset();
}    


Parameter *
SingleDomParamIter::operator()(void)
{
    // check if we still have parameters in the model
    // if not return 0, indicating we are done
    TaggedObject *theComponent = myIter();
    if (theComponent == 0)
	return 0;
    else {
	Parameter *result = (Parameter *)theComponent;
	return result;
    }
}

    
    
