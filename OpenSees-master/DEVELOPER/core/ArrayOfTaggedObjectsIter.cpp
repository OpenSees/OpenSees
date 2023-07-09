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
// $Date: 2000-09-15 08:23:30 $
// $Source: /usr/local/cvs/OpenSees/SRC/tagged/storage/ArrayOfTaggedObjectsIter.cpp,v $
                                                                        
                                                                        
// File: ~/OOP/tagged/storage/ArrayOfTaggedObjectsIter.C
//
// Written: fmk 
// Created: Fri Sep 20 15:27:47: 1996
// Revision: A
//
// Description: This file contains the method definitions for class 
// ArrayOfTaggedObjectsIter. ArrayOfTaggedObjectsIter is a class for iterating through the 
// elements of a domain. 

#include <ArrayOfTaggedObjectsIter.h>
#include <ArrayOfTaggedObjects.h>


// ArrayOfTaggedObjectsIter(SingleDomain &theDomain):
//	constructor that takes the model, just the basic iter

ArrayOfTaggedObjectsIter::ArrayOfTaggedObjectsIter(ArrayOfTaggedObjects &theComponents)
  :myComponents(theComponents), currIndex(0), numDone(0)
{
}


ArrayOfTaggedObjectsIter::~ArrayOfTaggedObjectsIter()
{
}    

void
ArrayOfTaggedObjectsIter::reset(void)
{
    currIndex = 0;
    numDone = 0;
}

TaggedObject *
ArrayOfTaggedObjectsIter::operator()(void)
{
  // check if we still have elements in the model
  // if not return 0, indicating we are done

  // have to remove this if delete from an iter
  // if (numDone >= myComponents.numComponents)
  //    return 0;

  // search through domains ele list till we find the next element
  while ((currIndex <= myComponents.positionLastEntry) 
	 && (myComponents.theComponents[currIndex] == 0))
      currIndex++;

  // if not at the end of the list return the element
  // NOTE: a BAD type cast is needed here - Maybe change domain
  if (currIndex < myComponents.sizeComponentArray) {
      TaggedObject *dc= myComponents.theComponents[currIndex];
      numDone++; currIndex++;
      return(dc);
  }
  return (0);
}

    
    
