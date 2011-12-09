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
// $Source: /usr/local/cvs/OpenSees/SRC/tagged/storage/ArrayOfTaggedObjectsIter.h,v $
                                                                        
                                                                        
// File: ~/tagged/storage/ArrayComponentIter
//
// Written: fmk 
// Created: Fri Sep 20 15:27:47: 1996
// Revision: A
//
// Description: This file contains the class definition for ArrayOfTaggedObjectsIter.
// ArrayOfTaggedObjectsIter is an iter for returning the DomainComponents of
// an object of type ArrayOfComponents.

#ifndef ArrayOfTaggedObjectsIter_h
#define ArrayOfTaggedObjectsIter_h

#include <TaggedObjectIter.h>

class ArrayOfTaggedObjects;

class ArrayOfTaggedObjectsIter: public TaggedObjectIter
{
  public:
    ArrayOfTaggedObjectsIter(ArrayOfTaggedObjects &theComponents);
    virtual ~ArrayOfTaggedObjectsIter();
    
    virtual void reset(void);
    virtual TaggedObject *operator()(void);
    
  private:
    ArrayOfTaggedObjects &myComponents;
    int currIndex;
    int numDone;
};

#endif





