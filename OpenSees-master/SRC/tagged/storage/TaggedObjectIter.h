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
// $Source: /usr/local/cvs/OpenSees/SRC/tagged/storage/TaggedObjectIter.h,v $
                                                                        
                                                                        
// File: ~/tagged/storage/TaggedObjectIter.h
//
// Written: fmk 
// Created: Fri Sep 20 15:27:47: 1996
// Revision: A
//
// Description: This file contains the class definition for TaggedObjectIter.
// TaggedObjectIter is an abstract base class. An TaggedObjectIter is an 
// iter for returning the Components of an object of class TaggedObjectStorage.
// It must be written for each subclass of TaggedObjectStorage (this is done 
// for efficiency reasons), hence the abstract base class.

#ifndef TaggedObjectIter_h
#define TaggedObjectIter_h

class TaggedObject;

class TaggedObjectIter
{
  public:
    TaggedObjectIter() {};
    virtual ~TaggedObjectIter() {};

    virtual void reset(void) =0;
    virtual TaggedObject *operator()(void) =0;

  protected:
    
  private:

};

#endif





