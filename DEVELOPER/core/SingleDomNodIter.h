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
// $Source: /usr/local/cvs/OpenSees/SRC/domain/domain/single/SingleDomNodIter.h,v $
                                                                        
                                                                        
// File: ~/domain/domain/single/SingleDomNodIter.h
//
// Written: fmk 
// Created: Fri Sep 20 15:27:47: 1996
// Revision: A
//
// Description: This file contains the class definition for SingleDomNodIter.
// SingleDomNodIter is an iter for returning the nodes of an object of class
// SingleDomain. SingleDomNodIters must be written for each subclass of 
// SingleDomain where the storage of the nodes changes.

#ifndef SingleDomNodIter_h
#define SingleDomNodIter_h

#include <NodeIter.h>
class TaggedObjectStorage;
class TaggedObjectIter;

class SingleDomNodIter : public NodeIter
{
  public:
    SingleDomNodIter(TaggedObjectStorage *theStorage);
    virtual ~SingleDomNodIter();
    
    virtual void reset(void);
    virtual Node *operator()(void);
    
  private:
    TaggedObjectIter &myIter;
};

#endif

