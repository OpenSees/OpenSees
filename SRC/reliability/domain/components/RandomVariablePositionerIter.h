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
// $Source: /usr/local/cvs/OpenSees/SRC/reliability/domain/components/RandomVariablePositionerIter.h,v $

// Description: This file contains the class definition for RandomVariableIter.
// RandomVariablePositionerIter is an iter for returning the elements of an object of class
// SingleDomain. RandomVariablePositionerIters must be written for each subclass of 
// SingleDomain, wherin the elements are stored differently.

#ifndef RandomVariablePositionerIter_h
#define RandomVariablePositionerIter_h

class TaggedObjectStorage;
class TaggedObjectIter;

class RandomVariablePositioner;

class RandomVariablePositionerIter
{
  public:
    RandomVariablePositionerIter(TaggedObjectStorage *theStorage);
    virtual ~RandomVariablePositionerIter();

    virtual void reset(void);
    virtual RandomVariablePositioner *operator()(void);
    
  private:
    TaggedObjectIter &myIter;
};

#endif





