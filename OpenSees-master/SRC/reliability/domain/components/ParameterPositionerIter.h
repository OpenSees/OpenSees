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
// $Source: /usr/local/cvs/OpenSees/SRC/reliability/domain/components/ParameterPositionerIter.h,v $

// Description: This file contains the class definition for ParameterIter.
// ParameterPositionerIter is an iter for returning the elements of an object of class
// SingleDomain. ParameterPositionerIters must be written for each subclass of 
// SingleDomain, wherin the elements are stored differently.

#ifndef ParameterPositionerIter_h
#define ParameterPositionerIter_h

class TaggedObjectStorage;
class TaggedObjectIter;

class ParameterPositioner;

class ParameterPositionerIter
{
  public:
    ParameterPositionerIter(TaggedObjectStorage *theStorage);
    virtual ~ParameterPositionerIter();

    virtual void reset(void);
    virtual ParameterPositioner *operator()(void);
    
  private:
    TaggedObjectIter &myIter;
};

#endif





