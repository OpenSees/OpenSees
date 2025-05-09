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
                                                                        
// $Revision: 1.0 $
// $Date: 2025-05-29$
// $Source: /usr/local/cvs/OpenSees/SRC/domain/domain/single/SingleDomEQ_Iter.h,v $
                                                                        
                                                                        
// File: ~/domain/domain/single/SingleDomEQ_Iter.h
//
// Written: Yuli Huang (yulee@berkeley.edu)
// Created: 05/2020
// Revision: A
//
// Description: This file contains the class definition for 
// SingleDomEQ_Iter. SingleDomEQ_Iter is an iter for returning 
// the equation constraints  of an object of class SingleDomain. 
// SingleDomEQ_Iters must be written for each subclass of SingleDomain.
// where the storage of the EQ_Constraints changes.

#ifndef SingleDomEQ_Iter_h
#define SingleDomEQ_Iter_h

#include <EQ_ConstraintIter.h>

class TaggedObjectStorage;
class TaggedObjectIter;

class SingleDomEQ_Iter: public EQ_ConstraintIter
{
  public:
    SingleDomEQ_Iter(TaggedObjectStorage *theStorage);
    virtual ~SingleDomEQ_Iter();
    
    virtual void reset(void);
    virtual EQ_Constraint *operator()(void);
    
  private:
    TaggedObjectIter &myIter;
};

#endif

