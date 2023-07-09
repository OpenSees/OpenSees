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
// $Date: 2012-08-22 12:07:25 $
// $Source: /usr/local/cvs/OpenSees/SRC/domain/domain/single/SingleDomPC_Iter.h,v $
                                                                        
                                                                        
// File: ~/domain/domain/single/SingleDomPC_Iter.h
//
// Written: Minjie Zhu
// Created: Aug 22, 2012
// Revision: A
//
// Description: This file contains the class definition for 
// SingleDomPC_Iter. SingleDomPC_Iter is an iter for returning 
// the single point constraints  of an object of class SingleDomain. 
// SingleDomPC_Iters must be written for each subclass of SingleDomain 
// where the stoarge of the Pressure_Constraints changes.

#ifndef SingleDomPC_Iter_h
#define SingleDomPC_Iter_h

#include <Pressure_ConstraintIter.h>

class TaggedObjectStorage;
class TaggedObjectIter;


class SingleDomPC_Iter: public Pressure_ConstraintIter
{
  public:
    SingleDomPC_Iter(TaggedObjectStorage *theStorage);
    virtual ~SingleDomPC_Iter();
    
    virtual void reset(void);
    virtual Pressure_Constraint *operator()(void);
    
  private:
    TaggedObjectIter &myIter;
};

#endif

