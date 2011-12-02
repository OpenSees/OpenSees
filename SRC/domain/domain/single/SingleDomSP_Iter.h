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
// $Source: /usr/local/cvs/OpenSees/SRC/domain/domain/single/SingleDomSP_Iter.h,v $
                                                                        
                                                                        
// File: ~/domain/domain/single/SingleDomSP_Iter.h
//
// Written: fmk 
// Created: Fri Sep 20 15:27:47: 1996
// Revision: A
//
// Description: This file contains the class definition for 
// SingleDomSP_Iter. SingleDomSP_Iter is an iter for returning 
// the single point constraints  of an object of class SingleDomain. 
// SingleDomSP_Iters must be written for each subclass of SingleDomain 
// where the stoarge of the SP_Constraints changes.

#ifndef SingleDomSP_Iter_h
#define SingleDomSP_Iter_h

#include <SP_ConstraintIter.h>

class TaggedObjectStorage;
class TaggedObjectIter;


class SingleDomSP_Iter: public SP_ConstraintIter
{
  public:
    SingleDomSP_Iter(TaggedObjectStorage *theStorage);
    virtual ~SingleDomSP_Iter();
    
    virtual void reset(void);
    virtual SP_Constraint *operator()(void);
    
  private:
    TaggedObjectIter &myIter;
};

#endif

