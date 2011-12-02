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
// $Source: /usr/local/cvs/OpenSees/SRC/domain/domain/single/SingleDomAllSP_Iter.h,v $
                                                                        
                                                                        
// File: ~/domain/domain/single/SingleDomAllSP_Iter.h
//
// Written: fmk 
// Created: Fri Sep 20 15:27:47: 1996
// Revision: A
//
// Description: This file contains the class definition for 
// SingleDomAllSP_Iter. SingleDomAllSP_Iter is an iter for returning 
// all the single point constraints  of an object of class Domain, including
// those of the LoadPAtterns.
// SingleDomAllSP_Iters must be written for each subclass of SingleDomAllain 
// where the stoarge of the SP_Constraints changes.

#ifndef SingleDomAllSP_Iter_h
#define SingleDomAllSP_Iter_h

#include <SP_ConstraintIter.h>

class TaggedObjectStorage;
class TaggedObjectIter;
class Domain;
class LoadPatternIter;
class LoadPatternSPIter;
class LoadPattern;

class SingleDomAllSP_Iter: public SP_ConstraintIter
{
  public:
    SingleDomAllSP_Iter(Domain &theDomain);
    virtual ~SingleDomAllSP_Iter();
    
    virtual void reset(void);
    virtual SP_Constraint *operator()(void);
    
  private:
    Domain                *theDomain;
    bool                   doneDomainSPs;
    SP_ConstraintIter      *theDomainSPs;
    LoadPatternIter        *theLoadPatterns;
    LoadPattern            *currentLoadPattern;
    SP_ConstraintIter      *theLoadPatternSPs;
};

#endif
