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
// $Date: 2000-09-15 08:23:19 $
// $Source: /usr/local/cvs/OpenSees/SRC/domain/pattern/LoadPatternIter.h,v $
                                                                        
                                                                        
// File: ~/domain/loadcase/LoadPatternIter.h
//
// Written: fmk 
// Created: Fri Sep 20 15:27:47: 1996
// Revision: A
//
// Description: This file contains the class definition for LoadPatternIter.
// LoadPatternIter is an abstract base class. An LoadPatternIter is an 
// iter for returning the LoadPatterns of an object of class  LoadCasse. 
// LoadPatternIters must be written for each subclass of LoadCase.

#ifndef LoadPatternIter_h
#define LoadPatternIter_h

class LoadPattern;
class TaggedObjectStorage;
class TaggedObjectIter;

class LoadPatternIter
{
  public:
    LoadPatternIter(TaggedObjectStorage *theStorage);
    virtual ~LoadPatternIter();

    virtual LoadPattern *operator()(void);
    virtual void reset(void);

  protected:
    
  private:
    TaggedObjectIter &myIter;
    
};

#endif

