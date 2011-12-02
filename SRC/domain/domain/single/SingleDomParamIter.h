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
// $Date: 2006-12-13 18:17:37 $
// $Source: /usr/local/cvs/OpenSees/SRC/domain/domain/single/SingleDomParamIter.h,v $

// Description: This file contains the class definition for SingleDomParamIter.
// SingleDomParamIter is an iter for returning the parameters of an object of class
// SingleDomain. SingleDomParamIters must be written for each subclass of 
// SingleDomain, wherin the parameters are stored differently.

#ifndef SingleDomParamIter_h
#define SingleDomParamIter_h

#include <ParameterIter.h>

class TaggedObjectStorage;
class TaggedObjectIter;

class SingleDomParamIter: public ParameterIter
{
  public:
    SingleDomParamIter(TaggedObjectStorage *theStorage);
    virtual ~SingleDomParamIter();

    virtual void reset(void);
    virtual Parameter *operator()(void);
    
  private:
    TaggedObjectIter &myIter;
};

#endif





