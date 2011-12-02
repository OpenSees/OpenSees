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
// $Date: 2008-05-15 21:13:19 $
// $Source: /usr/local/cvs/OpenSees/SRC/reliability/domain/filter/FilterIter.h,v $

#ifndef FilterIter_h
#define FilterIter_h

class TaggedObjectStorage;
class TaggedObjectIter;

class Filter;

class FilterIter
{
  public:
    FilterIter(TaggedObjectStorage *theStorage);
    virtual ~FilterIter();

    virtual void reset(void);
    virtual Filter *operator()(void);
    
  private:
    TaggedObjectIter &myIter;
};

#endif





