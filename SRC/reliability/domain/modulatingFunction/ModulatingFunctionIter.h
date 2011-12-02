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
// $Date: 2008-05-15 21:13:32 $
// $Source: /usr/local/cvs/OpenSees/SRC/reliability/domain/modulatingFunction/ModulatingFunctionIter.h,v $

#ifndef ModulatingFunctionIter_h
#define ModulatingFunctionIter_h

class TaggedObjectStorage;
class TaggedObjectIter;

class ModulatingFunction;

class ModulatingFunctionIter
{
  public:
    ModulatingFunctionIter(TaggedObjectStorage *theStorage);
    virtual ~ModulatingFunctionIter();

    virtual void reset(void);
    virtual ModulatingFunction *operator()(void);
    
  private:
    TaggedObjectIter &myIter;
};

#endif





