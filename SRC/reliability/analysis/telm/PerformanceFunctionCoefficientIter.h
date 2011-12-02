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
// $Date: 2008-02-29 19:43:53 $
// $Source: /usr/local/cvs/OpenSees/SRC/reliability/analysis/telm/PerformanceFunctionCoefficientIter.h,v $
                                                                        
#ifndef PerformanceFunctionCoefficientIter_h
#define PerformanceFunctionCoefficientIter_h

//#include <NodeIter.h>
class TaggedObjectStorage;
class TaggedObjectIter;
class PerformanceFunctionCoeff;

class PerformanceFunctionCoefficientIter 
{
  public:
    PerformanceFunctionCoefficientIter(TaggedObjectStorage *theStorage);
    ~PerformanceFunctionCoefficientIter();
    
    void reset(void);
    PerformanceFunctionCoeff *operator()(void);
    
  private:
    TaggedObjectIter &myIter;
};

#endif

