/* ****************************************************************** **
**    OpenSees - Open System for Earthquake Engineering Simulation    **
**          Pacific Earthquake Engineering Research Center            **
**                                                                    **
**                                                                    **
** (C) Copyright 2001, The Regents of the University of California    **
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
**   Optimization module developed by:                                **
**   Quan Gu  (qgu@ucsd.edu)                                          **
**   Joel Conte (jpconte@ucsd.edu)                                    **
**   Philip Gill (pgill@ucsd.edu)                                     **
** ****************************************************************** */


//
// Written by Quan Gu (qgu@ucsd.edu)    March 2010
//



#ifndef DesignVariablePositionerIter_h
#define DesignVariablePositionerIter_h

class TaggedObjectStorage;
class TaggedObjectIter;

class DesignVariablePositioner;

class DesignVariablePositionerIter
{
  public:
    DesignVariablePositionerIter(TaggedObjectStorage *theStorage);
    virtual ~DesignVariablePositionerIter();

    virtual void reset(void);
    virtual DesignVariablePositioner *operator()(void);
    
  private:
    TaggedObjectIter &myIter;
};

#endif





