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
                                                                        
// $Revision: 1.2 $
// $Date: 2005-11-28 21:23:50 $
// $Source: /usr/local/cvs/OpenSees/SRC/analysis/model/DOF_GrpIter.h,v $
                                                                        
// Written: fmk 
// Created: Fri Sep 20 15:27:47: 1996
// Revision: A
//
// Description: This file contains the class definition for DOF_GrpIter.
// DOF_GrpIter is an abstract base class. A DOF_GrpIter is an iter 
// for returning the DOF_Groups of an object of class AnalysisModel. 
// DOF_GrpIters must be written for each subclass of AnalysisModel.


#ifndef DOF_GrpIter_h
#define DOF_GrpIter_h

class DOF_Group;
class TaggedObjectStorage;
class TaggedObjectIter;

class DOF_GrpIter 
{
  public:
    DOF_GrpIter();
    DOF_GrpIter(TaggedObjectStorage *);
    virtual ~DOF_GrpIter();

    virtual void reset(void);    
    virtual DOF_Group *operator()(void);

  protected:
    
  private:
    TaggedObjectIter *myIter;
};

#endif

