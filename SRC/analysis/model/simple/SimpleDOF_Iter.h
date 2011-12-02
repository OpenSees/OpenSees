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
// $Date: 2000-09-15 08:23:17 $
// $Source: /usr/local/cvs/OpenSees/SRC/analysis/model/simple/SimpleDOF_Iter.h,v $
                                                                        
                                                                        
// File: ~/analysis/model/simple/SimpleDOF_Iter.h
//
// Written: fmk 
// Created: Fri Sep 20 15:27:47: 1996
// Revision: A
//
// Description: This file contains the class definition for SimpleDOF_Iter.
// SimpleDOF_Iter is an iter for returning the DOF_Groups of an object of class
// SimpleAnalysisModel. SimpleDOF_Iters must be written for each subclass of 
// SimpleAnalusisModel, wherin the elements are stored differently.

#ifndef SimpleDOF_Iter_h
#define SimpleDOF_Iter_h

#include <DOF_GrpIter.h>

class AnalysisModel;

class SimpleDOF_Iter: public DOF_GrpIter
{
  public:
    SimpleDOF_Iter(AnalysisModel &theModel);
    virtual ~SimpleDOF_Iter();
    
    virtual void reset(void);
    virtual DOF_Group *operator()(void);
    
  private:
    AnalysisModel &myModel;
    int currIndex;
    int numDone;
};

#endif





