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
// $Source: /usr/local/cvs/OpenSees/SRC/analysis/model/simple/SimpleFE_Iter.h,v $
                                                                        
                                                                        
// File: ~/analysis/model/simple/SimpleFE_Iter.h
//
// Written: fmk 
// Created: Fri Sep 20 15:27:47: 1996
// Revision: A
//
// Description: This file contains the class definition for SimpleFE_Iter.
// SimpleFE_Iter is an iter for returning the elements of an object of class
// AnalysisModel. SimpleFE_Iters must be written for each subclass of 
// SimpleAnalusis, wherin the elements are stored differently.

#ifndef SimpleFE_Iter_h
#define SimpleFE_Iter_h

#include <FE_EleIter.h>

class AnalysisModel;

class SimpleFE_Iter: public FE_EleIter
{
  public:
    SimpleFE_Iter(AnalysisModel &theModel);
    virtual ~SimpleFE_Iter();
    
    virtual void reset(void);
    virtual FE_Element *operator()(void);
    
  private:
    AnalysisModel &myModel;
    int currIndex;
    int numDone;
};

#endif





