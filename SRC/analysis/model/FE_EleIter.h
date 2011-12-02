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
// $Source: /usr/local/cvs/OpenSees/SRC/analysis/model/FE_EleIter.h,v $
                                                                        
                                                                        
// File: ~/analysis/model/FE_EleIter.h
//
// Written: fmk 
// Created: Fri Sep 20 15:27:47: 1996
// Revision: A
//
// Description: This file contains the class definition for FE_EleIter.
// FE_EleIter is an abstract base class. A FE_EleIter is an iter 
// for returning the FE_Elements of an object of class AnalysisModel. 
// FE_EleIters must be written for each subclass of AnalysisModel.


#ifndef FE_EleIter_h
#define FE_EleIter_h

class FE_Element;

class FE_EleIter 
{
  public:
    FE_EleIter() {};
    virtual ~FE_EleIter() {};
    
    virtual FE_Element *operator()(void) =0;

  protected:
    
  private:

};

#endif

