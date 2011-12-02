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
// $Date: 2000-09-15 08:23:16 $
// $Source: /usr/local/cvs/OpenSees/SRC/analysis/analysis/TransientAnalysis.h,v $
                                                                        
                                                                        
#ifndef TransientAnalysis_h
#define TransientAnalysis_h

// File: ~/analysis/analysis/TransientAnalysis.h
// 
// Written: fmk 
// Created: 11/96
// Revision: A
//
// Description: This file contains the class definition for TransientAnalysis.
// TransientAnalysis is a subclass of AnalysisMethod, it is used to perform a 
// dynamic analysis on the FE\_Model. The class itself is an abstract base
// class. 
//
// What: "@(#) TransientAnalysis.h, revA"

#include <Analysis.h>

class TransientAnalysis: public Analysis
{
  public:
    TransientAnalysis(Domain &theDomain);
    virtual int analyze(int numSteps, double dT) = 0;
    
    virtual ~TransientAnalysis();
    
  protected:
    
  private:
};

#endif

