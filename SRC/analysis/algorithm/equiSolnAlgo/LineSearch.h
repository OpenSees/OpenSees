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
** ****************************************************************** */
                                                                        
// $Revision: 1.3 $
// $Date: 2003-02-14 23:00:42 $
// $Source: /usr/local/cvs/OpenSees/SRC/analysis/algorithm/equiSolnAlgo/LineSearch.h,v $
                                                                        
#ifndef LineSearch_h
#define LineSearch_h

// Written: fmk 
// Created: 11/01

// Description: This file contains the class definition for 
// LineSearch. LineSearch is an abstract base class, 
// i.e. no objects of it's type can be created.  Its subclasses seek
// to find a better solution to R(U)=0 than the solution Ui-1 + delta Ui
// would give, typically Ui = Ui-1 + factor * delta Ui.
// 
// What: "@(#)LineSearch.h, revA"

#include <MovableObject.h>
#include <OPS_Globals.h>

class SolutionAlgorithm;
class IncrementalIntegrator;
class LinearSOE;
//class OPS_Stream; //Jeremic@ucdavis.edu taken out since there is an include<iOPS_Stream.h> in LineSearch.h 

class LineSearch: public MovableObject
{
  public:
    LineSearch(int classTag);
    virtual ~LineSearch();

    // virtual functions
    virtual int newStep(LinearSOE &theSOE) =0;
    virtual int search(double s0, 
		       double s1, 
		       LinearSOE &theSOE, 
		       IncrementalIntegrator &theIntegrator) =0;
    virtual void Print(OPS_Stream &s, int flag =0) =0;
  protected:
    
  private:
};

#endif


