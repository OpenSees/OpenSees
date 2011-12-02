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
                                                                        
// $Revision: 1.3 $
// $Date: 2001-06-14 05:26:47 $
// $Source: /usr/local/cvs/OpenSees/SRC/convergenceTest/ConvergenceTest.h,v $
                                                                        
                                                                        
#ifndef ConvergenceTest_h
#define ConvergenceTest_h

// File: ~/convergenceTest/ConvergenceTest.h
//
// Written: fmk 
// Date: 09/98
// Revised:
//
// Purpose: This file contains the class definition for ConvergenceTest,
// which is an abstract class. Objects of concrete subclasses can be used 
// to test the convergence of an algorithm. 

#include <MovableObject.h>
#include <Vector.h>
#include <bool.h>

class EquiSolnAlgo;


class ConvergenceTest: public MovableObject
{
  public:
    // constructors and destructor
    ConvergenceTest(int classTag);	
    virtual ~ConvergenceTest();

    virtual ConvergenceTest *getCopy( int iterations ) = 0 ;

    virtual int setEquiSolnAlgo(EquiSolnAlgo &theAlgorithm) =0;
    virtual int start(void) =0;
    virtual int test(void) = 0;
    
    virtual int getNumTests(void) =0;    
    virtual int getMaxNumTests(void) =0;        
    virtual double getRatioNumToMax(void) =0;            
    virtual const Vector &getNorms(void) =0;
    
    
  protected:

  private:
};


#endif

