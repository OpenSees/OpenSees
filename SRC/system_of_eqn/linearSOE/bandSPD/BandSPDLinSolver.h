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
// $Date: 2009-05-11 20:55:31 $
// $Source: /usr/local/cvs/OpenSees/SRC/system_of_eqn/linearSOE/bandSPD/BandSPDLinSolver.h,v $
                                                                        
                                                                        
#ifndef BandSPDLinSolver_h
#define BandSPDLinSolver_h

// Written: fmk 
// Created: Tue Sep 26 16:27:47: 1996
//
// Description: This file contains the class definition for BandSPDLinSolver.
// BandSPDLinSolver is an abstract base class and thus no objects of it's type
// can be instantiated. It has pure virtual functions which must be
// implemented in it's derived classes.  Instances of BandSPDLinSolver 
// are used to solve a system of equations of type BandSPDLinSOE.
//
// What: "@(#) BandSPDLinSolver.h, revA"

#include <LinearSOESolver.h>
class BandSPDLinSOE;

class BandSPDLinSolver : public LinearSOESolver
{
  public:
    BandSPDLinSolver(int classTag);    
    virtual ~BandSPDLinSolver();

    virtual int solve(void) = 0;
    virtual int setLinearSOE(BandSPDLinSOE &theSOE);
    
  protected:
    BandSPDLinSOE *theSOE;

  private:

};

#endif

