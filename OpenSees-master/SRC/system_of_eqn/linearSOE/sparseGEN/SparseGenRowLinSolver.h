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
                                                                        
#ifndef SparseGenRowLinSolver_h
#define SparseGenRowLinSolver_h

// $Revision: 1.1 $
// $Date: 2005-04-08 02:38:18 $
// $Source: /usr/local/cvs/OpenSees/SRC/system_of_eqn/linearSOE/sparseGEN/SparseGenRowLinSolver.h,v $

// Written: fmk 
// Created: 04/05
// Revision: A
//
// Description: This file contains the class definition for SparseGenRowLinSolver.
// SparseGenRowLinSolver is an abstract base class and thus no objects of it's type
// can be instantiated. It has pure virtual functions which must be
// implemented in it's derived classes.  Instances of SparseGenRowLinSolver 
// are used to solve a system of equations of type SparseGenRowLinSOE.
//
// What: "@(#) SparseGenRowLinSolver.h, revA"

#include <LinearSOESolver.h>
class SparseGenRowLinSOE;

class SparseGenRowLinSolver : public LinearSOESolver
{
  public:
    SparseGenRowLinSolver(int classTag);    
    virtual ~SparseGenRowLinSolver();

    virtual int setLinearSOE(SparseGenRowLinSOE &theSOE);
    
  protected:
    SparseGenRowLinSOE *theSOE;

  private:

};

#endif

