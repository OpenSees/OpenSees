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
                                                                        
#ifndef DistributedSparseGenRowLinSolver_h
#define DistributedSparseGenRowLinSolver_h

// $Revision: 1.1 $
// $Date: 2005-12-06 22:20:18 $
// $Source: /usr/local/cvs/OpenSees/SRC/system_of_eqn/linearSOE/sparseGEN/DistributedSparseGenRowLinSolver.h,v $

// Written: fmk 
// Created: 04/05
// Revision: A
//
// Description: This file contains the class definition for DistributedSparseGenRowLinSolver.
// DistributedSparseGenRowLinSolver is an abstract base class and thus no objects of it's type
// can be instantiated. It has pure virtual functions which must be
// implemented in it's derived classes.  Instances of DistributedSparseGenRowLinSolver 
// are used to solve a system of equations of type DistributedSparseGenRowLinSOE.
//
// What: "@(#) DistributedSparseGenRowLinSolver.h, revA"

#include <LinearSOESolver.h>
class DistributedSparseGenRowLinSOE;

class DistributedSparseGenRowLinSolver : public LinearSOESolver
{
  public:
    DistributedSparseGenRowLinSolver(int classTag);    
    virtual ~DistributedSparseGenRowLinSolver();

    virtual int setLinearSOE(DistributedSparseGenRowLinSOE &theSOE);
    
  protected:
    DistributedSparseGenRowLinSOE *theSOE;

  private:

};

#endif

