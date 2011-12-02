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
                                                                        
// $Revision: 1.1 $
// $Date: 2005-05-25 23:33:46 $
// $Source: /usr/local/cvs/OpenSees/SRC/system_of_eqn/linearSOE/diagonal/DistributedDiagonalSolver.h,v $

// Written: fmk 
// Created: 05/05
//
// Description: This file contains the class definition for DistributedDiagonalSolver.
// DistributedDiagonalSolver is a concrete class for solving a system stored using
// a DistributedDiagonalSOE. The solve() method is overwritten to allow subclassing.
//
// What: "@(#) DistributedDiagonalSolver.h, revA"

#ifndef DistributedDiagonalSolver_h
#define DistributedDiagonalSolver_h

#include <LinearSOESolver.h>
class DistributedDiagonalSOE;

class DistributedDiagonalSolver : public LinearSOESolver
{
  public:
    DistributedDiagonalSolver(int classTag);    
    DistributedDiagonalSolver(double minDiagTol=1.0e-18);    
    virtual ~DistributedDiagonalSolver();

    virtual int solve(void);
    virtual int setSize(void);
    virtual int setLinearSOE(DistributedDiagonalSOE &theSOE);
    int sendSelf(int commitTag, Channel &theChannel);
    int recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker);    
    
  protected:
    DistributedDiagonalSOE *theSOE;

  private:
    double minDiagTol;
};

#endif

