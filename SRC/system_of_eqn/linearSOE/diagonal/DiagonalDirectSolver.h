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
// $Date: 2005-01-27 22:22:50 $
// $Source: /usr/local/cvs/OpenSees/SRC/system_of_eqn/linearSOE/diagonal/DiagonalDirectSolver.h,v $

// Written: fmk 
// Created: Jan 2005
// Revision: A
//
// Description: This file contains the class definition for 
// DiagonalDirectSolver. DiagonalDirectSolver is a subclass 
// of LinearSOESOlver. It solves diagonal system directly!

// What: "@(#) DiagonalDirectSolver.h, revA"

#ifndef DiagonalDirectSolver_h
#define DiagonalDirectSolver_h

#include <DiagonalSolver.h>
class DiagonalSOE;

class DiagonalDirectSolver : public DiagonalSolver
{
  public:
    DiagonalDirectSolver(double tol=1.0e-18);    
    virtual ~DiagonalDirectSolver();

    virtual int solve(void);        
    virtual int setSize(void);    
    double getDeterminant(void);
    
    virtual int setDiagonalSOE(DiagonalSOE &theSOE);

    int sendSelf(int commitTag, Channel &theChannel);
    int recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker);
    
  protected:
    double minDiagTol;
    
  private:

};


#endif

