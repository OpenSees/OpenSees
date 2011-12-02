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
// $Date: 2000-09-15 08:23:30 $
// $Source: /usr/local/cvs/OpenSees/SRC/system_of_eqn/linearSOE/umfGEN/UmfpackGenLinSolver.h,v $
                                                                        
                                                                        
// File: ~/system_of_eqn/linearSOE/umfGEN/UmfpackGenLinSolver.h
//
// Written: fmk 
// Created: 11/98
// Revision: A
//
// Description: This file contains the class definition for 
// UmfpackGenLinSolver. It solves the UmfpackGenLinSOEobject by calling
// UMFPACK2.2.1 routines.
//
// What: "@(#) UmfpackGenLinSolver.h, revA"

#ifndef UmfpackGenLinSolver_h
#define UmfpackGenLinSolver_h

#include <LinearSOESolver.h>

class UmfpackGenLinSOE;

class UmfpackGenLinSolver : public LinearSOESolver
{
  public:
    UmfpackGenLinSolver();     
    ~UmfpackGenLinSolver();

    int solve(void);
    int setSize(void);

    int setLinearSOE(UmfpackGenLinSOE &theSOE);
    
    int sendSelf(int commitTag, Channel &theChannel);
    int recvSelf(int commitTag, Channel &theChannel, 
		 FEM_ObjectBroker &theBroker);    
    
  protected:

  private:
    int icntl[20];
    int keep[20];
    double cntl[10];
    int info[40];
    double rinfo[20];
    
    int *copyIndex;
    int lIndex;
    double *work;
    
    UmfpackGenLinSOE *theSOE;
};

#endif

