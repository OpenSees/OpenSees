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
// $Date: 2000-09-15 08:23:29 $
// $Source: /usr/local/cvs/OpenSees/SRC/system_of_eqn/linearSOE/bandSPD/BandSPDLinThreadSolver.h,v $
                                                                        
                                                                        
// File: ~/system_of_eqn/linearSOE/bandSPD/BandSPDLinThreadSolver.h
//
// Written: fmk 
// Created: Mar 1998
// Revision: A
//
// Description: This file contains the class definition for 
// BandSPDLinThreadSolver. It solves the BandSPDLinSOE in parallel
// using solaris threads.
//
// What: "@(#) BandSPDLinThreadSolver.h, revA"

#ifndef BandSPDLinThreadSolver_h
#define BandSPDLinThreadSolver_h

#include <BandSPDLinSolver.h>

class BandSPDLinThreadSolver : public BandSPDLinSolver
{
  public:
    BandSPDLinThreadSolver();    
    BandSPDLinThreadSolver(int numProcessors, int blockSize);        
    ~BandSPDLinThreadSolver();

    int solve(void);
    int setSize(void);
    
    int sendSelf(Channel &theChannel, FEM_ObjectBroker &theBroker);
    int recvSelf(Channel &theChannel, FEM_ObjectBroker &theBroker);    
    
  protected:

  private:
    int NP;
    int running;
    int blockSize;
};

#endif


