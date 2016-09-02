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
                                                                        
// $Revision: 1.0 $
// $Date: 2014-09-06 13:53:05 $
// $Source: /usr/local/cvs/OpenSees/SRC/system_of_eqn/linearSOE/sparseGEN/PFEMUnifiedSolver.h,v $

// Written: Minjie Zhu
// Created: September 2014
                                                                        
                                                                        
#ifndef PFEMUnifiedSolver_h
#define PFEMUnifiedSolver_h

//
// Description: This file contains the class definition for PFEMUnifiedSolver.
// A PFEMUnifiedSolver object can be constructed to solve a PFEMUnifiedLinSOE
// object. It obtains the solution by making calls on the
// The PFEMUnifiedSolver uses Fractional Step Method to solve PFEM equations. 
//
// What: "@(#) PFEMUnifiedSolver.h, revA"

#include <LinearSOESolver.h>

extern "C" {
#include <cs.h>
}

class PFEMUnifiedLinSOE;

class PFEMUnifiedSolver : public LinearSOESolver
{
public:
    PFEMUnifiedSolver();
    virtual ~PFEMUnifiedSolver();

    int solve();
    int setSize();
    int setLinearSOE(PFEMUnifiedLinSOE& theSOE);

    int sendSelf(int commitTag, Channel &theChannel);
    int recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker);  

private:
    
    PFEMUnifiedLinSOE* theSOE;

    css* Msym;
    csn* Mnum;
};

#endif

