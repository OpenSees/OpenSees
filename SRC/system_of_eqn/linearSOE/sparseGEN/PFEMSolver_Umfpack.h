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
// $Date: 2012-09-17 10:51:44 $
// $Source: /usr/local/cvs/OpenSees/SRC/system_of_eqn/linearSOE/sparseGEN/PFEMSolver_Umfpack.h,v $
                                                                        
                                                                        
#ifndef PFEMSolverUmfpack_h
#define PFEMSolver_Umfpack_h

// File: ~/system_of_eqn/linearSOE/sparseGEN/PFEMSolver_Umfpack.h
//
// Written: Minjie 
// Created: Sep 17 2012
//
// Description: This file contains the class definition for PFEMSolver_Umfpack.
// A PFEMSolver_Umfpack object can be constructed to solve a PFEMLinSOE
// object. It obtains the solution by making calls on the
// The PFEMSolver_Umfpack uses Fractional Step Method to solve PFEM equations. 
//
// What: "@(#) PFEMSolver_Umfpack.h, revA"

#include <PFEMSolver.h>
extern "C" {
#include <cs.h>
}
#include "../../../../OTHER/UMFPACK/umfpack.h"

class PFEMLinSOE;

class PFEMSolver_Umfpack : public PFEMSolver
{
public:
    PFEMSolver_Umfpack(double tol=1e-4, int niter=100, int p=0);
    virtual ~PFEMSolver_Umfpack();

    int solve();
    int setSize();
    virtual int setLinearSOE(PFEMLinSOE& theSOE);

    int sendSelf(int commitTag, Channel &theChannel);
    int recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker);  

private:
    
    void *Symbolic;
    double Control[UMFPACK_CONTROL], Info[UMFPACK_INFO];
    PFEMLinSOE* theSOE;
    int print;
    double ptol;
    int pmaxiter;
};

#endif

