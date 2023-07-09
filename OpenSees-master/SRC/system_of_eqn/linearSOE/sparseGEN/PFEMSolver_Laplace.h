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


#ifndef PFEMSolver_Laplace_h
#define PFEMSolver_Laplace_h

// File: ~/system_of_eqn/linearSOE/sparseGEN/PFEMSolver_Laplace.h
//
// Written: Minjie
// Created: Sep 17 2012
//
// Description: This file contains the class definition for PFEMSolver_Laplace.
// A PFEMSolver_Laplace object can be constructed to solve a PFEMLinSOE
// object. It obtains the solution by making calls on the
// The PFEMSolver_Laplace uses Fractional Step Method to solve PFEM equations.
//
// What: "@(#) PFEMSolver_Laplace.h, revA"

#include <PFEMSolver.h>
extern "C" {
#include <cs.h>
}
#include "../../../../OTHER/UMFPACK/umfpack.h"

class PFEMLinSOE;

class PFEMSolver_Laplace : public PFEMSolver
{
public:
    PFEMSolver_Laplace(bool once);
    virtual ~PFEMSolver_Laplace();

    int solve();
    int setSize();
    virtual int setLinearSOE(PFEMLinSOE& theSOE);

    int sendSelf(int commitTag, Channel &theChannel);
    int recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker);

private:

    void *MSym, *MNum, *LSym, *LNum;
    double Control[UMFPACK_CONTROL], Info[UMFPACK_INFO];
    PFEMLinSOE* theSOE;
    bool numonce, factored;
};

#endif
