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
// $Source: /usr/local/cvs/OpenSees/SRC/system_of_eqn/linearSOE/sparseGEN/PFEMCompressibleSolver.h,v $
                                                                        
                                                                        
#ifndef PFEMCompressibleSolver_h
#define PFEMCompressibleSolver_h

// File: ~/system_of_eqn/linearSOE/sparseGEN/PFEMCompressibleSolver.h
//
// Written: Minjie 
// Created: Sep 17 2012
//
// Description: This file contains the class definition for PFEMCompressibleSolver.
// A PFEMCompressibleSolver object can be constructed to solve a PFEMCompressibleLinSOE
// object. It obtains the solution by making calls on the
// The PFEMCompressibleSolver uses Fractional Step Method to solve PFEM equations. 
//
// What: "@(#) PFEMCompressibleSolver.h, revA"

#include <LinearSOESolver.h>
#include <Vector.h>
extern "C" {
#include <cs.h>
}
#include "../../../../OTHER/UMFPACK/umfpack.h"


class PFEMCompressibleLinSOE;

class PFEMCompressibleSolver : public LinearSOESolver
{
public:
    PFEMCompressibleSolver();
    virtual ~PFEMCompressibleSolver();

    int solve();
    int setSize();
    virtual int setLinearSOE(PFEMCompressibleLinSOE& theSOE);

    int sendSelf(int commitTag, Channel &theChannel);
    int recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker);  

private:
    
    PFEMCompressibleLinSOE* theSOE;
    double Control[UMFPACK_CONTROL], Info[UMFPACK_INFO];
};

#endif

