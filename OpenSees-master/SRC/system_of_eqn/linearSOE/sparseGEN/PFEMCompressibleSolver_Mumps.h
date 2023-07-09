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
// $Source: /usr/local/cvs/OpenSees/SRC/system_of_eqn/linearSOE/sparseGEN/PFEMCompressibleSolver_Mumps.h,v $
                                                                        
                                                                        
#ifndef PFEMCompressibleSolver_Mumps_h
#define PFEMCompressibleSolver_Mumps_h

// File: ~/system_of_eqn/linearSOE/sparseGEN/PFEMCompressibleSolver_Mumps.h
//
// Written: Minjie 
// Created: Sep 17 2012
//
// Description: This file contains the class definition for PFEMCompressibleSolver_Mumps.
// A PFEMCompressibleSolver_Mumps object can be constructed to solve a PFEMCompressibleLinSOE
// object. It obtains the solution by making calls on the
// The PFEMCompressibleSolver_Mumps uses quasi-incompressible method to solve PFEM equations. 
//
// What: "@(#) PFEMCompressibleSolver_Mumps.h, revA"

#include <PFEMCompressibleSolver.h>
#include <Vector.h>
extern "C" {
#include <cs.h>
}
#include <dmumps_c.h>

class PFEMCompressibleLinSOE;

class PFEMCompressibleSolver_Mumps : public PFEMCompressibleSolver
{
public:
    PFEMCompressibleSolver_Mumps(int r, int e, int s);
    virtual ~PFEMCompressibleSolver_Mumps();

    int solve();
    int setSize();
    int setLinearSOE(PFEMCompressibleLinSOE& theSOE);

    int sendSelf(int commitTag, Channel &theChannel);
    int recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker);  

private:
    
    PFEMCompressibleLinSOE* theSOE;

    DMUMPS_STRUC_C sid;
    int myid;

    static const int JOB_INIT = -1;
    static const int JOB_END = -2;
    static const int JOB_ANALYSIS = 1;
    static const int JOB_FACTORIZATION = 2;
    static const int JOB_SOLUTION = 3;
    static const int JOB_SOLVE = 6;
    static const int USE_COMM_WORLD = -987654;
    void ICNTL(DMUMPS_STRUC_C& id, int I, int val);

    int relax, err;
};

#else
class PFEMCompressibleSolver_Mumps{};

#endif

