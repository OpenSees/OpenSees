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
// $Source: /usr/local/cvs/OpenSees/SRC/system_of_eqn/linearSOE/sparseGEN/PFEMSolver_Mumps.h,v $
                                                                        
                                                                        
#ifndef PFEMSolver_Mumps_h
#define PFEMSolver_Mumps_h

// File: ~/system_of_eqn/linearSOE/sparseGEN/PFEMSolver_Mumps.h
//
// Written: Minjie 
// Created: Sep 17 2012
//
// Description: This file contains the class definition for PFEMSolver_Mumps.
// A PFEMSolver_Mumps object can be constructed to solve a PFEMLinSOE
// object. It obtains the solution by making calls on the
// The PFEMSolver_Mumps uses Fractional Step Method to solve PFEM equations. 
//
// What: "@(#) PFEMSolver_Mumps.h, revA"

#include <PFEMSolver.h>
#include <dmumps_c.h>

class PFEMLinSOE;

class PFEMSolver_Mumps : public PFEMSolver
{
public:
    PFEMSolver_Mumps(int r, int e, int h, int s);
    virtual ~PFEMSolver_Mumps();

    virtual int solve();
    virtual int setSize();
    virtual int setLinearSOE(PFEMLinSOE& theSOE);

    int sendSelf(int commitTag, Channel &theChannel);
    int recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker);  

private:
    
    PFEMLinSOE* theSOE;

    DMUMPS_STRUC_C sid, pid;
    int myid;

    static const int JOB_INIT = -1;
    static const int JOB_END = -2;
    static const int JOB_ANALYSIS = 1;
    static const int JOB_FACTORIZATION = 2;
    static const int JOB_SOLUTION = 3;
    static const int JOB_SOLVE = 6;
    static const int USE_COMM_WORLD = -987654;
    void ICNTL(DMUMPS_STRUC_C& id, int I, int val);

    int relax, err, host, sym;
};

#endif
