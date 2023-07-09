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
// $Date: 2015-3-18 11:02:54 $

// Written: Minjie Zhu
// Created: March 2015
                                                                        
                                                                        
#ifndef PFEMUnifiedSolver_Hybrid_h
#define PFEMUnifiedSolver_Hybrid_h

//
// Description: This file contains the class definition for PFEMUnifiedSolver_Hybrid.
// A PFEMUnifiedSolver_Hybrid object can be constructed to solve a PFEMUnifiedLinSOE
// object. It obtains the solution by making calls on the
// The PFEMUnifiedSolver_Hybrid uses Fractional Step Method to solve PFEM equations. 
//
// What: "@(#) PFEMUnifiedSolver_Hybrid.h, revA"

#include <LinearSOESolver.h>
#include <Vector.h>
#include <dmumps_c.h>
#include <cuda.h>
#include <cusolverDn.h>
#include <vector>

class PFEMGeneralLinSOE;

class PFEMUnifiedSolver_Hybrid : public LinearSOESolver
{
public:
    PFEMUnifiedSolver_Hybrid(int r, int e, int h, int s);
    virtual ~PFEMUnifiedSolver_Hybrid();

    int solve();
    int setSize();
    int setLinearSOE(PFEMGeneralLinSOE& theSOE);

    int sendSelf(int commitTag, Channel &theChannel);
    int recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker);  

private:

    typedef std::vector<int> Index;
    typedef std::vector<double> Value;

private:
    PFEMGeneralLinSOE* theSOE;
    
    DMUMPS_STRUC_C id;
    int myid;
    cusolverDnHandle_t handle;

    static const int JOB_INIT = -1;
    static const int JOB_END = -2;
    static const int JOB_ANALYSIS = 1;
    static const int JOB_FACTORIZATION = 2;
    static const int JOB_SOLUTION = 3;
    static const int USE_COMM_WORLD = -987654;
    void ICNTL(int I, int val);

    int relax, err, host, sym;
};

#endif

