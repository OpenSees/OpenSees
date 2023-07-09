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
                                                                        
// $Revision$
// $Date$
                                                                        
#ifndef PFEMSolver_Mumps_h
#define PFEMSolver_Mumps_h

//
// Written: Minjie 
//
#include <PFEMSolver.h>
#ifdef _MUMPS
#include <dmumps_c.h>
#endif

class PFEMLinSOE;

class PFEMSolver_Mumps : public PFEMSolver
{
public:
    PFEMSolver_Mumps(int r, int e, int add, int s, int p=0,
		     double tol=1e-4, int niter=100, double bittol=1e-6);
    virtual ~PFEMSolver_Mumps();

    virtual int solve();
    virtual int setSize();
    virtual int setLinearSOE(PFEMLinSOE& theSOE);

    int sendSelf(int commitTag, Channel &theChannel);
    int recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker);  

private:
    
    PFEMLinSOE* theSOE;
#ifdef _MUMPS
    DMUMPS_STRUC_C sid;

    static const int JOB_INIT = -1;
    static const int JOB_END = -2;
    static const int JOB_ANALYSIS = 1;
    static const int JOB_FACTORIZATION = 2;
    static const int JOB_SOLUTION = 3;
    static const int JOB_SOLVE = 6;
    static const int USE_COMM_WORLD = -987654;
    void ICNTL(DMUMPS_STRUC_C& id, int I, int val);
#endif

    int relax, err, add, sym, print;
    double ptol, Bitol;
    int pmaxiter;
};

#endif
