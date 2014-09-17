// $Revision: 1.0 $
// $Date: 2005/12/06 17:58:00 $
// $Source: /usr/local/cvs/OpenSees/SRC/system_of_eqn/linearSOE/diagonal/MPIDiagonalSolver.h,v $


/* 
 * @author: George petropoulos <gnp>
 *
 * @Description: This file contains the class definition for MPIDiagonalSolver. 
 *               MPIDiagonalSolver is a concrete class for solving a system stored using
 *               a MPIDiagonalSOE. The solve() method is overwritten to allow subclassing.
 *               
 * @Date: 12/05
 *
 * Copyright: ALL RIGHTS RESERVED BY AUTHOR
 *
 */


#ifndef MPIDiagonalSolver_h
#define MPIDiagonalSolver_h

#include <mpi.h>
#include <LinearSOESolver.h>
#include <ID.h>

class MPIDiagonalSOE;

class MPIDiagonalSolver : public LinearSOESolver
{
  public:
    MPIDiagonalSolver(int classTag);    
    MPIDiagonalSolver(double minDiagTol=1.0e-18);    
    virtual ~MPIDiagonalSolver();

    virtual int solve(void);
    virtual int setSize(void);
    virtual int setLinearSOE(MPIDiagonalSOE &theSOE);
    int sendSelf(int commitTag, Channel &theChannel);
    int recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker);    

    //    void intersectionsAB(ID& arrayA, int* arrayB, int sizeA, int sizeB, double* A, double* sharedA, double* B, double* sharedB);
    void intersectionsAB(ID& arrayA, int* arrayB, int sizeA, int sizeB, double* A, double* sharedA, double* B, double* sharedB, int* storage, int neighbor_pid);

  protected:
    MPIDiagonalSOE *theSOE;

  private:
    double minDiagTol;
    bool notSet;
};

#endif

