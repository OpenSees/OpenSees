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
                                                                        
// $Revision: 1.1 $
// $Date: 2005-12-06 22:21:03 $
// $Source: /usr/local/cvs/OpenSees/SRC/system_of_eqn/linearSOE/sparseGEN/DistributedSuperLU.h,v $
                                                                        
                                                                        
#ifndef DistributedSuperLU_h
#define DistributedSuperLU_h

// Written: fmk 
//
// Description: This file contains the class definition for DistributedSuperLU.
// A DistributedSuperLU object can be constructed to solve a SparseGenColLinSOE
// object. It obtains the solution by making calls on the
// the SuperLU library developed at UC Berkeley by Prof. James Demmel, 
// Xiaoye S. Li and John R. Gilbert.
// The SuperLU library contains a set of subroutines to solve a sparse
// linear system  $AX=B$. It uses Gaussian elimination with partial
// pivoting (GEPP). The columns of A may be preordered before
// factorization; the preordering for sparsity is completely separate
// from the factorization and a number of ordering schemes are provided. 
//
// What: "@(#) DistributedSuperLU.h, revA"


#include <SparseGenColLinSolver.h>

class DistributedSuperLU : public SparseGenColLinSolver
{
  public:
    DistributedSuperLU(int npRow, int npCol);
    DistributedSuperLU();
    ~DistributedSuperLU();

    int solve(void);
    int setSize(void);

    virtual int setProcessID(int domainTag);
    virtual int setChannels(int numChannels, Channel **theChannels);
    int sendSelf(int commitTag, Channel &theChannel);
    int recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker);    
    
  protected:

  private:

    bool gridInit;
    int npRow, npCol;

    int processID;
    int numChannels;
    Channel **theChannels;

    double *b;
    int *rowA;
};

#endif

