// File: ~/system_of_eqn/pardiso/PARDISOLinSolver.h
//
// Written: M. Salehi opensees.net@gmail.com
// website : http://opensees.net
// Created: 02/19
// Revision: A
//
// Description: This file contains the class definition for 
// PARDISOLinSolver. It solves the Sparse General SOE by calling
// some "C" functions. The solver used here is generalized sparse
// solver. The user can choose three different ordering schema.
//
// What: "@(#) PARDISOLinSolver.h, revA"


#ifndef PARDISOSymLinSolver_h
#define PARDISOSymLinSolver_h

#include <LinearSOESolver.h>
#include <PARDISOSymLinSOE.h>

//https://software.intel.com/en-us/articles/pardiso-parameter-table

class PARDISOSymLinSolver : public LinearSOESolver
{
  public:
	  PARDISOSymLinSolver();
    ~PARDISOSymLinSolver();

    int solve(void);
    int setSize(void);

    int setLinearSOE(PARDISOSymLinSOE &theSOE);
	
    int sendSelf(int cTag, Channel &theChannel);
    int recvSelf(int cTag, 
		 Channel &theChannel, 
		 FEM_ObjectBroker &theBroker);
  protected:

  private:
	  PARDISOSymLinSOE *theSOE;
};

#endif

