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
                                                                        
// $Revision: 1.3 $
// $Date: 2005-03-17 20:47:09 $
// $Source: /usr/local/cvs/OpenSees/SRC/system_of_eqn/linearSOE/sparseGEN/SuperLU.h,v $
                                                                        
                                                                        
#ifndef SuperLU_h
#define SuperLU_h

// File: ~/system_of_eqn/linearSOE/sparseGEN/SuperLU.h
//
// Written: fmk 
// Created: 11/96
//
// Description: This file contains the class definition for SuperLU.
// A SuperLU object can be constructed to solve a SparseGenColLinSOE
// object. It obtains the solution by making calls on the
// the SuperLU library developed at UC Berkeley by Prof. James Demmel, 
// Xiaoye S. Li and John R. Gilbert.
// The SuperLU library contains a set of subroutines to solve a sparse
// linear system  $AX=B$. It uses Gaussian elimination with partial
// pivoting (GEPP). The columns of A may be preordered before
// factorization; the preordering for sparsity is completely separate
// from the factorization and a number of ordering schemes are provided. 
//
// What: "@(#) SuperLU.h, revA"

#include <SparseGenColLinSolver.h>
#include <slu_ddefs.h>
#include <supermatrix.h>

class SuperLU : public SparseGenColLinSolver
{
  public:
    SuperLU(int permSpec = 0, 
	    double drop_tol = 0.0, 
	    int panelSize = 6, 
	    int relax = 6,
	    char symmetric = 'N');
    ~SuperLU();

    int solve(void);
    int setSize(void);

    int sendSelf(int commitTag, Channel &theChannel);
    int recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker);    
    
  protected:

  private:
    SuperMatrix A,L,U,B,AC;
    int *perm_r;
    int *perm_c;
    int *etree;
    int sizePerm;
    int relax, permSpec, panelSize;
    double drop_tol;
    char symmetric;
    superlu_options_t options;
    SuperLUStat_t stat;
};

#endif

