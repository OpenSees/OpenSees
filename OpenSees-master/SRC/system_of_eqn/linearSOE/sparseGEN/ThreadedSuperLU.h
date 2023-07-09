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
                                                                        
// $Revision: 1.2 $
// $Date: 2002-01-25 20:27:11 $
// $Source: /usr/local/cvs/OpenSees/SRC/system_of_eqn/linearSOE/sparseGEN/ThreadedSuperLU.h,v $
                                                                        
                                                                        
#ifndef ThreadedSuperLU_h
#define ThreadedSuperLU_h

// File: ~/system_of_eqn/linearSOE/sparseGEN/ThreadedSuperLU.h
//
// Written: fmk 
// Created: 11/96
// Revision: A
//
// Description: This file contains the class interface for ThreadedSuperLU.
// This is a class that uses the threads version of SuperLU
//
// What: "@(#) ThreadedSuperLU.h, revA"

#include <SparseGenColLinSolver.h>
#include <supermatrix.h>
#include <pdsp_defs.h>

class ThreadedSuperLU : public SparseGenColLinSolver
{
  public:
    ThreadedSuperLU(int numThreads = 2,
		    int permSpec = 0, int panelSize = 6, 
		    int relax = 6, double thresh = 0.0);     

    ~ThreadedSuperLU();

    int solve(void);
    int setSize(void);

    int sendSelf(int commitTag, Channel &theChannel);
    int recvSelf(int commitTag, Channel &theChannel, 
		 FEM_ObjectBroker &theBroker);        
    
  protected:

  private:
    SuperMatrix A,L,U,B,AC;
    int *perm_r;
    int *perm_c;
    int *etree;
    int sizePerm;
    int relax, permSpec, panelSize;
    float thresh;

    int numThreads;
    pdgstrf_options_t pdgstrf_options;
    yes_no_t refact, usepr;
    fact_t fact;
    trans_t trans;
    void *work;
    int lwork;
    Gstat_t gStat;
};

#endif


