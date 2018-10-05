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
// $Date: 2002-06-08 16:17:26 $
// $Source: /usr/local/cvs/OpenSees/SRC/system_of_eqn/linearSOE/itpack/ItpackLinSolver.h,v $

#ifndef ItpackLinSolver_h
#define ItpackLinSolver_h

// Written: MHS
// Created: Sept 2001
//
// Description: This file contains the class definition for ItpackLinSolver.
// ItpackLinSolver is a concrete subclass of LinearSOE. It stores full
// unsymmetric linear system of equations using 1d arrays in Fortran style

#include <LinearSOESolver.h>

// Adaptive methods
#define ItpackJCG         1
#define ItpackJSI         2
#define ItpackSOR         3
#define ItpackSSORCG      4
#define ItpackSSORSI      5
#define ItpackRSCG        6
#define ItpackRSSI        7

// Textbook methods
#define ItpackJ           8
#define ItpackGS          9
#define ItpackSORFixed   10
#define ItpackSSORFixed  11
#define ItpackRS         12

class ItpackLinSOE;

class ItpackLinSolver : public LinearSOESolver
{
 public:
  ItpackLinSolver(int method, int maxIter = 100, double omega = 1.0);
  ItpackLinSolver();
  virtual ~ItpackLinSolver();
  
  int solve(void);
  int setSize(void);
  int setLinearSOE(ItpackLinSOE &theSOE);
  
  int sendSelf(int commitTag, Channel &theChannel);
  int recvSelf(int commitTag, Channel &theChannel, 
	       FEM_ObjectBroker &theBroker);
  
 private:
  // Sparse Ax=b represntation
  ItpackLinSOE *theSOE;
  
  int *IA;
  int *JA;
  
  // Size of system, i.e., the number of equations
  int n;
  
  // Parameter arrays sent to ITPACK subroutines
  int iparm[12];
  double rparm[12];
  
  // Workspace arrays sent to ITPACK subroutines
  int *iwksp;
  double *wksp;
  
  // Length of wksp array
  int nwksp;
  
  // Maximum number of iterations
  int maxIter;

  // Integer indicating which method to use
  int method;

  // Parameter for SOR and SSOR fixed omega methods
  double omega;
};

#endif
