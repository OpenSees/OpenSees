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
                                                                       
// Written: fmk 
// Created: 11/14

// this class is based on the S5 solver, written to work with Cula
// S4 library (S5 not available for apple at time of writing!) 

#ifndef _CULASPARSESOLVER4
#define _CULASPARSESOLVER4

#include <SparseGenRowLinSolver.h>
#include <cula_sparse.h>
class SparseGenRowLinSOE;


class CulaSparseSolverS4 : public SparseGenRowLinSolver
{
 public:
  CulaSparseSolverS4(double tolerance,
		     int maxNumIterations,
		     int preCond,
		     int solver);
  
  ~CulaSparseSolverS4(void);
  
  int solve(void);
  int setSize(void);
  
  int sendSelf(int commitTag, Channel &theChannel);
  int recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker);    
  
  int setLinearSOE(SparseGenRowLinSOE &theLinearSOE);
  
 private:
  int n;// order of matrix
  int nnz;// Number of no-zero ones in A;
  
  double tolerance;
  int maxNumIterations;
  int preCond;			
  int solver;		      
  
  SparseGenRowLinSOE *theSOE;
  
  // status returned by each and every cula routine
  culaStatus status;
  
  // configuration structures
  culaIterativeConfig config;
  
  // character buffer used for results and error messages
  char buf[512];
  
  // result structure
  culaIterativeResult result;
  
  // structures for the different solvers - all only containing an integer!
  culaCgOptions cg;
  culaBicgOptions biCG;
  culaBicgstabOptions biCGstab;
  culaBicgstablOptions biCGstabl;
  culaGmresOptions gmRes;
  culaMinresOptions minRes;
  
  // structures for the different preconditioners
  culaJacobiOptions jacobi;
  culaBlockjacobiOptions blockJacobi;
  culaIlu0Options ilu0;   
};

#endif
