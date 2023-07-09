/* ****************************************************************** **
**    OpenSees - Open System for Earthquake Engineering Simulation    **
**          Pacific Earthquake Engineering Research Center            **
**                                                                    **
**                                                                    **
** (C) Copyright 2001, The Regents of the University of California    **
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
// $Date: 2008-12-03 23:44:54 $
// $Source: /usr/local/cvs/OpenSees/SRC/analysis/algorithm/equiSolnAlgo/accelerator/MillerAccelerator.h,v $

// Written: MHS
// Created: April 2002

// Description: This file contains the class definition for 
// MillerAccelerator. 

#ifndef MillerAccelerator_h
#define MillerAccelerator_h

#include <Accelerator.h>
#include <IncrementalIntegrator.h>

class MillerAccelerator: public Accelerator
{
 public:
  MillerAccelerator(int maxDim = 3, double tol = 1.0e-2,
		    int tangent = CURRENT_TANGENT);
  virtual ~MillerAccelerator();
  
  int newStep(LinearSOE &theSOE);
  int accelerate(Vector &v, LinearSOE &theSOE, 
		 IncrementalIntegrator &theIntegrator);
  int updateTangent(IncrementalIntegrator &theIntegrator);
  bool updateTangent(void);
  
  int getTangent(void) {return theTangent;}

  void Print(OPS_Stream &s, int flag=0);
  
  int sendSelf(int commitTag, Channel &theChannel);
  int recvSelf(int commitTag, Channel &theChannel, 
	       FEM_ObjectBroker &theBroker);

 protected:
  
 private:
  // Iteration count
  int iteration;
  
  // Number of equations
  int numEqns;
  
  // Current dimension of the subspace
  int dimension;

  // Maximum number of vectors to use for acceleration
  int maxDimension;
  
  // Tolerance for dropping vectors
  double tolerance;
  
  // Workspace array -- length N*(2*MVEC+2)
  double *work;
  
  // Storage for f(x) := J^{-1}R(x)
  double *fData;
  
  int theTangent;
};

#endif
