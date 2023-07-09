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
// $Date: 2007-10-26 03:56:45 $
// $Source: /usr/local/cvs/OpenSees/SRC/analysis/algorithm/equiSolnAlgo/accelerator/SecantAccelerator2.h,v $

// Written: MHS
// Created: April 2002

// Description: This file contains the class definition for 
// SecantAccelerator2. 

#ifndef SecantAccelerator2_h
#define SecantAccelerator2_h

#include <Accelerator.h>
#include <IncrementalIntegrator.h>

class SecantAccelerator2: public Accelerator
{
 public:
  SecantAccelerator2(int maxIter, int tangent);
  SecantAccelerator2(int maxIter, int tangent, double r1, double r2);
  virtual ~SecantAccelerator2();
  
  int newStep(LinearSOE &theSOE);
  int accelerate(Vector &v, LinearSOE &theSOE, 
		 IncrementalIntegrator &theIntegrator);
  int updateTangent(IncrementalIntegrator &theIntegrator);
  
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
  
  // "Cut-out" factors
  double R1;
  double R2;
  
  // Correction and RHS vectors from last iteration
  Vector *vOld;
  Vector *rOld;
  
  int maxIterations;
  int theTangent;
  bool cutOut;
};

#endif
