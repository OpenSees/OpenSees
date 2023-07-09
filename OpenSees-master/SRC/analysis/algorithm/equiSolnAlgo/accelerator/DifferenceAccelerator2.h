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

// $Revision$
// $Date$
// $Source$

// Written: MHS
// Created: April 2002

// Description: This file contains the class definition for 
// DifferenceAccelerator2. 

#ifndef DifferenceAccelerator2_h
#define DifferenceAccelerator2_h

#include <Accelerator.h>
#include <IncrementalIntegrator.h>

class DifferenceAccelerator2 : public Accelerator
{
 public:
  DifferenceAccelerator2(int maxDim = 3, int tangent = CURRENT_TANGENT);
  virtual ~DifferenceAccelerator2();
  
  int newStep(LinearSOE &theSOE);
  int accelerate(Vector &v, LinearSOE &theSOE, 
		 IncrementalIntegrator &theIntegrator);
  int updateTangent(IncrementalIntegrator &theIntegrator);
  
  void Print(OPS_Stream &s, int flag=0);
  
  int sendSelf(int commitTag, Channel &theChannel);
  int recvSelf(int commitTag, Channel &theChannel, 
	       FEM_ObjectBroker &theBroker);

 protected:
  
 private:
  // Current dimension of Krylov subspace
  int dimension;
  
  // Size information
  int numEqns;
  int maxDimension;
  
  // Storage for update vectors
  Vector **v;
  
  // Storage for subspace vectors
  Vector **Av;
  
  // Array data sent to LAPACK subroutine
  double *AvData;
  double *rData;
  double *work;
  
  // Length of work array
  int lwork;

  // Which tangent to form at restart
  int theTangent;
};

#endif
