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
// $Source: /usr/local/cvs/OpenSees/SRC/analysis/algorithm/equiSolnAlgo/accelerator/RaphsonAccelerator.h,v $

// Written: MHS
// Created: April 2002

// Description: This file contains the class definition for 
// RaphsonAccelerator. 

#ifndef RaphsonAccelerator_h
#define RaphsonAccelerator_h

#include "Accelerator.h"
#include <IncrementalIntegrator.h>

class RaphsonAccelerator: public Accelerator
{
 public:
  RaphsonAccelerator(int tangent = CURRENT_TANGENT);
  virtual ~RaphsonAccelerator();
  
  int newStep(LinearSOE &theSOE);
  int accelerate(Vector &v, LinearSOE &theSOE, 
		 IncrementalIntegrator &theIntegrator);
  int updateTangent(IncrementalIntegrator &theIntegrator);
  bool updateTangent(void) {return true;}

  int getTangent(void) {return theTangent;}

  void Print(OPS_Stream &s, int flag=0);
  
  int sendSelf(int commitTag, Channel &theChannel);
  int recvSelf(int commitTag, Channel &theChannel, 
	       FEM_ObjectBroker &theBroker);

 protected:
  
 private:
  // Flag indicating which tangent to form
  int theTangent;

  // Total number of iterations for the time step
  int totalIter;
};

#endif
