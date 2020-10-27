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
// $Source: /usr/local/cvs/OpenSees/SRC/analysis/algorithm/equiSolnAlgo/accelerator/MillerAccelerator.cpp,v $
                                                                        
// Written: MHS
// Created: April 2002

// Description: This file contains the class implementation for 
// MillerAccelerator. 

#include <MillerAccelerator.h>

#include <Vector.h>
#include <LinearSOE.h>

MillerAccelerator::MillerAccelerator(int max, double tol, int tangent)
  :Accelerator(ACCELERATOR_TAGS_Miller),
   iteration(0), numEqns(0), dimension(0), maxDimension(max), tolerance(tol),
   work(0), fData(0), theTangent(tangent)
{
  if (maxDimension < 1)
    maxDimension = 1;

  if (maxDimension > 10)
    maxDimension = 10;
}

MillerAccelerator::~MillerAccelerator()
{
  if (work != 0)
    delete [] work;

  if (fData != 0)
    delete [] fData;
}

int 
MillerAccelerator::newStep(LinearSOE &theSOE)
{
  int newNumEqns = theSOE.getNumEqn();

  if (newNumEqns != numEqns) {
    if (fData != 0) {
      delete [] fData;
      fData = 0;
    }
    if (work != 0) {
      delete [] work;
      work = 0;
    }
    numEqns = newNumEqns;
  }

  if (fData == 0)
    fData = new double [numEqns];

  // Make sure max dim <= num equations
  if (maxDimension > numEqns)
    maxDimension = numEqns;

  // Length of work array -- N*(2*MVEC+2)
  int lwork = numEqns*(2*maxDimension+2);
  
  if (work == 0)
    work = new double [lwork];

  // Reset iteration counter
  iteration = 1;
  dimension = (theTangent == CURRENT_TANGENT) ? maxDimension : 0;

  return 0;
}

#ifdef _WIN32

extern "C" int NACCEL(int *n, int *itr, int *mvec,
			       double *tol, double *u, double *f);

#define naccel_ NACCEL

#else

extern "C" int naccel_(int *n, int *itr, int *mvec,
		       double *tol, double *u, double *f);

#endif

int
MillerAccelerator::accelerate(Vector &vStar, LinearSOE &theSOE, 
			      IncrementalIntegrator &theIntegrator)
{
  // Vector for fData
  Vector fVec(fData, numEqns);
  fVec = vStar;

   // Vector size
  int N = numEqns;

  // Newton iteration count
  int ITR = iteration;

  // Maximum number of vectors used for acceleration
  int MVEC = maxDimension;

  // Tolerance for dropping linearly dependent vectors
  double TOL = tolerance;

  // Work space
  double *U = work;

  // Accelerated correction on output
  double *F = fData;

  // Call the naccel subroutine
  naccel_(&N, &ITR, &MVEC, &TOL, U, F);

  vStar = fVec;

  iteration++;
  dimension++;

  return 0; 
}

int
MillerAccelerator::updateTangent(IncrementalIntegrator &theIntegrator)
{
  /*
  if (dimension >= maxDimension) {
    dimension = 0;
    if (theTangent != NO_TANGENT) {
      iteration = 1; // reset Newton iteration if tangent is formed
      theIntegrator.formTangent(theTangent);
      //opserr << "MillerAccelerator::updateTangent() tangent formed" << endln;
      return 1;
    }
    else
      return 0;
  }
  else
    return 0;
  */

  if (dimension < maxDimension)
    return 0;

  switch (theTangent) {
  case CURRENT_TANGENT:
    iteration = 1;
    dimension = 0;
    theIntegrator.formTangent(CURRENT_TANGENT);
    return 1;
    break;
  case INITIAL_TANGENT:
    dimension = 0;
    theIntegrator.formTangent(INITIAL_TANGENT);
    return 0;
    break;
  case NO_TANGENT:
    dimension = 0;
    return 0;
    break;
  default:
    return 0;
  }
}

bool
MillerAccelerator::updateTangent(void)
{
  if (dimension >= maxDimension) {
    dimension = 0;
    iteration = 1; // reset Newton iteration if tangent is formed
    return true;
  }
  else
    return false;
}

void
MillerAccelerator::Print(OPS_Stream &s, int flag)
{
  s << "MillerAccelerator" << endln;
  s << "\tMax subspace dimension: " << maxDimension << endln;
  s << "\tTolerance: " << tolerance << endln;
}

int
MillerAccelerator::sendSelf(int commitTag, Channel &theChannel)
{
  return -1;
}

int
MillerAccelerator::recvSelf(int commitTag, Channel &theChannel, 
			    FEM_ObjectBroker &theBroker)
{
  return -1;
}
