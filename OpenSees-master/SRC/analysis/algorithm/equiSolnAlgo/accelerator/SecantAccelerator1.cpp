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
                                                                        
// $Revision: 1.2 $
// $Date: 2008-08-26 18:26:07 $
// $Source: /usr/local/cvs/OpenSees/SRC/analysis/algorithm/equiSolnAlgo/accelerator/SecantAccelerator1.cpp,v $
                                                                        
// Written: MHS
// Created: April 2002

// Description: This file contains the class implementation for 
// SecantAccelerator1. 

#include <SecantAccelerator1.h>

#include <Vector.h>
#include <LinearSOE.h>
#include <IncrementalIntegrator.h>

SecantAccelerator1::SecantAccelerator1(int maxIter, int tangent)
  :Accelerator(ACCELERATOR_TAGS_Secant),
  iteration(0), numEqns(0), R1(0.0), R2(0.0),
   vOld(0), rOld(0), maxIterations(maxIter), theTangent(tangent),
   cutOut(false)
{

}

SecantAccelerator1::SecantAccelerator1(int maxIter, int tangent,
				       double r1, double r2)
  :Accelerator(ACCELERATOR_TAGS_Secant),
   iteration(0), numEqns(0), R1(r1), R2(r2),
   vOld(0), rOld(0), maxIterations(maxIter), theTangent(tangent),
   cutOut(true)
{

}

SecantAccelerator1::~SecantAccelerator1()
{
  if (vOld != 0)
    delete vOld;

  if (rOld != 0)
    delete rOld;
}

int 
SecantAccelerator1::newStep(LinearSOE &theSOE)
{
  int newNumEqns = theSOE.getNumEqn();

  if (vOld != 0 && vOld->Size() != newNumEqns) {
    delete vOld;
    vOld = 0;
  }

  if (rOld != 0 && rOld->Size() != newNumEqns) {
    delete rOld;
    rOld = 0;
  }

  numEqns = newNumEqns;

  if (vOld == 0)
    vOld = new Vector(numEqns);

  if (rOld == 0)
    rOld = new Vector(numEqns);

  // Reset iteration counter
  iteration = 0;

  return 0;
}

int
SecantAccelerator1::accelerate(Vector &vStar, LinearSOE &theSOE, 
			       IncrementalIntegrator &theIntegrator)
{
  // Current right hand side
  const Vector &rNew  = theSOE.getB();

  // No acceleration on first iteration
  if (iteration == 0) {
    // do nothing
  }
  else {
    // Store gamma in rOld ... \gamma = R_i - R_{i-1}
    rOld->addVector(-1.0, rNew, 1.0);
    
    double den = 1.0 / ((*vOld)^(*rOld));
    double C   = ((*vOld)^rNew) * den;
    double A   = 1.0-C;

    // Check "cut-out" criteria
    if (cutOut && (A > R1 || A < 1.0/R1)) {
      // do nothing
      //opserr << "SecantAccelerator1::accelerate() -- cut out, A = " << A << endln;
    }
    else
      vStar *= A;
  }

  // Store old values for next iteration
  *rOld = rNew;
  *vOld = vStar;

  iteration++;

  return 0; 
}

int
SecantAccelerator1::updateTangent(IncrementalIntegrator &theIntegrator)
{
  if (iteration > maxIterations) {
    //opserr << "SecantAccelerator1::updateTangent() tangent formed" << endln;
    iteration = 0;
    if (theTangent != NO_TANGENT) {
      theIntegrator.formTangent(theTangent);
      return 1;
    }
    else
      return 0;
  }
  else
    return 0;
}

void
SecantAccelerator1::Print(OPS_Stream &s, int flag)
{
  s << "SecantAccelerator1" << endln;
  s << "\tIterations till reform tangent: " << maxIterations << endln;
  if (cutOut)
    s << "\tCut-out factors -- R1: " << R1 << ", R2: " << R2 << endln;
  else
    s << "\tNo cut-outs" << endln;
}

int
SecantAccelerator1::sendSelf(int commitTag, Channel &theChannel)
{
  return -1;
}

int
SecantAccelerator1::recvSelf(int commitTag, Channel &theChannel, 
			     FEM_ObjectBroker &theBroker)
{
  return -1;
}
