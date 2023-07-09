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
// $Date: 2008-09-16 18:15:42 $
// $Source: /usr/local/cvs/OpenSees/SRC/analysis/algorithm/equiSolnAlgo/accelerator/PeriodicAccelerator.cpp,v $
                                                                        
// Written: MHS
// Created: April 2002

// Description: This file contains the class implementation for 
// PeriodicAccelerator. 

#include <PeriodicAccelerator.h>

#include <Vector.h>
#include <LinearSOE.h>
#include <IncrementalIntegrator.h>

#include <ID.h>
#include <Channel.h>

PeriodicAccelerator::PeriodicAccelerator(int iter, int tangent)
  :Accelerator(ACCELERATOR_TAGS_Periodic),
   iteration(0), totalIter(0), maxIter(iter), theTangent(tangent)
{
  if (maxIter < 1)
    maxIter = 1;
}

PeriodicAccelerator::~PeriodicAccelerator()
{

}

int 
PeriodicAccelerator::newStep(LinearSOE &theSOE)
{
  totalIter = 0;

  // Reset iteration counter
  iteration = (theTangent == CURRENT_TANGENT) ? maxIter : 0;

  return 0;
}

int
PeriodicAccelerator::accelerate(Vector &vStar, LinearSOE &theSOE,
				IncrementalIntegrator &theIntegrator)
{
  iteration++;
  totalIter++;

  return 0; 
}

int
PeriodicAccelerator::updateTangent(IncrementalIntegrator &theIntegrator)
{
  /*
  if (theTangent == NO_TANGENT)
    return 0;

  else if (theTangent == SECOND_TANGENT) {
    if (totalIter == maxIter) {
      theIntegrator.formTangent(CURRENT_TANGENT);
      return 1;
    }
    else
      return 0;
  }

  else { // CURRENT_TANGENT or INITIAL_TANGENT
    if (iteration >= maxIter) {
      iteration = 0;
      theIntegrator.formTangent(theTangent);
      if (theTangent == CURRENT_TANGENT)
	return 1;
      else
	return 0;
    }
    else
      return 0;
  }
  */

  if (iteration < maxIter)
    return 0;

  switch (theTangent) {
  case CURRENT_TANGENT:
    iteration = 0;
    theIntegrator.formTangent(CURRENT_TANGENT);
    return 1;
    break;
  case INITIAL_TANGENT:
    iteration = 0;
    theIntegrator.formTangent(INITIAL_TANGENT);
    return 0;
    break;
  case NO_TANGENT:
    iteration = 0;
    return 0;
    break;
  default:
    return 0;
  }
}

bool
PeriodicAccelerator::updateTangent(void)
{
  if (iteration > maxIter) {
    iteration = 0;
    return true;
  }
  else 
    return false;
}

void
PeriodicAccelerator::Print(OPS_Stream &s, int flag)
{
  s << "PeriodicAccelerator" << endln;
  s << "\tIterations till restart: " << maxIter << endln;
}

int
PeriodicAccelerator::sendSelf(int commitTag, Channel &theChannel)
{
  static ID data(2);
  data(0) = theTangent;
  data(1) = maxIter;
  return theChannel.sendID(0, commitTag, data);
  
}

int
PeriodicAccelerator::recvSelf(int commitTag, Channel &theChannel, 
			    FEM_ObjectBroker &theBroker)
{
  static ID data(2);
  int res = theChannel.recvID(0, commitTag, data);
  theTangent = data(0);
  maxIter = data(1);
  return res;
}
