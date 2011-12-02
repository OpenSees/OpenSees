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
// $Source: /usr/local/cvs/OpenSees/SRC/analysis/algorithm/equiSolnAlgo/accelerator/RaphsonAccelerator.cpp,v $
                                                                        
// Written: MHS
// Created: April 2002

// Description: This file contains the class implementation for 
// RaphsonAccelerator. 

#include <RaphsonAccelerator.h>

#include <Vector.h>
#include <LinearSOE.h>
#include <IncrementalIntegrator.h>

RaphsonAccelerator::RaphsonAccelerator(int tangent)
  :Accelerator(ACCELERATOR_TAGS_Raphson), theTangent(tangent), totalIter(0)
{

}

RaphsonAccelerator::~RaphsonAccelerator()
{

}

int 
RaphsonAccelerator::newStep(LinearSOE &theSOE)
{
  totalIter = 0;

  return 0;
}

int
RaphsonAccelerator::accelerate(Vector &vStar, LinearSOE &theSOE,
			       IncrementalIntegrator &theIntegrator)
{
  totalIter++;

  return 0; 
}

int
RaphsonAccelerator::updateTangent(IncrementalIntegrator &theIntegrator)
{
  /*
  if (theTangent == NO_TANGENT)
    return 0;

  else if (theTangent == SECOND_TANGENT) {
    if (totalIter == 1) {
      theIntegrator.formTangent(CURRENT_TANGENT);
      return 1;
    }
    else 
      return 0;
  }

  else { // CURRENT_TANGENT or INITIAL_TANGENT
    theIntegrator.formTangent(theTangent);
    return 1;
  }
  */

  switch (theTangent) {
  case CURRENT_TANGENT:
    theIntegrator.formTangent(CURRENT_TANGENT);
    return 1;
    break;
  case INITIAL_TANGENT:
    theIntegrator.formTangent(INITIAL_TANGENT);
    return 0;
    break;
  default:
    return 0;
    break;
  }
}

void
RaphsonAccelerator::Print(OPS_Stream &s, int flag)
{
  s << "RaphsonAccelerator" << endln;
}

int
RaphsonAccelerator::sendSelf(int commitTag, Channel &theChannel)
{
  return -1;
}

int
RaphsonAccelerator::recvSelf(int commitTag, Channel &theChannel, 
			    FEM_ObjectBroker &theBroker)
{
  return -1;
}
