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
// $Date: 2005-06-20 21:35:18 $
// $Source: /usr/local/cvs/OpenSees/SRC/system_of_eqn/linearSOE/diagonal/DistributedDiagonalSolver.cpp,v $

// Written: fmk 
// Created: 05/05
//
// What: "@(#) DistributedDiagonalSolver.C, revA"

#include <DistributedDiagonalSolver.h>
#include <DistributedDiagonalSOE.h>
#include <Channel.h>

DistributedDiagonalSolver::DistributedDiagonalSolver(int classTag)
:LinearSOESolver(classTag),
 theSOE(0), minDiagTol(0.0)
{
  
}    

DistributedDiagonalSolver::DistributedDiagonalSolver(double tol)
:LinearSOESolver(SOLVER_TAGS_DistributedDiagonalSolver),
 theSOE(0), minDiagTol(tol)
{

}    

DistributedDiagonalSolver::~DistributedDiagonalSolver()    
{

}    

int 
DistributedDiagonalSolver::setLinearSOE(DistributedDiagonalSOE &theProfileSPDSOE)
{
  
  theSOE = &theProfileSPDSOE;
  return 0;
}

int 
DistributedDiagonalSolver::setSize(void)
{
  return 0;
}


int 
DistributedDiagonalSolver::solve(void)
{
  int size = theSOE->size;
  int processID = theSOE->processID;

  Channel **theChannels = theSOE->theChannels;
  int numChannels = theSOE->numChannels;

  int numShared = theSOE->numShared;
  ID &myDOFs = theSOE->myDOFs;
  ID &myDOFsShared = theSOE->myDOFsShared;
  
  double *X = theSOE->X;
  double *B = theSOE->B;
  double *A = theSOE->A;
  double *dataShared = theSOE->dataShared;
  Vector *vectShared = theSOE->vectShared;

  //
  // first copy A & B contributions to sharedData
  //

  for (int i=0; i<numShared; i++)
    dataShared[i] = 0.0;

  // assuming numShared < size .. could do an if statement

  for (int i=0; i<numShared; i++) {
    int dof = myDOFsShared(i);
    int loc = myDOFs.getLocation(dof);
    if (loc >= 0) {
      dataShared[i] = A[loc];
      dataShared[i+numShared] = B[loc];
    }
  }

  //
  // use P0 to gather & send back out
  //

  if (numShared != 0) {
    if (processID != 0) {
      Channel *theChannel = theChannels[0];
      theChannel->sendVector(0, 0, *vectShared);
      theChannel->recvVector(0, 0, *vectShared);
    } 
    else {

      static Vector otherShared(1);
      otherShared.resize(2*numShared);
      for (int i=0; i<numChannels; i++) {
	Channel *theChannel = theChannels[i];
	theChannel->recvVector(0, 0, otherShared);
	*vectShared += otherShared;
      }
      for (int i=0; i<numChannels; i++) {
	Channel *theChannel = theChannels[i];
	theChannel->sendVector(0, 0, *vectShared);
      }
    }
  }
  
  
  //
  // set the corresponding A & B entries
  //
  
  
  for (int i=0; i<numShared; i++) {
    int dof = myDOFsShared(i);
    int loc = myDOFs.getLocation(dof);
    if (loc >= 0) {
      A[loc] = dataShared[i];
      B[loc] = dataShared[i+numShared];
    }
  }  

  //
  // now solve
  //
  
  for (int i=0; i<size; i++) {
    X[i] = B[i]/A[i];
  }

  return 0;
}

int
DistributedDiagonalSolver::sendSelf(int cTag,
			       Channel &theChannel)
{
  static Vector data(1);
  data(0) = minDiagTol;
  return theChannel.sendVector(0, cTag, data);
}


int 
DistributedDiagonalSolver::recvSelf(int cTag,
			       Channel &theChannel, 
			       FEM_ObjectBroker &theBroker)
{
  static Vector data(1);
  theChannel.recvVector(0, cTag, data);

  minDiagTol = data(0);
  return 0;
}

