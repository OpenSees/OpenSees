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
// $Date: 2005-02-09 19:58:47 $
// $Source: /usr/local/cvs/OpenSees/SRC/system_of_eqn/linearSOE/diagonal/DiagonalDirectSolver.cpp,v $

// Written: fmk 
// Created: Jan 2005
// Revision: A
//
// Description: This file contains the implementation for ProfileSPDLinSOESolver
// Description: This file contains the class definition for 
// DiagonalDirectSolver. DiagonalDirectSolver is a subclass 
// of LinearSOESOlver. It solves a DiagonalSOE object directly!

// What: "@(#) DiagonalDirectSolver.C, revA"

#include <DiagonalDirectSolver.h>
#include <DiagonalSOE.h>
#include <math.h>
#include <Channel.h>
#include <FEM_ObjectBroker.h>

void* OPS_DiagonalDirectSolver()
{
    DiagonalSolver *theSolver = new DiagonalDirectSolver();   
    return new DiagonalSOE(*theSolver);
}

DiagonalDirectSolver::DiagonalDirectSolver(double tol)
:DiagonalSolver(SOLVER_TAGS_DiagonalDirectSolver),
 minDiagTol(tol)
{

}

    
DiagonalDirectSolver::~DiagonalDirectSolver()
{

}

int
DiagonalDirectSolver::setSize(void)
{
  if (theSOE == 0) {
    opserr << "DiagonalDirectSolver::setSize()";
    opserr << " No system has been set!\n";
    return -1;
  }
  return 0;
}


int 
DiagonalDirectSolver::solve(void)
{

  // check for quick returns
  if (theSOE == 0) {
    opserr << "DiagonalDirectSolver::solve(void): ";
    opserr << " - No ProfileSPDSOE has been assigned\n";
    return -1;
  }
    
  if (theSOE->size == 0)
    return 0;
  
  // set some pointers
  double *Aptr = theSOE->A;
  double *Bptr = theSOE->B;
  double *Xptr = theSOE->X;
  int size = theSOE->size;

  if (theSOE->isAfactored == false)  {
    
    // FACTOR & SOLVE
    double invD;
    for (int i=0; i<size; i++) {
      
      double aii = *Aptr;

      // check that the diag > the tolerance specified
      if (aii == 0.0) {
	opserr << "DiagonalDirectSolver::solve() - ";
	opserr << " aii = 0 (i, aii): (" << i << ", " << aii << ")\n"; 
	return(-2);
      }
      if (fabs(aii) <= minDiagTol) {
	opserr << "DiagonalDirectSolver::solve() - ";
	opserr << " aii < minDiagTol (i, aii): (" << i;
	opserr << ", " << aii << ")\n"; 
	return(-2);
      }		

      // store the inverse 1/Aii in A; and solve for Xi
      invD = 1.0/aii; 
      *Xptr++ = invD * *Bptr++;
      *Aptr++ = invD;
    }

    theSOE->isAfactored = true;

  } else {

    // JUST SOLVE
    for (int i=0; i<size; i++) {
      *Xptr++ = *Aptr++ * *Bptr++;
    }
  }	
    
  return 0;
}

double
DiagonalDirectSolver::getDeterminant(void) 
{
  double determinant = 0.0;
  return determinant;
}

int 
DiagonalDirectSolver::setDiagonalSOE(DiagonalSOE &theNewSOE)
{
  if (theSOE != 0) {
    opserr << "DiagonalDirectSolver::setProfileSOE() - ";
    opserr << " has already been called \n";	
    return -1;
  }
  
  theSOE = &theNewSOE;
  return 0;
}


int
DiagonalDirectSolver::sendSelf(int cTag,
			       Channel &theChannel)
{
  static Vector data(1);
  data(0) = minDiagTol;
  return theChannel.sendVector(0, cTag, data);
}


int 
DiagonalDirectSolver::recvSelf(int cTag,
			       Channel &theChannel, 
			       FEM_ObjectBroker &theBroker)
{
  static Vector data(1);
  theChannel.recvVector(0, cTag, data);

  minDiagTol = data(0);
  return 0;
}


