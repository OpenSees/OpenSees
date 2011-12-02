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
                                                                        
// $Revision: 1.1 $
// $Date: 2006-03-15 00:26:02 $
// $Source: /usr/local/cvs/OpenSees/SRC/system_of_eqn/linearSOE/mumps/MumpsParallelSolver.h,v $
                                                                        
                                                                        
#ifndef MumpsParallelSolver_h
#define MumpsParallelSolver_h

// Written: fmk 
// Created: 02/06
//
// Description: This file contains the class definition for Mumps.
// A Mumps object can be constructed to solve a MumpsSOE. It obtains 
// the solution by making calls on the the Mumps library, which is based 
// on public domain software developed during the Esprit IV European project 
// PARASOL (1996-1999) by CERFACS, ENSEEIHT-IRIT and RAL. Since this first 
// public domain version in 1999, the developments are supported by CERFACS, 
// ENSEEIHT-IRIT, and INRIA Rhone-Alpes. 
// Main contributors to MUMPS are Patrick Amestoy, Iain Duff, Abdou Guermouche,
// Jacko Koster, Jean-Yves L'Excellent, and Stephane Pralet.

// Up-to-date copies of the MUMPS package can be obtained from the Web pages 
// http://www.enseeiht.fr/apo/MUMPS/ or http://graal.ens-lyon.fr/MUMPS

// What: "@(#) Mumps.h, revA"

#include <LinearSOESolver.h>
extern "C" {
#include <dmumps_c.h>
}

class MumpsParallelSOE;

class MumpsParallelSolver : public LinearSOESolver
{
  public:
  MumpsParallelSolver();
  MumpsParallelSolver(int MPI_COMM,
		      int ICNTL7 =7);

  virtual ~MumpsParallelSolver();
  
  int solve(void);
  int setSize(void);
  
  int sendSelf(int commitTag, Channel &theChannel);
  int recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker);    

  int setLinearSOE(MumpsParallelSOE &theSOE);

 protected:

 private:
  bool init;
  DMUMPS_STRUC_C id;
  MumpsParallelSOE *theMumpsSOE;

  int rank;
  int np;
};

#endif

