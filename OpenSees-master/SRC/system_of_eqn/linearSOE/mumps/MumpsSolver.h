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
                                                                        
// $Revision: 1.6 $
// $Date: 2008-04-15 07:15:29 $
// $Source: /usr/local/cvs/OpenSees/SRC/system_of_eqn/linearSOE/mumps/MumpsSolver.h,v $
                                                                        
                                                                        
#ifndef MumpsSolver_h
#define MumpsSolver_h

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

class MumpsSOE;

class MumpsSolver : public LinearSOESolver
{
  public:
  MumpsSolver(int ICNTL7=7, int ICNTL14=20);
	      
  virtual ~MumpsSolver();
  
  int solve(void);
  int setSize(void);
  
  int sendSelf(int commitTag, Channel &theChannel);
  int recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker);    

  int setLinearSOE(MumpsSOE &theSOE);

 protected:

 private:

  int initializeMumps(void);
  int solveAfterInitialization(void);

  DMUMPS_STRUC_C id;
  MumpsSOE *theMumpsSOE;
  bool init;
  int icntl14;
  int icntl7;
  bool needsSetSize;
};

#endif

