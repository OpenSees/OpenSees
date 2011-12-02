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
// $Date: 2001-11-19 22:44:00 $
// $Source: /usr/local/cvs/OpenSees/SRC/system_of_eqn/eigenSOE/SymBandEigenSolver.h,v $

// Written: MHS
// Created: Oct 2001
//
// Description: This file contains the class definition for
// SymBandEigenSolver, which computes the eigenvalues and eigenvectors
// of a symmetric banded matrix using the LAPACK subroutine dsbevx.

#ifndef SymBandEigenSolver_h
#define SymBandEigenSolver_h

#include <EigenSolver.h>
#include <SymBandEigenSOE.h>

class SymBandEigenSolver : public EigenSolver
{
 public:
  SymBandEigenSolver();    
  virtual ~SymBandEigenSolver();
  
  virtual int solve(void) {return this->solve(theSOE->size);}
  virtual int solve(int nModes);
  virtual int setSize(void);
  virtual int setEigenSOE(SymBandEigenSOE &theSOE);
  
  virtual const Vector &getEigenvector(int mode);
  virtual double getEigenvalue(int mode);
  
  int sendSelf(int commitTag, Channel &theChannel);
  int recvSelf(int commitTag, Channel &theChannel, 
	       FEM_ObjectBroker &theBroker);
  
 protected:
  
 private:
  SymBandEigenSOE *theSOE;
  int numModes;

  double *eigenvalue;
  double *eigenvector;
  Vector *eigenV;
  double *work;
};

#endif
