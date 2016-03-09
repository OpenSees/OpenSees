// $Revision: 1.2 $
// $Date: 2009-05-11 21:01:10 $
// $Source: /usr/local/cvs/OpenSees/SRC/system_of_eqn/eigenSOE/EigenSolver.h,v $

// Written: Jun Peng
// Created: Sat Feb. 6, 1999
// Revision: A
//
// Description: This file contains the class definition of EigenSOE.
// EigenSOE is a subclass of Solver.
// This is an abstract base class and thus no objects of it's type
// can be instantiated. Instances of EigenSolver are used to solve
// a EigenSOE. (perform eigen analysis)
//
// This class is inheritanted from the base class of Solver
// which was created by fmk (Frank).


#ifndef EigenSolver_h
#define EigenSolver_h

#include <MovableObject.h>
#include <Vector.h>

class EigenSOE;

class EigenSolver : public MovableObject
{
  public:
     EigenSolver(int classTag);
     virtual ~EigenSolver();

     virtual int solve(int numModes, bool generalized, bool findSmallest = true) =0;     
     virtual const Vector &getEigenvector(int mode) = 0;
     virtual double getEigenvalue(int mode) = 0;     

     virtual int setSize() = 0;
     
  private:

};

#endif
