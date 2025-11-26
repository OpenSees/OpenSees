// $Revision: 1.3 $
// $Date: 2009-05-14 22:45:57 $
// $Source: /usr/local/cvs/OpenSees/SRC/system_of_eqn/eigenSOE/EigenSOE.h,v $

// Written: Jun Peng
// Created: Sat Feb. 6, 1999
// Revision: A
//
// Description: This file contains the class definition of EigenSOE.
// EigenSOE is a subclass of SystemOfEqn.
// It has pure virtual functions which must be implemented in it's derived
// subclasses. To solve the genreal eigen value equations means that
// by the given K and M, find the corresponding eigen value and eigen
// vectors.
//
// This class is inheritanted from the base class of SystemOfEqn
// which was created by fmk (Frank).


#ifndef EigenSOE_h
#define EigenSOE_h

#ifndef _bool_h
#include <bool.h>
#endif

#include <MovableObject.h>

class EigenSolver;
class AnalysisModel;
class Graph;
class Matrix;
class Vector;
class ID;
class LinearSOE;

class EigenSOE : public MovableObject
{
  public:
     EigenSOE(EigenSolver &theSolver, int classTag);
     EigenSOE(int classTag);
     virtual ~EigenSOE();
     
     virtual int solve(int numModes, bool generalized, bool findSmallest = true);
     virtual int setLinks(AnalysisModel &theModel);    
     virtual int setLinearSOE(LinearSOE &theSOE) {return -1;};
     
     // pure virtual functions
     virtual int addA(const Matrix &, const ID &, double fact = 1.0) = 0;
     virtual int addM(const Matrix &, const ID &, double fact = 1.0) = 0;

     virtual int setSize(Graph &theGraph) = 0;
     virtual void zeroA(void) = 0;
     virtual void zeroM(void) = 0;

     // methods to get the eigenvectors and eigenvalues
     virtual const Vector &getEigenvector(int mode);
     virtual double getEigenvalue(int mode);          
     
  protected:
     virtual int setSolver(EigenSolver &newSolver);
     EigenSolver *getSolver(void);
     EigenSolver *theSolver;
     
  private:     
     
};

#endif



