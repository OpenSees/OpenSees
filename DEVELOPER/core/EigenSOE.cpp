// $Revision: 1.5 $
// $Date: 2009-05-14 22:45:57 $
// $Source: /usr/local/cvs/OpenSees/SRC/system_of_eqn/eigenSOE/EigenSOE.cpp,v $

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


#include <EigenSOE.h>
#include <EigenSolver.h>
#include <AnalysisModel.h>

EigenSOE::EigenSOE(EigenSolver &theEigenSolver, int classTag)
  :MovableObject(classTag), theSolver(&theEigenSolver)
{

}


EigenSOE::EigenSOE(int classTag)
  :MovableObject(classTag), theSolver(0)
{

}


EigenSOE::~EigenSOE()
{
  if (theSolver != 0)
    delete theSolver;
}

int 
EigenSOE::solve(int numModes, bool generalized, bool findSmallest)
{
  if (theSolver == 0)
    return -1;
  else
    return (theSolver->solve(numModes, generalized, findSmallest));
}

int
EigenSOE::setSolver(EigenSolver &newSolver)
{
    theSolver = &newSolver;
    return 0;
}

EigenSolver *
EigenSOE::getSolver(void)
{	
    return theSolver;
}

const Vector &
EigenSOE::getEigenvector(int mode) {
    return theSolver->getEigenvector(mode);
}

double 
EigenSOE::getEigenvalue(int mode) {
    return theSolver->getEigenvalue(mode);
}


int 
EigenSOE::setLinks(AnalysisModel &theModel)
{
  return 0;
}

