// File: ~/system_of_eqn/linearSOE/LawSolver/SymArpackSolver.h
//
// Written: Jun Peng
// Created: 2/1999
// Revision: A
//
// Description: This file contains the class definition for 
// SymArpackSolver. It solves the SymArpackSOE object by calling
// some "C" functions. The solver used here is generalized sparse
// solver. The user can choose three different ordering schema.
//
// What: "@(#) SymArpackSolver.h, revA"

#ifndef SymArpackSolver_h
#define SymArpackSolver_h

#include <EigenSolver.h>
#include <SymArpackSOE.h>

class SymArpackSOE;

class SymArpackSolver : public EigenSolver
{
  public:
    SymArpackSolver(int numE = 0);     
    virtual ~SymArpackSolver();

    virtual int solve(int numModes, bool generalized, bool findSmallest = true);    
    virtual int setSize(void);

    virtual int setEigenSOE(SymArpackSOE &theSOE); 
	
    virtual const Vector &getEigenvector(int mode);
    virtual double getEigenvalue(int mode);
    
    int sendSelf(int cTag, Channel &theChannel);
    int recvSelf(int cTag, Channel &theChannel, 
		 FEM_ObjectBroker &theBroker);
    
  protected:

  private:
    SymArpackSOE *theSOE;
    bool factored;

    int theNev;
    double *value;
    double *vector;
    Vector *eigenV;
    
    void myMv(int n, double *v, double *result);
    void myCopy(int n, double *v, double *result);
    int getNCV(int n, int nev);
    
};

#endif

