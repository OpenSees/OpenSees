// File: ~/system_of_eqn/linearSOE/LawSolver/SymSparseLinSOE.h
//
// Written: Jun Peng  (junpeng@stanford.edu)
//          Prof. Kincho H. Law
//          Stanford University
// Created: 12/98
// Revision: A
//
// Description: This file contains the class definition for 
// SymSparseLinSOE.h. It stores the sparse matrix A in a fashion
// that only store the none zero entries.
//
// What: "@(#) SymSparseLinSOE.h, revA"
//
// Almost all the information (Matrix A and Vector B) is stored as 
// global variables in the file "symbolic.h".


#ifndef SymSparseLinSOE_h
#define SymSparseLinSOE_h

#include <LinearSOE.h>
#include <Vector.h>

extern "C" {
   #include <FeStructs.h>
}

class SymSparseLinSolver;

class SymSparseLinSOE : public LinearSOE
{
  public:
    SymSparseLinSOE(SymSparseLinSolver &theSolver, int lSparse);        
    SymSparseLinSOE(int N, int NNZ, int *rowStartA, int *colA,
		    SymSparseLinSolver &theSolver, int lSparse);        

    ~SymSparseLinSOE();

    int getNumEqn(void) const;
    int setSize(Graph &theGraph);
    int addA(const Matrix &, const ID &, double fact = 1.0);
    int addB(const Vector &, const ID &, double fact = 1.0);    
    int setB(const Vector &, double fact = 1.0);        
    
    void zeroA(void);
    void zeroB(void);
    
    const Vector &getX(void);
    const Vector &getB(void);    
    double normRHS(void);

    void setX(int loc, double value);        
    void setX(const Vector &x);        
    int setSymSparseLinSolver(SymSparseLinSolver &newSolver);    

    int sendSelf(int commitTag, Channel &theChannel);
    int recvSelf(int commitTag, Channel &theChannel, 
		 FEM_ObjectBroker &theBroker);

    friend class SymSparseLinSolver;

  protected:
    
  private:
    int size;            // order of A
    int nnz;             // number of non-zeros in A
    double *B, *X;       // 1d arrays containing coefficients of B and X
    int *colA, *rowStartA;  //These are (ADJNCY, XADJ) pair.

    Vector *vectX;
    Vector *vectB;
    int Bsize;
    bool factored;

    int      LSPARSE;
    int      nblks;
    int      *xblk,  *invp;
    double   *diag, **penv;
    int      *rowblks;
    OFFDBLK  **begblk;
    OFFDBLK  *first;

};

#endif

