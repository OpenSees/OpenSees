// File: ~/system_of_eqn/linearSOE/LawSolver/SymArpackSOE.h
//
// Written: Jun Peng
// Created: 12/98
// Revision: A
//
// Description: This file contains the class definition for 
// SymArpackSOE.h. It stores the sparse matrix A in a fashion
// that only store the none zero.
//
// What: "@(#) SymArpackSOE.h, revA"
//
// Almost all the information (Matrix A and Vector B) is stored as 
// global variables in the file "symbolic.h".


#ifndef SymArpackSOE_h
#define SymArpackSOE_h

#include <EigenSOE.h>
#include <Vector.h>

extern "C" {
   #include <FeStructs.h>
}



class SymArpackSolver;
class AnalysisModel;

class SymArpackSOE : public EigenSOE
{
  public:
    SymArpackSOE(SymArpackSolver &theSolver, 
		 AnalysisModel &theModel,
		 double shift = 0.0);        

    virtual ~SymArpackSOE();

    virtual int getNumEqn(void) const;
    virtual int setSize(Graph &theGraph);
    
    virtual int addA(const Matrix &, const ID &, double fact = 1.0);
    virtual int addM(const Matrix &, const ID &, double fact = 1.0);    
      
    virtual void zeroA(void);
    virtual void zeroM(void);
    
    virtual double getShift(void);
    
    int sendSelf(int commitTag, Channel &theChannel);
    int recvSelf(int commitTag, Channel &theChannel, 
		 FEM_ObjectBroker &theBroker);

    friend class SymArpackSolver;


  protected:
    
  private:
    int size;               // order of A
    int nnz;                // number of non-zeros in A
    int *colA, *rowStartA;  //These are (ADJNCY, XADJ) pair.

    bool factored;
    double shift;
    AnalysisModel *theModel;

    int      nblks;
    int      *xblk,  *invp;
    double   *diag, **penv;
    int      *rowblks;
    OFFDBLK  **begblk;
    OFFDBLK  *first;
};

#endif

