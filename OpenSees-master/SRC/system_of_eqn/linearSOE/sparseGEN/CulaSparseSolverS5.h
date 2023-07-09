#ifndef _CULASPARSESOLVERS5
#define _CULASPARSESOLVERS5

#include <cula_sparse.h>
#include <SparseGenRowLinSolver.h>
#include <SparseGenRowLinSOE.h>
#include <SparseGenColLinSOE.h>
#include <SparseGenColLinSolver.h>



class CulaSparseSolverS5 : public SparseGenRowLinSolver
{
 public:
   CulaSparseSolverS5(void);
   CulaSparseSolverS5(double relTol,int maxInteration,int preCond,int solver,int single,int host);
   ~CulaSparseSolverS5(void);
   
   int solve(void);
   int setSize(void);
   
   int sendSelf(int commitTag, Channel &theChannel);
   int recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker);    
   
   int setLinearSOE(SparseGenRowLinSOE &theLinearSOE);
   //int setLinearSOE(SparseGenColLinSOE &theLinearSOE);
   
   
 private:
   int n;// order of matrix
   int nnz;// Number of no-zero ones in A;
   
   //double* A,X,B;//Matrix A, X, and B
   //int* rowA,colStartA;
   
   double relTol;
   int maxInteration;
   int preCond;			//fainv
   int solver;				//Bicg
   int single;
   int host;
   
   
   float *Xsingle;
   float *Bsingle;
   float *Asingle;
      
   SparseGenRowLinSOE *theSOE;
   //SparseGenColLinSOE *theSOE;
   
   // status returned by each and every cula routine
   culaSparseStatus status;
   
   // create library handle
   culaSparseHandle handle;
   culaSparsePlan plan;
   
   // configuration structures
   culaSparseConfig config;
   culaSparseCsrOptions CsrOpt;
   
   // character buffer used for results and error messages
   char buf[512];
   
   // result structure
   culaSparseResult result;
};

#endif
