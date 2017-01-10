#pragma once
#include <SparseGenRowLinSolver.h>
#include <SparseGenRowLinSOE.h>
#include <SparseGenColLinSOE.h>
#include <SparseGenColLinSolver.h>

typedef int (*CUSPSOLVE)(double* Aptr, double* Bptr, double* Xptr,int n,int nnz,int* rowPtr,int*  colInd,int maxInt,double relTol,int pre,int solv);

class CuSPSolver~ :
	public SparseGenRowLinSolver
{
public:
	CuSPSolver(void);
	CuSPSolver(int maxInt,double relTol,int pre,int solv);
	~CuSPSolver(void);
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




	double relTolerance;
	int maxInteration;
	int preCond;			
	int solver;				
	int single;
	int host;


	float *Xsingle;
	float *Bsingle;
	float *Asingle;

	int error;

	//CUSPSOLVE SolveFunc;
	



	


};

