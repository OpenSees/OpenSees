#include "CuSPSolver.h"
#include <Windows.h>

CuSPSolver::CuSPSolver(void):SparseGenRowLinSolver(SOLVER_TAGS_CuSP)
{
  single=0;
  error=0;
  
  HINSTANCE hDLL=LoadLibrary("CuSPSolver.dll");
  if (hDLL)
    {
      SolveFunc=(CUSPSOLVE)GetProcAddress(hDLL,"CuSPSolve");
      
      if (!SolveFunc)
	{
	  error=1;
	  return;
	}
    }else
    {
      error=2;
      return;
    }
  
  maxInteration=100000;
  relTolerance=1e-6;
  
  preCond=0;
  solver=0;
}



CuSPSolver::CuSPSolver(int maxInt,double relTol,int pre,int solv):SparseGenRowLinSolver(SOLVER_TAGS_CuSP)
{
  single=0;
  error=0;
  
  HINSTANCE hDLL=LoadLibrary("CuSPSolver.dll");
  if (hDLL)
    {
      SolveFunc=(CUSPSOLVE)GetProcAddress(hDLL,"CuSPSolve");
      
      if (!SolveFunc)
	{
	  error=1;
	  return;
	}
    }else
    {
      error=2;
      return;
    }
  
  maxInteration=maxInt;
  relTolerance=relTol;
  preCond=pre;
  solver=solv;
}


CuSPSolver::~CuSPSolver(void)
{
  if (single==1)
    {
      delete Xsingle;
      delete Bsingle;
      delete Asingle;
    }
}


int
CuSPSolver::setSize()
{
  int n=theSOE->size;
  int nnz=theSOE->nnz;
  
  if (single==1)
    {
      Xsingle=new float[n];
      Bsingle=new float[n];
      Asingle=new float[nnz];
    }
  
  return 0;
}

int
CuSPSolver::setLinearSOE(SparseGenRowLinSOE &theLinearSOE)
{
  theSOE = &theLinearSOE;
  return 0;
}


int
CuSPSolver::sendSelf(int cTAg, Channel &theChannel)
{
  // doing nothing
  return 0;
}


int
CuSPSolver::recvSelf(int cTag,
		     Channel &theChannel, FEM_ObjectBroker &theBroker)
{
  // nothing to do
  return 0;
}

int CuSPSolver::solve(void)
{
  if (error!=0)
    {
      switch (error)
	{
	case 1:opserr<<"DLL Load Error!\n";break;
	case 2:opserr<<"Function Load Error!The version of DLL maybe incorrect.\n";break;
	}
      return -1;
    }
  
  
  
  n = theSOE->size;
  nnz = theSOE->nnz;
  
  double *Xptr = theSOE->X;
  double *Bptr = theSOE->B;
  double *Aptr = theSOE->A;
  
  
  int *rowPtr=theSOE->rowStartA;
  int *colInd=theSOE->colA;
  
  
  int errorcode=SolveFunc(Aptr,Bptr,Xptr,n,nnz,rowPtr,colInd,maxInteration,relTolerance,preCond,solver);
  if (errorcode!=0)
    {
      switch (errorcode)
	{
	case -1:opserr<<"Wrong Preconditioner! Please check the TCL file.\n";break;
	case -2:opserr<<"Wrong Solver! Please check the TCL file.\n"; break;
	}
      return -1;
    }
  
  return 0;
}




