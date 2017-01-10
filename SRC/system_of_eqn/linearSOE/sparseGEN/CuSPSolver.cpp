// Written: fmk 
// modified by gourik@cdac.in (CDAC, Pune)
// nayakshweta19@gmail.com (ME student PICT,Pune)
//

#include "CuSPSolver.h"
#include <iostream>
using std::cout;

int cuspsolver(double *C, double *B, double *host_x,int size,int nnz,int *rowStartA, int *colA,int maxInteration,double relTolerance,int preCond,int solver);

CuSPSolver::CuSPSolver(void):SparseGenRowLinSolver(SOLVER_TAGS_CuSP)
{
  single=0;
  error=0;
 #ifdef _WIN32
  {
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
  }
 #endif
  maxInteration=100000;
  relTolerance=1e-6;
  
  preCond=0;
  solver=0;
}



CuSPSolver::CuSPSolver(int maxInt,double relTol,int pre,int solv):SparseGenRowLinSolver(SOLVER_TAGS_CuSP)
{
  single=0;
  error=0;
  #ifdef _WIN32 
	{  
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
	 }
   #endif
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
/*****************************************************************************/
int
CuSPSolver::setLinearSOE(SparseGenRowLinSOE &theLinearSOE)
{
  //std::cout<<"In set SOE";
  theSOE = &theLinearSOE;
  return 0;
}

/********************************************************************************/
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
 #ifdef _WIN32
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
	  
	 }  
  
  #endif
  n = theSOE->size;
  nnz = theSOE->nnz;
  
  double *Xptr = theSOE->X;
  double *Bptr = theSOE->B;
  double *Aptr = theSOE->A;
  
  
  int *rowPtr=theSOE->rowStartA;
  int *colInd=theSOE->colA;
/*std::cout<<"\nIN CUSPSOLVERL.CPP\t";/*<<maxInteration<<"\n";
std::cout<<"size\t"<<n<<"\n";
std::cout<<"non zero\t"<<nnz<<"\n";
std::cout<<"A\n";
for(int i=0;i<nnz;i++)
std::cout<<i<<"\t"<<Aptr[i]<<"\n";
std::cout<<"\n";
std::cout<<"rowPtr\n";
for(int i=0;i<n;i++)
std::cout<<rowPtr[i]<<"\t";
std::cout<<"\n";
std::cout<<"colInd\n";
for(int i=0;i<n;i++)
std::cout<<colInd[i]<<"\t";
std::cout<<"\n";*/
	#ifdef _WIN32
	{
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
	}
       #endif
theSOE->getA();
  int errorcode=cuspsolver(Aptr,Bptr,Xptr,n,nnz,rowPtr,colInd,maxInteration,relTolerance,preCond,solver);
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




