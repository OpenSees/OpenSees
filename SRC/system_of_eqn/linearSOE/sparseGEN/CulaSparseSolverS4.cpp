/* ****************************************************************** **
**    OpenSees - Open System for Earthquake Engineering Simulation    **
**          Pacific Earthquake Engineering Research Center            **
**                                                                    **
**                                                                    **
** (C) Copyright 1999, The Regents of the University of California    **
** All Rights Reserved.                                               **
**                                                                    **
** Commercial use of this program without express permission of the   **
** University of California, Berkeley, is strictly prohibited.  See   **
** file 'COPYRIGHT'  in main directory for information on usage and   **
** redistribution,  and for a DISCLAIMER OF ALL WARRANTIES.           **
**                                                                    **
** Developed by:                                                      **
**   Frank McKenna (fmckenna@ce.berkeley.edu)                         **
**   Gregory L. Fenves (fenves@ce.berkeley.edu)                       **
**   Filip C. Filippou (filippou@ce.berkeley.edu)                     **
**                                                                    **
** ****************************************************************** */
                                                                       
// Written: fmk 
// Created: 11/14

// this class is based on the S5 solver, written to work with Cula
// S4 library (S5 not available for apple at time of writing!) 

#include <CulaSparseSolverS4.h>
#include <SparseGenRowLinSOE.h>

static int numCulaS4Solvers = 0;

CulaSparseSolverS4::CulaSparseSolverS4(double tol,
				       int maxIterations,
				       int pCond,
				       int slver)
 :SparseGenRowLinSolver(SOLVER_TAGS_CulaSparseS4)
{

  if (numCulaS4Solvers == 0) {
    // initialize once
    numCulaS4Solvers++;
    status = culaSparseInitialize();
    if (status != culaNoError) {
	culaGetErrorInfoString( status, culaGetErrorInfo(), buf, sizeof(buf) );
	opserr << "CulaSparseSolverS4S4::CulaSparseSolverS4S4 : init Error: " << buf;
	solver = 10; // set to invalid so solver returns error
	return;
    }
  }
  
  tolerance = tol;
  maxNumIterations=maxIterations;
  preCond=pCond;
  solver=slver;

  status = culaIterativeConfigInit(&config);
  config.tolerance = tolerance;
  config.maxIterations = maxNumIterations;
  config.indexing = 0;
  //  config.debug = 1;

  switch (preCond)
    {
    case 0: 
      break;
    case 1: 
      status = culaJacobiOptionsInit(&jacobi);
      break;
    case 2: 
      status = culaBlockjacobiOptionsInit(&blockJacobi);
      break;
    case 3: 
      status = culaIlu0OptionsInit(&ilu0);   
      break;
    default:
      opserr << "WARNING invalid preconditioner - CulaSparseSolverS4S4 using blockJacobi preconditioner\n";
    }

  if (status != culaNoError) {
    culaGetErrorInfoString( status, culaGetErrorInfo(), buf, sizeof(buf) );
    opserr << "CulaSparseSolverS4S4::CulaSparseSolverS4S4 : init pre Error: " << buf;
    solver = 10; // set to invalid so solver returns error
    return;
  }

  
  switch (solver)
    {
    case 0: 
      status = culaCgOptionsInit(&cg);
      break;
    case 1: 
      status = culaBicgOptionsInit(&biCG);
      break;
    case 2: 
      status = culaBicgstabOptionsInit(&biCGstab);
      break;
    case 3: 
      status = culaBicgstablOptionsInit(&biCGstabl);
      break;
    case 4: 
      status = culaGmresOptionsInit(&gmRes);
      break;
    case 5: 
      status = culaMinresOptionsInit(&minRes);
      break;
    default:
      status=culaArgumentError;
      return;
    }
    
    if (status != culaNoError) {
      culaGetErrorInfoString( status, culaGetErrorInfo(), buf, sizeof(buf) );
      opserr << "CulaSparseSolverS4S4::CulaSparseSolverS4S4 : init Solver Error: " << buf;
      solver = 10; // set to invalid so solver returns error
      return;
    }
}


CulaSparseSolverS4::~CulaSparseSolverS4(void)
{

}


int
CulaSparseSolverS4::setSize()
{
  return 0;
}

int
CulaSparseSolverS4::setLinearSOE(SparseGenRowLinSOE &theLinearSOE)
{
  theSOE = &theLinearSOE;
  return 0;
}


int
CulaSparseSolverS4::sendSelf(int cTAg, Channel &theChannel)
{
  return 0;
}


int
CulaSparseSolverS4::recvSelf(int cTag,
			     Channel &theChannel, FEM_ObjectBroker &theBroker)
{
  return 0;
}

int CulaSparseSolverS4::solve(void)
{
  if (theSOE == 0) {
    opserr << "WARNING CulaSparse::solve(void)- ";
    opserr << " No LinearSOE object has been set!!!\n";
    return -1;
  }
  
  if (status != culaNoError) {
    return -1;
  }
  
  n = theSOE->size;
  nnz = theSOE->nnz;
  
  double *Xptr = theSOE->X;
  double *Bptr = theSOE->B;
  double *Aptr = theSOE->A;
  
  int *rowPtr=theSOE->rowStartA;
  int *colInd=theSOE->colA;
  
  switch (solver) 
    {
    case 0: 
      switch (preCond) 
	{ 
	case 0: 
	  status = culaDcsrCg(&config, &cg, NULL, n, nnz, Aptr, colInd, rowPtr, Xptr, Bptr, &result);
	  break;
        case 1: 
	  status = culaDcsrCgJacobi(&config, &cg, &jacobi, n, nnz, Aptr, colInd, rowPtr, Xptr, Bptr, &result);
	  break;
	case 2: 
	  status = culaDcsrCgBlockjacobi(&config, &cg, &blockJacobi, n, nnz, Aptr, colInd, rowPtr, Xptr, Bptr, &result);
	  break;
	case 3: 
	  status = culaDcsrCgIlu0(&config, &cg, &ilu0, n, nnz, Aptr, colInd, rowPtr, Xptr, Bptr, &result);
	  break;
	default:
	  status=culaArgumentError;
	}
      break;
    case 1: 
      switch (preCond) 
	{ 
	case 0: 
	  status = culaDcsrBicg(&config, &biCG, NULL, n, nnz, Aptr, colInd, rowPtr, Xptr, Bptr, &result);
	  break;
        case 1: 
	  status = culaDcsrBicgJacobi(&config, &biCG, &jacobi, n, nnz, Aptr, colInd, rowPtr, Xptr, Bptr, &result);
	  break;
	case 2: 
	  status = culaDcsrBicgBlockjacobi(&config, &biCG, &blockJacobi, n, nnz, Aptr, colInd, rowPtr, Xptr, Bptr, &result);
	  break;
	case 3: 
	  status = culaDcsrBicgIlu0(&config, &biCG, &ilu0, n, nnz, Aptr, colInd, rowPtr, Xptr, Bptr, &result);
	  break;
	default:
	  status=culaArgumentError;
	}
      break;
    case 2: 
      status = culaBicgstabOptionsInit(&biCGstab);
      switch (preCond) 
	{
 	case 0: 
	  status = culaDcsrBicgstab(&config, &biCGstab, NULL, n, nnz, Aptr, colInd, rowPtr, Xptr, Bptr, &result);
	  break;
        case 1: 
	  status = culaDcsrBicgstabJacobi(&config, &biCGstab, &jacobi, n, nnz, Aptr, colInd, rowPtr, Xptr, Bptr, &result);
	  break;
	case 2: 
	  status = culaDcsrBicgstabBlockjacobi(&config, &biCGstab, &blockJacobi, n, nnz, Aptr, colInd, rowPtr, Xptr, Bptr, &result);
	  break;
	case 3: 
	  status = culaDcsrBicgstabIlu0(&config, &biCGstab, &ilu0, n, nnz, Aptr, colInd, rowPtr, Xptr, Bptr, &result);
	  break;
	default:
	  status=culaArgumentError;
	}
      break;

    case 3: 
      switch (preCond) 
	{
 	case 0: 
	  status = culaDcsrBicgstabl(&config, &biCGstabl, NULL, n, nnz, Aptr, colInd, rowPtr, Xptr, Bptr, &result);
	  break;
        case 1: 
	  status = culaDcsrBicgstablJacobi(&config, &biCGstabl, &jacobi, n, nnz, Aptr, colInd, rowPtr, Xptr, Bptr, &result);
	  break;
	case 2: 
	  status = culaDcsrBicgstablBlockjacobi(&config, &biCGstabl, &blockJacobi, n, nnz, Aptr, colInd, rowPtr, Xptr, Bptr, &result);
	  break;
	case 3: 
	  status = culaDcsrBicgstablIlu0(&config, &biCGstabl, &ilu0, n, nnz, Aptr, colInd, rowPtr, Xptr, Bptr, &result);
	  break;
	default:
	  status=culaArgumentError;
	}
      break;

    case 4: 
      switch (preCond) 
	{
 	case 0: 
	  status = culaDcsrGmres(&config, &gmRes, NULL, n, nnz, Aptr, colInd, rowPtr, Xptr, Bptr, &result);
	  break;
        case 1: 
	  status = culaDcsrGmresJacobi(&config, &gmRes, &jacobi, n, nnz, Aptr, colInd, rowPtr, Xptr, Bptr, &result);
	  break;
	case 2: 
	  status = culaDcsrGmresBlockjacobi(&config, &gmRes, &blockJacobi, n, nnz, Aptr, colInd, rowPtr, Xptr, Bptr, &result);
	  break;
	case 3: 
	  status = culaDcsrGmresIlu0(&config, &gmRes, &ilu0, n, nnz, Aptr, colInd, rowPtr, Xptr, Bptr, &result);
	  break;
	default:
	  status=culaArgumentError;
	}
      break;
    case 5: 
      switch (preCond) 
	{
 	case 0: 
	  status = culaDcsrMinres(&config, &minRes, NULL, n, nnz, Aptr, colInd, rowPtr, Xptr, Bptr, &result);
	  break;
        case 1: 
	  status = culaDcsrMinresJacobi(&config, &minRes, &jacobi, n, nnz, Aptr, colInd, rowPtr, Xptr, Bptr, &result);
	  break;
	case 2: 
	  status = culaDcsrMinresBlockjacobi(&config, &minRes, &blockJacobi, n, nnz, Aptr, colInd, rowPtr, Xptr, Bptr, &result);
	  break;
	case 3: 
	  status = culaDcsrMinresIlu0(&config, &minRes, &ilu0, n, nnz, Aptr, colInd, rowPtr, Xptr, Bptr, &result);
	  break;
	default:
	  status=culaArgumentError;
	}
      break;
    default:
      status=culaArgumentError;
    }
  
  // see if solver failed 
  if (status != culaNoError) {
    culaGetErrorInfoString( status, culaGetErrorInfo(), buf, sizeof(buf) );
    opserr << "CulaSparseSolverS4S4::CulaSparseSolverS4S4 : init Error" << buf;
    return -1;
  }

  return 0;
}


