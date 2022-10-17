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

// $Revision: 1.13 $
// $Date: 2007-04-13 22:27:55 $
// $Source: /usr/local/cvs/OpenSees/SRC/analysis/algorithm/equiSolnAlgo/KrylovNewton.cpp,v $

// Written: MHS
// Created: June 2001
//
// Description: This file contains the class definition for 
// KrylovNewton.  KrylovNewton is a class which uses a Krylov
// subspace accelerator on the modified Newton method.
// The accelerator is described by Carlson and Miller in
// "Design and Application of a 1D GWMFE Code"
// from SIAM Journal of Scientific Computing (Vol. 19, No. 3,
// pp. 728-765, May 1998)

#include <KrylovNewton.h>
#include <AnalysisModel.h>
#include <StaticAnalysis.h>
#include <IncrementalIntegrator.h>
#include <LinearSOE.h>
#include <Channel.h>
#include <FEM_ObjectBroker.h>
#include <ConvergenceTest.h>
#include <Matrix.h>
#include <Vector.h>
#include <ID.h>

// Constructor
KrylovNewton::KrylovNewton(int theTangentToUse, int maxDim)
:EquiSolnAlgo(EquiALGORITHM_TAGS_KrylovNewton),
 tangent(theTangentToUse),
 v(0), Av(0), AvData(0), rData(0), work(0), lwork(0),
 numEqns(0), maxDimension(maxDim)
{
  if (maxDimension < 0)
    maxDimension = 0;
}

KrylovNewton::KrylovNewton(ConvergenceTest &theT, int theTangentToUse, int maxDim)
:EquiSolnAlgo(EquiALGORITHM_TAGS_KrylovNewton),
 tangent(theTangentToUse),
 v(0), Av(0), AvData(0), rData(0), work(0), lwork(0),
 numEqns(0), maxDimension(maxDim)
{
  if (maxDimension < 0)
    maxDimension = 0;
}

// Destructor
KrylovNewton::~KrylovNewton()
{
  if (v != 0) {
    for (int i = 0; i < maxDimension+1; i++)
      delete v[i];
    delete [] v;
  }

  if (Av != 0) {
    for (int i = 0; i < maxDimension+1; i++)
      delete Av[i];
    delete [] Av;
  }

  if (AvData != 0)
    delete [] AvData;

  if (rData != 0)
    delete [] rData;

  if (work != 0)
    delete [] work;
}

int 
KrylovNewton::solveCurrentStep(void)
{
  // set up some pointers and check they are valid
  // NOTE this could be taken away if we set Ptrs as protecetd in superclass
  AnalysisModel *theAnaModel = this->getAnalysisModelPtr();
  IncrementalIntegrator *theIntegrator = this->getIncrementalIntegratorPtr();
  LinearSOE *theSOE = this->getLinearSOEptr();
  
  if ((theAnaModel == 0) || (theIntegrator == 0) || (theSOE == 0)
      || (theTest == 0)){
    opserr << "WARNING KrylovNewton::solveCurrentStep() - setLinks() has";
    opserr << " not been called - or no ConvergenceTest has been set\n";
    return -5;
  }	

  // Get size information from SOE
  numEqns  = theSOE->getNumEqn();
  if (maxDimension > numEqns)
    maxDimension = numEqns;

  if (v == 0) {
    // Need to allocate an extra vector for "next" update
    v = new Vector*[maxDimension+1];
    for (int i = 0; i < maxDimension+1; i++)
      v[i] = new Vector(numEqns);
  }

  if (Av == 0) {
    Av = new Vector*[maxDimension+1];
    for (int i = 0; i < maxDimension+1; i++)
      Av[i] = new Vector(numEqns);
  }

  if (AvData == 0)
    AvData = new double [maxDimension*numEqns];

  if (rData == 0)
    // The LAPACK least squares subroutine overwrites the RHS vector
    // with the solution vector ... these vectors are not the same
    // size, so we need to use the max size
    rData = new double [(numEqns > maxDimension) ? numEqns : maxDimension];

  // Length of work vector should be >= 2*min(numEqns,maxDimension)
  // See dgels subroutine documentation
  lwork = 2 * ((numEqns < maxDimension) ? numEqns : maxDimension);
  
  if (work == 0)
    work = new double [lwork];

  // Evaluate system residual R(y_0)
  if (theIntegrator->formUnbalance() < 0) {
    opserr << "WARNING KrylovNewton::solveCurrentStep() -";
    opserr << "the Integrator failed in formUnbalance()\n";	
    return -2;
  }


  // set itself as the ConvergenceTest objects EquiSolnAlgo
  theTest->setEquiSolnAlgo(*this);
  if (theTest->start() < 0) {
    opserr << "KrylovNewton::solveCurrentStep() -";
    opserr << "the ConvergenceTest object failed in start()\n";
    return -3;
  }
  
  
  // Evaluate system Jacobian J = R'(y)|y_0
  if (theIntegrator->formTangent(tangent) < 0){
    opserr << "WARNING KrylovNewton::solveCurrentStep() -";
    opserr << "the Integrator failed in formTangent()\n";
    return -1;
  }    

  // Loop counter
  int k = 1;

  // Current dimension of Krylov subspace
  int dim = 0;

  int result = -1;

  do {

    // Clear the subspace if its dimension has exceeded max
    if (dim > maxDimension) {
      dim = 0;
      if (theIntegrator->formTangent(tangent) < 0){
	opserr << "WARNING KrylovNewton::solveCurrentStep() -";
	opserr << "the Integrator failed to produce new formTangent()\n";
	return -1;
      }
    }

    // Solve for residual f(y_k) = J^{-1} R(y_k)
    if (theSOE->solve() < 0) {
      opserr << "WARNING KrylovNewton::solveCurrentStep() -";
      opserr << "the LinearSysOfEqn failed in solve()\n";	
      return -3;
    }

    // Solve least squares A w_{k+1} = r_k
    if (this->leastSquares(dim) < 0) {
      opserr << "WARNING KrylovNewton::solveCurrentStep() -";
      opserr << "the Integrator failed in leastSquares()\n";
      return -1;
    }		    

    // Update system with v_k
    if (theIntegrator->update(*(v[dim])) < 0) {
      opserr << "WARNING KrylovNewton::solveCurrentStep() -";
      opserr << "the Integrator failed in update()\n";	
      return -4;
    }	

    // Evaluate system residual R(y_k)
    if (theIntegrator->formUnbalance() < 0) {
      opserr << "WARNING KrylovNewton::solveCurrentStep() -";
      opserr << "the Integrator failed in formUnbalance()\n";	
      return -2;
    }

    // Increase current dimension of Krylov subspace
    dim++;

    result = theTest->test();
    this->record(k++);

  } while (result == -1);
  
  if (result == -2) {
    opserr << "KrylovNewton::solveCurrentStep() -";
    opserr << "the ConvergenceTest object failed in test()\n";
    return -3;
  }
  
  // note - if positive result we are returning what the convergence
  // test returned which should be the number of iterations
  return result;
}

int
KrylovNewton::sendSelf(int cTag, Channel &theChannel)
{
  static ID data(2);
  data(0) = tangent;
  data(1) = maxDimension;
  if (theChannel.sendID(cTag, 0, data) < 0) {
    opserr << "KrylovNewton::sendSelf() - failed\n";
    return -1;
  }
  return 0;
}

int
KrylovNewton::recvSelf(int cTag, Channel &theChannel, 
					   FEM_ObjectBroker &theBroker)
{
  static ID data(2);
  if (theChannel.recvID(cTag, 0, data) <  0) {
    opserr << "KrylovNewton::recvSelf() - failed\n";
    return -1;
  }
  tangent = data(0);
  maxDimension = data(1);
  return 0;
}

void
KrylovNewton::Print(OPS_Stream &s, int flag)
{
  s << "KrylovNewton";
  s << "\n\tMax subspace dimension: " << maxDimension;
  s << "\n\tNumber of equations: " << numEqns << endln;
}

#ifdef _WIN32

extern "C" int  DGELS(char *T, int *M, int *N, int *NRHS,
		      double *A, int *LDA, double *B, int *LDB,
		      double *WORK, int *LWORK, int *INFO);

#else

extern "C" int dgels_(char *T, int *M, int *N, int *NRHS,
		      double *A, int *LDA, double *B, int *LDB,
		      double *WORK, int *LWORK, int *INFO);

#endif

int
KrylovNewton::leastSquares(int k)
{
  LinearSOE *theSOE = this->getLinearSOEptr();  
  const Vector &r = theSOE->getX();

  // v_{k+1} = w_{k+1} + q_{k+1}
  *(v[k])  = r;
  *(Av[k]) = r;

  // Subspace is empty
  if (k == 0)
    return 0;

  // Compute Av_k = f(y_{k-1}) - f(y_k) = r_{k-1} - r_k
  Av[k-1]->addVector(1.0, r, -1.0);

  int i,j;

  // Put subspace vectors into AvData
  Matrix A(AvData, numEqns, k);
  for (i = 0; i < k; i++) {
    Vector &Ai = *(Av[i]);
    for (j = 0; j < numEqns; j++)
      A(j,i) = Ai(j);
  }

  // Put residual vector into rData (need to save r for later!)
  Vector B(rData, numEqns);
  B = r;
  
  // No transpose
  //  char *trans = "N";
  char trans[] = "N";

  // The number of right hand side vectors
  int nrhs = 1;

  // Leading dimension of the right hand side vector
  int ldb = (numEqns > k) ? numEqns : k;

  // Subroutine error flag
  int info = 0;

  // Call the LAPACK least squares subroutine
#ifdef _WIN32
  DGELS(trans, &numEqns, &k, &nrhs, AvData, &numEqns, rData, &ldb, work, &lwork, &info);
#else
  dgels_(trans, &numEqns, &k, &nrhs, AvData, &numEqns, rData, &ldb, work, &lwork, &info);
#endif

  // Check for error returned by subroutine
  if (info < 0) {
    opserr << "WARNING KrylovNewton::leastSquares() - \n";
    opserr << "error code " << info << " returned by LAPACK dgels\n";
    return info;
  }
  
  // Compute the correction vector
  double cj;
  for (j = 0; j < k; j++) {
    
    // Solution to least squares is written to rData
    cj = rData[j];
    
    // Compute w_{k+1} = c_1 v_1 + ... + c_k v_k
    v[k]->addVector(1.0, *(v[j]), cj);
    
    // Compute least squares residual q_{k+1} = r_k - (c_1 Av_1 + ... + c_k Av_k)
    v[k]->addVector(1.0, *(Av[j]), -cj);
  }
  
  return 0;
}

