/* ****************************************************************** **
**    OpenSees - Open System for Earthquake Engineering Simulation    **
**          Pacific Earthquake Engineering Research Center            **
**                                                                    **
**                                                                    **
** (C) Copyright 2001, The Regents of the University of California    **
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
                                                                        
// $Revision$
// $Date$
// $Source$
                                                                        
// Written: MHS
// Created: April 2002

// Description: This file contains the class implementation for 
// KrylovAccelerator2. 

#include <KrylovAccelerator2.h>

#include <Matrix.h>
#include <Vector.h>
#include <LinearSOE.h>
#include <IncrementalIntegrator.h>

#include <ID.h>
#include <Channel.h>
#include <math.h>

KrylovAccelerator2::KrylovAccelerator2(int max, int tangent)
  :Accelerator(ACCELERATOR_TAGS_Krylov),
   dimension(0), numEqns(0), maxDimension(max),
   v(0), Av(0), AvData(0), rData(0), work(0), lwork(0), theTangent(tangent)
{
  if (maxDimension < 0)
    maxDimension = 0;
}

KrylovAccelerator2::~KrylovAccelerator2()
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
KrylovAccelerator2::newStep(LinearSOE &theSOE)
{
  int newNumEqns = theSOE.getNumEqn();

  if (numEqns != newNumEqns) {
    if (v != 0) {
      for (int i = 0; i < maxDimension+1; i++)
	delete v[i];
      delete [] v;
      v = 0;
    }
    
    if (Av != 0) {
      for (int i = 0; i < maxDimension+1; i++)
	delete Av[i];
      delete [] Av;
      Av = 0;
    }
    
    if (AvData != 0) {
      delete [] AvData;
      AvData = 0;
    }
    
    if (rData != 0) {
      delete [] rData;
      rData = 0;
    }
    
    if (work != 0) {
      delete [] work;
      work = 0;
    }
  }

  numEqns = newNumEqns;
  if (maxDimension > numEqns)
    maxDimension = numEqns;

  if (v == 0) {
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
  if (lwork < 1)
    lwork = 1;

  if (work == 0)
    work = new double [lwork];

  // Reset dimension of subspace
  dimension = 0;

  return 0;
}

#ifdef _WIN32

extern "C" int DGELS(char *T, unsigned int *SZ, int *M, int *N, int *NRHS,
			      double *A, int *LDA, double *B, int *LDB,
			      double *WORK, int *LWORK, int *INFO);

#else

extern "C" int dgels_(char *T, int *M, int *N, int *NRHS,
		      double *A, int *LDA, double *B, int *LDB,
		      double *WORK, int *LWORK, int *INFO);

#endif

int
KrylovAccelerator2::accelerate(Vector &vStar, LinearSOE &theSOE, 
			      IncrementalIntegrator &theIntegrator)
{
  const Vector &R = theSOE.getB();

  int k = dimension;

  // Store residual for differencing at next iteration
  *(Av[k]) = R;

  // If subspace is not empty
  if (dimension > 0) {

    // Compute Av_k = f(y_{k-1}) - f(y_k) = r_{k-1} - r_k
    Av[k-1]->addVector(1.0, R, -1.0);
    
    int i,j;
    
    // Put subspace vectors into AvData
    Matrix A(AvData, numEqns, k);
    for (i = 0; i < k; i++) {
      Vector &Ai = *(Av[i]);
      for (j = 0; j < numEqns; j++)
	A(j,i) = Ai(j);
    }

    for (i = 0; i < k; i++) {
      for (int j = i+1; j < k; j++) {
	double sum = 0.0;
	double sumi = 0.0;
	double sumj = 0.0;
	for (int ii = 0; ii < numEqns; ii++) {
	  sum += A(ii,i)*A(ii,j);
	  sumi += A(ii,i)*A(ii,i);
	  sumj += A(ii,j)*A(ii,j);
	}
	sumi = sqrt(sumi);
	sumj = sqrt(sumj);
	sum = sum/(sumi*sumj);
	//if (fabs(sum) > 0.99)
	  //opserr << sum << ' ' << i << ' ' << j << "   ";
      }
    }

    // Put residual vector into rData (need to save r for later!)
    Vector B(rData, numEqns);
    B = R;
    
    // No transpose
    char trans[] = "N";
    
    // The number of right hand side vectors
    int nrhs = 1;
    
    // Leading dimension of the right hand side vector
    int ldb = (numEqns > k) ? numEqns : k;
    
    // Subroutine error flag
    int info = 0;
    
    // Call the LAPACK least squares subroutine
#ifdef _WIN32
    unsigned int sizeC = 1;
    DGELS(trans, &sizeC, &numEqns, &k, &nrhs, AvData, &numEqns,
	  rData, &ldb, work, &lwork, &info);
#else
    //SUBROUTINE DGELS( TRANS, M, N, NRHS, A, LDA, B, LDB, WORK, LWORK,
    //		      $                  INFO )

    dgels_(trans, &numEqns, &k, &nrhs, AvData, &numEqns,
	   rData, &ldb, work, &lwork, &info);
#endif
    
    // Check for error returned by subroutine
    if (info < 0) {
      opserr << "WARNING KrylovAccelerator2::accelerate() - \n";
      opserr << "error code " << info << " returned by LAPACK dgels\n";
      return info;
    }
    
    Vector Q(numEqns);
    Q = R;

    // Compute the correction vector
    double cj;
    for (j = 0; j < k; j++) {
      
      // Solution to least squares is written to rData
      cj = rData[j];
      
      // Compute w_{k+1} = c_1 v_1 + ... + c_k v_k
      vStar.addVector(1.0, *(v[j]), cj);

      // Compute least squares residual
      // q_{k+1} = r_k - (c_1 Av_1 + ... + c_k Av_k)
      Q.addVector(1.0, *(Av[j]), -cj);
    }

    theSOE.setB(Q);
    //opserr << "Q: " << Q << endln;
  }

  theSOE.solve();
  vStar.addVector(1.0, theSOE.getX(), 1.0);

  // Put accelerated vector into storage for next iteration
  *(v[k]) = vStar;

  dimension++;

  return 0; 
}

int
KrylovAccelerator2::updateTangent(IncrementalIntegrator &theIntegrator)
{
  if (dimension > maxDimension) {
    dimension = 0;
    if (theTangent != NO_TANGENT) {
      //opserr << "KrylovAccelerator2::updateTangent() tangent formed" << endln;
      theIntegrator.formTangent(theTangent);
      return 1;
    }
    else
      return 0;
  }
  else
    return 0;
}

bool
KrylovAccelerator2::updateTangent(void)
{
  if (dimension > maxDimension) {
    dimension = 0;
    return true;
  }
  else
    return false;
}

void
KrylovAccelerator2::Print(OPS_Stream &s, int flag)
{
  s << "KrylovAccelerator2" << endln;
  s << "\tMax subspace dimension: " << maxDimension << endln;
}

int
KrylovAccelerator2::sendSelf(int commitTag, Channel &theChannel)
{
  return -1;
}

int
KrylovAccelerator2::recvSelf(int commitTag, Channel &theChannel, 
			     FEM_ObjectBroker &theBroker)
{
  return -1;
}
