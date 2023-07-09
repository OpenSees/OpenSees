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
                                                                        
// $Revision: 1.9 $
// $Date: 2009-05-20 17:31:36 $
// $Source: /usr/local/cvs/OpenSees/SRC/system_of_eqn/eigenSOE/SymBandEigenSolver.cpp,v $

// Written: MHS
// Created: Oct 2001
//
// Description: This file contains the class definition for
// SymBandEigenSolver, which computes the eigenvalues and eigenvectors
// of a symmetric banded matrix using the LAPACK subroutine dsbevx.

#include <SymBandEigenSolver.h>
#include <math.h>
#include <stdio.h>
#include <AnalysisModel.h>
#include <DOF_GrpIter.h>
#include <DOF_Group.h>
#include <FE_EleIter.h>
#include <FE_Element.h>
#include <Integrator.h>

SymBandEigenSolver::SymBandEigenSolver()
:EigenSolver(EigenSOLVER_TAGS_SymBandEigenSolver),
 theSOE(0), numModes(0), eigenvalue(0), eigenvector(0), eigenV(0)
{

}

SymBandEigenSolver::~SymBandEigenSolver()
{
  if (eigenvalue != 0)
    delete [] eigenvalue;
  if (eigenvector != 0)
    delete [] eigenvector;
  if (eigenV != 0)
    delete eigenV;
}

#ifdef _WIN32

extern "C" int DSBEVX(char *jobz, char *range, char *uplo, int *n, int *kd,
			       double *ab, int *ldab, double *q, int *ldq,
			       double *vl, double *vu, int *il, int *iu, double *abstol,
			       int *m, double *w, double *z, int *ldz,
			       double *work, int *iwork, int *ifail, int *info);

#else

extern "C" int dsbevx_(char *jobz, char *range, char *uplo, int *n, int *kd,
		       double *ab, int *ldab, double *q, int *ldq,
		       double *vl, double *vu, int *il, int *iu, double *abstol,
		       int *m, double *w, double *z, int *ldz,
		       double *work, int *iwork, int *ifail, int *info);

#endif

int
SymBandEigenSolver::solve(int nModes, bool generalized, bool findSmallest)
{

  if (generalized == true) {
    opserr << "SymBandEigenSolver::solve() - only does standard problem\n";
    return -1;
  }

  if (theSOE == 0) {
    opserr << "SymBandEigenSolver::solve() -- no EigenSOE has been set yet\n";
    return -1;
  }
  
  // Set number of modes
  numModes = nModes;

  // Number of equations
  int n = theSOE->size;

  // Check for quick return
  if (numModes < 1) {
    numModes = 0;
    return 0;
  }

  // Simple check
  if (numModes > n)
    numModes = n;

  // Allocate storage for eigenvalues
  if (eigenvalue != 0)
    delete [] eigenvalue;
  eigenvalue = new double [n];

  // Real work array (see LAPACK dsbevx subroutine documentation)
  double *work = new double [7*n];

  // Integer work array (see LAPACK dsbevx subroutine documentation)
  int *iwork = new int [5*n];

  // Leading dimension of eigenvectors
  int ldz = n;

  // Allocate storage for eigenvectors
  if (eigenvector != 0)
    delete [] eigenvector;
  eigenvector = new double [ldz*numModes];

  // Number of superdiagonals
  int kd = theSOE->numSuperD;

  // Matrix data
  double *ab = theSOE->A;

  // Leading dimension of the matrix
  int ldab = kd + 1;

  // Leading dimension of q
  int ldq = n;

  // Orthogonal matrix used in reduction to tridiagonal form
  // (see LAPACK dsbevx subroutine documentation)
  double *q = new double [ldq*n];

  // Index ranges [1,numModes] of eigenpairs to compute
  int il = 1;
  int iu = numModes;

  // Compute eigenvalues and eigenvectors
  char *jobz = "V";

  // Selected eigenpairs are based on index range [il,iu]
  char *range = "I";

  // Upper triagle of matrix is stored
  char *uplo = "U";
  
  // Return value
  int *ifail = new int [n];
  int info = 0;

  // Number of eigenvalues returned
  int m = 0;

  // Not used
  double vl = 0.0;
  double vu = 1.0;

  // Not used ... I think!
  double abstol = -1.0;


  // if Mass matrix we make modifications to A:
  //         A -> M^(-1/2) A M^(-1/2)
  double *M = theSOE->M;
  double *A = theSOE->A;
  int numSuperD = theSOE->numSuperD;
  int size = n;
  if (M != 0) {
    int i,j;
    bool singular = false;
    // form M^(-1/2) and check for singular mass matrix
    for (int k=0; k<size; k++) {
      if (M[k] == 0.0) {
	singular = true;
	// alternative is to set as a small no ~ 1e-10 times smallest m(i,i) != 0.0
	opserr << "SymBandEigenSolver::solve() - M matrix singular\n";
	return -1;
      } else {
	M[k] = 1.0/sqrt(M[k]);
      }
    }

    // make modifications to A
    //   Aij -> Mi Aij Mj  (based on new M)
    for (i=0; i<size; i++) {
      double *AijPtr = A +(i+1)*(numSuperD+1) - 1;
      int minColRow = i - numSuperD;
      if (minColRow < 0) minColRow = 0;
      for (j=i; j>=minColRow; j--) {
	*AijPtr *= M[j]*M[i];
	AijPtr--;
      }
    }
  }

  // Call the LAPACK eigenvalue subroutine
#ifdef _WIN32
  unsigned int sizeC = 1;
  DSBEVX(jobz, range, uplo, &n, &kd, ab, &ldab,
	 q, &ldq, &vl, &vu, &il, &iu, &abstol, &m,
	 eigenvalue, eigenvector, &ldz, work, iwork, ifail, &info);
#else
  dsbevx_(jobz, range, uplo, &n, &kd, ab, &ldab,
	  q, &ldq, &vl, &vu, &il, &iu, &abstol, &m,
	  eigenvalue, eigenvector, &ldz, work, iwork, ifail, &info);
#endif

  delete [] q;
  delete [] work;
  delete [] iwork;
  delete [] ifail;

  if (info < 0) {
    opserr << "SymBandEigenSolver::solve() -- invalid argument number " << -info << " passed to LAPACK dsbevx\n";
    return info;
  }

  if (info > 0) {
    opserr << "SymBandEigenSolver::solve() -- LAPACK dsbevx returned error code " << info << endln;
    return -info;
  }

  if (m < numModes) {
    opserr << "SymBandEigenSolver::solve() -- LAPACK dsbevx only computed " << m << " eigenvalues, " <<
      numModes << "were requested\n";

    numModes = m;
  }

  theSOE->factored = true;

  // make modifications to the eigenvectors
  //   Eij -> Mi Eij  (based on new M)


  M = theSOE->M;
  if (M != 0) {
    
    for (int j=0; j<numModes; j++) {
      double *eigVectJptr = &eigenvector[j*ldz];
      double *MPtr = M;
      for (int i=0; i<size; i++) 
	*eigVectJptr++ *= *MPtr++;
    }
  }

  return 0;
}

int
SymBandEigenSolver::setEigenSOE(SymBandEigenSOE &theBandSOE)
{
  theSOE = &theBandSOE;

  return 0;
}

const Vector &
SymBandEigenSolver::getEigenvector(int mode)
{
  if (mode < 1 || mode > numModes) {
    opserr << "SymBandEigenSolver::getEigenVector() -- mode " << mode << " is out of range (1 - "
	   << numModes << ")\n";

    eigenV->Zero();
    return *eigenV;  
  }
  
  int size = theSOE->size;

  int index = (mode - 1) * size;
  
  Vector &vec = *eigenV;
  if (eigenvector != 0) {
    for (int i = 0; i < size; i++) {
      vec(i) = eigenvector[index++];
    }	
  }
  else {
    opserr << "SymBandEigenSolver::getEigenVector() -- eigenvectors not yet computed\n";
    eigenV->Zero();
  }      
  
  return *eigenV;  
}

double
SymBandEigenSolver::getEigenvalue(int mode)
{
  if (mode < 1 || mode > numModes) {
    opserr << "SymBandEigenSolver::getEigenvalue() -- mode " << mode << " is out of range (1 - "
	   << numModes << ")\n";

    return 0.0;
  }
  
  if (eigenvalue != 0)
    return eigenvalue[mode-1];
  else {
    opserr << "SymBandEigenSolver::getEigenvalue() -- eigenvalues not yet computed\n";
    return 0.0;
  }      
}

int
SymBandEigenSolver::setSize()
{
  int size = theSOE->size;    

  if (eigenV == 0 || eigenV->Size() != size) {
    if (eigenV != 0)
      delete eigenV;
    
    eigenV = new Vector(size);
    if (eigenV == 0 || eigenV->Size() != size) {
      opserr << "SymBandEigenSolver::ssetSize() -- ran out of memory for eigenvector of size " << size << endln;
      return -2;	    
    }
  }
  
  return 0;
}

int    
SymBandEigenSolver::sendSelf(int commitTag, Channel &theChannel)
{
  return 0;
}

int
SymBandEigenSolver::recvSelf(int commitTag, Channel &theChannel, 
			     FEM_ObjectBroker &theBroker)
{
  return 0;
}
