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
** Reliability module developed by:                                   **
**   Terje Haukaas (haukaas@ce.berkeley.edu)                          **
**   Armen Der Kiureghian (adk@ce.berkeley.edu)                       **
**                                                                    **
** ****************************************************************** */
                                                                        
// $Revision: 1.12 $
// $Date: 2008-05-11 19:52:54 $
// $Source: /usr/local/cvs/OpenSees/SRC/reliability/analysis/misc/MatrixOperations.cpp,v $


//
// Written by Terje Haukaas (haukaas@ce.berkeley.edu)
//

#include <MatrixOperations.h>
#include <Vector.h>
#include <Matrix.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>


MatrixOperations::MatrixOperations(Matrix passedMatrix)
{
	int rows = passedMatrix.noRows();
	int cols = passedMatrix.noCols();

	theMatrix = new Matrix(rows,cols);
	(*theMatrix) = passedMatrix;

	theLowerCholesky = new Matrix(rows,cols);
	theInverseLowerCholesky = new Matrix(rows,cols);
	theInverse = new Matrix(rows,cols);
	theTranspose = new Matrix(rows,cols);
	theSquareRoot = new Matrix(rows,cols);

	theMatrixNorm = 0;
	theTrace = 0; 

	// Set number of modes
	numModes = passedMatrix.noRows();
	eigenV = 0;
	eigenvalue = 0;
	eigenvector = 0;
}


MatrixOperations::~MatrixOperations()
{
	delete theMatrix;
	delete theLowerCholesky;
	delete theInverseLowerCholesky;
	delete theInverse;
	delete theTranspose;
	delete theSquareRoot;
    
	if (eigenvalue != 0)
      delete [] eigenvalue;
    if (eigenvector != 0)
      delete [] eigenvector;
    if (eigenV != 0)
      delete eigenV;
}



int
MatrixOperations::setMatrix(Matrix passedMatrix)
{
	int rows = passedMatrix.noRows();
	int cols = passedMatrix.noCols();

	// delete
	delete theMatrix;
	delete theLowerCholesky;
	delete theInverseLowerCholesky;
	delete theInverse;
	delete theTranspose;
	delete theSquareRoot;

	// reallocate
	theMatrix = new Matrix(rows,cols);
	(*theMatrix) = passedMatrix;
	theLowerCholesky = new Matrix(rows,cols);
	theInverseLowerCholesky = new Matrix(rows,cols);
	theInverse = new Matrix(rows,cols);
	theTranspose = new Matrix(rows,cols);
	theSquareRoot = new Matrix(rows,cols);

	return 0;
}


Matrix 
MatrixOperations::getMatrix()
{
	return (*theMatrix);
}


Matrix 
MatrixOperations::getLowerCholesky()
{
	if (theLowerCholesky == 0) {
		opserr << "MatrixOperations::getLowerCholesky() - this" << endln
			<< " matrix has not been computed." << endln;
		return (*theMatrix);
	}

	return (*theLowerCholesky);
}


Matrix 
MatrixOperations::getInverseLowerCholesky()
{
	if (theInverseLowerCholesky == 0) {
		opserr << "MatrixOperations::getInverseLowerCholesky() - this" << endln
			<< " matrix has not been computed." << endln;
		return (*theMatrix);
	}

	return (*theInverseLowerCholesky);
}


Matrix 
MatrixOperations::getInverse()
{
	if (theInverse == 0) {
		opserr << "MatrixOperations::getInverse() - this" << endln
			<< " matrix has not been computed." << endln;
		return (*theMatrix);
	}
	
	return (*theInverse);
}


Matrix 
MatrixOperations::getTranspose()
{
	if (theTranspose == 0) {
		opserr << "MatrixOperations::getTranspose() - this" << endln
			<< " matrix has not been computed." << endln;
		return (*theMatrix);
	}
	
	return (*theTranspose);
}


Matrix 
MatrixOperations::getSquareRoot()
{
	if (theSquareRoot == 0) {
		opserr << "MatrixOperations::getSquareRoot() - this" << endln
			<< " matrix has not been computed." << endln;
		return (*theMatrix);
	}
	
	return (*theSquareRoot);
}


double
MatrixOperations::getMatrixNorm()
{
	return theMatrixNorm;
}


double
MatrixOperations::getTrace()
{
	return theTrace;
}


int
MatrixOperations::computeLowerCholesky()
{
	Matrix passedMatrix = (*theMatrix);

	// This algorithm is more or less like given in Professor Der Kiureghians
	// handout/brief in CE193 on Cholesky decomposition.

	// Should first check that the matrix is quadratic, etc. 

	int sizeOfPassedMatrix = passedMatrix.noCols();
	int i, j, k;
	double sumOfLambda_i_k_squared = 0;
	double sumOfLambda_i_kLambda_j_k = 0;
	Matrix lambda( sizeOfPassedMatrix , sizeOfPassedMatrix );

	for ( i=0 ; i<sizeOfPassedMatrix ; i++ ) 
	{
		for ( j=0 ; j<sizeOfPassedMatrix ; j++ )
		{
			sumOfLambda_i_k_squared = 0;
			sumOfLambda_i_kLambda_j_k = 0;
			lambda(i,j) = 0.0;
			for ( k=0 ; k<i ; k++ )
			{
				sumOfLambda_i_k_squared += pow ( lambda(i,k) , 2 );
			}
			for ( k=0 ; k<j ; k++ )
			{
				sumOfLambda_i_kLambda_j_k += lambda(i,k) * lambda(j,k);
			}
			if ( i == j )
			{
				if ( ( passedMatrix(i,j) - sumOfLambda_i_k_squared ) < 1.0e-8) {
					opserr << "WARNING: MatrixOperations::computeLowerCholesky()" << endln
						<< " ... matrix may be close to singular. " << endln;
				}
				lambda(i,j) = sqrt ( passedMatrix(i,j) - sumOfLambda_i_k_squared );
			}
			if ( i > j )
			{
				if ( fabs(lambda(j,j)) < 1.0e-8) {
					opserr << "WARNING: MatrixOperations::computeLowerCholesky()" << endln
						<< " ... matrix may be close to singular. " << endln;
				}
				lambda(i,j) = ( passedMatrix(i,j) - sumOfLambda_i_kLambda_j_k ) / lambda(j,j);
			}
			if ( i < j )
			{
				lambda(i,j) = 0.0;
			}
		}
	}

	(*theLowerCholesky) = lambda;

	return 0;
}




int
MatrixOperations::computeInverseLowerCholesky()
{
	Matrix passedMatrix = (*theMatrix);

	// This algorithm is more or less like given in Professor Der Kiureghians
	// handout/brief in CE193 on Cholesky decomposition.

	// Should first check that the matrix is quadratic, etc. 

	int sizeOfPassedMatrix = passedMatrix.noCols();
	
	this->computeLowerCholesky();
	Matrix lambda = this->getLowerCholesky();

	Matrix gamma( sizeOfPassedMatrix , sizeOfPassedMatrix );
	int i, j, k;
	double sumOfLambda_i_kGamma_k_j = 0;
	for ( i=0 ; i<sizeOfPassedMatrix ; i++ ) 
	{
		for ( j=0 ; j<sizeOfPassedMatrix ; j++ )
		{
			sumOfLambda_i_kGamma_k_j = 0;
			gamma(i,j) = 0.0;
			for ( k=j ; k<i ; k++ )
			{
				sumOfLambda_i_kGamma_k_j += lambda(i,k) * gamma(k,j);
			}
			if ( i == j )
			{
				gamma(i,j) = 1 / lambda(i,i);
			}
			if ( i > j )
			{
				if ( fabs(lambda(i,i)) < 1.0e-8) {
					opserr << "WARNING: MatrixOperations::computeInverseLowerCholesky()" << endln
						<< " ... matrix may be close to singular. " << endln;
				}
				gamma(i,j) = - sumOfLambda_i_kGamma_k_j / lambda(i,i);
			}
			if ( i < j )
			{
				gamma(i,j) = 0.0;
			}
		}
	}

	(*theInverseLowerCholesky) = gamma;

	return 0;
}


int
MatrixOperations::computeCholeskyAndItsInverse()
{
	Matrix &passedMatrix = (*theMatrix);

	// This algorithm is more or less like given in Professor Der Kiureghians
	// handout/brief in CE193 on Cholesky decomposition.

	// Should first check that the matrix is quadratic, etc. 

	int sizeOfPassedMatrix = passedMatrix.noCols();
	double sumOfLambda_i_k_squared = 0.0;
	double sumOfLambda_i_kLambda_j_k = 0.0;
	double sumOfLambda_i_kGamma_k_j = 0.0;
	Matrix lambda( sizeOfPassedMatrix , sizeOfPassedMatrix );
	Matrix gamma( sizeOfPassedMatrix , sizeOfPassedMatrix );
	int i, j, k;

	// Lower Cholesky
	for ( i=0 ; i<sizeOfPassedMatrix ; i++ ) 
	{
		for ( j=0 ; j<sizeOfPassedMatrix ; j++ )
		{
			sumOfLambda_i_k_squared = 0.0;
			sumOfLambda_i_kLambda_j_k = 0.0;
			lambda(i,j) = 0.0;
			for ( k=0 ; k<i ; k++ )
			{
				sumOfLambda_i_k_squared += pow ( lambda(i,k) , 2.0 );
			}
			for ( k=0 ; k<j ; k++ )
			{
				sumOfLambda_i_kLambda_j_k += lambda(i,k) * lambda(j,k);
			}
			if ( i == j )
			{
				if ( ( passedMatrix(i,j) - sumOfLambda_i_k_squared ) < 1.0e-8) {
					opserr << "WARNING: MatrixOperations::computeCholeskyAndItsInverse()" << endln
						<< " ... matrix may be close to singular. " << endln;
				}
				lambda(i,j) = sqrt ( passedMatrix(i,j) - sumOfLambda_i_k_squared );
			}
			if ( i > j )
			{
				if ( fabs(lambda(j,j)) < 1.0e-8) {
					opserr << "WARNING: MatrixOperations::computeCholeskyAndItsInverse()" << endln
						<< " ... matrix may be close to singular. " << endln;
				}
				lambda(i,j) = ( passedMatrix(i,j) - sumOfLambda_i_kLambda_j_k ) / lambda(j,j);
			} 
			if ( i < j )
			{
				lambda(i,j) = 0.0;
			}
		}
	}

	// Inverse lower Cholesky
	for ( i=0 ; i<sizeOfPassedMatrix ; i++ ) 
	{
		for ( j=0 ; j<sizeOfPassedMatrix ; j++ )
		{
			sumOfLambda_i_kGamma_k_j = 0;
			gamma(i,j) = 0.0;
			for ( k=j ; k<i ; k++ )
			{
				sumOfLambda_i_kGamma_k_j += lambda(i,k) * gamma(k,j);
			}
			if ( i == j )
			{
				gamma(i,j) = 1 / lambda(i,i);
			}
			if ( i > j )
			{
				if ( fabs(lambda(i,i)) < 1.0e-8) {
					opserr << "WARNING: MatrixOperations::computeCholeskyAndItsInverse()" << endln
						<< " ... matrix may be close to singular. " << endln;
				}
				gamma(i,j) = - sumOfLambda_i_kGamma_k_j / lambda(i,i);
			}
			if ( i < j )
			{
				gamma(i,j) = 0.0;
			}
		}
	}

	(*theLowerCholesky) = lambda;
	(*theInverseLowerCholesky) = gamma;

	return 0;
}






int
MatrixOperations::computeTranspose()
{
	Matrix &A = (*theMatrix);

	int sizeOfA = A.noCols();

	Matrix B(sizeOfA,sizeOfA);

	for (int i=0; i<sizeOfA; i++) {
		for (int j=0; j<sizeOfA; j++) {
			B(i,j) = A(j,i);
		}
	}

	(*theTranspose) = B;

	return 0;
}





int
MatrixOperations::computeSquareRoot()
{
	Matrix &A = (*theMatrix);

	int sizeOfA = A.noCols();

	Matrix B(sizeOfA,sizeOfA);

	for (int i=0; i<sizeOfA; i++) {
		for (int j=0; j<sizeOfA; j++) {
			B(i,j) = sqrt( A(i,j) );
		}
	}

	(*theSquareRoot) = B;

	return 0;
}



int
MatrixOperations::computeInverse()
{

	Matrix &A = (*theMatrix);

	// Return the inverse matrix B such that A*B=I
	int sizeOfA = A.noCols();
	Matrix B ( sizeOfA, sizeOfA );
	Matrix AB ( sizeOfA, 2*sizeOfA );
	int i, j, k;

	// Set up the AB matrix
	for ( i=0 ; i<sizeOfA ; i++ )
	{
		for ( j=0 ; j<(sizeOfA*2) ; j++ )
		{
			if ( j < sizeOfA )
			{
				AB(i,j) = A(i,j);
			}
			else
			{
				if ( j == (sizeOfA+i) )
				{
					AB(i,j) = 1.0;
				}
				else
				{
					AB(i,j) = 0.0;
				}
			}
		}
	}

	// The Gauss-Jordan method
	double pivot;
	double ABii;

	for ( k=0 ; k<sizeOfA ; k++ )
	{
		for ( i=k ; i<sizeOfA ; i++ )
		{
			ABii = AB(i,i);
			pivot = AB(i,k);
			for ( j=k ; j<(sizeOfA*2) ; j++ )
			{
				if ( i == k )
				{
					AB(i,j) = AB(i,j) / ABii;
				}
				else
				{
					AB(i,j) = AB(i,j) - pivot * AB(k,j);
				}
			}
		}
	}
	for ( k=1 ; k<sizeOfA ; k++ )
	{
		for ( i=k ; i<sizeOfA ; i++ )
		{
			pivot = AB( (sizeOfA-i-1), (sizeOfA-k) );
			for ( j=(sizeOfA-k) ; j<(sizeOfA*2) ; j++ )
			{
					AB((sizeOfA-i-1),(j)) = 
						AB((sizeOfA-i-1),(j)) - 
						pivot * AB((sizeOfA-k),(j));
			}
		}
	}


	// Collect the B matrix
	for ( i=0 ; i<sizeOfA ; i++ )
	{
		for ( j=sizeOfA ; j<(sizeOfA*2) ; j++ )
		{
			B(i,(j-sizeOfA)) = AB(i,j);
		}
	}

	(*theInverse) = B;
	
	return 0;
}



int
MatrixOperations::computeMatrixNorm()
{
	Matrix &passedMatrix = (*theMatrix);

	int numberOfColumns = passedMatrix.noCols();
	int numberOfRows = passedMatrix.noRows();
	double sum = 0;
	for ( int i=0 ; i<numberOfRows ; i++ )
	{
		for ( int j=0 ; j<numberOfColumns ; j++ )
		{
			sum += passedMatrix(i,j) * passedMatrix(i,j);
		}
	}

	theMatrixNorm = sqrt(sum);

	return 0;
}


int
MatrixOperations::computeTrace()
{
	Matrix &passedMatrix = (*theMatrix);

	int numberOfColumns = passedMatrix.noCols();
	int numberOfRows = passedMatrix.noRows();
	
	if (numberOfColumns != numberOfRows) {
		opserr << "MatrixOperations::computeTrace() - can not" << endln
			<< " compute the trace of a non-quadratic matrix." << endln;
		return -1;
	}

	double product = 1.0;
	for ( int i=0 ; i<numberOfRows ; i++ )
	{
		product = product * passedMatrix(i,i);
	}

	theTrace = product;

	return 0;
}





/// === Quan Apr. 2007 ==

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

int MatrixOperations::performEigenAnalysis(int pBeginMod, int pEndMod)
{

  beginMod=pBeginMod;
  endMod = pEndMod;
 
  Matrix * A =  theMatrix;

//  opserr<< "hessian:"<< * theMatrix<<endln;

  if (A ==0) {opserr<<"error, Hessian does not exist in MatrixOperator::performEigenAnalysis !"<<endln; exit(-1);}

  // Number of equations
  int n = A->noRows();
  if (n != A->noCols()) {
	opserr<<"MatrixOperations::performEigenAnalysis wrong. m!=n"<<endln;
	exit(-1);
  }

  numModes = endMod - beginMod +1; //Q
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
  setSizeOfEigenVector(ldz); // alloc memory here.

  // Allocate storage for eigenvectors
  if (eigenvector != 0)
    delete [] eigenvector;
  eigenvector = new double [ldz*numModes];

  // Number of superdiagonals
  int kd = n-1;

  // Matrix data

  
  //double *ab =  theSOE->A;  //
  double *ab = new double [n*n];
  int i;
  for (i=0; i<n*n; i++)
	  ab[i]=0.0;

  

  for (i=0; i<n; i++)
	  for (int j=i; j<n; j++)
		  ab[kd+i-j+n*j]= (*A)(i,j);
  
// opserr<<"Hessian:"<<(*A)<<endln;
 //opserr<<"ab:"<<endln;
 /* for (i=0; i<n; i++)
	  for (int j=0; j<n; j++)
         opserr<< ab[n*i+j]<<endln;*/
  // Leading dimension of the matrix
  int ldab = kd + 1;

  // Leading dimension of q
  int ldq = n;

  // Orthogonal matrix used in reduction to tridiagonal form
  // (see LAPACK dsbevx subroutine documentation)
  double *q = new double [ldq*n];

  // Index ranges [1,numModes] of eigenpairs to compute
  int il = beginMod;
  int iu = endMod;

  // Compute eigenvalues and eigenvectors
  char jobz[3] = "V";

  // Selected eigenpairs are based on index range [il,iu]
  char range[3] = "I";

  // Upper triagle of matrix is stored
  char uplo[3] = "U";
  
  // Return value
  int *ifail = new int [n];
  int info = 0;

  // Number of eigenvalues returned
  int m = 0;

  // Not used
  double vl = 0.0;
  double vu = 1.0;

  // Not used 
  double abstol = -1.0;


  
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
    opserr << "Hessian::performEigenAnalysis() -- invalid argument number " << -info << " passed to LAPACK dsbevx\n";
    return info;
  }

  if (info > 0) {
    opserr << "Hessian::performEigenAnalysis() -- LAPACK dsbevx returned error code " << info << endln;
    return -info;
  }

  if (m < numModes) {
    opserr << "Hessian::performEigenAnalysis() -- LAPACK dsbevx only computed " << m << " eigenvalues, " <<
      numModes << "were requested\n";

    numModes = m;
  }

   

  
  return 0;

}



const Vector &
MatrixOperations::getEigenvector(int mode)
{
  if (mode < beginMod || mode > endMod) {
    opserr << "Hessian::getEigenVector() -- mode " << mode << " is out of range (1 - "
	   << numModes << ")\n";

    eigenV->Zero();
    return *eigenV;  
  }
  
  int size = theMatrix->noRows();

  int index = (mode - beginMod) * size;
  
  Vector &vec = *eigenV;
  if (eigenvector != 0) {
    for (int i = 0; i < size; i++) {
      vec(i) = eigenvector[index++];
    }	
  }
  else {
    opserr << "Hessian::getEigenVector() -- eigenvectors not yet computed\n";
    eigenV->Zero();
  }      
  
  return *eigenV;  
}

double
MatrixOperations::getEigenvalue(int mode)
{
  if (mode < beginMod || mode > endMod) {
    opserr << "Hessian::getEigenvalue() -- mode " << mode << " is out of range (1 - "
	   << numModes << ")\n";

    return 0.0;
  }
  
  if (eigenvalue != 0)
    return eigenvalue[mode-beginMod];
  else {
    opserr << "Hessian::getEigenvalue() -- eigenvalues not yet computed\n";
    return 0.0;
  }      
}

int
MatrixOperations::setSizeOfEigenVector(int newSize)
{
    

  if (eigenV == 0 || eigenV->Size() != newSize) {
    if (eigenV != 0)
      delete eigenV;
    
    eigenV = new Vector(newSize);
    if (eigenV == 0 || eigenV->Size() != newSize) {
      opserr << "Hessian::ssetSize() -- ran out of memory for eigenvector of size " << newSize << endln;
      return -2;	    
    }
  }
  
  return 0;
}
