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
                                                                        
// $Revision: 1.10 $
// $Date: 2007-11-01 17:57:41 $
// $Source: /usr/local/cvs/OpenSees/SRC/reliability/analysis/misc/MatrixOperations.cpp,v $


//
// Written by Terje Haukaas (haukaas@ce.berkeley.edu)
//

#include <MatrixOperations.h>
#include <Vector.h>
#include <Matrix.h>
#include <math.h>
#include <string.h>


MatrixOperations::MatrixOperations(const Matrix &passedMatrix)
{
	int rows = passedMatrix.noRows();
	int cols = passedMatrix.noCols();

	theMatrix = new Matrix(rows,cols);
	(*theMatrix) = passedMatrix;

	theLowerCholesky = new Matrix(rows,cols);
	theInverseLowerCholesky = new Matrix(rows,cols);
}


MatrixOperations::~MatrixOperations()
{
	delete theMatrix;
	delete theLowerCholesky;
	delete theInverseLowerCholesky;
}



int
MatrixOperations::setMatrix(const Matrix &passedMatrix)
{
	int rows = passedMatrix.noRows();
	int cols = passedMatrix.noCols();

	// delete
	delete theMatrix;
	delete theLowerCholesky;
	delete theInverseLowerCholesky;

	// reallocate
	theMatrix = new Matrix(rows,cols);
	(*theMatrix) = passedMatrix;
	theLowerCholesky = new Matrix(rows,cols);
	theInverseLowerCholesky = new Matrix(rows,cols);

	return 0;
}


const Matrix& 
MatrixOperations::getMatrix()
{
	return (*theMatrix);
}


const Matrix& 
MatrixOperations::getLowerCholesky()
{
	if (theLowerCholesky == 0) {
		opserr << "MatrixOperations::getLowerCholesky() - this" << endln
			<< " matrix has not been computed." << endln;
		return (*theMatrix);
	}

	return (*theLowerCholesky);
}


const Matrix&
MatrixOperations::getInverseLowerCholesky()
{
	if (theInverseLowerCholesky == 0) {
		opserr << "MatrixOperations::getInverseLowerCholesky() - this" << endln
			<< " matrix has not been computed." << endln;
		return (*theMatrix);
	}

	return (*theInverseLowerCholesky);
}

int
MatrixOperations::computeLowerCholesky()
{
	const Matrix &passedMatrix = (*theMatrix);

	// This algorithm is more or less like given in Professor Der Kiureghians
	// handout/brief in CE193 on Cholesky decomposition.

	// Should first check that the matrix is quadratic, etc. 

	int sizeOfPassedMatrix = passedMatrix.noCols();
	int i, j, k;
	double sumOfLambda_i_k_squared = 0;
	double sumOfLambda_i_kLambda_j_k = 0;
	//Matrix lambda( sizeOfPassedMatrix , sizeOfPassedMatrix );
	Matrix &lambda = *theLowerCholesky;

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

	//(*theLowerCholesky) = lambda;

	return 0;
}




int
MatrixOperations::computeInverseLowerCholesky()
{
	const Matrix &passedMatrix = (*theMatrix);

	// This algorithm is more or less like given in Professor Der Kiureghians
	// handout/brief in CE193 on Cholesky decomposition.

	// Should first check that the matrix is quadratic, etc. 

	int sizeOfPassedMatrix = passedMatrix.noCols();
	
	this->computeLowerCholesky();
	const Matrix &lambda = this->getLowerCholesky();

	//Matrix gamma( sizeOfPassedMatrix , sizeOfPassedMatrix );
	Matrix &gamma = *theInverseLowerCholesky;
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

	//(*theInverseLowerCholesky) = gamma;

	return 0;
}


int
MatrixOperations::computeCholeskyAndItsInverse()
{
	const Matrix &passedMatrix = (*theMatrix);

	// This algorithm is more or less like given in Professor Der Kiureghians
	// handout/brief in CE193 on Cholesky decomposition.

	// Should first check that the matrix is quadratic, etc. 

	int sizeOfPassedMatrix = passedMatrix.noCols();
	double sumOfLambda_i_k_squared = 0.0;
	double sumOfLambda_i_kLambda_j_k = 0.0;
	double sumOfLambda_i_kGamma_k_j = 0.0;
	//Matrix lambda( sizeOfPassedMatrix , sizeOfPassedMatrix );
	//Matrix gamma( sizeOfPassedMatrix , sizeOfPassedMatrix );
	Matrix &lambda = *theLowerCholesky;
	Matrix &gamma = *theInverseLowerCholesky;
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

	//(*theLowerCholesky) = lambda;
	//(*theInverseLowerCholesky) = gamma;

	return 0;
}
