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
                                                                        
// $Revision: 1.5 $
// $Date: 2003-03-04 00:39:26 $
// $Source: /usr/local/cvs/OpenSees/SRC/reliability/analysis/misc/MatrixOperations.cpp,v $


//
// Written by Terje Haukaas (haukaas@ce.berkeley.edu)
//

#include <MatrixOperations.h>
#include <Vector.h>
#include <Matrix.h>
#include <math.h>
#include <string.h>


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

}


MatrixOperations::~MatrixOperations()
{
	delete theMatrix;
	delete theLowerCholesky;
	delete theInverseLowerCholesky;
	delete theInverse;
	delete theTranspose;
	delete theSquareRoot;
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

	// Return the invers matrix B such that A*B=I
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





/*
double
MatrixOperations::vectorDotProduct(Vector vector1, Vector vector2)
{
	int sizeOfVector1 = vector1.Size();
	int sizeOfVector2 = vector2.Size();
	double sum = 0;

	if ( sizeOfVector1 != sizeOfVector2 )
	{
		opserr << "vectorDotProduct can not be performed" << endln;
	}
	else
	{
		for ( int i=0 ; i<sizeOfVector1 ; i++ )
		{
			sum += vector1(i) * vector2(i);
		}
	}

	return sum;
}
*/


/*
Vector 
MatrixOperations::matrixTimesVector(Matrix passedMatrix, Vector passedVector)
{
	int numberOfMatrixRows = passedMatrix.noRows();
	int sizeOfPassedVector = passedVector.Size();
	Vector result(numberOfMatrixRows);
	int i, j;
	double sum;

	for ( i=0 ; i<numberOfMatrixRows ; i++ )
	{
		sum = 0.0;
		for ( j=0 ; j<sizeOfPassedVector ; j++ )
		{
			sum += passedMatrix(i,j) * passedVector(j);
		}
		result(i) = sum;
	}

	return result;
}
*/
