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
                                                                        
// $Revision: 1.9 $
// $Date: 2008-05-11 19:52:54 $
// $Source: /usr/local/cvs/OpenSees/SRC/reliability/analysis/misc/MatrixOperations.h,v $


//
// Written by Terje Haukaas (haukaas@ce.berkeley.edu)
//

#ifndef MatrixOperations_h
#define MatrixOperations_h

#include <Vector.h>
#include <Matrix.h>

class MatrixOperations
{

public:
	MatrixOperations(Matrix passedMatrix);
	~MatrixOperations();
	
	int setMatrix(Matrix passedMatrix);

	int computeLowerCholesky();
	int computeInverseLowerCholesky();
	int computeCholeskyAndItsInverse();
	int computeMatrixNorm();
	int computeInverse();
	int computeTranspose();
	int computeSquareRoot();
	int computeTrace();
	

	Matrix getMatrix();
	Matrix getLowerCholesky();
	Matrix getInverseLowerCholesky();
	Matrix getInverse();
	Matrix getTranspose();
	Matrix getSquareRoot();
	double getMatrixNorm();
	double getTrace();

	// -- eigen analysis--
	int performEigenAnalysis (int beginMod, int endMod);
    int setSizeOfEigenVector(int newSize);
    const Vector & getEigenvector(int mode);

    double getEigenvalue(int mode);


protected:

private:
	Matrix *theMatrix;
	Matrix *theLowerCholesky;
	Matrix *theInverseLowerCholesky;
	Matrix *theInverse;
	Matrix *theTranspose;
	Matrix *theSquareRoot;
	double theMatrixNorm;
	double theTrace; 

	// eigen analysis -- Quan
	int endMod;
	int beginMod;
    int numModes;

	double *eigenvalue;
	double *eigenvector;
	Vector *eigenV;

};

#endif

