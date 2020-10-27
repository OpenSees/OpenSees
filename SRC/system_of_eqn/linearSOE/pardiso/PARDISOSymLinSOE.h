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

// $Revision: 1.2 $
// $Date: 2006-01-10 00:42:51 $
// $Source: /usr/local/cvs/OpenSees/SRC/system_of_eqn/linearSOE/pardiso/PARDISOSymLinSOE.h,v $
//
// Written: M. Salehi opensees.net@gmail.com
// website : http://opensees.net
// Created: 02/19
// Revision: A
//
//
// Description: This file contains the class definition for PARDISOSymLinSOE
// PARDISOSymLinSOE is a subclass of LinearSOE. It stores the matrix equation
// Ax=b using the sparse row-compacted storage scheme for storing the 
// matrix A. 
//
// What: "@(#) PARDISOSymLinSOE.h, revA"

#ifndef PARDISOSymLinSOE_h
#define PARDISOSymLinSOE_h

#include <LinearSOE.h>
#include <Vector.h>

class PARDISOSymLinSolver;

class PARDISOSymLinSOE : public LinearSOE
{
public:
	PARDISOSymLinSOE(PARDISOSymLinSolver &theSolver);

	~PARDISOSymLinSOE();

	int getNumEqn(void) const;
	int setSize(Graph &theGraph);
	int addA(const Matrix &, const ID &, double fact = 1.0);
	int addB(const Vector &, const ID &, double fact = 1.0);
	int setB(const Vector &, double fact = 1.0);

	void zeroA(void);
	void zeroB(void);

	const Vector &getX(void);
	const Vector &getB(void);
	double normRHS(void);

	void setX(int loc, double value);
	void setX(const Vector &x);
	int setPARDISOSymLinSolver(PARDISOSymLinSolver &newSolver);

	int sendSelf(int commitTag, Channel &theChannel);
	int recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker);
	friend class PARDISOSymLinSolver;
protected:

private:
	int size;            // order of A
	int nnz;             // number of non-zeros in A
	double *A, *B, *X;   // 1d arrays containing coefficients of A, B and X
	int *colA, *rowStartA; // int arrays containing info about coeficientss in A
	Vector *vectX;
	Vector *vectB;
	int Asize, Bsize;    // size of the 1d array holding A
	bool factored;
	int getRow(int index, int* rowStartA, int nrow) {
		for (int i = 0; i < nrow; i++)
		{
			if (index >= rowStartA[i] - 1 && index < rowStartA[i + 1] - 1)
				return i;
		}
		return -1;
	}

	int getSymSize(int* genColA, int* genRowStartA, int nnz, int nrow, int* nsymCol, int* nsymRowStartA)
	{
		int currentRow = -1;
		int nskippedCol = 0;
		(*nsymRowStartA) = 0;
		for (int i = 0; i < nnz; i++)
		{
			int col = genColA[i] - 1;
			int row = getRow(i, genRowStartA, nrow);
			if (row > col)
			{
				nskippedCol++;
			}
			else
			{
				if (row != currentRow) // new row
				{
					(*nsymRowStartA)++;
					currentRow = row;
				}
			}
		}

		(*nsymRowStartA)++;
		*nsymCol = nnz - nskippedCol;

		return 0;
	};

	int getSym(int* genColA, int* genRowStartA, int nnz, int nrow, int* _symColA,
		int* _symRowStartA, int nsymCol, int nsymRowStartA)
	{
		

		int currentRow = -1;
		int skipped = 0;
		int ncola = 0;
		int nrowstarta = 0;
		for (int i = 0; i < nnz; i++)
		{
			int col = genColA[i] - 1;
			int row = getRow(i, genRowStartA, nrow);
			if (row > col)
			{
				skipped++;
			}
			else
			{
				_symColA[ncola++] = col + 1;
				if (row != currentRow) // new row
				{
					_symRowStartA[nrowstarta++] = i - skipped + 1;
					currentRow = row;
				}
			}
		}

		_symRowStartA[nrowstarta++] = nsymCol + 1;

		/*opserr << "genColA :";
		for (int i = 0; i < nnz; i++)
			opserr << genColA[i] << " ";
		opserr << endln;

		opserr << "genRowStartA :";
		for (int i = 0; i < nrow + 1; i++)
			opserr << genRowStartA[i] << " ";
		opserr << endln;

		opserr << "symColA " << nsymCol << ":";
		for (int i = 0; i < nsymCol; i++)
			opserr << _symColA[i] << " ";
		opserr << endln;

		opserr << "symRowStartA " << nsymRowStartA << ":";
		for (int i = 0; i < nsymRowStartA; i++)
			opserr << _symRowStartA[i] << " ";
		opserr << endln;*/

		return 0;
	};
};


#endif

