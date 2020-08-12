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
                                                                        
// $Revision: 1.12 $
// $Date: 2007/07/16 22:57:03 $
// $Source: /usr/local/cvs/OpenSees/SRC/matrix/Matrix.h,v $
                                                                        
                                                                        
#ifndef Matrix_h
#define Matrix_h 

// Written: fmk 
// Created: 11/96
// Revision: A
//
// Description: This file contains the class definition for Matrix.
// Matrix is a concrete class implementing the matrix abstraction.
// Matrix class is used to provide the abstraction for the most
// general type of matrix, that of an unsymmetric full matrix.
//
// What: "@(#) Matrix.h, revA"

#include <OPS_Globals.h>

class Vector;
class ID;
class Message;

#define MATRIX_VERY_LARGE_VALUE 1.0e213

class Matrix
{
  public:
    // constructors and destructor
    Matrix();	
    Matrix(int nrows, int ncols);
    Matrix(double *data, int nrows, int ncols);    
    Matrix(const Matrix &M);    
#ifdef USE_CXX11
    Matrix( Matrix &&M);    
#endif
    ~Matrix();

    // utility methods
    int setData(double *newData, int nRows, int nCols);
    inline int noRows() const;
    inline int noCols() const;
    void Zero(void);
    int resize(int numRow, int numCol);
    Vector diagonal() const;
    
    int  Assemble(const Matrix &,const ID &rows, const ID &cols, 
		  double fact = 1.0);  
    
    int Solve(const Vector &V, Vector &res) const;
    int Solve(const Matrix &M, Matrix &res) const;
    int Invert(Matrix &res) const;

    int addMatrix(double factThis, const Matrix &other, double factOther);
    int addMatrixTranspose(double factThis, const Matrix &other, double factOther);
    int addMatrixProduct(double factThis, const Matrix &A, const Matrix &B, double factOther); // AB
    int addMatrixTransposeProduct(double factThis, const Matrix &A, const Matrix &B, double factOther); // A'B
    int addMatrixTripleProduct(double factThis, const Matrix &A, const Matrix &B, double factOther); // A'BA
    int addMatrixTripleProduct(double factThis, const Matrix &A, const Matrix &B, const Matrix &C, double otherFact); //A'BC
#if _DLL
	inline double* GetData() { return this->data; }
	void Print() {
		opserr << "[ ";
		for (int i = 0; i < this->numRows; i++)
		{
			for (int j = 0; j < this->numCols - 1; j++)
			{
				opserr << this->operator()(i, j) << ", ";
			}
			opserr << this->operator()(i, this->numCols - 1) << " ";
			opserr << ";" << endln;
		}
		opserr << "] " << endln;
	}
#endif
    // overloaded operators 
    inline double &operator()(int row, int col);
    inline double operator()(int row, int col) const;
    Matrix operator()(const ID &rows, const ID & cols) const;
    
    Matrix &operator=(const Matrix &M);

#ifdef USE_CXX11
    Matrix &operator=(Matrix &&M);
#endif
    
    // matrix operations which will preserve the derived type and
    // which can be implemented efficiently without many constructor calls.

    // matrix-scalar operations
    Matrix &operator+=(double fact);
    Matrix &operator-=(double fact);
    Matrix &operator*=(double fact);
    Matrix &operator/=(double fact); 

    // matrix operations which generate a new Matrix. They are not the
    // most efficient to use, as constructors must be called twice. They
    // however are usefull for matlab like expressions involving Matrices.

    // matrix-scalar operations
    Matrix operator+(double fact) const;
    Matrix operator-(double fact) const;
    Matrix operator*(double fact) const;
    Matrix operator/(double fact) const;
    
    // matrix-vector operations
    Vector operator*(const Vector &V) const;
    Vector operator^(const Vector &V) const;    

    
    // matrix-matrix operations
    Matrix operator+(const Matrix &M) const;
    Matrix operator-(const Matrix &M) const;
    Matrix operator*(const Matrix &M) const;
//     Matrix operator/(const Matrix &M) const;    
    Matrix operator^(const Matrix &M) const;
    Matrix &operator+=(const Matrix &M);
    Matrix &operator-=(const Matrix &M);

    // methods to read/write to/from the matrix
    void Output(OPS_Stream &s) const;
    //    void Input(istream &s);
    
    // methods added by Remo
    int  Assemble(const Matrix &V, int init_row, int init_col, double fact = 1.0);
    int  Assemble(const Vector &V, int init_row, int init_col, double fact = 1.0);
    int  AssembleTranspose(const Matrix &V, int init_row, int init_col, double fact = 1.0);
    int  AssembleTranspose(const Vector &V, int init_row, int init_col, double fact = 1.0);
    int  Extract(const Matrix &V, int init_row, int init_col, double fact = 1.0);

    int Eigen3(const Matrix &M);

    friend OPS_Stream &operator<<(OPS_Stream &s, const Matrix &M);
    //    friend istream &operator>>(istream &s, Matrix &M);    
    friend Matrix operator*(double a, const Matrix &M);
    
    
    friend class Vector;    
    friend class Message;
    friend class UDP_Socket;
    friend class TCP_Socket;
    friend class TCP_SocketSSL;
    friend class TCP_SocketNoDelay;
    friend class MPI_Channel;
    friend class MySqlDatastore;
    friend class BerkeleyDbDatastore;

  protected:

  private:
    static double MATRIX_NOT_VALID_ENTRY;
    static double *matrixWork;
    static int *intWork;
    static int sizeDoubleWork;
    static int sizeIntWork;

    int numRows;
    int numCols;
    int dataSize;
    double *data;
    int fromFree;
};


/********* INLINED MATRIX FUNCTIONS ***********/
inline int 
Matrix::noRows() const 
{
  return numRows;
}

inline int 
Matrix::noCols() const 
{
  return numCols;
}


inline double &
Matrix::operator()(int row, int col)
{ 
#ifdef _G3DEBUG
  if ((row < 0) || (row >= numRows)) {
    opserr << "Matrix::operator() - row " << row << " our of range [0, " <<  numRows-1 << endln;
    return data[0];
  } else if ((col < 0) || (col >= numCols)) {
    opserr << "Matrix::operator() - row " << col << " our of range [0, " <<  numCols-1 << endln;
    return MATRIX_NOT_VALID_ENTRY;
  }
#endif
  return data[col*numRows + row];
}


inline double 
Matrix::operator()(int row, int col) const
{ 
#ifdef _G3DEBUG
  if ((row < 0) || (row >= numRows)) {
    opserr << "Matrix::operator() - row " << row << " our of range [0, " <<  numRows-1 << endln;
    return data[0];
  } else if ((col < 0) || (col >= numCols)) {
    opserr << "Matrix::operator() - row " << col << " our of range [0, " <<  numCols-1 << endln;
    return MATRIX_NOT_VALID_ENTRY;
  }
#endif
  return data[col*numRows + row];
}

#endif




