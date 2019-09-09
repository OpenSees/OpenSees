#include "stdafx.h"
#include "MatrixWrapper.h"
#include "Matrix.h"

using namespace OpenSees;

OpenSees::MatrixWrapper::MatrixWrapper(array<double, 2>^ data)
{
	int nRows = data->GetLength(0);
	int nCols = data->GetLength(1);
	double* _data = new double[nRows*nCols];
	int c = 0;
	for (int i = 0; i < nCols; i++)
		for (int j = 0; j < nRows; j++)
		_data[c++] = data[j,i];
	_Matrix = new Matrix(_data, nRows, nCols);
}

OpenSees::MatrixWrapper::MatrixWrapper()
{
	_Matrix = 0;
}

OpenSees::MatrixWrapper::MatrixWrapper(int row, int col)
{
	_Matrix = new Matrix(row,col);
}

MatrixWrapper::~MatrixWrapper()
{
	if (_Matrix != 0)
		delete _Matrix;
}


