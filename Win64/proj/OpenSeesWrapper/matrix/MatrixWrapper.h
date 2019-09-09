#pragma once
#include "VectorWrapper.h"
#include "Matrix.h"
using namespace System;
namespace OpenSees {
	public ref class MatrixWrapper
	{
	public:
		MatrixWrapper(array<double, 2>^ data);
		MatrixWrapper(int row, int col);
		MatrixWrapper();
		void Zero() { _Matrix->Zero(); }
		int noRows() { return _Matrix->noRows(); }
		int noCols() { return _Matrix->noCols(); }

		VectorWrapper^ Solve(VectorWrapper^ v)
		{
			int size = v->Size();
			VectorWrapper^ ans = gcnew VectorWrapper(size);
			_Matrix->Solve(*v->_Vector, *ans->_Vector);
			return ans;
		}

		MatrixWrapper^ Solve(MatrixWrapper^ v)
		{
			MatrixWrapper^ ans = gcnew MatrixWrapper(v->noRows(),v->noCols());
			_Matrix->Solve(*v->_Matrix, *ans->_Matrix);
			return ans;
		}

		MatrixWrapper^ Invert()
		{
			MatrixWrapper^ ans = gcnew MatrixWrapper(_Matrix->noRows(), _Matrix->noCols());
			_Matrix->Invert(*ans->_Matrix);
			return ans;
		}

		double Get(int row,int col) {
			return _Matrix->operator()(row, col);
		};

		void Set(int row, int col, double value)
		{
			_Matrix->operator()(row, col) = value;
		};

		int SetData(array<double,2>^ data, double selfFact, double otherFact) {
			if (this->noRows() != data->GetLength(0) || this->noCols() != data->GetLength(1)) return -1;
			for (int i = 0; i < data->GetLength(0); i++)
				for (int j = 0; j < data->GetLength(1); j++)
				{
					if (selfFact == 0)
					{
						if (otherFact == 1)
							this->_Matrix->operator()(i,j) = data[i,j];
						else
							this->_Matrix->operator()(i, j) = data[i,j] * otherFact;
					}
					else
					{
						if (otherFact == 1)
							this->_Matrix->operator()(i, j) = this->_Matrix->operator()(i, j)*selfFact + data[i, j];
						else
							this->_Matrix->operator()(i, j) = this->_Matrix->operator()(i, j)*selfFact + data[i, j] * otherFact;
					}
				}
			return 0;
		}

		property double default[int, int]
		{
			double get(int row,int col) { return _Matrix->operator()(row,col); }
			void set(int row,int col, double value) { _Matrix->operator()(row,col) = value; }
		}

		~MatrixWrapper();

		static array<double, 2>^ GetArray(MatrixWrapper^ mat) {
			int ncol = mat->noCols();
			int nrow = mat->noRows();
			array<double, 2>^ _mat = gcnew array<double, 2>(nrow, ncol);
			for (int i = 0; i < nrow; i++)
			{
				for (int j = 0; j < ncol; j++)
				{
					_mat[i, j] = mat[i, j];
				}
			}
			return _mat;
		};

		void Add(double factThis, MatrixWrapper^ other, double factOther) {
			this->_Matrix->addMatrix(factThis, *other->_Matrix, factOther);
			return;
		}


		void Print() {
			this->_Matrix->Print();
		}

	internal:
		Matrix * _Matrix;

		//inline static const Matrix& GetMatrix(MatrixWrapper^ mat) {
		//	// TODO :: test
		//	Matrix& v = *mat->_Matrix;
		//	return v;
		//};
		MatrixWrapper(Matrix* mat) {
			_Matrix = mat;
		};

		static MatrixWrapper^ GetMatrixWrapper(Matrix* mat) {
			// TODO :: test
			if (mat == 0) return nullptr;
			MatrixWrapper^ wvec = gcnew MatrixWrapper();
			wvec->_Matrix = mat;
			return wvec;
		};

		static MatrixWrapper^ GetMatrixWrapper(const Matrix &mat) {
			int ncol = mat.noCols();
			int nrow = mat.noRows();
			MatrixWrapper^ _mat = gcnew MatrixWrapper(nrow,ncol);
			for (int i = 0; i < nrow; i++)
			{
				for (int j = 0; j < ncol; j++)
				{
					_mat[i,j] = mat(i,j);
				}
			}
			return _mat;
		};

		static array<double,2>^ GetArray(const Matrix &mat) {
			int ncol = mat.noCols();
			int nrow = mat.noRows();
			array<double, 2>^ _mat = gcnew array<double, 2>(nrow, ncol);
			for (int i = 0; i < nrow; i++)
			{
				for (int j = 0; j < ncol; j++)
				{
					_mat[i, j] = mat(i, j);
				}
			}
			return _mat;
		};

		static array<double, 2>^ GetArray(Matrix *mat) {
			int ncol = mat->noCols();
			int nrow = mat->noRows();
			array<double, 2>^ _mat = gcnew array<double, 2>(nrow, ncol);
			for (int i = 0; i < nrow; i++)
			{
				for (int j = 0; j < ncol; j++)
				{
					_mat[i, j] = mat->operator()(i, j);
				}
			}
			return _mat;
		};


	private:

	};
}
