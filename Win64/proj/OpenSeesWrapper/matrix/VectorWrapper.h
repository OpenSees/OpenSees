#pragma once
#include <Vector.h>

namespace OpenSees {

	public ref class VectorWrapper
	{
	public:
		VectorWrapper(array<double>^ data);
		VectorWrapper();
		VectorWrapper(int size);
		array<double>^ ToArray() {
			int size = _Vector->Size();
			array<double>^ data = gcnew array<double>(size);
			for (int i = 0; i < size; i++)
				data[i] = _Vector->operator[](i);
			return data;
		}

		int Size()
		{
			return _Vector->Size();
		}

		double Get(int index) {
			return _Vector->operator[](index);
		};


		void Set(int index, double value)
		{
			_Vector->operator()(index) = value;
		};

		property double default[int]
		{
			double get(int index) { return _Vector->operator[](index); }
			void set(int index, double value) { _Vector->operator()(index) = value; }
		}

			void Zero()
		{
			_Vector->Zero();
		}
		~VectorWrapper();

		void Add(double factThis, VectorWrapper^ other, double factOther) {
			this->_Vector->addVector(factThis, *other->_Vector, factOther);
			return;
		}

		int SetData(array<double>^ data, double selfFact, double otherFact) {
			if (this->Size() != data->Length) return -1;
			for (int i = 0; i < data->Length; i++)
			{
				if (selfFact == 0)
				{
					if (otherFact == 1)
						this->_Vector->operator()(i) = data[i];
					else
						this->_Vector->operator()(i) = data[i] * otherFact;
				}
				else
				{
					if (otherFact == 1)
						this->_Vector->operator()(i) = this->_Vector->operator()(i)*selfFact + data[i];
					else
						this->_Vector->operator()(i) = this->_Vector->operator()(i)*selfFact + data[i] * otherFact;
				}
			}
				
			return 0;
		}

		static VectorWrapper^ operator += (VectorWrapper^ vec1, VectorWrapper^ vec2) {
			if (vec1->Size() != vec2->Size())
			{
				System::Console::WriteLine("incompatible vectors ");
				return vec1;
			}

			(*vec1->_Vector) += *vec2->_Vector;
			return vec1;
		}

		static VectorWrapper^ operator -= (VectorWrapper^ vec1, VectorWrapper^ vec2) {
			if (vec1->Size() != vec2->Size())
			{
				System::Console::WriteLine("incompatible vectors ");
				return vec1;
			}

			(*vec1->_Vector) -= *vec2->_Vector;
			return vec1;
		}

		static VectorWrapper^ operator * (VectorWrapper^ vec1, double fact) {
			
			(*vec1->_Vector) = (*vec1->_Vector)*fact;
			return vec1;
		}

		static array<double>^ GetArray(VectorWrapper^ vec) {
			int sz = vec->Size();
			array<double>^ _vec = gcnew array<double>(sz);
			for (int i = 0; i < sz; i++)
			{
				_vec[i] = vec[i];
			}
			return _vec;
		};

		void Print() {
			this->_Vector->Print();
		}

	internal:
		Vector * _Vector;

		VectorWrapper(Vector * vec) {
			_Vector = vec;
		};

		VectorWrapper(double* data, int size) {
			_Vector = new Vector(data, size);
		};

		static VectorWrapper^ GetVectorWrapper(Vector* vec) {
			if (vec == 0) return nullptr;
			int length = vec->Size();
			if (length < 1)return nullptr;

			VectorWrapper^ wvec = gcnew VectorWrapper(length);
			for (int i = 0; i < length; i++)
				wvec[i] = vec->operator()(i);
			return wvec;
		};

		static VectorWrapper^ GetVectorWrapper(const Vector &vec) {
			int sz = vec.Size();
			VectorWrapper^ v = gcnew VectorWrapper(sz);
			for (int i = 0; i < sz; i++)
			{
				v[i] = vec(i);
			}
			return v;
		};

		static array<double>^ GetArray(const Vector &vec) {
			int sz = vec.Size();
			array<double>^ _vec = gcnew array<double>(sz);
			for (int i = 0; i < sz; i++)
			{
				_vec[i] = vec(i);
			}
			return _vec;
		};

		static array<double>^ GetArray(Vector *vec) {
			int sz = vec->Size();
			array<double>^ _vec = gcnew array<double>(sz);
			for (int i = 0; i < sz; i++)
			{
				_vec[i] = vec->operator()(i);
			}
			return _vec;
		};

		/*static const Vector& GetVector(VectorWrapper^ vec) {
			Vector& v = *vec->_Vector;
			return v;
		};*/

	private:

	};
}