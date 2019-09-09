#include "stdafx.h"
#include "Vector.h"
#include "VectorWrapper.h"

using namespace OpenSees;
OpenSees::VectorWrapper::VectorWrapper(array<double>^ data)
{
	int size = data->Length;
	double* _data = new double[size];
	for (int i = 0; i < size; i++)
		_data[i] = data[i];
	_Vector = new Vector(_data, size);
}

VectorWrapper::VectorWrapper(int size)
{
	_Vector = new Vector(size);
}

VectorWrapper::VectorWrapper()
{
	_Vector = 0;
}

VectorWrapper::~VectorWrapper()
{
	if (_Vector != 0)
		delete _Vector;
}


//double VectorWrapper::Get(int x)
//{
//	return _Vector->operator[](x);
//}
//
//void VectorWrapper::Set(int x, double value)
//{
//	(*_Vector)[x] = value;
//}



