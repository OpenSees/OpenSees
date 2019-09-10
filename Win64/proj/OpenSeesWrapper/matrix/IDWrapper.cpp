#include "stdafx.h"
#include "IDWrapper.h"

using namespace OpenSees;
IDWrapper::IDWrapper(array<int>^ data)
{
	int size = data->Length;
	int* _data = new int[size];
	for (int i = 0; i < size; i++)
		_data[i] = data[i];
	_ID = new ID(_data, size);
}
IDWrapper::IDWrapper(int size)
{
	_ID = new ID(size);
}

IDWrapper::IDWrapper()
{
	_ID = 0;
}

IDWrapper::~IDWrapper()
{
	if (_ID != 0)
		delete _ID;
}





