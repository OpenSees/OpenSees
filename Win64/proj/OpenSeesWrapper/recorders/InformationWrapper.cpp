#include "stdafx.h"

#include "InformationWrapper.h"

using namespace OpenSees::Recorders;

InformationWrapper::InformationWrapper() {
	_Information = new Information();
}

InformationWrapper::InformationWrapper(int val) {
	_Information = new Information(val);
}

InformationWrapper::InformationWrapper(double val) {
	_Information = new Information(val);
}

InformationWrapper::InformationWrapper(IDWrapper^ val) {
	_Information = new Information(*val->_ID);
}

InformationWrapper::InformationWrapper(VectorWrapper^ val) {
	_Information = new Information(*val->_Vector);
}

InformationWrapper::InformationWrapper(MatrixWrapper^ val) {
	_Information = new Information(*val->_Matrix);
}

InformationWrapper::InformationWrapper(IDWrapper^ val1, VectorWrapper^ val2) {
	_Information = new Information(*val1->_ID, *val2->_Vector);
}





