#include "stdafx.h"
#include "ArpackSOEWrapper.h"
#include <ArpackSOE.h>

using namespace OpenSees::Systems::Eigens;
ArpackSOEWrapper::ArpackSOEWrapper(double shift)
{
	this->_EigenSOE = new ArpackSOE(shift);
}

ArpackSOEWrapper::ArpackSOEWrapper()
{
	this->_EigenSOE = new ArpackSOE();
}

ArpackSOEWrapper::~ArpackSOEWrapper()
{
	if (_EigenSOE != 0)
		delete _EigenSOE;
}

