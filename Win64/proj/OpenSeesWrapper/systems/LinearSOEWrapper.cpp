#include "stdafx.h"
#include "LinearSOEWrapper.h"

using namespace OpenSees::Systems::Linears;
LinearSOEWrapper::LinearSOEWrapper()
{
	
}

LinearSOEWrapper::~LinearSOEWrapper()
{
	if (_LinearSOE != 0)
		delete _LinearSOE;
}

