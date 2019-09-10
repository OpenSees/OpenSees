#include "stdafx.h"
#include "EigenSOEWrapper.h"

using namespace OpenSees::Systems::Eigens;
EigenSOEWrapper::EigenSOEWrapper()
{

}

EigenSOEWrapper::~EigenSOEWrapper()
{
	if (_EigenSOE != 0)
		delete _EigenSOE;
}

