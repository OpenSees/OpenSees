#include "stdafx.h"
#include "GraphNumbererWrapper.h"

using namespace OpenSees::GraphNumberers;
GraphNumbererWrapper::GraphNumbererWrapper()
{
	
}

GraphNumbererWrapper::~GraphNumbererWrapper()
{
	if (_GraphNumberer != 0)
		delete _GraphNumberer;
}

RCMWrapper::RCMWrapper(bool GPS)
{
	_GraphNumberer = new RCM(GPS);
}

RCMWrapper::~RCMWrapper()
{
	if (_GraphNumberer != 0)
		delete _GraphNumberer;
}



