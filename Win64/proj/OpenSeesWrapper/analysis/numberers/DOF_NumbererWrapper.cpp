#include "stdafx.h"
#include <DOF_Numberer.h>
#include <GraphNumberer.h>
#include "DOF_NumbererWrapper.h"

using namespace OpenSees;
using namespace OpenSees::Numberers;
DOF_NumbererWrapper::DOF_NumbererWrapper(GraphNumbererWrapper^ graphNumberer)
{
	if (graphNumberer != nullptr)
	{
		GraphNumberer* _graphNumberer = graphNumberer->_GraphNumberer;
		_DOFNumberer = new DOF_Numberer(*_graphNumberer);
	}
}

DOF_NumbererWrapper::~DOF_NumbererWrapper()
{
	if (_DOFNumberer != 0)
		delete _DOFNumberer;
}


PlainNumbererWrapper::PlainNumbererWrapper()
	:DOF_NumbererWrapper(nullptr)
{
	_DOFNumberer = new PlainNumberer();
}

PlainNumbererWrapper::~PlainNumbererWrapper()
{
	if (_DOFNumberer != 0)
		delete _DOFNumberer;
}

