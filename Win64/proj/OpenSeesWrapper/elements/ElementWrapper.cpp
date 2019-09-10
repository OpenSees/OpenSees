#include "stdafx.h"
#include "ElementWrapper.h"
#include <InformationProxy.h>

using namespace OpenSees;
using namespace OpenSees::Elements;
ElementWrapper::ElementWrapper()
{
	//_TaggedObject = (TaggedObject *) this->_Element;
}

ElementWrapper::~ElementWrapper()
{
	if (_Element != 0)
		delete _Element;
}






