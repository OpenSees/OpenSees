#include "stdafx.h"
#include "TaggedObjectWrapper.h"

using namespace OpenSees;
TaggedObjectWrapper::TaggedObjectWrapper()
{
	
}

TaggedObjectWrapper::~TaggedObjectWrapper()
{
	if (_TaggedObject != 0)
		delete _TaggedObject;
}

int TaggedObjectWrapper::GetTag()
{
	return _TaggedObject->getTag();
}
