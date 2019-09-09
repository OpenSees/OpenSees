#pragma once
#include <OPS_Globals.h>
#include <TaggedObject.h>

namespace OpenSees {
	public ref class TaggedObjectWrapper abstract
	{
	public:
		TaggedObjectWrapper();
		~TaggedObjectWrapper();
		int GetTag();
		void PrintSelf(int flag)
		{
			if (_TaggedObject != 0)
				_TaggedObject->Print(opserr, flag);
		}

	internal:
		TaggedObject * _TaggedObject;
	private:

	};
}