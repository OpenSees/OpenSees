#pragma once
#include <OPS_Globals.h>
#include <TaggedObject.h>
#include <Domain.h>

namespace OpenSees {
	public ref class TaggedObjectWrapper abstract
	{
	public:
		TaggedObjectWrapper();
		~TaggedObjectWrapper();
		virtual int GetTag();
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