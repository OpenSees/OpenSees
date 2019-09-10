#pragma once

#include "../../taggeds/TaggedObjectWrapper.h"
#include "../../actors/IMovableObjectWrapper.h"

namespace OpenSees {

		public ref class DomainComponentWrapper abstract :TaggedObjectWrapper, IMovableObjectWrapper
		{
		public:
			DomainComponentWrapper();
			DomainComponentWrapper(int tag);
			~DomainComponentWrapper();
		internal:
			
		private:

		};
	
}
