#pragma once
#include "../taggeds/TaggedObjectWrapper.h"
#include "../actors/IMovableObjectWrapper.h"

namespace OpenSees {
	namespace Materials {
		public ref class MaterialWrapper abstract : TaggedObjectWrapper, IMovableObjectWrapper
		{
		public:
			MaterialWrapper();
			~MaterialWrapper() {

			};
		private:

		};
	}
}
