#pragma once
#include <Element.h>
#include "../domains/components/DomainComponentWrapper.h"

using namespace OpenSees;
using namespace System;

namespace OpenSees {
	namespace Elements {

		public ref class BaseElementWrapper abstract : DomainComponentWrapper
		{
		public:
			BaseElementWrapper() {};
			~BaseElementWrapper() {
				if (_Element != 0)
					delete _Element;
			};
		internal:
			Element* _Element;
		private:

		};
	}
}
