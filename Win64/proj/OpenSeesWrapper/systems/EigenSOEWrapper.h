#pragma once
#include "../actors/IMovableObjectWrapper.h"
#include <EigenSOE.h>

namespace OpenSees {
	namespace Systems {
		namespace Eigens {
			public ref class EigenSOEWrapper abstract : IMovableObjectWrapper
			{
			public:
				EigenSOEWrapper();
				~EigenSOEWrapper();
			internal:
				EigenSOE * _EigenSOE;
			private:

			};
		}
	}
}
