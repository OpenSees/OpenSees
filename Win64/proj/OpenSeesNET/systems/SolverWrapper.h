#pragma once
#include "../actors/IMovableObjectWrapper.h"

namespace OpenSees {
	namespace Systems {
		namespace Linears {
			public ref class SolverWrapper abstract : IMovableObjectWrapper
			{
			public:
				SolverWrapper();
				~SolverWrapper();
			internal:

			private:

			};
		}
	}
}