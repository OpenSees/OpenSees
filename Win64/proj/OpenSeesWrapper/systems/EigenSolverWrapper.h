#pragma once
#include "../actors/IMovableObjectWrapper.h"
#include <EigenSolver.h>

namespace OpenSees {
	namespace Systems {
		namespace Eigens {
			public ref class EigenSolverWrapper abstract : IMovableObjectWrapper
			{
			public:
				EigenSolverWrapper();
				~EigenSolverWrapper();
			internal:

			private:

			};
		}
	}
}
