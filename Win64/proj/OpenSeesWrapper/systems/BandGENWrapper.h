#pragma once

#include <BandGenLinLapackSolver.h>
#include <BandGenLinSOE.h>
#include <BandGenLinSolver.h>
#include <DistributedBandGenLinSOE.h>

#include "../actors/IMovableObjectWrapper.h"
#include "LinearSOESolverWrapper.h"
#include "LinearSOEWrapper.h"

namespace OpenSees {
	namespace Systems {
		namespace Linears {

			public ref class BandGenLinSolverWrapper abstract : LinearSOESolverWrapper
			{
			public:
				BandGenLinSolverWrapper();
				~BandGenLinSolverWrapper() {
					if (_BandGenLinSolver != 0)
						delete _BandGenLinSolver;
				};
			internal:
				BandGenLinSolver * _BandGenLinSolver;
			private:

			};

			public ref class BandGenLinLapackSolverWrapper : BandGenLinSolverWrapper
			{
			public:
				BandGenLinLapackSolverWrapper();
				~BandGenLinLapackSolverWrapper() {
					if (_BandGenLinSolver != 0)
						delete _BandGenLinSolver;
				};
			internal:
				
			private:

			};

			public ref class BandGenLinSOEWrapper : LinearSOEWrapper
			{
			public:
				BandGenLinSOEWrapper() {};
				BandGenLinSOEWrapper(BandGenLinSolverWrapper^ theSolver);
				~BandGenLinSOEWrapper() {
					if (_LinearSOE != 0)
						delete _LinearSOE;
				};
			internal:

			private:

			};

			public ref class DistributedBandGenLinSOEWrapper : BandGenLinSOEWrapper
			{
			public:
				DistributedBandGenLinSOEWrapper(BandGenLinSolverWrapper^ theSolver);
				~DistributedBandGenLinSOEWrapper() {
					if (_LinearSOE != 0)
						delete _LinearSOE;
				};
			internal:

			private:

			};
		}
	}
}
