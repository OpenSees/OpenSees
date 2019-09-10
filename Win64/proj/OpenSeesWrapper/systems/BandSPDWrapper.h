#pragma once

#include <BandSPDLinLapackSolver.h>
#include <BandSPDLinSOE.h>
#include <BandSPDLinSolver.h>
#include <DistributedBandSPDLinSOE.h>

#include "../actors/IMovableObjectWrapper.h"
#include "LinearSOESolverWrapper.h"
#include "LinearSOEWrapper.h"

namespace OpenSees {
	namespace Systems {
		namespace Linears {

			public ref class BandSPDLinSolverWrapper abstract : LinearSOESolverWrapper
			{
			public:
				BandSPDLinSolverWrapper();
				~BandSPDLinSolverWrapper() {
					if (_BandSPDLinSolver != 0)
						delete _BandSPDLinSolver;
				};
			internal:
				BandSPDLinSolver * _BandSPDLinSolver;
			private:

			};

			public ref class BandSPDLinLapackSolverWrapper : BandSPDLinSolverWrapper
			{
			public:
				BandSPDLinLapackSolverWrapper();
				~BandSPDLinLapackSolverWrapper() {
					if (_BandSPDLinSolver != 0)
						delete _BandSPDLinSolver;
				};
			internal:
				
			private:

			};

			public ref class BandSPDLinSOEWrapper : LinearSOEWrapper
			{
			public:
				BandSPDLinSOEWrapper() {};
				BandSPDLinSOEWrapper(BandSPDLinSolverWrapper^ theSolver);
				~BandSPDLinSOEWrapper() {
					if (_LinearSOE != 0)
						delete _LinearSOE;
				};
			internal:

			private:

			};

			public ref class DistributedBandSPDLinSOEWrapper : BandSPDLinSOEWrapper
			{
			public:
				DistributedBandSPDLinSOEWrapper(BandSPDLinSolverWrapper^ theSolver);
				~DistributedBandSPDLinSOEWrapper() {
					if (_LinearSOE != 0)
						delete _LinearSOE;
				};
			internal:

			private:

			};
		}
	}
}
