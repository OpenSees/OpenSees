#pragma once

#include <DiagonalDirectSolver.h>
#include <DiagonalSOE.h>
#include <DiagonalSolver.h>
#include <DistributedDiagonalSOE.h>
#include <DistributedDiagonalSolver.h>

#include "../actors/IMovableObjectWrapper.h"
#include "LinearSOESolverWrapper.h"
#include "LinearSOEWrapper.h"

namespace OpenSees {
	namespace Systems {
		namespace Linears {

			public ref class DiagonalSolverWrapper abstract : LinearSOESolverWrapper
			{
			public:
				DiagonalSolverWrapper() {};
				~DiagonalSolverWrapper() {
					if (_DiagonalSolver != 0)
						delete _DiagonalSolver;
				};
			internal:
				DiagonalSolver * _DiagonalSolver;
			private:

			};

			public ref class DiagonalDirectSolverWrapper : DiagonalSolverWrapper
			{
			public:
				DiagonalDirectSolverWrapper();
				~DiagonalDirectSolverWrapper() {
					if (_DiagonalSolver != 0)
						delete _DiagonalSolver;
				};
			internal:

			private:

			};

			public ref class DiagonalSOEWrapper : LinearSOEWrapper
			{
			public:
				DiagonalSOEWrapper() {};
				DiagonalSOEWrapper(DiagonalSolverWrapper^ theSolver);
				~DiagonalSOEWrapper() {
					if (_LinearSOE != 0)
						delete _LinearSOE;
				};
			internal:

			private:

			};

			

			public ref class DistributedDiagonalSolverWrapper abstract : LinearSOESolverWrapper
			{
			public:
				DistributedDiagonalSolverWrapper();
				~DistributedDiagonalSolverWrapper() {
					if (_DistributedDiagonalSolver != 0)
						delete _DistributedDiagonalSolver;
				};
			internal:
				DistributedDiagonalSolver * _DistributedDiagonalSolver;
			private:

			};

			

			public ref class DistributedDiagonalSOEWrapper : DiagonalSOEWrapper
			{
			public:
				DistributedDiagonalSOEWrapper(DistributedDiagonalSolverWrapper^ theSolver);
				~DistributedDiagonalSOEWrapper() {
					if (_LinearSOE != 0)
						delete _LinearSOE;
				};
			internal:

			private:

			};
		}
	}
}
