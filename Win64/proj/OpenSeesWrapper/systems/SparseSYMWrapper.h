#pragma once


#include <SymSparseLinSOE.h>
#include <SymSparseLinSolver.h>

#include "LinearSOESolverWrapper.h"
#include "LinearSOEWrapper.h"

using namespace OpenSees;
namespace OpenSees {
	namespace Systems {
		namespace Linears {

			public ref class SymSparseLinSolverWrapper : LinearSOESolverWrapper
			{
			public:
				SymSparseLinSolverWrapper();
				~SymSparseLinSolverWrapper() {
					if (_SymSparseLinSolver != 0)
						delete _SymSparseLinSolver;
				};
			internal:
				SymSparseLinSolver * _SymSparseLinSolver;
			private:

			};

			public ref class SymSparseLinSOEWrapper : LinearSOEWrapper {
			public :
				SymSparseLinSOEWrapper(SymSparseLinSolverWrapper^ solver, int lSparse);
				~SymSparseLinSOEWrapper() {
					if (_LinearSOE != 0)
						delete _LinearSOE;
				}
			};
		}
	}
}
