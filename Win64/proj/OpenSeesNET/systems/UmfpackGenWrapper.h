#pragma once


#include <UmfpackGenLinSOE.h>
#include <UmfpackGenLinSolver.h>

#include "LinearSOESolverWrapper.h"
#include "LinearSOEWrapper.h"

using namespace OpenSees;
namespace OpenSees {
	namespace Systems {
		namespace Linears {

			public ref class UmfpackGenLinSolverWrapper : LinearSOESolverWrapper
			{
			public:
				UmfpackGenLinSolverWrapper();
				~UmfpackGenLinSolverWrapper() {
					if (_UmfpackGenLinSolver != 0)
						delete _UmfpackGenLinSolver;
				};
			internal:
				UmfpackGenLinSolver * _UmfpackGenLinSolver;
			private:

			};

			public ref class UmfpackGenLinSOEWrapper : LinearSOEWrapper {
			public :
				UmfpackGenLinSOEWrapper(UmfpackGenLinSolverWrapper^ solver);
				~UmfpackGenLinSOEWrapper() {
					if (_LinearSOE != 0)
						delete _LinearSOE;
				}
			};
		}
	}
}
