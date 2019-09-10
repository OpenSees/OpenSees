#pragma once


#include <SparseGenColLinSolver.h>
#include <SparseGenColLinSOE.h>

#include <SuperLU.h>

#if _PARDISO
#include <PARDISOGenLinSolver.h>
#include <PARDISOGenLinSOE.h>
#include <PARDISOSymLinSolver.h>
#include <PARDISOSymLinSOE.h>
#endif

#include <SparseGenRowLinSolver.h>
#include <SparseGenRowLinSOE.h>

#include "../actors/IMovableObjectWrapper.h"
#include "../analysis/models/AnalysisModelWrapper.h"

#include "LinearSOESolverWrapper.h"
#include "LinearSOEWrapper.h"

using namespace OpenSees;
namespace OpenSees {
	namespace Systems {
		namespace Linears {

			public ref class SparseGenColLinSolverWrapper abstract : LinearSOESolverWrapper
			{
			public:
				SparseGenColLinSolverWrapper() {};
				~SparseGenColLinSolverWrapper() {
					if (_SparseGenColLinSolver != 0)
						delete _SparseGenColLinSolver;
				};
			internal:
				SparseGenColLinSolver * _SparseGenColLinSolver;
			private:

			};

			public ref class SuperLUWrapper : SparseGenColLinSolverWrapper
			{
			public:
				SuperLUWrapper();
				SuperLUWrapper(int permSpec,
					double drop_tol,
					int panelSize,
					int relax,
					char symmetric) ;
				~SuperLUWrapper() {
					if (_SparseGenColLinSolver != 0)
						delete _SparseGenColLinSolver;
				};
			internal:
				
			private:


			};
			public ref class SparseGenRowLinSolverWrapper abstract : LinearSOESolverWrapper
			{
			public:
				SparseGenRowLinSolverWrapper() {};
				~SparseGenRowLinSolverWrapper() {
					if (_SparseGenRowLinSolver != 0)
						delete _SparseGenRowLinSolver;
				};
			internal:
				SparseGenRowLinSolver * _SparseGenRowLinSolver;
			private:

			};

			public ref class SparseGenColLinSOEWrapper : LinearSOEWrapper {
			public :
				SparseGenColLinSOEWrapper(SparseGenColLinSolverWrapper^ solver);
				~SparseGenColLinSOEWrapper() {
					if (_LinearSOE != 0)
						delete _LinearSOE;
				}
			};

			public ref class SparseGenRowLinSOEWrapper : LinearSOEWrapper {
			public:
				SparseGenRowLinSOEWrapper(SparseGenRowLinSolverWrapper^ solver);
				~SparseGenRowLinSOEWrapper() {
					if (_LinearSOE != 0)
						delete _LinearSOE;
				}
			};

#if _PARDISO
			public ref class PARDISOGenLinSolverWrapper : LinearSOESolverWrapper
			{
			public:
				PARDISOGenLinSolverWrapper();
				~PARDISOGenLinSolverWrapper() {
					if (_PARDISOGenLinSolver != 0)
						delete _PARDISOGenLinSolver;
				};
			internal:
				PARDISOGenLinSolver * _PARDISOGenLinSolver;
			private:

			};

			public ref class PARDISOGenLinSOEWrapper : LinearSOEWrapper {
			public:
				PARDISOGenLinSOEWrapper(PARDISOGenLinSolverWrapper^ solver);
				~PARDISOGenLinSOEWrapper() {
					if (_LinearSOE != 0)
						delete _LinearSOE;
				}
			};

			public ref class PARDISOSymLinSolverWrapper : LinearSOESolverWrapper
			{
			public:
				PARDISOSymLinSolverWrapper();
				~PARDISOSymLinSolverWrapper() {
					if (_PARDISOSymLinSolver != 0)
						delete _PARDISOSymLinSolver;
				};
			internal:
				PARDISOSymLinSolver * _PARDISOSymLinSolver;
			private:

			};

			public ref class PARDISOSymLinSOEWrapper : LinearSOEWrapper {
			public:
				PARDISOSymLinSOEWrapper(PARDISOSymLinSolverWrapper^ solver);
				~PARDISOSymLinSOEWrapper() {
					if (_LinearSOE != 0)
						delete _LinearSOE;
				}
			};
#endif

		}
	}
}
