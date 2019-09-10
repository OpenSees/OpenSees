#pragma once

#include <FullGenEigenSOE.h>
#include <FullGenEigenSolver.h>
#include <FullGenLinLapackSolver.h>
#include <FullGenLinSOE.h>
#include <FullGenLinSolver.h>

#include "../actors/IMovableObjectWrapper.h"
#include "../analysis/models/AnalysisModelWrapper.h"

#include "LinearSOESolverWrapper.h"
#include "LinearSOEWrapper.h"
#include "EigenSolverWrapper.h"
#include "EigenSOEWrapper.h"

using namespace OpenSees;
namespace OpenSees {
	namespace Systems {
		namespace Eigens {
			public ref class FullGenEigenSolverWrapper : EigenSolverWrapper
			{
			public:
				FullGenEigenSolverWrapper();
				~FullGenEigenSolverWrapper() {
					if (_FullGenEigenSolver != 0)
						delete _FullGenEigenSolver;
				};
			internal:
				FullGenEigenSolver * _FullGenEigenSolver;
			private:

			};

			public ref class FullGenEigenSOEWrapper : EigenSOEWrapper
			{
			public:
				FullGenEigenSOEWrapper(FullGenEigenSolverWrapper^ theSolver, AnalysisModelWrapper^ model);
				~FullGenEigenSOEWrapper() {
					if (_EigenSOE != 0)
						delete _EigenSOE;
				};
			internal:
				
			private:

			};
		}
		namespace Linears {

			
			public ref class FullGenLinSolverWrapper abstract : LinearSOESolverWrapper
			{
			public:
				FullGenLinSolverWrapper() {};
				~FullGenLinSolverWrapper() {
					if (_FullGenLinSolver != 0)
						delete _FullGenLinSolver;
				};
			internal:
				FullGenLinSolver * _FullGenLinSolver;
			private:

			};

			public ref class FullGenLinLapackSolverWrapper : FullGenLinSolverWrapper
			{
			public:
				FullGenLinLapackSolverWrapper();
				~FullGenLinLapackSolverWrapper() {
					if (_FullGenLinSolver != 0)
						delete _FullGenLinSolver;
				};
			internal:
			private:

			};

			public ref class FullGenLinSOEWrapper : LinearSOEWrapper
			{
			public:
				FullGenLinSOEWrapper() {};
				FullGenLinSOEWrapper(FullGenLinSolverWrapper^ theSolver);
				FullGenLinSOEWrapper(int N,FullGenLinSolverWrapper^ theSolver);
				~FullGenLinSOEWrapper() {
					if (_LinearSOE != 0)
						delete _LinearSOE;
				};
			internal:

			private:

			};
		
		}
	}
}
