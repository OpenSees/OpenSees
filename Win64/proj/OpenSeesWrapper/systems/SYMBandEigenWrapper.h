#pragma once


#include <SymBandEigenSOE.h>
#include <SymBandEigenSolver.h>

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

			public ref class SymBandEigenSolverWrapper : EigenSolverWrapper
			{
			public:
				SymBandEigenSolverWrapper();
				~SymBandEigenSolverWrapper() {
					if (_SymBandEigenSolver != 0)
						delete _SymBandEigenSolver;
				};
			internal:
				SymBandEigenSolver * _SymBandEigenSolver;
			private:

			};

			public ref class SymBandEigenSOEWrapper : EigenSOEWrapper {
			public :
				SymBandEigenSOEWrapper(SymBandEigenSolverWrapper^ solver, AnalysisModelWrapper^ model);
				~SymBandEigenSOEWrapper() {
					if (_EigenSOE != 0)
						delete _EigenSOE;
				}
			};
		}
	}
}
