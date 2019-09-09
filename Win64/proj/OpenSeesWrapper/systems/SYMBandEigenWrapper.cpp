#include "stdafx.h"
#include "SYMBandEigenWrapper.h"


using namespace OpenSees;
using namespace OpenSees::Systems::Eigens;

SymBandEigenSolverWrapper::SymBandEigenSolverWrapper() {
	_SymBandEigenSolver = new SymBandEigenSolver();
}

SymBandEigenSOEWrapper::SymBandEigenSOEWrapper(SymBandEigenSolverWrapper^ solver, AnalysisModelWrapper^ model) {
	_EigenSOE = new SymBandEigenSOE(*solver->_SymBandEigenSolver, *model->_AnalysisModel);
}

