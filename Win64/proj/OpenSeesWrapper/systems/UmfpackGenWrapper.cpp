#include "stdafx.h"
#include "UmfpackGenWrapper.h"

using namespace OpenSees::Systems::Linears;

UmfpackGenLinSolverWrapper::UmfpackGenLinSolverWrapper() {
	_UmfpackGenLinSolver = new UmfpackGenLinSolver();
}

UmfpackGenLinSOEWrapper::UmfpackGenLinSOEWrapper(UmfpackGenLinSolverWrapper^ solver) {
	_LinearSOE = new UmfpackGenLinSOE(*solver->_UmfpackGenLinSolver);
}
