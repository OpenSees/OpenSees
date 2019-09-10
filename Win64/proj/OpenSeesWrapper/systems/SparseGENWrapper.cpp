#include "stdafx.h"
#include "SparseGENWrapper.h"

using namespace OpenSees::Systems::Linears;

SuperLUWrapper::SuperLUWrapper(int permSpec,
	double drop_tol,
	int panelSize,
	int relax,
	char symmetric) {
	_SparseGenColLinSolver = new SuperLU(permSpec, drop_tol, panelSize, relax, symmetric);
}

SuperLUWrapper::SuperLUWrapper() {
	_SparseGenColLinSolver = new SuperLU();
}

SparseGenColLinSOEWrapper::SparseGenColLinSOEWrapper(SparseGenColLinSolverWrapper^ solver) {
	_LinearSOE = new SparseGenColLinSOE(*solver->_SparseGenColLinSolver);
}

SparseGenRowLinSOEWrapper::SparseGenRowLinSOEWrapper(SparseGenRowLinSolverWrapper^ solver) {
	_LinearSOE = new SparseGenRowLinSOE(*solver->_SparseGenRowLinSolver);
}

#if _PARDISO
PARDISOGenLinSolverWrapper::PARDISOGenLinSolverWrapper() {
	_PARDISOGenLinSolver = new PARDISOGenLinSolver();
}

PARDISOGenLinSOEWrapper::PARDISOGenLinSOEWrapper(PARDISOGenLinSolverWrapper^ solver) {
	_LinearSOE = new PARDISOGenLinSOE(*solver->_PARDISOGenLinSolver);
}

PARDISOSymLinSolverWrapper::PARDISOSymLinSolverWrapper() {
	_PARDISOSymLinSolver = new PARDISOSymLinSolver();
}

PARDISOSymLinSOEWrapper::PARDISOSymLinSOEWrapper(PARDISOSymLinSolverWrapper^ solver) {
	_LinearSOE = new PARDISOSymLinSOE(*solver->_PARDISOSymLinSolver);
}
#endif
