#include "stdafx.h"
#include "SparseSymWrapper.h"

using namespace OpenSees::Systems::Linears;

SymSparseLinSolverWrapper::SymSparseLinSolverWrapper() {
	_SymSparseLinSolver = new SymSparseLinSolver();
}

SymSparseLinSOEWrapper::SymSparseLinSOEWrapper(SymSparseLinSolverWrapper^ solver, int lSparse) {
	_LinearSOE = new SymSparseLinSOE(*solver->_SymSparseLinSolver,lSparse);
}
