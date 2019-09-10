#include "stdafx.h"
#include "DiagonalWrapper.h"

using namespace OpenSees::Systems::Linears;
DiagonalDirectSolverWrapper::DiagonalDirectSolverWrapper() {
	_DiagonalSolver = new DiagonalDirectSolver();
}

DiagonalSOEWrapper::DiagonalSOEWrapper(DiagonalSolverWrapper^ theSolver) {
	_LinearSOE = new DiagonalSOE(*theSolver->_DiagonalSolver);
}

DistributedDiagonalSolverWrapper::DistributedDiagonalSolverWrapper() {
	_DistributedDiagonalSolver = new DistributedDiagonalSolver();
}

DistributedDiagonalSOEWrapper::DistributedDiagonalSOEWrapper(DistributedDiagonalSolverWrapper^ theSolver){
	_LinearSOE = new DistributedDiagonalSOE(*theSolver->_DistributedDiagonalSolver);
}