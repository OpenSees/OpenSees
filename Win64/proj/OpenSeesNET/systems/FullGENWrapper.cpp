#include "stdafx.h"
#include "FullGENWrapper.h"

using namespace OpenSees::Systems::Linears; 
using namespace OpenSees::Systems::Eigens;

FullGenEigenSolverWrapper::FullGenEigenSolverWrapper() {
	_FullGenEigenSolver = new FullGenEigenSolver();
}

FullGenEigenSOEWrapper::FullGenEigenSOEWrapper(FullGenEigenSolverWrapper^ theSolver, AnalysisModelWrapper^ model) {
	_EigenSOE = new FullGenEigenSOE(*theSolver->_FullGenEigenSolver, *model->_AnalysisModel);
}


FullGenLinLapackSolverWrapper::FullGenLinLapackSolverWrapper() {
	_FullGenLinSolver = new FullGenLinLapackSolver();
}

FullGenLinSOEWrapper::FullGenLinSOEWrapper(FullGenLinSolverWrapper^ theSolver){
	_LinearSOE = new FullGenLinSOE(*theSolver->_FullGenLinSolver);
}

FullGenLinSOEWrapper::FullGenLinSOEWrapper(int N,FullGenLinSolverWrapper^ theSolver) {
	_LinearSOE = new FullGenLinSOE(N, *theSolver->_FullGenLinSolver);
}