#include "stdafx.h"
#include "BandSPDWrapper.h"

using namespace OpenSees::Systems::Linears;
BandSPDLinSolverWrapper::BandSPDLinSolverWrapper() {
}

BandSPDLinLapackSolverWrapper::BandSPDLinLapackSolverWrapper() {
	_BandSPDLinSolver = new BandSPDLinLapackSolver();
}

BandSPDLinSOEWrapper::BandSPDLinSOEWrapper(BandSPDLinSolverWrapper^ theSolver) {
	_LinearSOE = new BandSPDLinSOE(*theSolver->_BandSPDLinSolver);
}

DistributedBandSPDLinSOEWrapper::DistributedBandSPDLinSOEWrapper(BandSPDLinSolverWrapper^ theSolver){
	_LinearSOE = new DistributedBandSPDLinSOE(*theSolver->_BandSPDLinSolver);
}