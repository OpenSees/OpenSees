#include "stdafx.h"
#include "BandGENWrapper.h"

using namespace OpenSees::Systems::Linears;
BandGenLinSolverWrapper::BandGenLinSolverWrapper()
{

}

BandGenLinLapackSolverWrapper::BandGenLinLapackSolverWrapper() {
	_BandGenLinSolver = new BandGenLinLapackSolver();
}

BandGenLinSOEWrapper::BandGenLinSOEWrapper(BandGenLinSolverWrapper^ theSolver) {
	_LinearSOE = new BandGenLinSOE(*theSolver->_BandGenLinSolver);
}

DistributedBandGenLinSOEWrapper::DistributedBandGenLinSOEWrapper(BandGenLinSolverWrapper^ theSolver){
	_LinearSOE = new DistributedBandGenLinSOE(*theSolver->_BandGenLinSolver);
}