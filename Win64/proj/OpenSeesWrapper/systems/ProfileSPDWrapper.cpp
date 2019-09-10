#include "stdafx.h"
#include "ProfileSPDWrapper.h"

using namespace OpenSees::Systems::Linears; 

SProfileSPDLinSolverWrapper::SProfileSPDLinSolverWrapper()
{
	_SProfileSPDLinSolver = new SProfileSPDLinSolver();
}

ProfileSPDLinDirectSolverWrapper::ProfileSPDLinDirectSolverWrapper()
{
	_ProfileSPDLinSolver = new ProfileSPDLinDirectSolver();
}

ProfileSPDLinSubstrSolverWrapper::ProfileSPDLinSubstrSolverWrapper()
{
	_ProfileSPDLinSolver = new ProfileSPDLinSubstrSolver();
}

ProfileSPDLinSOEWrapper::ProfileSPDLinSOEWrapper(ProfileSPDLinSolverWrapper^ solver)
{
	_LinearSOE = new ProfileSPDLinSOE(*solver->_ProfileSPDLinSolver);
}

DistributedProfileSPDLinSOEWrapper::DistributedProfileSPDLinSOEWrapper(ProfileSPDLinSolverWrapper^ solver)
{
	_LinearSOE = new DistributedProfileSPDLinSOE(*solver->_ProfileSPDLinSolver);
}


SProfileSPDLinSOEWrapper::SProfileSPDLinSOEWrapper(SProfileSPDLinSolverWrapper^ solver)
{
	_LinearSOE = new SProfileSPDLinSOE(*solver->_SProfileSPDLinSolver);
}