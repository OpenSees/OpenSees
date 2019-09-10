#include "stdafx.h"
#include "AnalysisWrapper.h"

using namespace OpenSees;
using namespace OpenSees::Analysis;
using namespace OpenSees::Algorithms;
using namespace OpenSees::Handlers;
using namespace OpenSees::Integrators::Static;
using namespace OpenSees::Integrators::Transient;
using namespace OpenSees::Numberers;
using namespace OpenSees::ConvergenceTests;

AnalysisWrapper::AnalysisWrapper()
{
	
}

AnalysisWrapper::~AnalysisWrapper()
{
	
}


StaticAnalysisWrapper::StaticAnalysisWrapper(DomainWrapper^ theDomain,	ConstraintHandlerWrapper^ theHandler,	DOF_NumbererWrapper^ theNumberer,	AnalysisModelWrapper^ theModel,	EquiSolnAlgoWrapper^ theSolnAlgo,	LinearSOEWrapper^ theSOE,	StaticIntegratorWrapper^ theIntegrator,	ConvergenceTestWrapper^ theTest)
{
	ConvergenceTest* _convergenceTest = 0;
	if (theTest != nullptr)
		_convergenceTest = theTest->_ConvergenceTest;

	_StaticAnalysis = new StaticAnalysis(
		*theDomain->_Domain,
		*theHandler->_ConstraintHandler,
		*theNumberer->_DOFNumberer,
		*theModel->_AnalysisModel,
		*theSolnAlgo->_EquiSolnAlgo,
		*theSOE->_LinearSOE,
		*theIntegrator->_StaticIntegrator,
		_convergenceTest);
}
int 
StaticAnalysisWrapper::Analyze(int steps)
{
	return _StaticAnalysis->analyze(steps);
}

StaticAnalysisWrapper::~StaticAnalysisWrapper()
{
	if (_StaticAnalysis != 0)
		delete _StaticAnalysis;//
}


TransientAnalysisWrapper::TransientAnalysisWrapper()
{

}

TransientAnalysisWrapper::~TransientAnalysisWrapper()
{

}



DirectIntegrationAnalysisWrapper::DirectIntegrationAnalysisWrapper(DomainWrapper ^ theDomain, ConstraintHandlerWrapper ^ theHandler,	DOF_NumbererWrapper ^ theNumberer, AnalysisModelWrapper ^ theModel, EquiSolnAlgoWrapper ^ theSolnAlgo,	LinearSOEWrapper ^ theSOE, TransientIntegratorWrapper ^ theIntegrator, ConvergenceTestWrapper ^ theTest)
{
	ConvergenceTest* _convergenceTest = 0;
	if (theTest != nullptr)
		_convergenceTest = theTest->_ConvergenceTest;

	_DirectIntegrationAnalysis = new DirectIntegrationAnalysis(
		*theDomain->_Domain,
		*theHandler->_ConstraintHandler,
		*theNumberer->_DOFNumberer,
		*theModel->_AnalysisModel,
		*theSolnAlgo->_EquiSolnAlgo,
		*theSOE->_LinearSOE,
		*theIntegrator->_TransientIntegrator,
		_convergenceTest);
}
int 
DirectIntegrationAnalysisWrapper::Analyze(int steps, double dt)
{
	return _DirectIntegrationAnalysis->analyze(steps, dt);
}

DirectIntegrationAnalysisWrapper::~DirectIntegrationAnalysisWrapper()
{
	if (_DirectIntegrationAnalysis != 0)
		delete _DirectIntegrationAnalysis;
}


PFEMAnalysisAnalysisWrapper::PFEMAnalysisAnalysisWrapper(DomainWrapper ^ theDomain, ConstraintHandlerWrapper ^ theHandler, DOF_NumbererWrapper ^ theNumberer, AnalysisModelWrapper ^ theModel, EquiSolnAlgoWrapper ^ theSolnAlgo, LinearSOEWrapper ^ theSOE, TransientIntegratorWrapper ^ theIntegrator, ConvergenceTestWrapper ^ theTest, double max, double min, double g, double r)
{
	ConvergenceTest* _convergenceTest = 0;
	if (theTest != nullptr)
		_convergenceTest = theTest->_ConvergenceTest;

	_PFEMAnalysis = new PFEMAnalysis(
		*theDomain->_Domain,
		*theHandler->_ConstraintHandler,
		*theNumberer->_DOFNumberer,
		*theModel->_AnalysisModel,
		*theSolnAlgo->_EquiSolnAlgo,
		*theSOE->_LinearSOE,
		*theIntegrator->_TransientIntegrator,
		_convergenceTest, max, min, g, r);
}
int
PFEMAnalysisAnalysisWrapper::Analyze()
{
	return _PFEMAnalysis->analyze();
}

PFEMAnalysisAnalysisWrapper::~PFEMAnalysisAnalysisWrapper()
{
	if (_PFEMAnalysis != 0)
		delete _PFEMAnalysis;
}


VariableTimeStepDirectIntegrationAnalysisWrapper::VariableTimeStepDirectIntegrationAnalysisWrapper(DomainWrapper ^ theDomain, ConstraintHandlerWrapper ^ theHandler, DOF_NumbererWrapper ^ theNumberer, AnalysisModelWrapper ^ theModel, EquiSolnAlgoWrapper ^ theSolnAlgo, LinearSOEWrapper ^ theSOE, TransientIntegratorWrapper ^ theIntegrator, ConvergenceTestWrapper ^ theTest)
{
	ConvergenceTest* _convergenceTest = 0;
	if (theTest != nullptr)
		_convergenceTest = theTest->_ConvergenceTest;

	_VariableTimeStepDirectIntegrationAnalysis = new VariableTimeStepDirectIntegrationAnalysis(
		*theDomain->_Domain,
		*theHandler->_ConstraintHandler,
		*theNumberer->_DOFNumberer,
		*theModel->_AnalysisModel,
		*theSolnAlgo->_EquiSolnAlgo,
		*theSOE->_LinearSOE,
		*theIntegrator->_TransientIntegrator,
		_convergenceTest);
}
int
VariableTimeStepDirectIntegrationAnalysisWrapper::Analyze(int numSteps, double dT, double dtMin, double dtMax, int Jd)
{
	return _VariableTimeStepDirectIntegrationAnalysis->analyze(numSteps, dT, dtMin, dtMax, Jd);
}

VariableTimeStepDirectIntegrationAnalysisWrapper::~VariableTimeStepDirectIntegrationAnalysisWrapper()
{
	if (_VariableTimeStepDirectIntegrationAnalysis != 0)
		delete _VariableTimeStepDirectIntegrationAnalysis;
}