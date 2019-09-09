#include "stdafx.h"
#include "StaticIntegratorWrapper.h"

using namespace OpenSees;
using namespace OpenSees::Integrators;
using namespace OpenSees::Integrators::Static;
StaticIntegratorWrapper::StaticIntegratorWrapper()
{
	
}

StaticIntegratorWrapper::~StaticIntegratorWrapper()
{
	
}

LoadControlWrapper::LoadControlWrapper(double deltaLambda, int numIncr, double minLambda, double maxlambda)
{
	_StaticIntegrator = new LoadControl(deltaLambda, numIncr, minLambda, maxlambda);
}

LoadControlWrapper::~LoadControlWrapper()
{
	if (_StaticIntegrator != 0)
		delete _StaticIntegrator;
}

DisplacementControlWrapper::DisplacementControlWrapper(int node, int dof, double increment, DomainWrapper^ domain, int numIncr, double min, double max)
{
	_StaticIntegrator = new DisplacementControl(node, dof, increment, domain->_Domain, numIncr, min, max);
}

DisplacementControlWrapper::~DisplacementControlWrapper()
{
	if (_StaticIntegrator != 0)
		delete _StaticIntegrator;
}

ArcLengthWrapper::ArcLengthWrapper(double arcLength, double alpha)
{
	_StaticIntegrator = new ArcLength(arcLength, alpha);
}

ArcLengthWrapper::~ArcLengthWrapper()
{
	if (_StaticIntegrator != 0)
		delete _StaticIntegrator;
}

ArcLength1Wrapper::ArcLength1Wrapper(double arcLength, double alpha)
{
	_StaticIntegrator = new ArcLength1(arcLength, alpha);
}

ArcLength1Wrapper::~ArcLength1Wrapper()
{
	if (_StaticIntegrator != 0)
		delete _StaticIntegrator;
}


HSConstraintWrapper::HSConstraintWrapper(double arcLength, double psi_u, double psi_f, double u_ref)
{
	_StaticIntegrator = new HSConstraint(arcLength, psi_u, psi_f, u_ref);
}

HSConstraintWrapper::~HSConstraintWrapper()
{
	if (_StaticIntegrator != 0)
		delete _StaticIntegrator;
}

LoadPathWrapper::LoadPathWrapper()
{
	_StaticIntegrator = new LoadPath();
}

LoadPathWrapper::LoadPathWrapper(VectorWrapper^ theLoadPath)
{
	_StaticIntegrator = new LoadPath(*theLoadPath->_Vector);
}

LoadPathWrapper::~LoadPathWrapper()
{
	if (_StaticIntegrator != 0)
		delete _StaticIntegrator;
}

MinUnbalDispNormWrapper::MinUnbalDispNormWrapper(double lambda1, int specNumIterStep, double dlambda1min, double dlambda1max, int signFirstStepMethod)
{
	_StaticIntegrator = new MinUnbalDispNorm(lambda1, specNumIterStep, dlambda1min, dlambda1max, signFirstStepMethod);
}

MinUnbalDispNormWrapper::~MinUnbalDispNormWrapper()
{
	if (_StaticIntegrator != 0)
		delete _StaticIntegrator;
}

QuadraticMethodWrapper::QuadraticMethodWrapper(double arcLength, int type)
{
	_StaticIntegrator = new EQPath(arcLength, type);
}

QuadraticMethodWrapper::~QuadraticMethodWrapper()
{
	if (_StaticIntegrator != 0)
		delete _StaticIntegrator;
}