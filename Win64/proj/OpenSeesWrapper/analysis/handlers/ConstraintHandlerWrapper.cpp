#include "stdafx.h"
#include "ConstraintHandlerWrapper.h"

using namespace OpenSees;
using namespace OpenSees::Handlers;
ConstraintHandlerWrapper::ConstraintHandlerWrapper()
{
}

ConstraintHandlerWrapper::~ConstraintHandlerWrapper()
{
}

PlainHandlerWrapper::PlainHandlerWrapper()
{
	_ConstraintHandler = new PlainHandler();
}

PlainHandlerWrapper::~PlainHandlerWrapper()
{
	if (_ConstraintHandler != 0)
		delete _ConstraintHandler;
}

TransformationConstraintHandlerWrapper::TransformationConstraintHandlerWrapper()
{
	_ConstraintHandler = new TransformationConstraintHandler();
}

TransformationConstraintHandlerWrapper::~TransformationConstraintHandlerWrapper()
{
	if (_ConstraintHandler != 0)
		delete _ConstraintHandler;
}

LagrangeConstraintHandlerWrapper::LagrangeConstraintHandlerWrapper(double alphaSP, double alphaMP)
{
	_ConstraintHandler = new LagrangeConstraintHandler(alphaSP, alphaMP);
}

LagrangeConstraintHandlerWrapper::~LagrangeConstraintHandlerWrapper()
{
	if (_ConstraintHandler != 0)
		delete _ConstraintHandler;
}

PenaltyConstraintHandlerWrapper::PenaltyConstraintHandlerWrapper(double alphaSP, double alphaMP)
{
	_ConstraintHandler = new PenaltyConstraintHandler(alphaSP, alphaMP);
}

PenaltyConstraintHandlerWrapper::~PenaltyConstraintHandlerWrapper()
{
	if (_ConstraintHandler != 0)
		delete _ConstraintHandler;
}