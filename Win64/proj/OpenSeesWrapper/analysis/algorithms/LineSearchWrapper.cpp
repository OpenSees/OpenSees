#include "stdafx.h"
#include <InitialInterpolatedLineSearch.h>
#include <BisectionLineSearch.h>
#include <RegulaFalsiLineSearch.h>
#include <SecantLineSearch.h>

#include "LineSearchWrapper.h"

using namespace OpenSees;
using namespace OpenSees::Algorithms::LineSearchs;

LineSearchWrapper::LineSearchWrapper()
{
	
}

OpenSees::Algorithms::LineSearchs::LineSearchWrapper::~LineSearchWrapper()
{
	
}

BisectionLineSearchWrapper::BisectionLineSearchWrapper()
{
	_LineSearch = new BisectionLineSearch();
}

BisectionLineSearchWrapper::BisectionLineSearchWrapper(double tolerance, int maxIter, double minEta, double maxEta, int flag)
{
	_LineSearch = new BisectionLineSearch(tolerance, maxIter, minEta, maxEta, flag);
}

BisectionLineSearchWrapper::~BisectionLineSearchWrapper()
{
	if (_LineSearch != 0)
		delete _LineSearch;
}


InitialInterpolatedLineSearchWrapper::InitialInterpolatedLineSearchWrapper()
{
	_LineSearch = new InitialInterpolatedLineSearch();
}

InitialInterpolatedLineSearchWrapper::InitialInterpolatedLineSearchWrapper(double tolerance, int maxIter, double minEta, double maxEta, int flag)
{
	_LineSearch = new InitialInterpolatedLineSearch(tolerance, maxIter, minEta, maxEta, flag);
}

InitialInterpolatedLineSearchWrapper::~InitialInterpolatedLineSearchWrapper()
{
	if (_LineSearch != 0)
		delete _LineSearch;
}

SecantLineSearchWrapper::SecantLineSearchWrapper()
{
	_LineSearch = new SecantLineSearch();
}

SecantLineSearchWrapper::SecantLineSearchWrapper(double tolerance, int maxIter, double minEta, double maxEta, int flag)
{
	_LineSearch = new SecantLineSearch(tolerance, maxIter, minEta, maxEta, flag);
}

SecantLineSearchWrapper::~SecantLineSearchWrapper()
{
	if (_LineSearch != 0)
		delete _LineSearch;
}


RegulaFalsiLineSearchWrapper::RegulaFalsiLineSearchWrapper()
{
	_LineSearch = new RegulaFalsiLineSearch();
}

RegulaFalsiLineSearchWrapper::RegulaFalsiLineSearchWrapper(double tolerance, int maxIter, double minEta, double maxEta, int flag)
{
	_LineSearch = new RegulaFalsiLineSearch(tolerance, maxIter, minEta, maxEta, flag);
}

RegulaFalsiLineSearchWrapper::~RegulaFalsiLineSearchWrapper()
{
	if (_LineSearch != 0)
		delete _LineSearch;
}