#include "stdafx.h"
#include <Linear.h>
#include <AcceleratedNewton.h>
#include <BFGS.h>
#include <Broyden.h>
#include <KrylovNewton.h>
#include <ModifiedNewton.h>
#include <NewtonHallM.h>
#include <NewtonLineSearch.h>
#include <PeriodicNewton.h>
#include <NewtonRaphson.h>
#include "AlgorithmsWrapper.h"

using namespace OpenSees;
using namespace OpenSees::Algorithms;

SolutionAlgorithmWrapper::SolutionAlgorithmWrapper()
{
}

EquiSolnAlgoWrapper::EquiSolnAlgoWrapper() {
}

EquiSolnAlgoWrapper::~EquiSolnAlgoWrapper()
{
	
}

LinearWrapper::LinearWrapper() {
	_EquiSolnAlgo = new Linear(0, 0);
}

LinearWrapper::LinearWrapper(int theTangent, int factorOnce) {
	_EquiSolnAlgo = new Linear(theTangent, factorOnce);
}

LinearWrapper::~LinearWrapper()
{
	if (_EquiSolnAlgo != 0)
		delete _EquiSolnAlgo;
}

NewtonRaphsonWrapper::NewtonRaphsonWrapper()
{
	_EquiSolnAlgo = new NewtonRaphson();
}

NewtonRaphsonWrapper::NewtonRaphsonWrapper(int theTangentToUse, double iFact, double cFact)
{
	_EquiSolnAlgo = new NewtonRaphson(theTangentToUse, iFact, cFact);
}

NewtonRaphsonWrapper::~NewtonRaphsonWrapper()
{
	if (_EquiSolnAlgo != 0)
		delete _EquiSolnAlgo;
}

AcceleratedNewtonWrapper::AcceleratedNewtonWrapper()
{
	_EquiSolnAlgo = new AcceleratedNewton();
}

AcceleratedNewtonWrapper::AcceleratedNewtonWrapper(int tangent)
{
	_EquiSolnAlgo = new AcceleratedNewton(tangent);
}

AcceleratedNewtonWrapper::AcceleratedNewtonWrapper(ConvergenceTestWrapper ^ theTest, AcceleratorWrapper ^ theAccel, int tangent)
{
	_EquiSolnAlgo = new AcceleratedNewton(*theTest->_ConvergenceTest, theAccel->_Accelerator, tangent);
}

AcceleratedNewtonWrapper::~AcceleratedNewtonWrapper()
{
	if (_EquiSolnAlgo != 0)
		delete _EquiSolnAlgo;
}

BFGSWrapper::BFGSWrapper()
{
	_EquiSolnAlgo = new BFGS();
}

BFGSWrapper::BFGSWrapper( int tangent, int n)
{
	_EquiSolnAlgo = new BFGS(tangent, n);
}

BFGSWrapper::BFGSWrapper(ConvergenceTestWrapper ^ theTest, int tangent, int n)
{
	_EquiSolnAlgo = new BFGS(*theTest->_ConvergenceTest, tangent, n);
}

BFGSWrapper::~BFGSWrapper()
{
	if (_EquiSolnAlgo != 0)
		delete _EquiSolnAlgo;
}

BroydenWrapper::BroydenWrapper()
{
	_EquiSolnAlgo = new Broyden();
}

BroydenWrapper::BroydenWrapper(int tangent, int n)
{
	_EquiSolnAlgo = new Broyden(tangent, n);
}

BroydenWrapper::BroydenWrapper(ConvergenceTestWrapper ^ theTest, int tangent, int n)
{
	_EquiSolnAlgo = new Broyden(*theTest->_ConvergenceTest, tangent, n);
}

BroydenWrapper::~BroydenWrapper()
{
	if (_EquiSolnAlgo != 0)
		delete _EquiSolnAlgo;
}

KrylovNewtonWrapper::KrylovNewtonWrapper()
{
	_EquiSolnAlgo = new KrylovNewton();
}

KrylovNewtonWrapper::KrylovNewtonWrapper(int tangent, int maxDim)
{
	_EquiSolnAlgo = new KrylovNewton(tangent, maxDim);
}

KrylovNewtonWrapper::KrylovNewtonWrapper(ConvergenceTestWrapper ^ theTest, int tangent, int maxDim)
{
	_EquiSolnAlgo = new KrylovNewton(*theTest->_ConvergenceTest, tangent, maxDim);
}

KrylovNewtonWrapper::~KrylovNewtonWrapper()
{
	if (_EquiSolnAlgo != 0)
		delete _EquiSolnAlgo;
}

ModifiedNewtonWrapper::ModifiedNewtonWrapper(int tangent)
{
	_EquiSolnAlgo = new ModifiedNewton(tangent);
}

OpenSees::Algorithms::ModifiedNewtonWrapper::ModifiedNewtonWrapper(int tangent, double iFact, double cFact)
{
	_EquiSolnAlgo = new ModifiedNewton(tangent, iFact, cFact);
}

OpenSees::Algorithms::ModifiedNewtonWrapper::ModifiedNewtonWrapper(ConvergenceTestWrapper ^ theTest, int tangent, double iFact, double cFact)
{
	_EquiSolnAlgo = new ModifiedNewton(*theTest->_ConvergenceTest, tangent, iFact, cFact);
}

ModifiedNewtonWrapper::~ModifiedNewtonWrapper()
{
	if (_EquiSolnAlgo != 0)
		delete _EquiSolnAlgo;
}

NewtonHallMWrapper::NewtonHallMWrapper()
{
	_EquiSolnAlgo = new NewtonHallM();
}

OpenSees::Algorithms::NewtonHallMWrapper::NewtonHallMWrapper(double initFactor, int method, double alpha, double c)
{
	_EquiSolnAlgo = new NewtonHallM(initFactor, method, alpha, c);
}

NewtonHallMWrapper::~NewtonHallMWrapper()
{
	if (_EquiSolnAlgo != 0)
		delete _EquiSolnAlgo;
}

NewtonLineSearchWrapper::NewtonLineSearchWrapper()
{
	_EquiSolnAlgo = new NewtonLineSearch();
}

NewtonLineSearchWrapper::NewtonLineSearchWrapper(ConvergenceTestWrapper ^ theTest, LineSearchWrapper ^ theLineSearch)
{
	_EquiSolnAlgo = new NewtonLineSearch(*theTest->_ConvergenceTest, theLineSearch->_LineSearch);
}

NewtonLineSearchWrapper::~NewtonLineSearchWrapper()
{
	if (_EquiSolnAlgo != 0)
		delete _EquiSolnAlgo;
}

PeriodicNewtonWrapper::PeriodicNewtonWrapper()
{
	_EquiSolnAlgo = new PeriodicNewton();
}

PeriodicNewtonWrapper::PeriodicNewtonWrapper(int tangent, int maxCount)
{
	_EquiSolnAlgo = new PeriodicNewton(tangent, maxCount);
}

PeriodicNewtonWrapper::PeriodicNewtonWrapper(ConvergenceTestWrapper ^ theTest, int tangent, int maxCount)
{
	_EquiSolnAlgo = new PeriodicNewton(*theTest->_ConvergenceTest, tangent, maxCount);
}

PeriodicNewtonWrapper::~PeriodicNewtonWrapper()
{
	if (_EquiSolnAlgo != 0)
		delete _EquiSolnAlgo;
}
