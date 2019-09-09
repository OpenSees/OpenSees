#include "stdafx.h"
#include <DifferenceAccelerator.h>
#include <DifferenceAccelerator2.h>
#include <KrylovAccelerator.h>
#include <KrylovAccelerator2.h>
#include <MillerAccelerator.h>
#include <PeriodicAccelerator.h>
#include <RaphsonAccelerator.h>
#include <SecantAccelerator1.h>
#include <SecantAccelerator2.h>
#include <SecantAccelerator3.h>
#include "AcceleratorWrapper.h"

using namespace OpenSees;
using namespace OpenSees::Algorithms::Accelerators;

AcceleratorWrapper::AcceleratorWrapper() {
}

AcceleratorWrapper::~AcceleratorWrapper()
{
	
}

RaphsonAcceleratorWrapper::~RaphsonAcceleratorWrapper()
{
	if (_Accelerator != 0)
		delete _Accelerator;
}

//DifferenceAcceleratorWrapper::~DifferenceAcceleratorWrapper()
//{
//	if (_Accelerator != 0)
//		delete _Accelerator;
//}

DifferenceAccelerator2Wrapper::~DifferenceAccelerator2Wrapper()
{
	if (_Accelerator != 0)
		delete _Accelerator;
}

KrylovAcceleratorWrapper::~KrylovAcceleratorWrapper()
{
	if (_Accelerator != 0)
		delete _Accelerator;
}

KrylovAccelerator2Wrapper::~KrylovAccelerator2Wrapper()
{
	if (_Accelerator != 0)
		delete _Accelerator;
}

MillerAcceleratorWrapper::~MillerAcceleratorWrapper()
{
	if (_Accelerator != 0)
		delete _Accelerator;
}

PeriodicAcceleratorWrapper::~PeriodicAcceleratorWrapper()
{
	if (_Accelerator != 0)
		delete _Accelerator;
}


SecantAccelerator1Wrapper::~SecantAccelerator1Wrapper()
{
	if (_Accelerator != 0)
		delete _Accelerator;
}

SecantAccelerator2Wrapper::~SecantAccelerator2Wrapper()
{
	if (_Accelerator != 0)
		delete _Accelerator;
}

SecantAccelerator3Wrapper::~SecantAccelerator3Wrapper()
{
	if (_Accelerator != 0)
		delete _Accelerator;
}


//DifferenceAcceleratorWrapper::DifferenceAcceleratorWrapper(int maxDim, int tangent) {
//	//_Accelerator = new DifferenceAccelerator(maxDim, tangent);
//}

DifferenceAccelerator2Wrapper::DifferenceAccelerator2Wrapper(int maxDim, int tangent) {
	_Accelerator = new DifferenceAccelerator2(maxDim, tangent);
}

KrylovAcceleratorWrapper::KrylovAcceleratorWrapper(int maxDim, int tangent) {
	_Accelerator = new KrylovAccelerator(maxDim, tangent);
}

KrylovAccelerator2Wrapper::KrylovAccelerator2Wrapper(int maxDim, int tangent) {
	_Accelerator = new KrylovAccelerator2(maxDim, tangent);
}

MillerAcceleratorWrapper::MillerAcceleratorWrapper(int maxDim, double tol,	int tangent) {
	//_Accelerator = new MillerAccelerator(maxDim, tol, tangent);
}

PeriodicAcceleratorWrapper::PeriodicAcceleratorWrapper(int iter, int tangent) {
	_Accelerator = new PeriodicAccelerator(iter, tangent);
}

RaphsonAcceleratorWrapper::RaphsonAcceleratorWrapper( int tangent) {
	_Accelerator = new RaphsonAccelerator( tangent);
}

SecantAccelerator1Wrapper::SecantAccelerator1Wrapper(int maxIter, int tangent) {
	_Accelerator = new SecantAccelerator1(maxIter,tangent);
}

SecantAccelerator2Wrapper::SecantAccelerator2Wrapper(int maxIter, int tangent) {
	_Accelerator = new SecantAccelerator2(maxIter, tangent);
}

SecantAccelerator3Wrapper::SecantAccelerator3Wrapper(int maxIter, int tangent) {
	_Accelerator = new SecantAccelerator3(maxIter, tangent);
}
