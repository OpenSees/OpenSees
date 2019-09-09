#include "stdafx.h"
#include "ConvergenceTestWrapper.h"

using namespace OpenSees;
using namespace OpenSees::ConvergenceTests;

ConvergenceTestWrapper::ConvergenceTestWrapper()
{
	
}

ConvergenceTestWrapper::~ConvergenceTestWrapper()
{
	if (_ConvergenceTest != 0)
		delete _ConvergenceTest;
}


CTestEnergyIncrWrapper::CTestEnergyIncrWrapper(double tol, int maxNumIter, int printFlag, int normType)
{
	_ConvergenceTest = new CTestEnergyIncr(tol, maxNumIter, printFlag, normType);
}

CTestEnergyIncrWrapper::~CTestEnergyIncrWrapper() {
	if (_ConvergenceTest != 0)
		delete _ConvergenceTest;
}


CTestFixedNumIterWrapper::CTestFixedNumIterWrapper(int maxNumIter, int printFlag, int normType)
{
	_ConvergenceTest = new CTestFixedNumIter(maxNumIter, printFlag, normType);
}

CTestFixedNumIterWrapper::~CTestFixedNumIterWrapper() {
	if (_ConvergenceTest != 0)
		delete _ConvergenceTest;
}


CTestNormDispIncrWrapper::CTestNormDispIncrWrapper(double tol, int maxNumIter, int printFlag, int normType, double maxTol)
{
	_ConvergenceTest = new CTestNormDispIncr(tol, maxNumIter, printFlag);
}

CTestNormDispIncrWrapper::~CTestNormDispIncrWrapper() {
	if (_ConvergenceTest != 0)
		delete _ConvergenceTest;
}

CTestNormUnbalanceWrapper::CTestNormUnbalanceWrapper(double tol, int maxNumIter, int printFlag, int normType, int maxincr, double maxTol)
{
	_ConvergenceTest = new CTestNormUnbalance(tol, maxNumIter, printFlag, normType, maxincr, maxTol);
}

CTestNormUnbalanceWrapper::~CTestNormUnbalanceWrapper() {
	if (_ConvergenceTest != 0)
		delete _ConvergenceTest;
}


CTestPFEMWrapper::CTestPFEMWrapper(double tv, double tp, double tv2, double tp2, double tvrel, double tprel,
	int maxNumIter, int maxincr,
	int printFlag, int normType)
{
	_ConvergenceTest = new CTestPFEM(tv, tp, tv2, tp2, tvrel, tprel, maxNumIter, maxincr, printFlag, normType);
}

CTestPFEMWrapper::~CTestPFEMWrapper() {
	if (_ConvergenceTest != 0)
		delete _ConvergenceTest;
}


CTestRelativeEnergyIncrWrapper::CTestRelativeEnergyIncrWrapper(double tol, int maxNumIter, int printFlag, int normType)
{
	_ConvergenceTest = new CTestRelativeEnergyIncr(tol, maxNumIter, printFlag, normType);
}

CTestRelativeEnergyIncrWrapper::~CTestRelativeEnergyIncrWrapper() {
	if (_ConvergenceTest != 0)
		delete _ConvergenceTest;
}


CTestRelativeNormDispIncrWrapper::CTestRelativeNormDispIncrWrapper(double tol, int maxNumIter, int printFlag, int normType)
{
	_ConvergenceTest = new CTestRelativeNormDispIncr(tol, maxNumIter, printFlag, normType);
}

CTestRelativeNormDispIncrWrapper::~CTestRelativeNormDispIncrWrapper() {
	if (_ConvergenceTest != 0)
		delete _ConvergenceTest;
}


CTestRelativeNormUnbalanceWrapper::CTestRelativeNormUnbalanceWrapper(double tol, int maxNumIter, int printFlag, int normType)
{
	_ConvergenceTest = new CTestRelativeNormUnbalance(tol, maxNumIter, printFlag, normType);
}

CTestRelativeNormUnbalanceWrapper::~CTestRelativeNormUnbalanceWrapper() {
	if (_ConvergenceTest != 0)
		delete _ConvergenceTest;
}


CTestRelativeTotalNormDispIncrWrapper::CTestRelativeTotalNormDispIncrWrapper(double tol, int maxNumIter, int printFlag, int normType)
{
	_ConvergenceTest = new CTestRelativeTotalNormDispIncr(tol, maxNumIter, printFlag, normType);
}

CTestRelativeTotalNormDispIncrWrapper::~CTestRelativeTotalNormDispIncrWrapper() {
	if (_ConvergenceTest != 0)
		delete _ConvergenceTest;
}


NormDispAndUnbalanceWrapper::NormDispAndUnbalanceWrapper(double tolDisp,
	double tolUnbalance,
	int maxNumIter,
	int printFlag,
	int normType, int maxincr)
{
	_ConvergenceTest = new NormDispAndUnbalance(tolDisp, tolUnbalance, maxNumIter, printFlag, normType, maxincr);
}

NormDispAndUnbalanceWrapper::~NormDispAndUnbalanceWrapper() {
	if (_ConvergenceTest != 0)
		delete _ConvergenceTest;
}



NormDispOrUnbalanceWrapper::NormDispOrUnbalanceWrapper(double tolDisp,
	double tolUnbalance,
	int maxNumIter,
	int printFlag,
	int normType, int maxincr)
{
	_ConvergenceTest = new NormDispOrUnbalance(tolDisp, tolUnbalance, maxNumIter, printFlag, normType, maxincr);
}

NormDispOrUnbalanceWrapper::~NormDispOrUnbalanceWrapper() {
	if (_ConvergenceTest != 0)
		delete _ConvergenceTest;
}




